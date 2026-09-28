#!/usr/bin/env python3

import argparse
import csv
import os
import re
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

PIPELINE_VERSION = "1.0.0"

# ------------------------------------------------------------
# HXB2 coordinates: 0-based, end-exclusive
# ------------------------------------------------------------
HXB2_LEN = 9719
LTR5 = (0, 634)
INTERNAL = (634, 9085)
LTR3 = (9085, 9719)
ANCHOR = 500

# Circular-architecture alignment thresholds.
MIN_IDENTITY = 0.75
MIN_ALN_BP = 50
STRONG_ANCHOR = 100
PUTATIVE_ANCHOR = 25
BOUNDARY_FLANK = 20

# Putative integration-site matching thresholds.
MIN_MATCH_BP = 100
MIN_MATCH_IDENTITY = 0.98
MIN_SHORTER_COVERAGE = 0.95

CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def raw_read(fasta_id):
    """Convert a FASTA ID such as 107_5_<READ> back to the READ in masterfile.csv."""
    fasta_id = re.sub(r"/[01]$", "", fasta_id)
    return re.sub(r"^\d+_\d+(?:_[A-Z][A-Z0-9]*)?_", "", fasta_id)


def overlap(a1, a2, b1, b2):
    return max(0, min(a2, b2) - max(a1, b1))


# ------------------------------------------------------------
# Read minimap2 PAF output.
# ------------------------------------------------------------
def read_paf(path):
    hits = defaultdict(list)
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        return hits

    with open(path) as f:
        for line in f:
            p = line.rstrip().split("\t")
            if len(p) < 12:
                continue

            tags = {}
            for item in p[12:]:
                tag = item.split(":", 2)
                if len(tag) == 3:
                    tags[tag[0]] = tag[2]

            aln_len = int(p[10])
            hits[p[0]].append({
                "qlen": int(p[1]),
                "qstart": int(p[2]),
                "qend": int(p[3]),
                "strand": p[4],
                "tname": p[5],
                "tlen": int(p[6]),
                "tstart": int(p[7]),
                "tend": int(p[8]),
                "aln_len": aln_len,
                "identity": int(p[9]) / aln_len if aln_len else 0,
                "tags": tags,
            })
    return hits


def good(hit):
    return hit["identity"] >= MIN_IDENTITY and hit["aln_len"] >= MIN_ALN_BP


def region_bp(hit, region):
    return overlap(hit["tstart"], hit["tend"], *region)


def crosses(hit, boundary, flank=BOUNDARY_FLANK):
    """Require aligned bases on both sides of a reference boundary."""
    cigar = hit["tags"].get("cg")
    if not cigar:
        return hit["tstart"] <= boundary - flank and hit["tend"] >= boundary + flank

    pos = hit["tstart"]
    left = right = 0

    for n, op in CIGAR.findall(cigar):
        n = int(n)
        if op in "M=X":
            left += overlap(pos, pos + n, boundary - flank, boundary)
            right += overlap(pos, pos + n, boundary, boundary + flank)
            pos += n
        elif op in "DN":
            pos += n

    return left >= flank and right >= flank


# ------------------------------------------------------------
# Report the active version and thresholds.
# ------------------------------------------------------------
def show_settings(_args):
    print(f"SMRTcap non-flanked HIV pipeline v{PIPELINE_VERSION}")
    print("Circularization thresholds:")
    print(f"  minimum identity: {100 * MIN_IDENTITY:.0f}%")
    print(f"  minimum alignment: {MIN_ALN_BP} bp")
    print(f"  strong anchor: {STRONG_ANCHOR} bp")
    print(f"  putative anchor: {PUTATIVE_ANCHOR} bp")
    print(f"  boundary flank: {BOUNDARY_FLANK} bp")
    print("Putative integration matching thresholds:")
    print(f"  minimum identity: {100 * MIN_MATCH_IDENTITY:.0f}%")
    print(f"  minimum shorter-fragment coverage: {100 * MIN_SHORTER_COVERAGE:.0f}%")
    print(f"  minimum alignment: {MIN_MATCH_BP} bp")


# ------------------------------------------------------------
# Build HXB2-derived references for 1-LTR and 2-LTR architecture.
# ------------------------------------------------------------
def make_refs(args):
    records = list(SeqIO.parse(args.hxb2, "fasta"))
    if len(records) != 1:
        raise ValueError("HXB2 FASTA must contain one sequence")

    hxb2 = str(records[0].seq).upper().replace("-", "")
    if len(hxb2) != HXB2_LEN:
        raise ValueError(f"Expected HXB2 length {HXB2_LEN}; found {len(hxb2)}")

    l5 = hxb2[slice(*LTR5)]
    l3 = hxb2[slice(*LTR3)]
    v5 = hxb2[INTERNAL[0]:INTERNAL[0] + ANCHOR]
    v3 = hxb2[INTERNAL[1] - ANCHOR:INTERNAL[1]]

    refs = [
        SeqRecord(Seq(v3 + l3 + l5 + v5), id="TWO_LTR", description=""),
        SeqRecord(Seq(v3 + l5 + v5), id="ONE_LTR_5", description=""),
        SeqRecord(Seq(v3 + l3 + v5), id="ONE_LTR_3", description=""),
    ]

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    SeqIO.write(refs, args.output, "fasta")


# ------------------------------------------------------------
# Read STRAND and integration metadata from masterfile.csv.
# ------------------------------------------------------------
def read_master(path):
    strands = {}
    integration_metadata = defaultdict(list)

    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))

    if not rows:
        return strands, integration_metadata

    fields = rows[0].keys()
    chromosome_col = "chromosome" if "chromosome" in fields else "CHROMOSOME" if "CHROMOSOME" in fields else None
    if chromosome_col is None:
        raise ValueError("chromosome column not found in masterfile.csv")

    for row in rows:
        participant = row.get("participant_id", "")
        read = raw_read(row.get("READ", ""))
        if not read:
            continue

        key = (participant, read)

        strand = row.get("STRAND", "").strip().lower()
        if strand in {"plus", "minus"}:
            old = strands.get(key)
            strands[key] = strand if not old or old == strand else "ambiguous"

        integration_metadata[key].append((
            row.get("CLONE_ID2", ""),
            row.get(chromosome_col, ""),
            row.get("INTEGRATION_SITE", ""),
        ))

    return strands, integration_metadata


# ------------------------------------------------------------
# Forward-normalize sequences using STRAND from masterfile.csv.
# ------------------------------------------------------------
def normalize(args):
    strands, _ = read_master(args.master)
    records = []

    for record in SeqIO.parse(args.sequences, "fasta"):
        strand = strands.get((args.participant_id, raw_read(record.id)), "")
        if strand not in {"plus", "minus"}:
            raise ValueError(f"No unambiguous STRAND found for {record.id}")

        seq = record.seq.reverse_complement() if strand == "minus" else record.seq
        records.append(SeqRecord(seq, id=record.id, description=""))

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    SeqIO.write(records, args.output, "fasta")


# ------------------------------------------------------------
# Split normalized all-HIV reads into normalized flanked and
# normalized non-flanked FASTAs using the known non-flanked IDs.
# ------------------------------------------------------------
def split_hiv(args):
    nonflanked_ids = {r.id for r in SeqIO.parse(args.nonflanked_ids, "fasta")}
    all_hiv = list(SeqIO.parse(args.all_hiv, "fasta"))
    all_ids = {r.id for r in all_hiv}

    missing = nonflanked_ids - all_ids
    if missing:
        raise ValueError(f"{len(missing)} non-flanked reads are absent from {args.all_hiv}")

    flanked = [r for r in all_hiv if r.id not in nonflanked_ids]
    nonflanked = [r for r in all_hiv if r.id in nonflanked_ids]

    SeqIO.write(flanked, args.flanked_output, "fasta")
    SeqIO.write(nonflanked, args.nonflanked_output, "fasta")


# ------------------------------------------------------------
# Circularization architecture.
# ------------------------------------------------------------
def ltr_only(hits):
    hits = [h for h in hits if good(h)]
    if not hits:
        return False

    ltr = max((max(region_bp(h, LTR5), region_bp(h, LTR3)) for h in hits), default=0)
    internal = max((region_bp(h, INTERNAL) for h in hits), default=0)
    return ltr >= PUTATIVE_ANCHOR and internal < PUTATIVE_ANCHOR


def circle_call(hits, is_ltr_only):
    two_strong = two_putative = one_strong = one_putative = weak_wrap = False
    ltr_len = LTR3[1] - LTR3[0]

    two_l3_end = ANCHOR + ltr_len
    two_l5_end = ANCHOR + 2 * ltr_len
    two_v5 = (two_l5_end, two_l5_end + ANCHOR)
    one_ltr_end = ANCHOR + ltr_len
    one_v5 = (one_ltr_end, one_ltr_end + ANCHOR)

    for hit in hits:
        if not good(hit):
            continue

        if hit["tname"] == "TWO_LTR":
            v3 = region_bp(hit, (0, ANCHOR))
            v5 = region_bp(hit, two_v5)
            junction = crosses(hit, two_l3_end)

            if junction and crosses(hit, ANCHOR) and crosses(hit, two_l5_end) \
                    and v3 >= STRONG_ANCHOR and v5 >= STRONG_ANCHOR:
                two_strong = True
            elif junction and (v3 >= PUTATIVE_ANCHOR or v5 >= PUTATIVE_ANCHOR):
                two_putative = True
            elif junction:
                weak_wrap = True

        elif hit["tname"] in {"ONE_LTR_5", "ONE_LTR_3"}:
            v3 = region_bp(hit, (0, ANCHOR))
            v5 = region_bp(hit, one_v5)
            architecture = crosses(hit, ANCHOR) and crosses(hit, one_ltr_end)

            if architecture and v3 >= STRONG_ANCHOR and v5 >= STRONG_ANCHOR:
                one_strong = True
            elif architecture and v3 >= PUTATIVE_ANCHOR and v5 >= PUTATIVE_ANCHOR:
                one_putative = True
            elif architecture:
                weak_wrap = True

    if two_strong:
        return "2-LTR circle — strong"
    if two_putative and (one_strong or one_putative):
        return "Circular HIV DNA — 1-LTR/2-LTR indeterminate"
    if two_putative:
        return "2-LTR circle — putative"
    if one_strong:
        return "1-LTR circle — strong"
    if one_putative:
        return "1-LTR circle — putative"
    if is_ltr_only:
        return "LTR 5'/3' ambiguous — no evidence of circularity"
    if weak_wrap:
        return "Circular HIV DNA — 1-LTR/2-LTR indeterminate"
    return "No circle-specific architecture detected"


# ------------------------------------------------------------
# Putative integration-site matching.
# ------------------------------------------------------------
def shorter_coverage(hit):
    """Fraction of the shorter HIV fragment covered by the alignment."""
    if hit["qlen"] <= hit["tlen"]:
        return (hit["qend"] - hit["qstart"]) / hit["qlen"]
    return (hit["tend"] - hit["tstart"]) / hit["tlen"]


def integration_assignment(hits, integration_metadata, participant_id):
    """Resolve qualifying forward-orientation matches to one integration clone/site."""
    qualifying = [
        hit for hit in hits
        if hit["strand"] == "+"
        and hit["aln_len"] >= MIN_MATCH_BP
        and hit["identity"] >= MIN_MATCH_IDENTITY
        and shorter_coverage(hit) >= MIN_SHORTER_COVERAGE
    ]

    if not qualifying:
        return "Origin indeterminate — no qualifying flanked-read match", "", "", "", "", "", "", ""

    # Keep the best qualifying alignment to each flanked READ.
    best = {}
    for hit in qualifying:
        read = raw_read(hit["tname"])
        score = (shorter_coverage(hit), hit["identity"], hit["aln_len"])
        if read not in best or score > best[read][0]:
            best[read] = (score, hit)

    matched_reads = sorted(best)
    identities = []
    coverages = []
    alignment_bp = []
    metadata = []

    for read in matched_reads:
        hit = best[read][1]
        identities.append(f"{100 * hit['identity']:.2f}")
        coverages.append(f"{100 * shorter_coverage(hit):.2f}")
        alignment_bp.append(str(hit["aln_len"]))

        values = set(integration_metadata.get((participant_id, read), []))
        if len(values) != 1:
            return (
                "Putative integration ambiguous — conflicting flanked-read metadata",
                ";".join(matched_reads),
                ";".join(identities),
                ";".join(coverages),
                ";".join(alignment_bp),
                "", "", "",
            )
        metadata.append(next(iter(values)))

    metrics = (
        ";".join(matched_reads),
        ";".join(identities),
        ";".join(coverages),
        ";".join(alignment_bp),
    )

    if len(metadata) == 1:
        clone, chromosome, site = metadata[0]
        return "Putatively integrated provirus", *metrics, clone, chromosome, site

    clones = {x[0] for x in metadata}
    chromosomes = {x[1] for x in metadata}
    sites = {x[2] for x in metadata}

    if len(clones) == 1 and "" not in clones and len(chromosomes) == 1 and len(sites) == 1:
        return (
            "Putatively integrated provirus",
            *metrics,
            next(iter(clones)),
            next(iter(chromosomes)),
            next(iter(sites)),
        )

    return "Putative integration ambiguous — multiple CLONE_ID2s", *metrics, "", "", ""


# ------------------------------------------------------------
# Classify forward-normalized non-flanked reads.
# ------------------------------------------------------------
def classify(args):
    linear = read_paf(args.linear_paf)
    circle = read_paf(args.circle_paf)
    flanked = read_paf(args.flanked_paf)
    _, integration_metadata = read_master(args.master)

    fields = [
        "participant_id",
        "sequence_id",
        "circle_annotation",
        "putative_integration_status",
        "matched_flanked_reads",
        "match_percent_identity",
        "match_shorter_fragment_coverage",
        "match_alignment_bp",
        "CLONE_ID_NEW",
        "CHROMOSOME_NEW",
        "INTEGRATION_SITE_NEW",
    ]

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    with open(args.output, "w", newline="") as f:
        out = csv.DictWriter(f, fieldnames=fields)
        out.writeheader()

        for record in SeqIO.parse(args.sequences, "fasta"):
            q = record.id
            circle_annotation = circle_call(circle[q], ltr_only(linear[q]))

            if circle_annotation == "No circle-specific architecture detected":
                status, matched, identity, coverage, aln_bp, clone, chromosome, site = integration_assignment(
                    flanked[q], integration_metadata, args.participant_id
                )
            else:
                status = ""
                matched = identity = coverage = aln_bp = ""
                clone = chromosome = site = ""

            out.writerow({
                "participant_id": args.participant_id,
                "sequence_id": raw_read(q),
                "circle_annotation": circle_annotation,
                "putative_integration_status": status,
                "matched_flanked_reads": matched,
                "match_percent_identity": identity,
                "match_shorter_fragment_coverage": coverage,
                "match_alignment_bp": aln_bp,
                "CLONE_ID_NEW": clone,
                "CHROMOSOME_NEW": chromosome,
                "INTEGRATION_SITE_NEW": site,
            })


# ------------------------------------------------------------
# Add non-flanked results to the combined masterfile.
# All derived circularization/matching columns stay blank for flanked reads.
# ------------------------------------------------------------
def update_master(args):
    with open(args.master, newline="") as f:
        rows = list(csv.DictReader(f))
        fields = list(rows[0].keys()) if rows else []

    chromosome_col = "chromosome" if "chromosome" in fields else "CHROMOSOME" if "CHROMOSOME" in fields else None
    if chromosome_col is None:
        raise ValueError("chromosome column not found in masterfile.csv")

    calls = {}
    if os.path.isfile(args.classification) and os.path.getsize(args.classification):
        with open(args.classification, newline="") as f:
            for row in csv.DictReader(f):
                calls[(row["participant_id"], raw_read(row["sequence_id"]))] = row

    new_fields = [
        "CIRCULARIZATION_ARCHITECTURE",
        "PUTATIVE_INTEGRATION_STATUS",
        "CHROMOSOME_NEW",
        "INTEGRATION_SITE_NEW",
        "CLONE_ID_NEW",
        "MATCHED_FLANKED_READS",
        "MATCH_PERCENT_IDENTITY",
        "MATCH_SHORTER_FRAGMENT_COVERAGE",
        "MATCH_ALIGNMENT_BP",
    ]

    for row in rows:
        # Derived fields are populated only for non-flanked reads.
        for field in new_fields:
            row[field] = ""

        if row.get(chromosome_col, "").strip().upper() != "HIV":
            continue

        call = calls.get((row.get("participant_id", ""), raw_read(row.get("READ", ""))))
        if not call:
            row["CIRCULARIZATION_ARCHITECTURE"] = "Not evaluated"
            row["PUTATIVE_INTEGRATION_STATUS"] = "Not evaluated"
            continue

        circle = call["circle_annotation"]
        row["CIRCULARIZATION_ARCHITECTURE"] = circle

        # This status is blank for circular/LTR-architecture reads.
        row["PUTATIVE_INTEGRATION_STATUS"] = call.get("putative_integration_status", "")

        row["MATCHED_FLANKED_READS"] = call.get("matched_flanked_reads", "")
        row["MATCH_PERCENT_IDENTITY"] = call.get("match_percent_identity", "")
        row["MATCH_SHORTER_FRAGMENT_COVERAGE"] = call.get("match_shorter_fragment_coverage", "")
        row["MATCH_ALIGNMENT_BP"] = call.get("match_alignment_bp", "")

        resolved = call.get("putative_integration_status") == "Putatively integrated provirus"

        if resolved:
            row["CHROMOSOME_NEW"] = call["CHROMOSOME_NEW"]
            row["INTEGRATION_SITE_NEW"] = call["INTEGRATION_SITE_NEW"]
            row["CLONE_ID_NEW"] = call["CLONE_ID_NEW"]

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    with open(args.output, "w", newline="") as f:
        out = csv.DictWriter(f, fieldnames=fields + [x for x in new_fields if x not in fields])
        out.writeheader()
        out.writerows(rows)


# ------------------------------------------------------------
# Command-line interface.
# ------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser()
    commands = parser.add_subparsers(dest="command", required=True)

    p = commands.add_parser("settings")
    p.set_defaults(func=show_settings)

    p = commands.add_parser("make-refs")
    p.add_argument("--hxb2", required=True)
    p.add_argument("--output", required=True)
    p.set_defaults(func=make_refs)

    p = commands.add_parser("normalize")
    p.add_argument("--participant-id", required=True)
    p.add_argument("--sequences", required=True)
    p.add_argument("--master", required=True)
    p.add_argument("--output", required=True)
    p.set_defaults(func=normalize)

    p = commands.add_parser("split-hiv")
    p.add_argument("--all-hiv", required=True)
    p.add_argument("--nonflanked-ids", required=True)
    p.add_argument("--flanked-output", required=True)
    p.add_argument("--nonflanked-output", required=True)
    p.set_defaults(func=split_hiv)

    p = commands.add_parser("classify")
    p.add_argument("--participant-id", required=True)
    p.add_argument("--sequences", required=True)
    p.add_argument("--linear-paf", required=True)
    p.add_argument("--circle-paf", required=True)
    p.add_argument("--flanked-paf", required=True)
    p.add_argument("--master", required=True)
    p.add_argument("--output", required=True)
    p.set_defaults(func=classify)

    p = commands.add_parser("update-master")
    p.add_argument("--master", required=True)
    p.add_argument("--classification", required=True)
    p.add_argument("--output", required=True)
    p.set_defaults(func=update_master)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
