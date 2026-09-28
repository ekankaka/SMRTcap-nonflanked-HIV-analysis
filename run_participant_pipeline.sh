#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

PID="${1:?Usage: bash run_participant_pipeline.sh PARTICIPANT_ID HXB2 WORK_DIR RESULTS_DIR}"
HXB2="${2:?Missing HXB2 FASTA}"
WORK_DIR="${3:?Missing WORK_DIR}"
RESULTS_DIR="${4:?Missing RESULTS_DIR}"

PREPARED_DIR="$WORK_DIR/prepared"
PARTICIPANT_WORK="$WORK_DIR/participants/$PID"
PARTICIPANT_RESULTS="$RESULTS_DIR/participants/$PID"

ALL_HIV="$PREPARED_DIR/${PID}_hiv.fasta"
NONFLANKED_RAW="$PREPARED_DIR/${PID}_hiv_non_flanked.fasta"
MASTER="$PREPARED_DIR/masterfile.csv"
CIRCLE_REFS="$WORK_DIR/circle_refs.fasta"

ALL_HIV_NORM="$PARTICIPANT_WORK/${PID}_hiv_forward_normalized.fasta"
NONFLANKED_NORM="$PARTICIPANT_WORK/${PID}_hiv_non_flanked_forward_normalized.fasta"
FLANKED_NORM="$PARTICIPANT_WORK/${PID}_hiv_flanked_forward_normalized.fasta"

[[ -s "$ALL_HIV" ]] || { echo "Missing participant FASTA: $ALL_HIV" >&2; exit 1; }
[[ -s "$MASTER" ]] || { echo "Missing masterfile: $MASTER" >&2; exit 1; }
[[ -s "$HXB2" ]] || { echo "Missing HXB2 FASTA: $HXB2" >&2; exit 1; }
[[ -s "$CIRCLE_REFS" ]] || { echo "Missing circle reference: $CIRCLE_REFS" >&2; exit 1; }

mkdir -p "$PARTICIPANT_WORK" "$PARTICIPANT_RESULTS"

# 1. Forward-normalize every HIV read using STRAND from masterfile.csv:
#    plus = keep as-is; minus = reverse-complement.
python3 "$SCRIPT_DIR/classify_nonflanked_hiv.py" normalize \
    --participant-id "$PID" \
    --sequences "$ALL_HIV" \
    --master "$MASTER" \
    --output "$ALL_HIV_NORM"

# Stop if this participant has no non-flanked HIV reads.
if [[ ! -s "$NONFLANKED_RAW" ]] || ! grep -q '^>' "$NONFLANKED_RAW"; then
    echo "No non-flanked HIV reads for participant $PID; classification skipped."
    exit 0
fi

# 2. Split the normalized reads into flanked and non-flanked FASTAs.
python3 "$SCRIPT_DIR/classify_nonflanked_hiv.py" split-hiv \
    --all-hiv "$ALL_HIV_NORM" \
    --nonflanked-ids "$NONFLANKED_RAW" \
    --flanked-output "$FLANKED_NORM" \
    --nonflanked-output "$NONFLANKED_NORM"

# 3. Evaluate circular architecture in normalized non-flanked reads.
minimap2 -x asm20 -c --cs=long --secondary=yes -N 20 "$HXB2" "$NONFLANKED_NORM" \
    > "$PARTICIPANT_WORK/${PID}_vs_HXB2.paf"

minimap2 -x asm20 -c --cs=long --secondary=yes -N 20 "$CIRCLE_REFS" "$NONFLANKED_NORM" \
    > "$PARTICIPANT_WORK/${PID}_vs_circle_refs.paf"

# 4. Match normalized non-flanked reads to normalized flanked reads.
if [[ -s "$FLANKED_NORM" ]] && grep -q '^>' "$FLANKED_NORM"; then
    minimap2 -x asm20 -c --cs=long --secondary=yes -N 50 "$FLANKED_NORM" "$NONFLANKED_NORM" \
        > "$PARTICIPANT_WORK/${PID}_vs_flanked_HIV.paf"
else
    : > "$PARTICIPANT_WORK/${PID}_vs_flanked_HIV.paf"
fi

# 5. Assign a putative integration site only when no circle-specific
#    architecture is detected and qualifying matches resolve consistently.
python3 "$SCRIPT_DIR/classify_nonflanked_hiv.py" classify \
    --participant-id "$PID" \
    --sequences "$NONFLANKED_NORM" \
    --linear-paf "$PARTICIPANT_WORK/${PID}_vs_HXB2.paf" \
    --circle-paf "$PARTICIPANT_WORK/${PID}_vs_circle_refs.paf" \
    --flanked-paf "$PARTICIPANT_WORK/${PID}_vs_flanked_HIV.paf" \
    --master "$MASTER" \
    --output "$PARTICIPANT_RESULTS/${PID}_nonflanked_HIV_classification.csv"

echo "Finished participant $PID"
