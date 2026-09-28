#!/usr/bin/env python3

from pathlib import Path
import re
import sys
import pandas as pd

if len(sys.argv) not in {3, 4}:
    raise SystemExit(
        "Usage: python3 combine_csv.py FINAL_RESULTS_DIR PREPARED_DIR [SAMPLE_MAP_CSV]"
    )

final_results = Path(sys.argv[1])
prepared_dir = Path(sys.argv[2])
sample_map_path = Path(sys.argv[3]) if len(sys.argv) == 4 else None
master_path = prepared_dir / "masterfile.csv"

if not final_results.is_dir():
    raise SystemExit(f"Missing final-results folder: {final_results}")

prepared_dir.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------
# 1. Find sample-specific *.annotated.csv files.
#    sample_id is the filename before ".annotated.csv".
# ------------------------------------------------------------
files = sorted(final_results.glob("*/*.annotated.csv"))
if not files:
    raise SystemExit("No *.annotated.csv files found in sample subfolders.")

sample_files = {}
for file in files:
    sample_id = file.name[:-len(".annotated.csv")]
    if not sample_id:
        raise ValueError(f"Could not determine sample_id from {file}")
    if re.search(r"[\s,;|]", sample_id):
        raise ValueError(f"sample_id cannot contain whitespace, ',', ';' or '|': {sample_id}")
    if sample_id in sample_files:
        raise ValueError(f"Duplicate sample_id found: {sample_id}")
    sample_files[sample_id] = file

# ------------------------------------------------------------
# 2. Optional sample-to-participant mapping.
#    Without a mapping file, each sample is treated as one participant.
# ------------------------------------------------------------
if sample_map_path is not None:
    if not sample_map_path.is_file():
        raise SystemExit(f"Missing sample mapping file: {sample_map_path}")

    sample_map = pd.read_csv(sample_map_path, dtype=str, keep_default_na=False)
    required = ["sample_id", "participant_id"]
    missing = [x for x in required if x not in sample_map.columns]
    if missing:
        raise ValueError(
            f"Missing column(s) in {sample_map_path}: {', '.join(missing)}"
        )

    sample_map["sample_id"] = sample_map["sample_id"].str.strip()
    sample_map["participant_id"] = sample_map["participant_id"].str.strip()

    if (sample_map["sample_id"] == "").any() or (sample_map["participant_id"] == "").any():
        raise ValueError("sample_id and participant_id cannot be blank in the mapping file")
    if sample_map["sample_id"].duplicated().any():
        duplicates = sorted(sample_map.loc[sample_map["sample_id"].duplicated(), "sample_id"].unique())
        raise ValueError(f"Duplicate sample_id(s) in mapping file: {', '.join(duplicates)}")
    if sample_map["sample_id"].str.contains(r"[\s,;|]", regex=True).any():
        raise ValueError("sample_id cannot contain whitespace, ',', ';' or '|' in the mapping file")
    if sample_map["participant_id"].str.contains(r"[,\\/|]", regex=True).any():
        raise ValueError("participant_id cannot contain ',', '/', '\\', or '|'")

    mapping = dict(zip(sample_map["sample_id"], sample_map["participant_id"]))
    missing_samples = sorted(set(sample_files) - set(mapping))
    if missing_samples:
        raise ValueError(
            "Sample(s) missing from mapping file: " + ", ".join(missing_samples)
        )

    extra_samples = sorted(set(mapping) - set(sample_files))
    if extra_samples:
        print("Warning: mapping rows not present in input: " + ", ".join(extra_samples))
else:
    mapping = {sample_id: sample_id for sample_id in sample_files}

# ------------------------------------------------------------
# 3. Combine annotated CSV files.
# ------------------------------------------------------------
tables = []
for sample_id, file in sample_files.items():
    df = pd.read_csv(file, dtype=str, keep_default_na=False)

    required = ["READ", "STRAND", "HIV_SEQ", "INTEGRATION_SITE"]
    missing = [x for x in required if x not in df.columns]
    if missing:
        raise ValueError(f"Missing column(s) in {file}: {', '.join(missing)}")

    chromosome_col = (
        "chromosome" if "chromosome" in df.columns
        else "CHROMOSOME" if "CHROMOSOME" in df.columns
        else None
    )
    if chromosome_col is None:
        raise ValueError(f"Missing chromosome column in {file}")

    df["STRAND"] = df["STRAND"].str.strip().str.lower()
    if not df["STRAND"].isin(["plus", "minus"]).all():
        bad = sorted(df.loc[~df["STRAND"].isin(["plus", "minus"]), "STRAND"].unique())
        raise ValueError(f"Unexpected STRAND value(s) in {file}: {bad}")

    if (df["READ"].str.strip() == "").any():
        raise ValueError(f"Blank READ value found in {file}")
    if df["READ"].str.contains(r"[\s,;|]", regex=True).any():
        raise ValueError(f"READ values cannot contain whitespace, ',', ';' or '|' in {file}")

    df["sample_id"] = sample_id
    df["participant_id"] = mapping[sample_id]
    tables.append(df)

master = pd.concat(tables, ignore_index=True)
first = ["sample_id", "participant_id", "READ", "STRAND"]
master = master[first + [c for c in master.columns if c not in first]]

chromosome_col = "chromosome" if "chromosome" in master.columns else "CHROMOSOME"

# Participant-level clone identifier based on genomic integration coordinate.
# The source CLONE_ID is retained but is not used for clone resolution.
master["CLONE_ID2"] = (
    master[chromosome_col].str.strip() + "_" + master["INTEGRATION_SITE"].str.strip()
)

master.to_csv(master_path, index=False)

# ------------------------------------------------------------
# 4. Create participant FASTAs from HIV_SEQ.
#    FASTA IDs are sample_id|READ so reads remain unambiguous
#    when several samples belong to the same participant.
# ------------------------------------------------------------
def write_fasta(rows, path):
    with open(path, "w") as out:
        for _, row in rows.iterrows():
            seq = re.sub(r"\s+", "", row["HIV_SEQ"]).upper()
            if seq:
                out.write(f">{row['sample_id']}|{row['READ']}\n{seq}\n")

# Remove prepared participant FASTAs from an earlier run.
for file in prepared_dir.glob("*_hiv.fasta"):
    file.unlink()
for file in prepared_dir.glob("*_hiv_non_flanked.fasta"):
    file.unlink()

for participant_id, df in master.groupby("participant_id", sort=True):
    df = df[df["HIV_SEQ"].str.strip() != ""].copy()
    df = df.drop_duplicates(["sample_id", "READ"], keep="first")

    write_fasta(df, prepared_dir / f"{participant_id}_hiv.fasta")

    nonflanked = df[df[chromosome_col].str.strip().str.upper() == "HIV"]
    write_fasta(nonflanked, prepared_dir / f"{participant_id}_hiv_non_flanked.fasta")

print(f"Combined {len(tables)} sample CSV files into {master_path}")
if sample_map_path is None:
    print("No sample mapping supplied: each sample_id is treated as one participant")
else:
    print(f"Applied sample mapping: {sample_map_path}")
print(f"Participants: {master['participant_id'].nunique()}")
print(f"Non-flanked HIV reads (chromosome == HIV): {(master[chromosome_col].str.strip().str.upper() == 'HIV').sum()}")
print(f"Prepared participant FASTAs in {prepared_dir}")
