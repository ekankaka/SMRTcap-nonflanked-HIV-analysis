#!/usr/bin/env python3

from pathlib import Path
import re
import sys
import pandas as pd

if len(sys.argv) != 3:
    raise SystemExit("Usage: python3 combine_csv.py FINAL_RESULTS_DIR INPUT_DIR")

final_results = Path(sys.argv[1])
input_dir = Path(sys.argv[2])
master_path = input_dir / "masterfile.csv"

if not final_results.is_dir():
    raise SystemExit(f"Missing final-results folder: {final_results}")

input_dir.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------
# 1. Combine annotated CSV files.
# ------------------------------------------------------------
tables = []

for file in sorted(final_results.glob("*/*.annotated.csv")):
    match = re.fullmatch(r"(\d+)-(\d+)(?:-([^.]+))?\.annotated\.csv", file.name)
    if not match:
        print(f"Skipping unexpected filename: {file}")
        continue

    participant_id, visit_id, repeat_test = match.groups()
    sample_id = f"{participant_id}-{visit_id}"
    repeat_test = repeat_test or ""
    sample_run_id = sample_id if not repeat_test else f"{sample_id}-{repeat_test}"

    df = pd.read_csv(file, dtype=str, keep_default_na=False)

    required = ["READ", "STRAND", "HIV_SEQ", "INTEGRATION_SITE"]
    missing = [x for x in required if x not in df.columns]
    if missing:
        raise ValueError(f"Missing column(s) in {file}: {', '.join(missing)}")

    chromosome_col = "chromosome" if "chromosome" in df.columns else "CHROMOSOME" if "CHROMOSOME" in df.columns else None
    if chromosome_col is None:
        raise ValueError(f"Missing chromosome column in {file}")

    df["STRAND"] = df["STRAND"].str.strip().str.lower()
    if not df["STRAND"].isin(["plus", "minus"]).all():
        bad = sorted(df.loc[~df["STRAND"].isin(["plus", "minus"]), "STRAND"].unique())
        raise ValueError(f"Unexpected STRAND value(s) in {file}: {bad}")

    df["sample_id"] = sample_id
    df["sample_run_id"] = sample_run_id
    df["participant_id"] = participant_id
    df["visit_id"] = visit_id
    df["repeat_test"] = repeat_test
    tables.append(df)

if not tables:
    raise SystemExit("No annotated CSV files matched the expected filename pattern.")

master = pd.concat(tables, ignore_index=True)
first = ["sample_id", "sample_run_id", "participant_id", "visit_id", "repeat_test", "READ", "STRAND"]
master = master[first + [c for c in master.columns if c not in first]]

chromosome_col = "chromosome" if "chromosome" in master.columns else "CHROMOSOME"

# Participant-level clone identifier based on genomic integration coordinate.
# The source CLONE_ID is retained but is not used for clone resolution.
master["CLONE_ID2"] = (
    master[chromosome_col].str.strip() + "_" + master["INTEGRATION_SITE"].str.strip()
)

master.to_csv(master_path, index=False)

# ------------------------------------------------------------
# 2. Create participant FASTAs from HIV_SEQ.
#    chromosome == HIV defines a non-flanked read.
# ------------------------------------------------------------
def write_fasta(rows, path):
    with open(path, "w") as out:
        for _, row in rows.iterrows():
            seq = re.sub(r"\s+", "", row["HIV_SEQ"]).upper()
            if seq:
                prefix = row["sample_run_id"].replace("-", "_")
                out.write(f">{prefix}_{row['READ']}\n{seq}\n")

# Remove participant FASTAs from an earlier run.
for file in input_dir.glob("*_hiv.fasta"):
    file.unlink()
for file in input_dir.glob("*_hiv_non_flanked.fasta"):
    file.unlink()

for participant_id, df in master.groupby("participant_id", sort=True):
    df = df[df["HIV_SEQ"].str.strip() != ""].copy()
    df = df.drop_duplicates(["sample_run_id", "READ"], keep="first")

    write_fasta(df, input_dir / f"{participant_id}_hiv.fasta")

    nonflanked = df[df[chromosome_col].str.strip().str.upper() == "HIV"]
    write_fasta(nonflanked, input_dir / f"{participant_id}_hiv_non_flanked.fasta")

print(f"Combined {len(tables)} CSV files into {master_path}")
print(f"Rows from repeat tests: {(master['repeat_test'] != '').sum()}")
print(f"Non-flanked HIV reads (chromosome == HIV): {(master[chromosome_col].str.strip().str.upper() == 'HIV').sum()}")
print(f"Created participant HIV FASTAs in {input_dir}")
