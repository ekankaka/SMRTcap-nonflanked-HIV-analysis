#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

PID="${1:?Usage: bash run_participant_pipeline.sh PARTICIPANT_ID HXB2 INPUT_DIR OUTPUT_DIR}"
HXB2="${2:?Missing HXB2 FASTA}"
INPUT_DIR="${3:?Missing INPUT_DIR}"
OUTPUT_DIR="${4:?Missing OUTPUT_DIR}"

ALL_HIV="$INPUT_DIR/${PID}_hiv.fasta"
NONFLANKED_RAW="$INPUT_DIR/${PID}_hiv_non_flanked.fasta"
MASTER="$INPUT_DIR/masterfile.csv"
WORK="$OUTPUT_DIR/work/$PID"
OUT="$OUTPUT_DIR/$PID"
CIRCLE_REFS="$OUTPUT_DIR/work/circle_refs.fasta"

ALL_HIV_NORM="$OUT/${PID}_hiv_forward_normalized.fasta"
NONFLANKED_NORM="$OUT/${PID}_hiv_non_flanked_forward_normalized.fasta"
FLANKED_NORM="$WORK/${PID}_hiv_flanked_forward_normalized.fasta"

[[ -s "$ALL_HIV" ]] || { echo "Missing participant FASTA: $ALL_HIV" >&2; exit 1; }
[[ -s "$MASTER" ]] || { echo "Missing masterfile: $MASTER" >&2; exit 1; }
[[ -s "$HXB2" ]] || { echo "Missing HXB2 FASTA: $HXB2" >&2; exit 1; }
[[ -s "$CIRCLE_REFS" ]] || { echo "Missing circle reference: $CIRCLE_REFS" >&2; exit 1; }

mkdir -p "$WORK" "$OUT"

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
#    Both groups use the same forward-orientation convention.
python3 "$SCRIPT_DIR/classify_nonflanked_hiv.py" split-hiv \
    --all-hiv "$ALL_HIV_NORM" \
    --nonflanked-ids "$NONFLANKED_RAW" \
    --flanked-output "$FLANKED_NORM" \
    --nonflanked-output "$NONFLANKED_NORM"

# 3. Evaluate circular architecture in normalized non-flanked reads.
minimap2 -x asm20 -c --cs=long --secondary=yes -N 20 "$HXB2" "$NONFLANKED_NORM" \
    > "$WORK/${PID}_vs_HXB2.paf"

minimap2 -x asm20 -c --cs=long --secondary=yes -N 20 "$CIRCLE_REFS" "$NONFLANKED_NORM" \
    > "$WORK/${PID}_vs_circle_refs.paf"

# 4. Match normalized non-flanked reads to normalized flanked reads.
if [[ -s "$FLANKED_NORM" ]] && grep -q '^>' "$FLANKED_NORM"; then
    minimap2 -x asm20 -c --cs=long --secondary=yes -N 50 "$FLANKED_NORM" "$NONFLANKED_NORM" \
        > "$WORK/${PID}_vs_flanked_HIV.paf"
else
    : > "$WORK/${PID}_vs_flanked_HIV.paf"
fi

# 5. Assign a putative integration site only when no circle-specific
#    architecture is detected and qualifying matches resolve consistently.
python3 "$SCRIPT_DIR/classify_nonflanked_hiv.py" classify \
    --participant-id "$PID" \
    --sequences "$NONFLANKED_NORM" \
    --linear-paf "$WORK/${PID}_vs_HXB2.paf" \
    --circle-paf "$WORK/${PID}_vs_circle_refs.paf" \
    --flanked-paf "$WORK/${PID}_vs_flanked_HIV.paf" \
    --master "$MASTER" \
    --output "$OUT/${PID}_nonflanked_HIV_classification.csv"

echo "Finished participant $PID"
