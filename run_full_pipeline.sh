#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
    cat <<'TXT'
Usage:
  bash run_full_pipeline.sh FINAL_RESULTS_DIR HXB2_FASTA [SAMPLE_MAP_CSV] [RESULTS_DIR]

Arguments:
  FINAL_RESULTS_DIR   Folder containing sample subfolders with *.annotated.csv files
  HXB2_FASTA          HXB2 reference FASTA
  SAMPLE_MAP_CSV      Optional CSV with sample_id,participant_id columns. Use - for none.
  RESULTS_DIR         Final results folder (default: results)

Notes:
  Intermediate files are written automatically to work/.

Examples:
  bash run_full_pipeline.sh final_results HXB2.fasta
  bash run_full_pipeline.sh final_results HXB2.fasta sample_to_participant.csv
  bash run_full_pipeline.sh final_results HXB2.fasta - my_results
TXT
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
fi

if (( $# < 2 || $# > 4 )); then
    usage >&2
    exit 1
fi

FINAL_RESULTS_DIR="$1"
HXB2="$2"
SAMPLE_MAP="${3:-}"
[[ "$SAMPLE_MAP" == "-" ]] && SAMPLE_MAP=""
RESULTS_DIR="${4:-results}"

WORK_DIR="work"
PREPARED_DIR="$WORK_DIR/prepared"
PARTICIPANT_WORK_DIR="$WORK_DIR/participants"
PARTICIPANT_RESULTS_DIR="$RESULTS_DIR/participants"

MASTER="$PREPARED_DIR/masterfile.csv"
CIRCLE_REFS="$WORK_DIR/circle_refs.fasta"
SUMMARY="$RESULTS_DIR/participant_hiv_summary.csv"
ALL_CLASS="$RESULTS_DIR/all_participants_nonflanked_HIV_classification.csv"
FINAL_OUTPUT="$RESULTS_DIR/masterfile_with_putative_integration_sites.csv"

# ------------------------------------------------------------
# Check software and input files before starting.
# ------------------------------------------------------------
for command_name in python3 minimap2; do
    command -v "$command_name" >/dev/null 2>&1 || {
        echo "Missing required program: $command_name" >&2
        exit 1
    }
done

python3 -c 'import pandas, Bio' >/dev/null 2>&1 || {
    echo "Missing Python package(s). Install with: python3 -m pip install -r $SCRIPT_DIR/requirements.txt" >&2
    exit 1
}

[[ -d "$FINAL_RESULTS_DIR" ]] || { echo "Missing folder: $FINAL_RESULTS_DIR" >&2; exit 1; }
[[ -s "$HXB2" ]] || { echo "Missing HXB2 FASTA: $HXB2" >&2; exit 1; }
if [[ -n "$SAMPLE_MAP" && ! -s "$SAMPLE_MAP" ]]; then
    echo "Missing sample mapping CSV: $SAMPLE_MAP" >&2
    exit 1
fi

mkdir -p "$PREPARED_DIR" "$PARTICIPANT_WORK_DIR" "$PARTICIPANT_RESULTS_DIR"

# Print the version, thresholds, and paths used in this run.
echo ""
python3 "$SCRIPT_DIR/classify_nonflanked_hiv.py" settings

echo "Annotated CSV input: $FINAL_RESULTS_DIR"
echo "HXB2 reference:      $HXB2"
echo "Sample mapping:      ${SAMPLE_MAP:-none (sample_id = participant_id)}"
echo "Work folder:         $WORK_DIR"
echo "Results folder:      $RESULTS_DIR"

# ------------------------------------------------------------
# 1. Combine annotated CSVs and create participant FASTAs.
# ------------------------------------------------------------
if [[ -n "$SAMPLE_MAP" ]]; then
    python3 "$SCRIPT_DIR/combine_csv.py" "$FINAL_RESULTS_DIR" "$PREPARED_DIR" "$SAMPLE_MAP"
else
    python3 "$SCRIPT_DIR/combine_csv.py" "$FINAL_RESULTS_DIR" "$PREPARED_DIR"
fi

# ------------------------------------------------------------
# 2. Build HXB2-derived references for circular HIV architecture.
# ------------------------------------------------------------
python3 "$SCRIPT_DIR/classify_nonflanked_hiv.py" make-refs \
    --hxb2 "$HXB2" \
    --output "$CIRCLE_REFS"

# ------------------------------------------------------------
# 3. Run each participant.
# ------------------------------------------------------------
shopt -s nullglob
files=("$PREPARED_DIR"/*_hiv.fasta)
((${#files[@]})) || { echo "No participant HIV FASTAs were created." >&2; exit 1; }

rm -f "$ALL_CLASS"
printf "participant_id,total_hiv,flanked_hiv,non_flanked_hiv,classification_run,non_flanked_circular,putative_integration_assigned,putative_status_failed_a_matching_threshold,putative_status_ambiguous_multiple_clone_id2,putative_status_ambiguous_conflicting_metadata,putative_status_not_evaluated\n" > "$SUMMARY"

for file in "${files[@]}"; do
    pid=$(basename "$file" _hiv.fasta)
    nonfile="$PREPARED_DIR/${pid}_hiv_non_flanked.fasta"

    total=$(grep -c '^>' "$file" || true)
    non=$(grep -c '^>' "$nonfile" 2>/dev/null || true)
    flanked=$((total - non))

    echo ""
    echo "=== Participant $pid ==="
    bash "$SCRIPT_DIR/run_participant_pipeline.sh" "$pid" "$HXB2" "$WORK_DIR" "$RESULTS_DIR"

    circular=0
    assigned=0
    failed_threshold=0
    ambiguous_clone=0
    conflicting_metadata=0
    not_evaluated=0

    if (( non > 0 )); then
        class="$PARTICIPANT_RESULTS_DIR/$pid/${pid}_nonflanked_HIV_classification.csv"

        if [[ ! -s "$ALL_CLASS" ]]; then
            cat "$class" > "$ALL_CLASS"
        else
            tail -n +2 "$class" >> "$ALL_CLASS"
        fi

        read -r classified circular assigned failed_threshold ambiguous_clone conflicting_metadata < <(
            awk -F ',' '
                NR > 1 {
                    classified++
                    if ($3 ~ /^1-LTR/ || $3 ~ /^2-LTR/ || $3 ~ /^Circular HIV DNA/) circular++
                    if ($3 == "No circle-specific architecture detected") {
                        if ($4 == "Putatively integrated provirus") assigned++
                        else if ($4 ~ /^Origin indeterminate/) failed_threshold++
                        else if ($4 ~ /^Putative integration ambiguous - multiple CLONE_ID2s/) ambiguous_clone++
                        else if ($4 ~ /^Putative integration ambiguous - conflicting flanked-read metadata/) conflicting_metadata++
                    }
                }
                END {
                    print classified+0, circular+0, assigned+0, failed_threshold+0, ambiguous_clone+0, conflicting_metadata+0
                }
            ' "$class"
        )

        not_evaluated=$((non - classified))
        if (( not_evaluated < 0 )); then
            not_evaluated=0
        fi
        run=yes
    else
        run=no
    fi

    printf "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" \
        "$pid" "$total" "$flanked" "$non" "$run" "$circular" "$assigned" \
        "$failed_threshold" "$ambiguous_clone" "$conflicting_metadata" "$not_evaluated" >> "$SUMMARY"
done

# ------------------------------------------------------------
# 4. Add circularization, match evidence, and putative integration
#    results to the combined masterfile.
# ------------------------------------------------------------
if [[ ! -s "$ALL_CLASS" ]]; then
    printf "participant_id,sequence_id,circle_annotation,putative_integration_status,CHROMOSOME_NEW,INTEGRATION_SITE_NEW,CLONE_ID_NEW,matched_flanked_reads,match_percent_identity,match_shorter_fragment_coverage,match_alignment_bp\n" > "$ALL_CLASS"
fi

python3 "$SCRIPT_DIR/classify_nonflanked_hiv.py" update-master \
    --master "$MASTER" \
    --classification "$ALL_CLASS" \
    --output "$FINAL_OUTPUT"

echo ""
echo "Finished."
echo "Main output:         $FINAL_OUTPUT"
echo "Participant summary: $SUMMARY"
echo "Intermediate files:  $WORK_DIR"
