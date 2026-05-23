#!/bin/bash
#SBATCH --job-name=eval_hrd
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00

# ============================================================
# run_eval_hrd.sbatch
#
# Generates HRD performance report (scatter + confusion matrix)
# across all runs, with parsing QC audit.
#
# Always produces:
#   hrd_performance_report.html          ← all samples
#   hrd_performance_report_filtered.html ← blacklisted samples excluded
#                                           (only if configs/blacklist.txt is
#                                            non-empty after stripping comments)
#
# Usage (ad hoc):
#   sbatch --export=ALL,REPO_DIR="$(pwd)",CURRENT_RUN=20260514 \
#     slurm/run_eval_hrd.sbatch
#
# Launched automatically by orchestrate_pipeline.sh after Stage 4.
# CURRENT_RUN defaults to empty string (no highlight).
# BLACKLIST_FILE defaults to configs/blacklist.txt in REPO_DIR.
# ============================================================

set -euo pipefail

REPO_DIR="${REPO_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
CURRENT_RUN="${CURRENT_RUN:-}"
MIN_RUN="${MIN_RUN:-20260511}"
BLACKLIST_FILE="${BLACKLIST_FILE:-${REPO_DIR}/configs/blacklist.txt}"

TSO_SAMPLES_ROOT="/storage/scratch01/groups/co/cn_extra/alleleSpecific/tso_samples"
GIS_SCORES="${TSO_SAMPLES_ROOT}/MYRIAD_GIS_SCORES/gis_scores.csv"

mkdir -p "${REPO_DIR}/slurm/logs"

echo "============================================================"
echo "HRD Performance Report"
echo "DATE:        $(date)"
echo "JOB ID:      ${SLURM_JOB_ID:-NA}"
echo "CURRENT RUN: ${CURRENT_RUN:-<all runs, no highlight>}"
echo "BLACKLIST:   ${BLACKLIST_FILE}"
echo "============================================================"

CURRENT_RUN_ARG=()
if [[ -n "${CURRENT_RUN}" ]]; then
    CURRENT_RUN_ARG=(--current-run "${CURRENT_RUN}")
fi

# ── Full report (all samples) ─────────────────────────────────────────────
echo "Generating full report..."
conda run -n shapeit4 python "${REPO_DIR}/pipeline/eval_hrd_performance.py" \
    --tso-samples-root "${TSO_SAMPLES_ROOT}" \
    --gis-scores       "${GIS_SCORES}"       \
    --min-run          "${MIN_RUN}"           \
    "${CURRENT_RUN_ARG[@]}"
echo "Written: ${TSO_SAMPLES_ROOT}/hrd_performance_report.html"

# ── Filtered report (blacklist applied) ──────────────────────────────────
EXCLUDE_SAMPLES=()
if [[ -f "${BLACKLIST_FILE}" ]]; then
    while IFS= read -r line; do
        [[ "$line" =~ ^#.*$ || -z "$line" ]] && continue
        EXCLUDE_SAMPLES+=("$line")
    done < "${BLACKLIST_FILE}"
fi

if [[ ${#EXCLUDE_SAMPLES[@]} -gt 0 ]]; then
    echo "Generating filtered report (${#EXCLUDE_SAMPLES[@]} samples excluded: ${EXCLUDE_SAMPLES[*]})..."
    conda run -n shapeit4 python "${REPO_DIR}/pipeline/eval_hrd_performance.py" \
        --tso-samples-root "${TSO_SAMPLES_ROOT}" \
        --gis-scores       "${GIS_SCORES}"       \
        --min-run          "${MIN_RUN}"           \
        --output           "${TSO_SAMPLES_ROOT}/hrd_performance_report_filtered.html" \
        --exclude-samples  "${EXCLUDE_SAMPLES[@]}" \
        "${CURRENT_RUN_ARG[@]}"
    echo "Written: ${TSO_SAMPLES_ROOT}/hrd_performance_report_filtered.html"
else
    echo "Blacklist empty or not found — skipping filtered report."
fi

echo "============================================================"
echo "HRD performance report complete: $(date)"
echo "============================================================"
