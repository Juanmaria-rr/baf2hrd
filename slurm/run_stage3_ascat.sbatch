#!/bin/bash
#SBATCH --job-name=baf_stage3_ascat
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00

# ============================================================
# run_stage3_ascat.sbatch
#
# Runs Stage 3 only (allele-specific CN inference) on the
# 20260406 ASCAT-downsampled merged CSV using the per-DP
# calibrated thresholds in thresholds_analysis/thresholds_lookup.json.
#
# Prerequisite: Stage 1+2 have already been run for the 20260406
#   data — the merged CSV must exist at:
#   .../20260406/merged_segments_baf_metrics_ALL_with20bins_enhancedFeatures.csv
#
# Output: .../20260406/scarHRD_inputs_thresholds/
# ============================================================

set -euo pipefail

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
mkdir -p "${REPO_DIR}/logs"

echo "============================================================"
echo "Stage 3 — ASCAT-downsampled (date: 20260406)"
echo "DATE:    $(date)"
echo "JOB ID:  ${SLURM_JOB_ID:-NA}"
echo "============================================================"

CONFIG="${REPO_DIR}/configs/pipeline_config_20260406.yaml"

conda run -n shapeit4 python "${REPO_DIR}/pipeline/baf_to_hrd_pipeline.py" \
    --config "${CONFIG}" \
    --stages 3

echo "============================================================"
echo "Stage 3 complete: $(date)"
echo "============================================================"
