#!/bin/bash
#SBATCH --job-name=baf_hrd_pipeline
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00

# ============================================================
# run_pipeline.sbatch
#
# Runs baf_to_hrd_pipeline.py (Stages 1–3) sequentially.
# Submitted by orchestrate_pipeline.sh after all Stage 0
# (BAM parsing) jobs have completed.
#
# Can also be submitted manually for runs where all samples
# already have TSV files:
#   sbatch run_pipeline.sbatch
# ============================================================

set -euo pipefail

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
mkdir -p "${SCRIPT_DIR}/logs"

echo "============================================================"
echo "BAF-to-HRD Pipeline — Stages 1–3"
echo "DATE:    $(date)"
echo "JOB ID:  ${SLURM_JOB_ID:-NA}"
echo "============================================================"

# PIPELINE_CONFIG can be injected by orchestrate_pipeline.sh via --export.
# Falls back to the last manual config for standalone runs.
CONFIG="${PIPELINE_CONFIG:-${REPO_DIR}/configs/pipeline_config_20260421.yaml}"
echo "CONFIG:  ${CONFIG}"

conda run -n shapeit4 python "${REPO_DIR}/pipeline/baf_to_hrd_pipeline.py" \
    --config "${CONFIG}"

echo "============================================================"
echo "Pipeline complete: $(date)"
echo "============================================================"
