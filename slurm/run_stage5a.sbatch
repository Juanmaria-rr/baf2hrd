#!/bin/bash
#SBATCH --job-name=eval_ascat
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00

# ============================================================
# run_stage5a.sbatch
#
# Stage 5a: compare threshold-inferred HRD scores from
# ASCAT-downsampled samples against the ASCAT gold standard.
#
# Submit after run_stage4_ascat.sbatch completes:
#   sbatch --dependency=afterok:<stage4_jobid> run_stage5a.sbatch
# ============================================================

set -euo pipefail

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
mkdir -p "${REPO_DIR}/logs"

echo "============================================================"
echo "Stage 5a — ASCAT evaluation"
echo "DATE:    $(date)"
echo "JOB ID:  ${SLURM_JOB_ID:-NA}"
echo "============================================================"

conda run -n shapeit4 python "${REPO_DIR}/pipeline/evaluate_ascat_downsampled.py"

echo "============================================================"
echo "Stage 5a complete: $(date)"
echo "============================================================"
