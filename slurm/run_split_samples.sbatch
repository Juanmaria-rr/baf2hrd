#!/bin/bash
#SBATCH --job-name=split_samples
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:15:00

# ============================================================
# run_split_samples.sbatch
#
# Generates an 80/20 sample-level train/val split from the
# 20260406 merged CSV.  Outputs:
#   autoresearch/calibration_samples.txt
#   autoresearch/validation_samples.txt
#
# After this job completes, run:
#   sbatch --dependency=afterok:<this_jobid> run_recalibrate.sbatch
#   sbatch --dependency=afterok:<this_jobid> autoresearch/run_prepare.sbatch
# ============================================================

set -euo pipefail

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
mkdir -p "${REPO_DIR}/logs"

echo "============================================================"
echo "Sample-level train/val split"
echo "DATE:    $(date)"
echo "JOB ID:  ${SLURM_JOB_ID:-NA}"
echo "============================================================"

conda run -n shapeit4 python "${REPO_DIR}/pipeline/split_samples.py"

echo "============================================================"
echo "Split complete: $(date)"
echo "============================================================"
