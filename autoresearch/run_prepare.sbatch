#!/bin/bash
#SBATCH --job-name=autoresearch_prepare
#SBATCH --output=logs/prepare_%j.log
#SBATCH --error=logs/prepare_%j.err
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=1:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jmroldan@cnio.es

# ============================================================
# run_prepare.sbatch
#
# Run ONCE before autoresearch experiments begin.
# Builds autoresearch_baf_hrd/val_data.csv from the Stage 2 merged CSV
# (20260406/merged_segments_baf_metrics_ALL_with20bins_enhancedFeatures.csv).
#
# Only rows with nMajor / nMinor present (ASCAT-profiled samples) and
# dp_filter == 8 are kept.  Output has 20 corrected + 20 raw BAF bin
# columns, plus ground-truth class labels.
#
# Submit:
#   sbatch autoresearch_baf_hrd/run_prepare.sbatch
# ============================================================

set -euo pipefail

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
mkdir -p "${SCRIPT_DIR}/logs"
mkdir -p "${SCRIPT_DIR}/results"

echo "=== run_prepare.sbatch  job=${SLURM_JOB_ID}  $(date) ==="

conda run -n shapeit4 python "${SCRIPT_DIR}/prepare.py"

echo "=== DONE  $(date) ==="
