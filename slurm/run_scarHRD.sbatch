#!/bin/bash
#SBATCH --job-name=scarHRD
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err
#SBATCH --partition=main
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jmroldan@cnio.es

# ============================================================
# run_scarHRD.sbatch
#
# Computes HRD scores via scarHRD for all *scarHRD_input*.txt
# files produced by Stage 3.
#
# Required (inject via --export):
#   REPO_DIR   — repository root (baf_hrd_pipeline/)
#   BASE_DIR   — pipeline run output dir (tso_samples/YYYYMMDD/)
#
# Input:  $BASE_DIR/scarHRD_inputs_thresholds/
# Output: $BASE_DIR/outputs_scarHRD/ALL_scarHRD_scores.csv
#
# Submitted automatically by orchestrate_pipeline.sh after
# run_pipeline.sbatch completes. Can also be run manually:
#   sbatch --export=ALL,REPO_DIR="$(pwd)",BASE_DIR="/path/to/YYYYMMDD" \
#     slurm/run_scarHRD.sbatch
# ============================================================

set -euo pipefail

REPO_DIR="${REPO_DIR:?REPO_DIR not set — submit with --export=ALL,REPO_DIR=\"\$(pwd)\",...}"
BASE_DIR="${BASE_DIR:?BASE_DIR not set — submit with --export=ALL,...,BASE_DIR=\"/path/to/YYYYMMDD\"}"

INPUT_DIR="${BASE_DIR}/scarHRD_inputs_thresholds"
OUT_DIR="${BASE_DIR}/outputs_scarHRD"

mkdir -p "${REPO_DIR}/slurm/logs"
mkdir -p "$OUT_DIR"

echo "============================================================"
echo "scarHRD — HRD scoring"
echo "DATE:      $(date)"
echo "JOB ID:    ${SLURM_JOB_ID:-NA}"
echo "BASE_DIR:  $BASE_DIR"
echo "INPUT_DIR: $INPUT_DIR"
echo "OUT_DIR:   $OUT_DIR"
echo "============================================================"

conda run -n r43 Rscript "${REPO_DIR}/pipeline/run_scarHRD.R" \
    "$INPUT_DIR" \
    "$OUT_DIR"

echo "============================================================"
echo "scarHRD complete: $(date)"
echo "============================================================"
