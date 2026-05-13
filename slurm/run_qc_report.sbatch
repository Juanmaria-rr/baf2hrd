#!/bin/bash
#SBATCH --job-name=qc_report
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00

# ============================================================
# run_qc_report.sbatch
#
# Generates QC plots for the BAF-to-HRD pipeline and writes
# a self-contained HTML report to {base_dir}/qc_report.html.
#
# Usage:
#   sbatch --export=ALL,REPO_DIR="$(pwd)",CONFIG=configs/pipeline_config_YYYYMMDD.yaml \
#     slurm/run_qc_report.sbatch
#
# CONFIG is required; REPO_DIR defaults to the submit directory.
# ============================================================

set -euo pipefail

REPO_DIR="${REPO_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
CONFIG="${CONFIG:?CONFIG not set — submit with: sbatch --export=ALL,REPO_DIR=\"\$(pwd)\",CONFIG=configs/pipeline_config_YYYYMMDD.yaml ...}"

mkdir -p "${REPO_DIR}/slurm/logs"

echo "============================================================"
echo "QC Report — BAF-to-HRD Pipeline"
echo "DATE:    $(date)"
echo "JOB ID:  ${SLURM_JOB_ID:-NA}"
echo "CONFIG:  ${CONFIG}"
echo "============================================================"

conda run -n shapeit4 python "${REPO_DIR}/pipeline/qc_report.py" \
    --config "${CONFIG}"

echo "============================================================"
echo "QC report complete: $(date)"
echo "============================================================"
