#!/bin/bash
#SBATCH --job-name=calibrate_purity
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=2:00:00

# ============================================================
# run_calibrate_purity.sbatch
#
# Fase 2 — Recalibración de thresholds con corrección de pureza.
# Usa outer_mass_purity_corrected (generado en Fase 1) para
# derivar nuevos thresholds_lookup_purity.json.
#
# Usage:
#   cd /path/to/baf_hrd_pipeline
#   sbatch --export=ALL,REPO_DIR="$(pwd)" slurm/run_calibrate_purity.sbatch
# ============================================================

set -euo pipefail

REPO_DIR="${REPO_DIR:?REPO_DIR not set — submit with: sbatch --export=ALL,REPO_DIR=\"\$(pwd)\" ...}"
mkdir -p "${REPO_DIR}/slurm/logs"

BASE="/storage/scratch01/groups/co/cn_extra/alleleSpecific/tso_samples"
INPUT="${BASE}/20260406_purity/merged_segments_baf_metrics_ALL_purity_corrected.csv"
OUTDIR="${BASE}/20260406_purity/thresholds_analysis_purity"

echo "============================================================"
echo "Threshold calibration — purity-corrected outer_mass (Fase 2)"
echo "DATE:    $(date)"
echo "JOB ID:  ${SLURM_JOB_ID:-NA}"
echo "INPUT:   ${INPUT}"
echo "OUTDIR:  ${OUTDIR}"
echo "FEATURE: outer_mass_purity_corrected"
echo "============================================================"

conda run -n shapeit4 python "${REPO_DIR}/pipeline/calibrate_thresholds.py" \
    --input   "${INPUT}" \
    --feature "outer_mass_purity_corrected" \
    --outdir  "${OUTDIR}"

echo "============================================================"
echo "Done: $(date)"
echo "============================================================"
