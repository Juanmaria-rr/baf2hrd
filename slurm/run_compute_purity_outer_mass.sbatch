#!/bin/bash
#SBATCH --job-name=purity_outer_mass
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=1:00:00

# ============================================================
# run_compute_purity_outer_mass.sbatch
#
# Adds outer_mass_purity_corrected to the ASCAT calibration CSV
# using analytical bin-center remapping (no Stage 1 re-run needed).
# Output goes to tso_samples/20260406_purity/.
#
# Usage:
#   cd /path/to/baf_hrd_pipeline
#   sbatch --export=ALL,REPO_DIR="$(pwd)" slurm/run_compute_purity_outer_mass.sbatch
# ============================================================

set -euo pipefail

REPO_DIR="${REPO_DIR:?REPO_DIR not set — submit with: sbatch --export=ALL,REPO_DIR=\"\$(pwd)\" ...}"
mkdir -p "${REPO_DIR}/slurm/logs"

INPUT="/storage/scratch01/groups/co/cn_extra/alleleSpecific/tso_samples/20260406/merged_segments_baf_metrics_ALL_with20bins_enhancedFeatures.csv"
OUTPUT_DIR="/storage/scratch01/groups/co/cn_extra/alleleSpecific/tso_samples/20260406_purity"

echo "============================================================"
echo "Purity-corrected outer_mass — Fase 1 recalibración"
echo "DATE:   $(date)"
echo "JOB ID: ${SLURM_JOB_ID:-NA}"
echo "INPUT:  ${INPUT}"
echo "OUTDIR: ${OUTPUT_DIR}"
echo "============================================================"

conda run -n shapeit4 python "${REPO_DIR}/pipeline/compute_purity_corrected_outer_mass.py" \
    --input    "${INPUT}" \
    --output-dir "${OUTPUT_DIR}"

echo "============================================================"
echo "Done: $(date)"
echo "============================================================"
