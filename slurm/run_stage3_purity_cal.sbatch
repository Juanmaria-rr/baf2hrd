#!/bin/bash
#SBATCH --job-name=stage3_purity_cal
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err
#SBATCH --partition=main
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jmroldan@cnio.es

# ============================================================
# run_stage3_purity_cal.sbatch
#
# Re-runs Stage 3 with purity-calibrated thresholds + corrected
# outer_mass feature, then runs scarHRD inline.
#
# Required (inject via --export):
#   REPO_DIR   — repository root (baf_hrd_pipeline/)
#   CONFIG     — path to pipeline_config_YYYYMMDD_purity_cal.yaml
#   BASE_DIR   — pipeline run output dir (e.g. tso_samples/YYYYMMDD_purity)
#
# Example:
#   REPO="$(pwd)"
#   for DATE in 20260511 20260514 20260520; do
#     sbatch --export=ALL,\
#       REPO_DIR="$REPO",\
#       CONFIG="$REPO/configs/pipeline_config_${DATE}_purity_cal.yaml",\
#       BASE_DIR="/storage/scratch01/groups/co/cn_extra/alleleSpecific/tso_samples/${DATE}_purity" \
#       slurm/run_stage3_purity_cal.sbatch
#   done
# ============================================================

set -euo pipefail

REPO_DIR="${REPO_DIR:?REPO_DIR not set}"
CONFIG="${CONFIG:?CONFIG not set}"
BASE_DIR="${BASE_DIR:?BASE_DIR not set}"

mkdir -p "${REPO_DIR}/slurm/logs"

echo "============================================================"
echo "Stage 3 purity-calibrated re-run"
echo "DATE:     $(date)"
echo "JOB ID:   ${SLURM_JOB_ID:-NA}"
echo "CONFIG:   ${CONFIG}"
echo "BASE_DIR: ${BASE_DIR}"
echo "============================================================"

# --- Stage 3 ---
conda run -n shapeit4 python "${REPO_DIR}/pipeline/baf_to_hrd_pipeline.py" \
    --config "${CONFIG}"

echo "Stage 3 complete: $(date)"

# --- Stage 4: scarHRD ---
INPUT_DIR="${BASE_DIR}/scarHRD_inputs_thresholds"
OUT_DIR="${BASE_DIR}/outputs_scarHRD"

mkdir -p "${OUT_DIR}"

echo "------------------------------------------------------------"
echo "Running scarHRD"
echo "INPUT_DIR: ${INPUT_DIR}"
echo "OUT_DIR:   ${OUT_DIR}"
echo "------------------------------------------------------------"

conda run -n r43 Rscript "${REPO_DIR}/pipeline/run_scarHRD.R" \
    "${INPUT_DIR}" \
    "${OUT_DIR}"

echo "============================================================"
echo "All done: $(date)"
echo "============================================================"
