#!/bin/bash
#SBATCH --job-name=parse_mpileup
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00

# ============================================================
# run_parse_mpileup.sbatch
#
# Converts raw bcftools-mpileup TSV to _qc.tsv (Stage 0b).
# Called by orchestrate_pipeline.sh with
#   --dependency=afterok:<mpileup_job_id>
#
# Usage (called by orchestrator — do not submit manually):
#   sbatch run_parse_mpileup.sbatch <SAMPLE> <OUTDIR>
#
# Arguments:
#   $1  SAMPLE  — sample base name (e.g. JBLAB17005)
#   $2  OUTDIR  — directory containing the raw mpileup TSV
#                 and where _qc.tsv will be written
# ============================================================

set -euo pipefail

SAMPLE="${1:?Usage: sbatch run_parse_mpileup.sbatch SAMPLE OUTDIR}"
OUTDIR="${2:?Usage: sbatch run_parse_mpileup.sbatch SAMPLE OUTDIR}"

REPO_DIR="${REPO_DIR:?REPO_DIR not set — submit via orchestrate_pipeline.sh or with --export=ALL,REPO_DIR=...}"

RAW_TSV="${OUTDIR}/${SAMPLE}_allele_counts_snp6_wrapper_qualityMetrics_q30Pass.tsv"
QC_TSV="${OUTDIR}/${SAMPLE}_qc.tsv"

echo "============================================================"
echo "DATE:    $(date)"
echo "JOB ID:  ${SLURM_JOB_ID:-NA}"
echo "SAMPLE:  $SAMPLE"
echo "INPUT:   $RAW_TSV"
echo "OUTPUT:  $QC_TSV"
echo "============================================================"

[[ -f "$RAW_TSV" ]] || { echo "[ERROR] Raw mpileup TSV not found: $RAW_TSV"; exit 1; }

conda run -n shapeit4 python "${REPO_DIR}/pipeline/parse_mpileup.py" "$RAW_TSV" "$QC_TSV"

echo "[INFO] Finished at $(date)"
echo "============================================================"
