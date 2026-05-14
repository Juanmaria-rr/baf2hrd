#!/bin/bash
#SBATCH --job-name=index_bam
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=2:00:00

# ============================================================
# run_index_bam.sbatch
#
# Indexes a BAM with samtools index.  Called by
# orchestrate_pipeline.sh when a BAM is missing its .bai or
# the existing .bai is empty (0 bytes, corrupt transfer stub).
# Stage 0 mpileup jobs depend on this job when indexing is needed.
#
# Usage (called by orchestrator — do not submit manually):
#   sbatch run_index_bam.sbatch <BAM>
# ============================================================

set -euo pipefail

export PATH="/home/jmroldan/miniconda3/envs/shapeit4/bin:$PATH"

BAM="${1:?Usage: sbatch run_index_bam.sbatch BAM}"
BAI="${BAM}.bai"

mkdir -p logs

echo "============================================================"
echo "DATE:    $(date)"
echo "JOB ID:  ${SLURM_JOB_ID:-NA}"
echo "BAM:     $BAM"
echo "BAI:     $BAI"
echo "============================================================"

samtools --version | head -n 1

[[ -f "$BAM" ]] || { echo "[ERROR] BAM not found: $BAM"; exit 1; }

# Remove empty stub so samtools index can write cleanly
if [[ -f "$BAI" && ! -s "$BAI" ]]; then
    echo "[INFO] Removing empty stub: $BAI"
    rm -f "$BAI"
fi

echo "[INFO] Running samtools index (-@ ${SLURM_CPUS_PER_TASK})..."
samtools index -@ "${SLURM_CPUS_PER_TASK}" "$BAM"

echo "[INFO] Index written: $BAI  ($(du -h "$BAI" | cut -f1))"
echo "[INFO] Done at $(date)"
echo "============================================================"
