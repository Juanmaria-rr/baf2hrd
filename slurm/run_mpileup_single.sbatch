#!/bin/bash
#SBATCH --job-name=mpileup_single
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00

# ============================================================
# run_mpileup_single.sbatch
#
# Single-sample bcftools mpileup → raw allele-count TSV.
# Called by orchestrate_pipeline.sh for samples that do not
# yet have a _qc.tsv file.
#
# Usage (called by orchestrator — do not submit manually):
#   sbatch run_mpileup_single.sbatch <SAMPLE> <OUTDIR>
#
# Arguments:
#   $1  SAMPLE   — sample base name (e.g. JBLAB17005)
#   $2  OUTDIR   — directory to write output files into
#
# Output:
#   {OUTDIR}/{SAMPLE}_allele_counts_snp6_wrapper_qualityMetrics_q30Pass.tsv
#   {OUTDIR}/{SAMPLE}_bam_global_stats_qualityMetrics_q30Pass.txt
#
# The _qc.tsv is produced by the subsequent run_parse_mpileup.sbatch job
# submitted with --dependency=afterok:<this_job_id> by the orchestrator.
# ============================================================

set -euo pipefail

SAMPLE="${1:?Usage: sbatch run_mpileup_single.sbatch SAMPLE OUTDIR}"
OUTDIR="${2:?Usage: sbatch run_mpileup_single.sbatch SAMPLE OUTDIR}"

# ============================================================
# Fixed paths — edit to match your cluster layout
# ============================================================
BAM_DIR="/storage/scratch01/groups/co/cn_extra/alleleSpecific/tso_samples/rg_bam"
FASTA="/home/_groups/co/reference/hg19/hg19.fa"
BED="/storage/scratch01/groups/co/cn_extra/alleleSpecific/tso_samples/snp6_nochrM_fixed.bed"
MIN_BQ=30
MIN_MQ=30

# ============================================================
# Resolve output paths
# ============================================================
mkdir -p "$OUTDIR" logs

OUT_PILEUP="${OUTDIR}/${SAMPLE}_allele_counts_snp6_wrapper_qualityMetrics_q30Pass.tsv"
OUT_STATS="${OUTDIR}/${SAMPLE}_bam_global_stats_qualityMetrics_q30Pass.txt"

# ============================================================
# Detect BAM (EXP08A or EXP08B)
# ============================================================
BAM_A="${BAM_DIR}/${SAMPLE}-TSO500-EXP08A.mdup.rg.bam"
BAM_B="${BAM_DIR}/${SAMPLE}-TSO500-EXP08B.mdup.rg.bam"

if [[ -f "$BAM_A" ]]; then
  BAM="$BAM_A"
  EXP_USED="EXP08A"
elif [[ -f "$BAM_B" ]]; then
  BAM="$BAM_B"
  EXP_USED="EXP08B"
else
  echo "[ERROR] No BAM found for $SAMPLE"
  echo "  Tried: $BAM_A"
  echo "  Tried: $BAM_B"
  exit 1
fi

# ============================================================
# Reproducibility log
# ============================================================
echo "============================================================"
echo "DATE:       $(date)"
echo "JOB ID:     ${SLURM_JOB_ID:-NA}"
echo "SAMPLE:     $SAMPLE"
echo "EXPERIMENT: $EXP_USED"
echo "BAM:        $BAM"
echo "FASTA:      $FASTA"
echo "BED:        $BED"
echo "OUTDIR:     $OUTDIR"
echo "BQ/MQ:      ${MIN_BQ}/${MIN_MQ}"
echo "============================================================"

bcftools --version | head -n 1
samtools --version | head -n 1

[[ -f "$FASTA" ]] || { echo "[ERROR] FASTA missing: $FASTA"; exit 1; }
[[ -f "$BED"   ]] || { echo "[ERROR] BED missing: $BED";   exit 1; }

# ============================================================
# Step 1: mpileup
# ============================================================
echo "[INFO] Running bcftools mpileup..."

bcftools mpileup \
  -f "$FASTA" \
  -R "$BED" \
  --min-BQ "$MIN_BQ" \
  --min-MQ "$MIN_MQ" \
  -a FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/SP,INFO/MQBZ,INFO/BQBZ \
  -Ou "$BAM" \
| bcftools query \
    -f '%CHROM\t%POS\t%REF\t%ALT\t[%DP]\t[%AD]\t[%ADF]\t[%ADR]\t[%SP]\t%MQBZ\t%BQBZ\n' \
> "$OUT_PILEUP"

echo "[INFO] Pileup written: $OUT_PILEUP"

# ============================================================
# Step 2: BAM stats
# ============================================================
echo "[INFO] Running samtools stats..."
samtools stats -@ "${SLURM_CPUS_PER_TASK}" "$BAM" > "$OUT_STATS"
echo "[INFO] BAM stats written: $OUT_STATS"

echo "[INFO] Job finished at $(date)"
echo "============================================================"
