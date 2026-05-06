#!/bin/bash
#SBATCH --job-name=stage0_array
#SBATCH --output=logs/%x.%A_%a.out
#SBATCH --error=logs/%x.%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
# Array index and concurrency are set at submission time, e.g.:
#   sbatch --array=0-419%20 run_stage0_array.sbatch bam_list.txt /path/to/tsv_outdir

# ============================================================
# run_stage0_array.sbatch
#
# SLURM array job: for each BAM in bam_list.txt, resolves it
# (local path or s3://), runs mpileup + parse_mpileup to produce
# _qc.tsv, then deletes the local BAM copy.
#
# Concurrency is controlled by the %N suffix on --array, which
# limits peak BAM disk usage to N × ~7 GB.
#
# BAM sources supported (one per line in bam_list.txt):
#   Local path:  /storage/.../rg_bam_20260505_all/JBLAB001-TSO500-EXP08A.mdup.rg.bam
#   S3 URI:      s3://cnio-data-sharing/all_bams/JBLAB001-TSO500-EXP08A.mdup.rg.bam
#   S3 key:      all_bams/JBLAB001-TSO500-EXP08A.mdup.rg.bam
#
# For local BAMs the file is used in place and NOT deleted after
# processing (it was already transferred — user decides when to
# clean up).  For S3 sources the BAM is downloaded to a temp dir
# and deleted on exit.
#
# Usage:
#   # Generate bam_list.txt — local paths
#   find /storage/.../rg_bam_20260505_all -name "*.mdup.rg.bam" > bam_list.txt
#
#   # Generate bam_list.txt — from S3
#   aws s3 ls s3://cnio-data-sharing/all_bams/ --recursive \
#     | grep '\.bam$' | awk '{print $4}' > bam_list.txt
#
#   # Submit (20 concurrent tasks = ~140 GB peak for S3 source)
#   N=$(wc -l < bam_list.txt)
#   sbatch --array=0-$(( N - 1 ))%20 \
#     run_stage0_array.sbatch bam_list.txt /path/to/tsv_outdir
#
# Skip logic: if OUTDIR/SAMPLE_qc.tsv already exists the task
# exits immediately — safe to re-submit after partial failures.
# ============================================================

set -euo pipefail

BAM_LIST="${1:?Usage: sbatch --array=0-N%BATCH run_stage0_array.sbatch BAM_LIST OUTDIR}"
OUTDIR="${2:?Usage: sbatch --array=0-N%BATCH run_stage0_array.sbatch BAM_LIST OUTDIR}"

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

# ============================================================
# Fixed paths
# ============================================================
FASTA="/home/_groups/co/reference/hg19/hg19.fa"
BED="/storage/scratch01/groups/co/cn_extra/alleleSpecific/tso_samples/snp6_nochrM_fixed.bed"
S3_BUCKET="s3://cnio-data-sharing"
MIN_BQ=30
MIN_MQ=30

# ============================================================
# Resolve this task's BAM entry
# ============================================================
LINE_NUM=$(( SLURM_ARRAY_TASK_ID + 1 ))
BAM_ENTRY=$(sed -n "${LINE_NUM}p" "$BAM_LIST")

[[ -z "$BAM_ENTRY" ]] && { echo "[ERROR] No entry at line ${LINE_NUM} in ${BAM_LIST}"; exit 1; }

# Detect source type
if [[ "$BAM_ENTRY" == s3://* ]]; then
    SOURCE="s3"
    S3_URI="$BAM_ENTRY"
elif [[ "$BAM_ENTRY" == /* || "$BAM_ENTRY" == ./* ]]; then
    SOURCE="local"
    LOCAL_BAM="$BAM_ENTRY"
else
    # Bare S3 key (no s3:// prefix)
    SOURCE="s3"
    S3_URI="${S3_BUCKET}/${BAM_ENTRY}"
fi

BAM_BASENAME="$(basename "${S3_URI:-$LOCAL_BAM}")"
# Strip .bam, then -EXPxxY, then -TSO500 — preserves modifiers like -HRD
SAMPLE=$(echo "${BAM_BASENAME%.bam}" | sed -E 's/-EXP[0-9]+[A-Z]?//; s/-TSO500//')

# ============================================================
# Skip if _qc.tsv already exists
# ============================================================
mkdir -p "$OUTDIR" logs
QC_TSV="${OUTDIR}/${SAMPLE}_qc.tsv"

if [[ -f "$QC_TSV" ]]; then
    echo "[SKIP] $SAMPLE — _qc.tsv already exists"
    exit 0
fi

# ============================================================
# Reproducibility log
# ============================================================
echo "============================================================"
echo "DATE:    $(date)"
echo "JOB ID:  ${SLURM_JOB_ID:-NA}  (array ${SLURM_ARRAY_JOB_ID:-NA}[${SLURM_ARRAY_TASK_ID}])"
echo "SAMPLE:  $SAMPLE"
echo "SOURCE:  $SOURCE"
[[ "$SOURCE" == "s3"    ]] && echo "S3_URI:  $S3_URI"
[[ "$SOURCE" == "local" ]] && echo "BAM:     $LOCAL_BAM"
echo "OUTDIR:  $OUTDIR"
echo "============================================================"

bcftools --version | head -n 1
samtools --version | head -n 1

[[ -f "$FASTA" ]] || { echo "[ERROR] FASTA missing: $FASTA"; exit 1; }
[[ -f "$BED"   ]] || { echo "[ERROR] BED missing: $BED";   exit 1; }

# ============================================================
# For S3: download to temp dir, delete on exit
# For local: use in place, no cleanup needed
# ============================================================
if [[ "$SOURCE" == "s3" ]]; then
    aws --version
    BAM_TMPDIR="/storage/scratch01/groups/co/cn_extra/alleleSpecific/tso_samples/bam_tmp_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
    cleanup() { echo "[INFO] Removing temp BAM dir: $BAM_TMPDIR"; rm -rf "$BAM_TMPDIR"; }
    trap cleanup EXIT
    mkdir -p "$BAM_TMPDIR"

    LOCAL_BAM="${BAM_TMPDIR}/${BAM_BASENAME}"
    echo "[INFO] Downloading BAM from S3..."
    aws s3 cp "$S3_URI" "$LOCAL_BAM"
    aws s3 cp "${S3_URI}.bai" "${LOCAL_BAM}.bai" 2>/dev/null || samtools index "$LOCAL_BAM"
    echo "[INFO] BAM size: $(du -sh "$LOCAL_BAM" | cut -f1)"
else
    [[ -f "$LOCAL_BAM" ]] || { echo "[ERROR] BAM not found: $LOCAL_BAM"; exit 1; }
    # Ensure index exists
    [[ -f "${LOCAL_BAM}.bai" ]] || samtools index "$LOCAL_BAM"
fi

# ============================================================
# Step 1: mpileup
# ============================================================
OUT_PILEUP="${OUTDIR}/${SAMPLE}_allele_counts_snp6_wrapper_qualityMetrics_q30Pass.tsv"
OUT_STATS="${OUTDIR}/${SAMPLE}_bam_global_stats_qualityMetrics_q30Pass.txt"

echo "[INFO] Running bcftools mpileup..."
bcftools mpileup \
    -f "$FASTA" \
    -R "$BED" \
    --min-BQ "$MIN_BQ" \
    --min-MQ "$MIN_MQ" \
    -a FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/SP,INFO/MQBZ,INFO/BQBZ \
    -Ou "$LOCAL_BAM" \
| bcftools query \
    -f '%CHROM\t%POS\t%REF\t%ALT\t[%DP]\t[%AD]\t[%ADF]\t[%ADR]\t[%SP]\t%MQBZ\t%BQBZ\n' \
> "$OUT_PILEUP"
echo "[INFO] Pileup written: $OUT_PILEUP"

echo "[INFO] Running samtools stats..."
samtools stats -@ "${SLURM_CPUS_PER_TASK}" "$LOCAL_BAM" > "$OUT_STATS"
echo "[INFO] BAM stats written: $OUT_STATS"

# ============================================================
# Step 2: parse mpileup → _qc.tsv
# ============================================================
echo "[INFO] Parsing mpileup TSV..."
conda run -n shapeit4 python "${REPO_DIR}/pipeline/parse_mpileup.py" "$OUT_PILEUP" "$QC_TSV"
echo "[INFO] _qc.tsv written: $QC_TSV"

echo "[INFO] Job finished at $(date)"
echo "============================================================"
