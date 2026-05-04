#!/bin/bash
# ============================================================
# orchestrate_pipeline.sh
#
# Unified entry-point for the BAF-to-HRD pipeline.
# Accepts samples as either:
#   A) Already-parsed TSV files (_qc.tsv, ASCAT .txt.gz, etc.)
#      → submitted directly to the pipeline (Stage 1+)
#   B) BAM files
#      → Stage 0a (mpileup) + Stage 0b (parse_mpileup) are
#        submitted first; the pipeline starts as a SLURM
#        dependency so it only runs after all BAM parsing is done
#
# Skip logic: if a sample already has a _qc.tsv in TSV_DIR,
# Stage 0 is skipped for that sample.
#
# Usage:
#   bash orchestrate_pipeline.sh
#
# Inspect submitted jobs:
#   squeue -u $USER
# ============================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ============================================================
# USER CONFIGURATION — edit these paths before running
# ============================================================

# Directory that contains *.mdup.rg.bam files
BAM_DIR="/storage/scratch01/groups/co/cn_extra/alleleSpecific/tso_samples/rg_bam"

# Directory where _qc.tsv files live (already parsed) or will be written (from BAMs)
# The pipeline (Stage 1+) reads from this directory.
TSV_DIR="/storage/scratch01/groups/co/cn_extra/alleleSpecific/tso_samples/20260225_GoldStandard_mpileup_stats_q30/results"

# BAM filename suffix patterns to try (EXP08A takes priority over EXP08B)
BAM_SUFFIXES=("-TSO500-EXP08A.mdup.rg.bam" "-TSO500-EXP08B.mdup.rg.bam")

# Known TSV suffixes that the pipeline can read directly (order = priority)
KNOWN_TSV_SUFFIXES=(
    "_qc.tsv"
    ".tumour.tso500_sparse_downsample_corrected.txt.gz"
    ".tumour.tso500_sparse_downsample.txt"
)

# Pipeline sbatch to submit once all Stage 0 jobs complete.
# Must be a self-contained sbatch that runs baf_to_hrd_pipeline.py.
PIPELINE_SBATCH="${SCRIPT_DIR}/run_pipeline.sbatch"

# Maximum simultaneous Stage 0 array tasks (throttles mpileup memory)
MAX_PARALLEL_MPILEUP=10

# ============================================================
# Helper: derive sample name from BAM path
# ============================================================
sample_from_bam() {
    local bam="$1"
    local base
    base="$(basename "$bam")"
    for sfx in "${BAM_SUFFIXES[@]}"; do
        if [[ "$base" == *"$sfx" ]]; then
            echo "${base%$sfx}"
            return
        fi
    done
    # Fallback: strip from last -TSO500-
    echo "${base%-TSO500-*}"
}

# ============================================================
# Helper: check if a sample already has a ready TSV
# ============================================================
has_tsv() {
    local sample="$1"
    for sfx in "${KNOWN_TSV_SUFFIXES[@]}"; do
        if [[ -f "${TSV_DIR}/${sample}${sfx}" ]]; then
            return 0
        fi
    done
    return 1
}

# ============================================================
# Sanity checks
# ============================================================
[[ -d "$BAM_DIR" ]] || { echo "[ERROR] BAM_DIR not found: $BAM_DIR"; exit 1; }
[[ -d "$TSV_DIR" ]] || mkdir -p "$TSV_DIR"
mkdir -p "${SCRIPT_DIR}/logs"

if [[ ! -f "$PIPELINE_SBATCH" ]]; then
    echo "[WARN] Pipeline sbatch not found: $PIPELINE_SBATCH"
    echo "       Stage 0 jobs will be submitted but the pipeline will NOT be"
    echo "       auto-chained.  Create run_pipeline.sbatch first, then re-run."
    PIPELINE_SBATCH=""
fi

# ============================================================
# Scan BAM directory
# ============================================================
echo "============================================================"
echo "BAF-to-HRD Pipeline Orchestrator"
echo "DATE:     $(date)"
echo "BAM_DIR:  $BAM_DIR"
echo "TSV_DIR:  $TSV_DIR"
echo "============================================================"

stage0_jids=()        # parse_mpileup job IDs (one per BAM-only sample)
n_skip=0
n_bam=0
n_tsv_only=0

# Collect all BAM files (deduplicate per sample: prefer EXP08A over EXP08B)
declare -A bam_for_sample
for bam in "${BAM_DIR}"/*.mdup.rg.bam; do
    [[ -f "$bam" ]] || continue
    sample="$(sample_from_bam "$bam")"
    # Only record the first match per sample (EXP08A wins because the glob is sorted)
    if [[ -z "${bam_for_sample[$sample]+_}" ]]; then
        bam_for_sample["$sample"]="$bam"
    fi
done

if [[ ${#bam_for_sample[@]} -eq 0 ]]; then
    echo "[INFO] No BAM files found in $BAM_DIR — nothing to do for Stage 0."
    echo "       If all samples already have TSV files the pipeline can be"
    echo "       submitted directly:"
    echo "         sbatch $PIPELINE_SBATCH"
    exit 0
fi

echo ""
echo "Found ${#bam_for_sample[@]} unique sample(s) in BAM_DIR."
echo ""

for sample in $(echo "${!bam_for_sample[@]}" | tr ' ' '\n' | sort); do
    bam="${bam_for_sample[$sample]}"

    if has_tsv "$sample"; then
        echo "  [SKIP]   $sample  (TSV already exists)"
        (( n_skip++ )) || true
        (( n_tsv_only++ )) || true
        continue
    fi

    echo "  [STAGE0] $sample"
    (( n_bam++ )) || true

    # ── Stage 0a: mpileup ──────────────────────────────────────────────────
    jid_mpileup=$(sbatch \
        --parsable \
        "${SCRIPT_DIR}/run_mpileup_single.sbatch" \
        "$sample" \
        "$TSV_DIR")

    echo "           mpileup job → $jid_mpileup"

    # ── Stage 0b: parse mpileup → _qc.tsv ─────────────────────────────────
    jid_parse=$(sbatch \
        --parsable \
        --dependency=afterok:"$jid_mpileup" \
        "${SCRIPT_DIR}/run_parse_mpileup.sbatch" \
        "$sample" \
        "$TSV_DIR")

    echo "           parse job   → $jid_parse  (dep: $jid_mpileup)"
    stage0_jids+=("$jid_parse")
done

# ============================================================
# Summary
# ============================================================
echo ""
echo "------------------------------------------------------------"
echo "Stage 0 jobs submitted : $n_bam"
echo "Samples skipped (TSV)  : $n_skip"
echo "------------------------------------------------------------"

# ============================================================
# Submit pipeline (Stage 1+)
# ============================================================
if [[ -z "$PIPELINE_SBATCH" ]]; then
    echo "[WARN] No pipeline sbatch configured — skipping pipeline submission."
    echo "       Create run_pipeline.sbatch and re-run, or submit manually."
    exit 0
fi

if [[ ${#stage0_jids[@]} -gt 0 ]]; then
    dep_str="afterok:$(IFS=':'; echo "${stage0_jids[*]}")"
    echo ""
    echo "Submitting pipeline with dependency: $dep_str"
    pipeline_jid=$(sbatch --parsable --dependency="$dep_str" "$PIPELINE_SBATCH")
else
    echo ""
    echo "All samples already have TSVs — submitting pipeline immediately."
    pipeline_jid=$(sbatch --parsable "$PIPELINE_SBATCH")
fi

echo "Pipeline job → $pipeline_jid"
echo ""
echo "Monitor:  squeue -u \$USER"
echo "Logs:     ${SCRIPT_DIR}/logs/"
echo "============================================================"
