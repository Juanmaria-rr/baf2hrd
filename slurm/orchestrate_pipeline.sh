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
#   bash orchestrate_pipeline.sh [OPTIONS]
#
# Options:
#   --bam-dir DIR       Directory with *.bam files (default: tso_samples/rg_bam)
#   --tsv-dir DIR       Directory where _qc.tsv files are/will be written
#                       (default: 20260225_GoldStandard_mpileup_stats_q30/results)
#   --segments-dir DIR  Pre-existing per-sample segments directory; skips
#                       split_segments.py output writes (reuse a prior run's dir)
#
# Examples:
#   # Default (existing samples)
#   bash orchestrate_pipeline.sh
#
#   # New batch
#   bash orchestrate_pipeline.sh \
#     --bam-dir /storage/.../tso_samples/rg_bam_20260505_batch01 \
#     --tsv-dir /storage/.../tso_samples/20260505_batch01_mpileup/results
#
# Inspect submitted jobs:
#   squeue -u $USER
# ============================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ============================================================
# DEFAULTS — overridable via CLI arguments
# ============================================================

# Directory that contains *.bam files
BAM_DIR="/storage/scratch01/groups/co/cn_extra/alleleSpecific/tso_samples/rg_bam"

# Directory where _qc.tsv files live (already parsed) or will be written (from BAMs)
# The pipeline (Stage 1+) reads from this directory.
TSV_DIR="/storage/scratch01/groups/co/cn_extra/alleleSpecific/tso_samples/20260225_GoldStandard_mpileup_stats_q30/results"

# Per-sample segments directory. Empty = let split_segments.py create a dated dir.
# Pass --segments-dir to reuse an existing one (skips split_segments.py writes).
SEGMENTS_DIR=""

# Combined TSO500 segment file — all samples in one CSV.
# split_segments.py reads this and writes one CSV per sample to SEGMENTS_DIR.
COMBINED_SEGMENTS="/storage/scratch01/groups/co/cn_extra/alleleSpecific/tso_samples/unique_copyright_100kb_abs_segtable_free_purity_filtered_unique.csv"

# ============================================================
# Parse CLI arguments
# ============================================================
while [[ $# -gt 0 ]]; do
    case "$1" in
        --bam-dir)      BAM_DIR="$2";      shift 2 ;;
        --tsv-dir)      TSV_DIR="$2";      shift 2 ;;
        --segments-dir) SEGMENTS_DIR="$2"; shift 2 ;;
        *) echo "[ERROR] Unknown argument: $1"; exit 1 ;;
    esac
done

# Python interpreter with pandas available
PYTHON="conda run -n shapeit4 python"

# Path to split_segments.py
SPLIT_SEGMENTS_PY="${SCRIPT_DIR}/../pipeline/split_segments.py"

# BAM filename suffix patterns to try (EXP08A takes priority over EXP08B).
# Handles both legacy mdup.rg naming and plain .bam naming from new transfers.
# HRD variants listed before plain so sample_from_bam yields SAMPLE-HRD correctly.
BAM_SUFFIXES=(
    "-TSO500-EXP08A.mdup.rg.bam"
    "-TSO500-EXP08B.mdup.rg.bam"
    "-TSO500-HRD-EXP08A.bam"
    "-TSO500-HRD-EXP08B.bam"
    "-TSO500-EXP08A.bam"
    "-TSO500-EXP08B.bam"
)

# Known TSV suffixes that the pipeline can read directly (order = priority)
KNOWN_TSV_SUFFIXES=(
    "_qc.tsv"
    ".tumour.tso500_sparse_downsample_corrected.txt.gz"
    ".tumour.tso500_sparse_downsample.txt"
)

# Pipeline sbatch to submit once all Stage 0 jobs complete.
# Must be a self-contained sbatch that runs baf_to_hrd_pipeline.py.
PIPELINE_SBATCH="${SCRIPT_DIR}/run_pipeline.sbatch"

# scarHRD sbatch to submit once Stages 1-3 complete.
SCARHRD_SBATCH="${SCRIPT_DIR}/run_scarHRD.sbatch"

# Maximum simultaneous Stage 0 array tasks (throttles mpileup memory)
MAX_PARALLEL_MPILEUP=10

# ============================================================
# Helper: derive sample name from BAM path
# Non-HRD: JBLAB17012-TSO500-EXP08B.{mdup.rg.,}bam  → JBLAB17012
# HRD:     JBLAB17012-TSO500-HRD-EXP08B.bam           → JBLAB17012-HRD
# ============================================================
sample_from_bam() {
    local bam="$1"
    local base
    base="$(basename "$bam")"

    # New plain-naming: regex captures JBLAB_ID and optional -HRD
    if [[ "$base" =~ ^(JBLAB[^-]+-TSO500(-HRD)?)-EXP[0-9]+[AB]\.bam$ ]]; then
        # ${BASH_REMATCH[1]} = JBLAB17012-TSO500[-HRD]; strip trailing -TSO500[-HRD]
        local prefix="${BASH_REMATCH[1]}"
        local hrd="${BASH_REMATCH[2]}"   # "-HRD" or ""
        echo "${prefix%-TSO500*}${hrd}"
        return
    fi

    # Legacy mdup.rg naming: JBLAB17012-TSO500-EXP08B.mdup.rg.bam
    for sfx in "${BAM_SUFFIXES[@]}"; do
        if [[ "$base" == *"$sfx" ]]; then
            echo "${base%$sfx}"
            return
        fi
    done

    # Last resort
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
# Pre-processing: split combined segments → per-sample CSVs
# (skipped when --segments-dir is passed — reuses an existing dir)
# ============================================================
echo "============================================================"
echo "Pre-processing: split_segments"
echo "============================================================"

if [[ -n "$SEGMENTS_DIR" ]]; then
    [[ -d "$SEGMENTS_DIR" ]] || { echo "[ERROR] --segments-dir not found: $SEGMENTS_DIR"; exit 1; }
    echo "[INFO] Reusing existing segments dir (no writes): $SEGMENTS_DIR"
else
    if [[ ! -f "$COMBINED_SEGMENTS" ]]; then
        echo "[ERROR] Combined segments file not found: $COMBINED_SEGMENTS"
        exit 1
    fi

    SPLIT_ARGS=(--all --combined "$COMBINED_SEGMENTS")

    # Capture the raw_segments_dir printed by split_segments.py
    SPLIT_OUTPUT=$($PYTHON "$SPLIT_SEGMENTS_PY" "${SPLIT_ARGS[@]}" 2>&1)
    echo "$SPLIT_OUTPUT"

    # Extract the output directory from the last raw_segments_dir= line
    SEGMENTS_DIR=$(echo "$SPLIT_OUTPUT" | grep "^\\[split_segments\\] raw_segments_dir=" | tail -1 | cut -d= -f2)

    if [[ -z "$SEGMENTS_DIR" ]]; then
        echo "[ERROR] Could not determine SEGMENTS_DIR from split_segments.py output."
        exit 1
    fi
fi

echo "[INFO] raw_segments_dir resolved to: $SEGMENTS_DIR"
echo ""

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
for bam in "${BAM_DIR}"/*.bam; do
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
        --export=ALL,BAM_DIR="$BAM_DIR" \
        "${SCRIPT_DIR}/run_mpileup_single.sbatch" \
        "$sample" \
        "$TSV_DIR")

    echo "           mpileup job → $jid_mpileup"

    # ── Stage 0b: parse mpileup → _qc.tsv ─────────────────────────────────
    jid_parse=$(sbatch \
        --parsable \
        --dependency=afterok:"$jid_mpileup" \
        --export=ALL,REPO_DIR="$REPO_DIR" \
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
# Generate run config
# ============================================================
RUN_DATE="$(date +%Y%m%d)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
TSO_SAMPLES_ROOT="/storage/scratch01/groups/co/cn_extra/alleleSpecific/tso_samples"
RUN_CONFIG="${REPO_DIR}/configs/pipeline_config_${RUN_DATE}.yaml"

# Only write if it does not already exist (idempotent re-runs same day)
if [[ ! -f "$RUN_CONFIG" ]]; then
    TEMPLATE="${REPO_DIR}/configs/pipeline_config_template.yaml"
    sed \
        -e "s|date: \"YYYYMMDD\"|date: \"${RUN_DATE}\"|" \
        -e "s|raw_segments_dir:.*|raw_segments_dir: \"${SEGMENTS_DIR}\"|" \
        -e "s|tso_samples/RESULTS_DIR|tso_samples/$(basename "${TSV_DIR}")|" \
        -e "s|THRESHOLDS_DATE|20260406|" \
        "$TEMPLATE" > "$RUN_CONFIG"
    echo "[INFO] Generated run config: $RUN_CONFIG"
else
    echo "[INFO] Re-using existing run config: $RUN_CONFIG"
fi

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
    pipeline_jid=$(sbatch --parsable --dependency="$dep_str" \
        --export=ALL,REPO_DIR="$REPO_DIR",PIPELINE_CONFIG="$RUN_CONFIG" "$PIPELINE_SBATCH")
else
    echo ""
    echo "All samples already have TSVs — submitting pipeline immediately."
    pipeline_jid=$(sbatch --parsable \
        --export=ALL,REPO_DIR="$REPO_DIR",PIPELINE_CONFIG="$RUN_CONFIG" "$PIPELINE_SBATCH")
fi

echo "Pipeline job → $pipeline_jid"

# ============================================================
# Submit scarHRD (Stage 4) dependent on pipeline completion
# ============================================================
BASE_DIR="${TSO_SAMPLES_ROOT}/${RUN_DATE}"

scarhrd_jid=$(sbatch --parsable \
    --dependency=afterok:"$pipeline_jid" \
    --export=ALL,REPO_DIR="$REPO_DIR",BASE_DIR="$BASE_DIR" \
    "$SCARHRD_SBATCH")

echo "scarHRD job → $scarhrd_jid  (dep: $pipeline_jid)"
echo ""
echo "Monitor:  squeue -u \$USER"
echo "Logs:     ${SCRIPT_DIR}/logs/"
echo "============================================================"
