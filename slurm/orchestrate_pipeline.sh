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
#   --bam-dir DIR            Directory with *.bam files (default: tso_samples/rg_bam)
#   --tsv-dir DIR            Directory where _qc.tsv files are/will be written
#   --segments-dir DIR       Pre-existing per-sample segments directory (reuse a prior run's dir)
#   --purity-file FILE       Purity/ploidy TSV (columns: sample, purity, ploidy).
#                            When provided, a purity-corrected run is submitted IN PARALLEL
#                            with the base run (same Stage 0 dependency), outputting to
#                            tso_samples/YYYYMMDD_purity/. Uses purity-calibrated thresholds.
#   --purity-thresholds FILE Override default purity-calibrated thresholds JSON.
#                            Default: tso_samples/20260406_purity/thresholds_analysis_purity/thresholds_lookup.json
#
# Examples:
#   # Default (existing samples)
#   bash orchestrate_pipeline.sh
#
#   # New batch with purity correction
#   bash orchestrate_pipeline.sh \
#     --bam-dir /storage/.../tso_samples/all_bam_tso \
#     --tsv-dir /storage/.../tso_samples/YYYYMMDD_mpileup_results \
#     --purity-file /storage/.../tso_samples/purity_ploidy_pipeline.tsv
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

# Purity correction: if PURITY_FILE is set, a second pipeline run is submitted in
# parallel using purity-corrected BAF and purity-calibrated thresholds.
# Outputs go to tso_samples/YYYYMMDD_purity/ alongside the base run.
PURITY_FILE=""
TSO_SAMPLES_ROOT_DEFAULT="/storage/scratch01/groups/co/cn_extra/alleleSpecific/tso_samples"
PURITY_THRESHOLDS="${TSO_SAMPLES_ROOT_DEFAULT}/20260406_purity/thresholds_analysis_purity/thresholds_lookup.json"

# Combined TSO500 segment file — all samples in one CSV.
# split_segments.py reads this and writes one CSV per sample to SEGMENTS_DIR.
COMBINED_SEGMENTS="/storage/scratch01/groups/co/cn_extra/alleleSpecific/tso_samples/unique_copyright_100kb_abs_segtable_free_purity_filtered_unique.csv"

# ============================================================
# Parse CLI arguments
# ============================================================
while [[ $# -gt 0 ]]; do
    case "$1" in
        --bam-dir)           BAM_DIR="$2";           shift 2 ;;
        --tsv-dir)           TSV_DIR="$2";           shift 2 ;;
        --segments-dir)      SEGMENTS_DIR="$2";      shift 2 ;;
        --purity-file)       PURITY_FILE="$2";       shift 2 ;;
        --purity-thresholds) PURITY_THRESHOLDS="$2"; shift 2 ;;
        *) echo "[ERROR] Unknown argument: $1"; exit 1 ;;
    esac
done

# Repo root (needed by Stage 0 sbatch exports — must be set before the submission loop)
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

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

# QC report sbatch to submit once Stage 4 completes.
QC_SBATCH="${SCRIPT_DIR}/run_qc_report.sbatch"

# HRD performance report sbatch to submit once Stage 4 completes.
EVAL_HRD_SBATCH="${SCRIPT_DIR}/run_eval_hrd.sbatch"

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
declare -A index_jid_for_bam  # BAM path → index_bam job ID (only when indexing needed)
n_skip=0
n_bam=0
n_tsv_only=0
n_index=0

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

# ============================================================
# Pre-flight: BAI index validation
# Submit index_bam jobs for any BAM whose .bai is missing or
# empty (0 bytes — corrupt transfer stub).  Stage 0 mpileup
# jobs for those samples depend on their index job.
# ============================================================
echo "------------------------------------------------------------"
echo "Pre-flight: BAI validation"
echo "------------------------------------------------------------"

for sample in $(echo "${!bam_for_sample[@]}" | tr ' ' '\n' | sort); do
    bam="${bam_for_sample[$sample]}"
    bai="${bam}.bai"

    if [[ -s "$bai" ]]; then
        continue  # valid index
    fi

    if [[ ! -f "$bai" ]]; then
        reason="missing"
    else
        reason="empty (0 bytes)"
    fi

    echo "  [INDEX]  $sample — BAI $reason"
    jid_index=$(sbatch \
        --parsable \
        "${SCRIPT_DIR}/run_index_bam.sbatch" \
        "$bam")
    index_jid_for_bam["$bam"]="$jid_index"
    echo "           index job  → $jid_index"
    (( n_index++ )) || true
done

if [[ $n_index -eq 0 ]]; then
    echo "  [OK] All BAI files valid — no indexing needed."
fi
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
    # Add index dependency only for BAMs that needed (re-)indexing
    mpileup_dep_flag=""
    if [[ -n "${index_jid_for_bam[$bam]+_}" ]]; then
        mpileup_dep_flag="--dependency=afterok:${index_jid_for_bam[$bam]}"
    fi

    jid_mpileup=$(sbatch \
        --parsable \
        $mpileup_dep_flag \
        --export=ALL,BAM_DIR="$BAM_DIR" \
        "${SCRIPT_DIR}/run_mpileup_single.sbatch" \
        "$sample" \
        "$TSV_DIR")

    if [[ -n "$mpileup_dep_flag" ]]; then
        echo "           mpileup job → $jid_mpileup  (dep: ${index_jid_for_bam[$bam]})"
    else
        echo "           mpileup job → $jid_mpileup"
    fi

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
echo "BAI index jobs         : $n_index"
echo "Stage 0 jobs submitted : $n_bam"
echo "Samples skipped (TSV)  : $n_skip"
echo "------------------------------------------------------------"

# ============================================================
# Generate run config
# ============================================================
RUN_DATE="$(date +%Y%m%d)"
TSO_SAMPLES_ROOT="${TSO_SAMPLES_ROOT_DEFAULT}"
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

# Write run_info.json for pipeline traceability (read by log_processed_samples.py)
RUN_INFO="${REPO_DIR}/configs/run_info_${RUN_DATE}.json"
if [[ ! -f "$RUN_INFO" ]]; then
    GIT_COMMIT="$(git -C "$REPO_DIR" rev-parse HEAD 2>/dev/null || echo 'unknown')"
    printf '{\n  "run_date": "%s",\n  "git_commit": "%s",\n  "config_file": "%s"\n}\n' \
        "$RUN_DATE" "$GIT_COMMIT" "$RUN_CONFIG" > "$RUN_INFO"
    echo "[INFO] Written run_info: $RUN_INFO"
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
    dep_str=""
    pipeline_jid=$(sbatch --parsable \
        --export=ALL,REPO_DIR="$REPO_DIR",PIPELINE_CONFIG="$RUN_CONFIG" "$PIPELINE_SBATCH")
fi

echo "Pipeline job (base) → $pipeline_jid"

# ============================================================
# Purity-corrected run (parallel to base when --purity-file set)
# ============================================================
purity_scarhrd_jid=""

if [[ -n "$PURITY_FILE" ]]; then
    PURITY_CONFIG="${REPO_DIR}/configs/pipeline_config_${RUN_DATE}_purity.yaml"
    PURITY_BASE_DIR="${TSO_SAMPLES_ROOT}/${RUN_DATE}_purity"

    if [[ ! -f "$PURITY_CONFIG" ]]; then
        sed \
            -e "s|base_dir: null|base_dir: \"${PURITY_BASE_DIR}\"|" \
            -e "s|purity_file: null.*|purity_file: \"${PURITY_FILE}\"|" \
            -e "s|thresholds_file:.*thresholds_lookup\.json.*|thresholds_file: \"${PURITY_THRESHOLDS}\"|" \
            "$RUN_CONFIG" > "$PURITY_CONFIG"
        echo "[INFO] Generated purity config: $PURITY_CONFIG"
    else
        echo "[INFO] Re-using existing purity config: $PURITY_CONFIG"
    fi

    if [[ -n "$dep_str" ]]; then
        purity_pipeline_jid=$(sbatch --parsable --dependency="$dep_str" \
            --export=ALL,REPO_DIR="$REPO_DIR",PIPELINE_CONFIG="$PURITY_CONFIG" "$PIPELINE_SBATCH")
    else
        purity_pipeline_jid=$(sbatch --parsable \
            --export=ALL,REPO_DIR="$REPO_DIR",PIPELINE_CONFIG="$PURITY_CONFIG" "$PIPELINE_SBATCH")
    fi
    echo "Pipeline job (purity) → $purity_pipeline_jid  [parallel, dep: ${dep_str:-none}]"

    purity_scarhrd_jid=$(sbatch --parsable \
        --dependency=afterok:"$purity_pipeline_jid" \
        --export=ALL,REPO_DIR="$REPO_DIR",BASE_DIR="$PURITY_BASE_DIR" \
        "$SCARHRD_SBATCH")
    echo "scarHRD job  (purity) → $purity_scarhrd_jid  (dep: $purity_pipeline_jid)"
fi

# ============================================================
# Submit scarHRD (Stage 4) dependent on base pipeline completion
# ============================================================
BASE_DIR="${TSO_SAMPLES_ROOT}/${RUN_DATE}"

scarhrd_jid=$(sbatch --parsable \
    --dependency=afterok:"$pipeline_jid" \
    --export=ALL,REPO_DIR="$REPO_DIR",BASE_DIR="$BASE_DIR" \
    "$SCARHRD_SBATCH")

echo "scarHRD job  (base)   → $scarhrd_jid  (dep: $pipeline_jid)"

qc_jid=$(sbatch --parsable \
    --dependency=afterok:"$scarhrd_jid" \
    --export=ALL,REPO_DIR="$REPO_DIR",CONFIG="$RUN_CONFIG" \
    "$QC_SBATCH")

echo "QC report             → $qc_jid  (dep: $scarhrd_jid)"

# eval_hrd waits for both scarHRD jobs (base + purity if applicable)
if [[ -n "$purity_scarhrd_jid" ]]; then
    eval_dep="afterok:${scarhrd_jid}:${purity_scarhrd_jid}"
else
    eval_dep="afterok:${scarhrd_jid}"
fi

eval_jid=$(sbatch --parsable \
    --dependency="$eval_dep" \
    --export=ALL,REPO_DIR="$REPO_DIR",CURRENT_RUN="$RUN_DATE" \
    "$EVAL_HRD_SBATCH")

echo "HRD eval              → $eval_jid  (dep: ${eval_dep#afterok:})"
echo ""
echo "Monitor:  squeue -u \$USER"
echo "Logs:     ${SCRIPT_DIR}/logs/"
echo "============================================================"
