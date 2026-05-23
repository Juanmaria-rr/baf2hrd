#!/usr/bin/env python3
"""
BAF-to-HRD Pipeline — Integrated

Three sequential stages:
  1. Per-sample BAF analysis → canonical segment + BAF-hit files + manifest
  2. Consolidation          → merged global dataset with enhanced features
  3. Allele-specific CN inference → scarHRD inputs per sample × DP filter

Usage:
  Sequential (all samples):    python baf_to_hrd_pipeline.py
  Parallel / SLURM (1 sample): python baf_to_hrd_pipeline.py <sample_id> <long_name> <folder_path>
"""

# ============================================================
# IMPORTS
# ============================================================

import os
import gc
import glob
import sys
import fcntl
import json
import shutil
import traceback
from pathlib import Path
from datetime import datetime
from contextlib import contextmanager

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from scipy.stats import pearsonr

from pptx import Presentation
from pptx.util import Inches
from PIL import Image

Image.MAX_IMAGE_PIXELS = 300000000

# ============================================================
# CONFIGURATION DEFAULTS
# ============================================================
# These values are overridden at runtime by apply_config() when main() is
# called.  OPERATIONS_DIR is the exception: it is used for module-level
# imports below and cannot be changed via --config at runtime.

DATE          = "20260421"
STAGES_TO_RUN = {1, 2, 3}
BASE_DIR      = f"/storage/scratch01/groups/co/cn_extra/alleleSpecific/tso_samples/{DATE}"

PURITY_FILE          = None
DP_VALUES            = [0, 2, 4, 8]
N_BINS               = 20
# TYPES_OF_FILES labels each processing pass per sample.  One entry = one pass
# (output tagged as "original_tso").  BAF_SUFFIXES (below) controls which file
# on disk is actually loaded — the two are independent.
TYPES_OF_FILES       = ["_qc.tsv"]
OPERATIONS_DIR       = os.path.dirname(os.path.abspath(__file__))
RAW_SEGMENTS_DIR     = "/storage/scratch01/groups/co/cn_extra/alleleSpecific/tso_samples/20260416"
THRESHOLDS_FILE      = "/storage/scratch01/groups/co/cn_extra/alleleSpecific/tso_samples/20260406/thresholds_analysis/thresholds_lookup.json"
MAX_GAP_BP_FOR_MERGE = 1

# Supported BAF file suffixes — order determines priority in multi-suffix folders.
# TSO:   produced by SLURM mpileup → awk pipeline (_qc.tsv).
# ASCAT: produced by correct_REFALT_downsampled_flagProblematic.py
#        (.tumour.tso500_sparse_downsample_corrected.txt.gz)
#        Older ASCAT runs may have the uncompressed .txt variant.
BAF_SUFFIXES = [
    "_qc.tsv",
    ".tumour.tso500_sparse_downsample_corrected.txt.gz",
    ".tumour.tso500_sparse_downsample.txt",
]

# ASCAT segment files live in a separate directory from TSO segment files.
# Set to a directory containing *nMajor/nMinor* CSV files for ASCAT samples;
# leave None to search only RAW_SEGMENTS_DIR.
ASCAT_SEGMENTS_DIR   = None

FOLDER_PATHS_DEFAULT = [
    "/storage/scratch01/groups/co/cn_extra/alleleSpecific/tso_samples/20260225_GoldStandard_mpileup_stats_q30/results"
]

# Derived paths — recomputed by apply_config()
OUTPUT_DIR               = BASE_DIR
UNIFIED_LOG_FILE         = os.path.join(BASE_DIR, "pipeline_unified.log")
SEGMENTS_PER_SAMPLE_DIR  = os.path.join(BASE_DIR, "segments_per_sample")
BAF_HITS_PER_SAMPLE_DIR  = os.path.join(BASE_DIR, "bafs_hits_per_sample")
BAF_CORRECTION_PLOTS_DIR = os.path.join(BASE_DIR, "baf_correction_plots")
SCARHRD_INPUT_DIR        = os.path.join(BASE_DIR, "scarHRD_inputs")
SCARHRD_THRESHOLDS_DIR   = os.path.join(BASE_DIR, "scarHRD_inputs_thresholds")
SCRIPTS_BACKUP_DIR       = os.path.join(BASE_DIR, "scripts_backup")
MANIFEST_FILE            = os.path.join(BASE_DIR, "pipeline_manifest.csv")
MERGED_OUTPUT_FILE       = os.path.join(BASE_DIR, "merged_segments_baf_metrics_ALL_with20bins_enhancedFeatures.csv")
SUMMARY_OUTPUT_FILE      = os.path.join(BASE_DIR, "merge_processing_summary.csv")

# Parallel-mode state — set by apply_config()
PARALLEL_MODE        = False
SAMPLE_ID_FILTER     = "MAIN"
LONG_NAME_FILTER     = None
FOLDER_PATH_OVERRIDE = None
FOLDER_PATHS         = FOLDER_PATHS_DEFAULT


def _parse_args():
    import argparse
    p = argparse.ArgumentParser(
        description="BAF-to-HRD Pipeline (Stages 1–3)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python baf_to_hrd_pipeline.py --config configs/run_20260421.yaml\n"
            "  python baf_to_hrd_pipeline.py --config configs/run_20260421.yaml --stages 3\n"
            "  python baf_to_hrd_pipeline.py --config configs/run_20260421.yaml SAMPLE LONGNAME /path\n"
        ),
    )
    p.add_argument("--config", metavar="YAML",     help="Pipeline config file (YAML)")
    p.add_argument("--date",   metavar="YYYYMMDD", help="Run date — overrides config")
    p.add_argument("--stages", metavar="N[,N]",    help="Stages to run, e.g. 2,3 — overrides config")
    p.add_argument("sample_id",   nargs="?", help="[parallel mode] Sample ID")
    p.add_argument("long_name",   nargs="?", help="[parallel mode] Long sample name")
    p.add_argument("folder_path", nargs="?", help="[parallel mode] Folder path containing BAF files")
    return p.parse_args()


def apply_config(args):
    """Populate all configuration globals from CLI args and optional YAML config."""
    import yaml

    global DATE, STAGES_TO_RUN, BASE_DIR, PURITY_FILE, DP_VALUES, N_BINS
    global TYPES_OF_FILES, RAW_SEGMENTS_DIR, THRESHOLDS_FILE, MAX_GAP_BP_FOR_MERGE
    global BAF_SUFFIXES, ASCAT_SEGMENTS_DIR, FOLDER_PATHS_DEFAULT
    global OUTPUT_DIR, UNIFIED_LOG_FILE, SEGMENTS_PER_SAMPLE_DIR
    global BAF_HITS_PER_SAMPLE_DIR, BAF_CORRECTION_PLOTS_DIR, SCARHRD_INPUT_DIR
    global SCARHRD_THRESHOLDS_DIR, SCRIPTS_BACKUP_DIR, MANIFEST_FILE
    global MERGED_OUTPUT_FILE, SUMMARY_OUTPUT_FILE
    global PARALLEL_MODE, SAMPLE_ID_FILTER, LONG_NAME_FILTER, FOLDER_PATH_OVERRIDE, FOLDER_PATHS

    cfg = {}
    if args.config:
        with open(args.config) as f:
            cfg = yaml.safe_load(f) or {}

    def get(section, key, default):
        return cfg.get(section, {}).get(key, default)

    DATE = args.date or get("run", "date", DATE)

    stages_raw = args.stages or get("run", "stages", None)
    STAGES_TO_RUN = (
        {int(s) for s in str(stages_raw).split(",")} if stages_raw else {1, 2, 3}
    )

    _tso_root = get(
        "paths", "tso_samples_root",
        "/storage/scratch01/groups/co/cn_extra/alleleSpecific/tso_samples",
    )
    BASE_DIR = get("run", "base_dir", None) or f"{_tso_root}/{DATE}"

    PURITY_FILE          = get("input",      "purity_file",          PURITY_FILE)
    DP_VALUES            = get("processing", "dp_values",            DP_VALUES)
    N_BINS               = get("processing", "n_bins",               N_BINS)
    TYPES_OF_FILES       = get("input",      "types_of_files",       TYPES_OF_FILES)
    RAW_SEGMENTS_DIR     = get("input",      "raw_segments_dir",     RAW_SEGMENTS_DIR)
    THRESHOLDS_FILE      = get("input",      "thresholds_file",      THRESHOLDS_FILE)
    MAX_GAP_BP_FOR_MERGE = get("processing", "max_gap_bp_for_merge", MAX_GAP_BP_FOR_MERGE)
    BAF_SUFFIXES         = get("input",      "baf_suffixes",         BAF_SUFFIXES)
    ASCAT_SEGMENTS_DIR   = get("input",      "ascat_segments_dir",   ASCAT_SEGMENTS_DIR)
    FOLDER_PATHS_DEFAULT = get("input",      "folder_paths",         FOLDER_PATHS_DEFAULT)

    OUTPUT_DIR               = BASE_DIR
    UNIFIED_LOG_FILE         = os.path.join(BASE_DIR, "pipeline_unified.log")
    SEGMENTS_PER_SAMPLE_DIR  = os.path.join(BASE_DIR, "segments_per_sample")
    BAF_HITS_PER_SAMPLE_DIR  = os.path.join(BASE_DIR, "bafs_hits_per_sample")
    BAF_CORRECTION_PLOTS_DIR = os.path.join(BASE_DIR, "baf_correction_plots")
    SCARHRD_INPUT_DIR        = os.path.join(BASE_DIR, "scarHRD_inputs")
    SCARHRD_THRESHOLDS_DIR   = os.path.join(BASE_DIR, "scarHRD_inputs_thresholds")
    SCRIPTS_BACKUP_DIR       = os.path.join(BASE_DIR, "scripts_backup")
    MANIFEST_FILE            = os.path.join(BASE_DIR, "pipeline_manifest.csv")
    MERGED_OUTPUT_FILE       = os.path.join(
        BASE_DIR, "merged_segments_baf_metrics_ALL_with20bins_enhancedFeatures.csv"
    )
    SUMMARY_OUTPUT_FILE      = os.path.join(BASE_DIR, "merge_processing_summary.csv")

    PARALLEL_MODE = args.sample_id is not None
    if PARALLEL_MODE:
        SAMPLE_ID_FILTER     = args.sample_id
        LONG_NAME_FILTER     = args.long_name
        FOLDER_PATH_OVERRIDE = args.folder_path
        FOLDER_PATHS         = [FOLDER_PATH_OVERRIDE]
    else:
        SAMPLE_ID_FILTER     = "MAIN"
        LONG_NAME_FILTER     = None
        FOLDER_PATH_OVERRIDE = None
        FOLDER_PATHS         = FOLDER_PATHS_DEFAULT

# ============================================================
# IMPORTS FROM functions.py
# ============================================================

if OPERATIONS_DIR in sys.path:
    sys.path.remove(OPERATIONS_DIR)
sys.path.insert(0, OPERATIONS_DIR)

functions_path = os.path.join(OPERATIONS_DIR, "functions.py")
if not os.path.exists(functions_path):
    raise SystemExit(f"ERROR: {functions_path} does not exist!")

try:
    from functions import (
        save_fig,
        safe_auc,
        assert_columns,
        prepare_baf_origin,
        add_title,
        build_track_only_slide,
        build_table_slide,
        build_slide_from_figure,
        create_roc_pr_hist_figure,
        create_input_data_scarHRD_justOdd,
    )
    print("✓ Imported mandatory functions from functions.py")
except ImportError as e:
    raise SystemExit(f"✗ Import error (mandatory): {e}")

try:
    from functions import (
        add_enhanced_features_to_merged,
        add_length_normalized_features,
    )
    ENHANCED_FEATURES_AVAILABLE = True
    print("✓ Imported enhanced-feature functions from functions.py")
except ImportError:
    ENHANCED_FEATURES_AVAILABLE = False
    print("⚠️  Enhanced features not available — proceeding without them")

    def add_enhanced_features_to_merged(merged, baf_temp):
        return merged

    def add_length_normalized_features(merged):
        return merged



# ============================================================
# UNIFIED LOGGER
# ============================================================

class UnifiedLogger:
    """Thread-safe logger writing to a shared log file + selective console output."""

    def __init__(self, log_file, sample_id="MAIN", verbose_console=False):
        self.log_file = log_file
        self.sample_id = sample_id
        self.verbose_console = verbose_console
        self.terminal = sys.stdout

        os.makedirs(os.path.dirname(log_file), exist_ok=True)

        self.important_keywords = [
            "✓", "⚠️", "ERROR", "WARNING",
            "===", "Sample:", "COMPLETE", "Processing DP",
            "FOLDER", "LOADING", "SEGMENT-LEVEL", "scarHRD",
        ]

    @contextmanager
    def _locked_file(self):
        f = open(self.log_file, "a")
        try:
            fcntl.flock(f.fileno(), fcntl.LOCK_EX)
            yield f
        finally:
            fcntl.flock(f.fileno(), fcntl.LOCK_UN)
            f.close()

    def write(self, message):
        if not message.strip():
            return
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        formatted_msg = f"{timestamp} [{self.sample_id:20s}] {message}"
        with self._locked_file() as f:
            f.write(formatted_msg + "\n")
        if self.verbose_console:
            self.terminal.write(message)
            if not message.endswith("\n"):
                self.terminal.write("\n")
        else:
            if any(kw in message for kw in self.important_keywords):
                self.terminal.write(message)
                if not message.endswith("\n"):
                    self.terminal.write("\n")

    def flush(self):
        self.terminal.flush()

    def close(self):
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        final_msg = f"\n{'='*70}\n[{self.sample_id}] Ended: {timestamp}\n{'='*70}\n"
        with self._locked_file() as f:
            f.write(final_msg)


# ============================================================
# SCRIPT 1 HELPER FUNCTIONS
# ============================================================

def build_output_stem(sample_name, file_type, baf_folder):
    return f"{sample_name}_{file_type}_{baf_folder}"


def create_sample_name_dict_custom(folder_path, max_samples=None):
    """
    Discover BAF files in folder_path for all suffixes in BAF_SUFFIXES.

    For each file found, derive:
      - long_name: filename stem (before the suffix), used to construct the
                   BAF file path and to look up segments/purity.
      - short_name: first two '_'-separated parts of long_name, used as the
                   human-readable sample ID.

    Returns {short_name: long_name}.
    """
    all_files = []
    for suffix in BAF_SUFFIXES:
        all_files.extend(glob.glob(os.path.join(folder_path, f"*{suffix}")))
    all_files = sorted(set(all_files))

    sample_dict = {}
    processed_long = set()

    for file_path in all_files:
        if max_samples is not None and len(processed_long) >= max_samples:
            break
        filename = os.path.basename(file_path)
        long_name = None
        for suffix in BAF_SUFFIXES:
            if filename.endswith(suffix):
                long_name = filename[: -len(suffix)]
                break
        if long_name is None or long_name in processed_long:
            continue
        parts = long_name.split("_")
        short_name = "_".join(parts[:2])
        # Avoid short_name collisions from different long_names
        if short_name in sample_dict:
            short_name = long_name  # fall back to full long_name as key
        sample_dict[short_name] = long_name
        processed_long.add(long_name)

    print(f"Processed {len(sample_dict)} samples")
    return sample_dict


def load_purity_mappings(purity_file):
    purity_mapping = {}
    ploidy_mapping = {}

    print(f"{'='*70}")
    print("LOADING PURITY METADATA (OPTIONAL)")
    print(f"{'='*70}")

    if purity_file is None:
        print("ℹ️  No purity file specified - skipping purity correction")
        print(f"{'='*70}\n")
        return purity_mapping, ploidy_mapping

    if not os.path.exists(purity_file):
        print(f"⚠️  Purity file NOT FOUND: {purity_file}")
        print("   Proceeding WITHOUT purity correction")
        print(f"{'='*70}\n")
        return purity_mapping, ploidy_mapping

    try:
        meta_df = pd.read_csv(purity_file, sep="\t")
        if len(meta_df) == 0:
            print("⚠️  Purity file is EMPTY - skipping purity correction")
            print(f"{'='*70}\n")
            return purity_mapping, ploidy_mapping

        meta_dict = meta_df.set_index("sample")[["purity", "ploidy"]].to_dict("index")
        for long_name, values in meta_dict.items():
            purity_mapping[long_name] = values["purity"]
            ploidy_mapping[long_name] = values["ploidy"]
            short_name = long_name.split("_vs_")[0]
            purity_mapping[short_name] = values["purity"]
            ploidy_mapping[short_name] = values["ploidy"]
            # HRD variants share the same tumor — register both forms
            purity_mapping[long_name  + "-HRD"] = values["purity"]
            ploidy_mapping[long_name  + "-HRD"] = values["ploidy"]
            purity_mapping[short_name + "-HRD"] = values["purity"]
            ploidy_mapping[short_name + "-HRD"] = values["ploidy"]

        print(f"✓ Loaded purity for {len(purity_mapping)} entries (including short name variants)")
    except Exception as e:
        print(f"⚠️  ERROR loading purity file: {e}")
        print("   Proceeding WITHOUT purity correction")
        purity_mapping = {}
        ploidy_mapping = {}

    print(f"{'='*70}\n")
    return purity_mapping, ploidy_mapping


def standardize_segments(df):
    """
    Normalize segment table to a common schema.
    Guaranteed output columns: chromosome, start, end, total_cn, nMajor, nMinor, segment_format.
    """
    out = df.copy()

    rename_map = {
        "chr": "chromosome", "Chr": "chromosome",
        "startpos": "start", "endpos": "end",
        "segVal": "total_cn", "CN": "total_cn",
    }
    out.rename(columns=rename_map, inplace=True)

    for col in ["chromosome", "start", "end"]:
        if col not in out.columns:
            raise ValueError(f"Missing required segment column: {col}")

    out["chromosome"] = (
        out["chromosome"].astype(str)
        .str.replace("^chr", "", regex=True)
        .str.replace(r"\.0$", "", regex=True)
        .str.strip()
        .replace({"X": "23", "Y": "24"})
    )
    out["chromosome"] = pd.to_numeric(out["chromosome"], errors="coerce")
    out = out.dropna(subset=["chromosome"]).copy()
    out["chromosome"] = out["chromosome"].astype(int)

    out["start"] = pd.to_numeric(out["start"], errors="coerce")
    out["end"]   = pd.to_numeric(out["end"],   errors="coerce")
    out = out.dropna(subset=["start", "end"]).copy()
    out["start"] = out["start"].astype(int)
    out["end"]   = out["end"].astype(int)

    has_nmajor = "nMajor" in out.columns
    has_nminor = "nMinor" in out.columns
    has_total  = "total_cn" in out.columns

    if has_nmajor and has_nminor:
        out["nMajor"] = pd.to_numeric(out["nMajor"], errors="coerce")
        out["nMinor"] = pd.to_numeric(out["nMinor"], errors="coerce")
        computed_total = out["nMajor"] + out["nMinor"]
        if has_total:
            existing = pd.to_numeric(out["total_cn"], errors="coerce").round()
            mismatch = (existing - computed_total.round()).abs() > 0.5
            n_mismatch = mismatch.sum()
            if n_mismatch > 0:
                print(f"  ⚠️  {n_mismatch} rows: existing total_cn differs from nMajor+nMinor "
                      f"by >0.5 — replacing with nMajor+nMinor")
        out["total_cn"] = computed_total
        out["segment_format"] = "ascat"
    elif has_total:
        out["total_cn"]        = pd.to_numeric(out["total_cn"], errors="coerce")
        out["nMajor"]          = np.nan
        out["nMinor"]          = np.nan
        out["segment_format"]  = "segval_only"
    else:
        raise ValueError("Unknown segment format: missing both total_cn and nMajor/nMinor")

    out["total_cn"] = pd.to_numeric(out["total_cn"], errors="coerce")
    return out


def _find_segment_file(sample_short, search_dirs):
    """
    Return the first existing segment file for sample_short across search_dirs.

    Tries the standard TSO filename pattern first; falls back to a glob for any
    CSV file whose name starts with sample_short in each directory.

    Samples sequenced with the HRD panel carry a '-HRD' suffix in their name
    but share the same reference segment file as the base sample (e.g.
    'JBLAB17012-HRD' → segments of 'JBLAB17012'). If no file is found under
    the full name, the lookup is retried without the '-HRD' suffix.
    """
    lookup_names = [sample_short]
    if sample_short.endswith("-HRD"):
        lookup_names.append(sample_short[: -len("-HRD")])

    for name in lookup_names:
        for d in search_dirs:
            if d is None or not os.path.isdir(d):
                continue
            standard = os.path.join(
                d, f"{name}_100kb_abs_segtable_free_purity_filtered_unique.csv"
            )
            if os.path.exists(standard):
                if name != sample_short:
                    print(f"  ℹ️  {sample_short}: segments resolved via base name '{name}'")
                return standard
            matches = sorted(glob.glob(os.path.join(d, f"{name}_*.csv")))
            if len(matches) > 1:
                print(f"  ⚠️  Multiple segment files matched '{name}*.csv' in {d}; "
                      f"using first: {os.path.basename(matches[0])}")
            if matches:
                if name != sample_short:
                    print(f"  ℹ️  {sample_short}: segments resolved via base name '{name}'")
                return matches[0]
    return None


def load_segments_for_samples(samples_dict, raw_segments_dir):
    print("Loading segment files...")
    search_dirs = [raw_segments_dir, ASCAT_SEGMENTS_DIR]
    df_list = []

    for sample_short, long_name in samples_dict.items():
        path = _find_segment_file(sample_short, search_dirs)
        if path is None:
            print(f"  ⚠️  Segments not found for {sample_short} in any search dir")
            continue
        df = pd.read_csv(path)
        df = standardize_segments(df)
        df["SampleID"] = sample_short
        df["sample"]   = long_name
        df_list.append(df)

    if len(df_list) == 0:
        return pd.DataFrame()

    combined = pd.concat(df_list, ignore_index=True)
    print(f"✓ Loaded {len(combined):,} segment rows across {len(df_list)} sample files")
    return combined


def apply_purity_correction(baf_df, purity, baf_col="BAF_n"):
    df = baf_df.copy()

    if pd.isna(purity) or purity <= 0 or purity > 1:
        print(f"    ⚠️  WARNING: Invalid purity {purity}, using raw BAF (no correction)")
        df["BAF_raw"]       = df[baf_col]
        df["BAF_corrected"] = df[baf_col]
        df["purity_used"]   = np.nan
        return df

    print(f"    Applying purity correction (purity={purity:.3f})...")
    baf_raw = df[baf_col].values
    baf_corrected = np.where(
        baf_raw >= 0.5,
        0.5 + (baf_raw - 0.5) / purity,
        0.5 - (0.5 - baf_raw) / purity,
    )
    baf_corrected = np.clip(baf_corrected, 0, 1)

    df["BAF_raw"]       = baf_raw
    df["BAF_corrected"] = baf_corrected
    df["purity_used"]   = purity

    n_clamped_upper = ((baf_corrected == 1.0) & (baf_raw < 1.0)).sum()
    n_clamped_lower = ((baf_corrected == 0.0) & (baf_raw > 0.0)).sum()
    mean_shift = np.abs(baf_corrected - baf_raw).mean()

    print(f"      BAF range before: [{baf_raw.min():.3f}, {baf_raw.max():.3f}]")
    print(f"      BAF range after:  [{baf_corrected.min():.3f}, {baf_corrected.max():.3f}]")
    print(f"      Mean absolute shift: {mean_shift:.4f}")
    print(f"      Values clamped to 1.0: {n_clamped_upper:,}")
    print(f"      Values clamped to 0.0: {n_clamped_lower:,}")
    return df


def calculate_segment_baf_summary_metrics(baf_hits_df):
    if len(baf_hits_df) == 0:
        return pd.DataFrame(columns=[
            "segment_id", "n_positions",
            "baf_mean", "baf_median", "baf_sd", "baf_min", "baf_max", "baf_q25", "baf_q75", "baf_iqr",
            "baf_mean_raw", "baf_median_raw", "baf_sd_raw", "baf_min_raw", "baf_max_raw",
            "baf_q25_raw", "baf_q75_raw", "baf_iqr_raw",
        ])

    summary = baf_hits_df.groupby("segment_id").agg({
        "BAF_corrected": [
            ("n_positions", "count"),
            ("baf_mean",   "mean"),
            ("baf_median", "median"),
            ("baf_sd",     "std"),
            ("baf_min",    "min"),
            ("baf_max",    "max"),
            ("baf_q25",    lambda x: x.quantile(0.25)),
            ("baf_q75",    lambda x: x.quantile(0.75)),
        ],
        "BAF_raw": [
            ("baf_mean_raw",   "mean"),
            ("baf_median_raw", "median"),
            ("baf_sd_raw",     "std"),
            ("baf_min_raw",    "min"),
            ("baf_max_raw",    "max"),
            ("baf_q25_raw",    lambda x: x.quantile(0.25)),
            ("baf_q75_raw",    lambda x: x.quantile(0.75)),
        ],
    }).reset_index()

    summary.columns = [
        "segment_id",
        "n_positions",
        "baf_mean", "baf_median", "baf_sd", "baf_min", "baf_max", "baf_q25", "baf_q75",
        "baf_mean_raw", "baf_median_raw", "baf_sd_raw", "baf_min_raw", "baf_max_raw",
        "baf_q25_raw", "baf_q75_raw",
    ]
    summary["baf_iqr"]     = summary["baf_q75"]     - summary["baf_q25"]
    summary["baf_iqr_raw"] = summary["baf_q75_raw"] - summary["baf_q25_raw"]
    return summary


def map_baf_to_segments_FINE(
    segments_df, baf_df,
    seg_chr="chromosome", seg_start="start", seg_end="end",
    baf_chr="CHROM", baf_pos="POS",
    baf_val_corrected="BAF_corrected", baf_val_raw="BAF_raw",
    n_bins=20,
):
    segments = segments_df.copy()
    baf      = baf_df.copy()

    for _df, _col in [(segments, seg_chr), (baf, baf_chr)]:
        _df[_col] = (
            _df[_col].astype(str)
            .str.replace(r"^chr", "", regex=True)
            .str.replace(r"\.0$", "", regex=True)
            .str.strip()
        )

    segments[seg_start] = pd.to_numeric(segments[seg_start], errors="coerce")
    segments[seg_end]   = pd.to_numeric(segments[seg_end],   errors="coerce")
    baf[baf_pos]        = pd.to_numeric(baf[baf_pos],        errors="coerce")

    bin_edges  = np.linspace(0, 1, n_bins + 1)
    bin_labels = [f"{bin_edges[i]:.2f}-{bin_edges[i+1]:.2f}" for i in range(n_bins)]

    print(f"  Creating {n_bins} BAF bins (width {bin_edges[1] - bin_edges[0]:.3f}): "
          f"{bin_labels[0]} → {bin_labels[-1]}")

    for label in bin_labels:
        segments[f"count_{label}"]     = 0
        segments[f"ratio_{label}"]     = np.nan
        segments[f"count_raw_{label}"] = 0
        segments[f"ratio_raw_{label}"] = np.nan

    segments["n_baf_points"] = 0
    baf_hits_list = []

    for idx, seg in segments.iterrows():
        seg_baf = baf[
            (baf[baf_chr] == str(seg[seg_chr])) &
            (baf[baf_pos] >= seg[seg_start]) &
            (baf[baf_pos] <= seg[seg_end])
        ].copy()

        n_points = len(seg_baf)
        segments.at[idx, "n_baf_points"] = n_points

        if n_points > 0:
            seg_baf["segment_id"] = seg["segment_id"]
            baf_hits_list.append(seg_baf)

            baf_vals_corrected = seg_baf[baf_val_corrected].values
            baf_vals_raw       = seg_baf[baf_val_raw].values

            for i in range(n_bins):
                lower = bin_edges[i]
                upper = bin_edges[i + 1]
                label = bin_labels[i]
                if i == n_bins - 1:
                    count_corr = ((baf_vals_corrected >= lower) & (baf_vals_corrected <= upper)).sum()
                    count_raw  = ((baf_vals_raw >= lower)       & (baf_vals_raw <= upper)).sum()
                else:
                    count_corr = ((baf_vals_corrected >= lower) & (baf_vals_corrected < upper)).sum()
                    count_raw  = ((baf_vals_raw >= lower)       & (baf_vals_raw < upper)).sum()
                segments.at[idx, f"count_{label}"]     = count_corr
                segments.at[idx, f"count_raw_{label}"] = count_raw

    baf_hits = pd.concat(baf_hits_list, ignore_index=True) if baf_hits_list else pd.DataFrame()

    for label in bin_labels:
        n = segments["n_baf_points"].replace(0, np.nan)
        segments[f"ratio_{label}"]     = segments[f"count_{label}"]     / n
        segments[f"ratio_raw_{label}"] = segments[f"count_raw_{label}"] / n

    outer_bins = [
        label for label in bin_labels
        if float(label.split("-")[1]) <= 0.25 or float(label.split("-")[0]) >= 0.75
    ]
    inner_bins = [
        label for label in bin_labels
        if 0.25 <= float(label.split("-")[0]) < 0.75 and float(label.split("-")[1]) <= 0.75
    ]

    segments["outer_mass"]     = segments[[f"ratio_{label}"     for label in outer_bins]].sum(axis=1).fillna(0)
    segments["inner_mass"]     = segments[[f"ratio_{label}"     for label in inner_bins]].sum(axis=1).fillna(0)
    segments["outer_mass_raw"] = segments[[f"ratio_raw_{label}" for label in outer_bins]].sum(axis=1).fillna(0)
    segments["inner_mass_raw"] = segments[[f"ratio_raw_{label}" for label in inner_bins]].sum(axis=1).fillna(0)

    print(f"  ✓ Mapped {len(baf_hits):,} BAF positions to {len(segments):,} segments")
    print(f"  ✓ Segments with BAF data: {(segments['n_baf_points'] > 0).sum():,}")
    return segments, baf_hits


def normalize_baf_df(df):
    """
    Normalize a raw BAF dataframe to a consistent schema regardless of source.

    TSO format (_qc.tsv):
      - CHROM: 'chr1', 'chr2', …
      - Count columns: Ref_Reads / Alt_Reads
      - Filter column: 'problematic_reason' (NaN → keep)

    ASCAT corrected format (.tumour.tso500_sparse_downsample_corrected.txt.gz):
      - CHROM: '1', '2', … (no chr prefix)
      - Count columns: REFcounts / ALTcounts  → renamed to Ref_Reads / Alt_Reads
      - Filter column: 'is_complex' (True/False string → drop Trues)

    Older ASCAT format (.tumour.tso500_sparse_downsample.txt):
      - CHROM: '1', '2', … (no chr prefix)
      - Count columns: REFcounts / ALTcounts  → renamed to Ref_Reads / Alt_Reads
      - Filter column: 'problematic' (True/False string → drop Trues)

    Returns the dataframe with:
      - Unified count columns: Ref_Reads / Alt_Reads
      - Added boolean column '_ok' (True = keep row)
    """
    df = df.copy()

    # Unify count column names so BAF_n computation always uses Ref_Reads / Alt_Reads.
    if "REFcounts" in df.columns and "Ref_Reads" not in df.columns:
        df = df.rename(columns={"REFcounts": "Ref_Reads", "ALTcounts": "Alt_Reads"})

    # Unify the problematic-position filter flag → boolean _ok column.
    if "problematic_reason" in df.columns:          # TSO _qc.tsv
        ok = df["problematic_reason"].isna()
    elif "is_complex" in df.columns:                # ASCAT corrected .txt.gz
        ok = ~df["is_complex"].astype(str).str.lower().str.strip().eq("true")
    elif "problematic" in df.columns:               # older ASCAT .txt
        ok = ~df["problematic"].astype(str).str.lower().str.strip().eq("true")
    else:
        ok = pd.Series(True, index=df.index)

    df["_ok"] = ok
    return df


def detect_baf_suffix(folder_path, long_name):
    """Return the first BAF_SUFFIXES entry for which a file exists in folder_path."""
    for suffix in BAF_SUFFIXES:
        candidate = os.path.join(folder_path, f"{long_name}{suffix}")
        if os.path.exists(candidate):
            return suffix
    return BAF_SUFFIXES[0]  # fallback — caller will handle missing file


def compute_feature_importance(df, label_col="loh_signal"):
    features = ["outer_mass", "inner_mass", "total_cn", "n_baf_points", "segment_len_bp"]
    d = df.dropna(subset=[label_col]).copy()
    uni_results = []

    for feat in features:
        dd = d[[label_col, feat]].dropna()
        if len(dd) > 0 and len(dd[label_col].unique()) >= 2:
            auc_val = safe_auc(dd[label_col].values, dd[feat].values)
        else:
            auc_val = np.nan
        uni_results.append((feat, auc_val))

    uni_df = pd.DataFrame(uni_results, columns=["feature", "univariate_auc"])

    X = d[features].copy()
    y = d[label_col].astype(int).values
    mask = X.notna().all(axis=1)
    X = X.loc[mask].values
    y = y[mask.values]

    if len(np.unique(y)) >= 2 and len(y) >= 10:
        scaler = StandardScaler()
        Xs = scaler.fit_transform(X)
        clf = LogisticRegression(penalty="l2", solver="lbfgs", max_iter=2000)
        clf.fit(Xs, y)
        coef = clf.coef_[0]
        proba = clf.predict_proba(Xs)[:, 1]
        model_auc = safe_auc(y, proba)
    else:
        coef      = np.full(len(features), np.nan)
        model_auc = np.nan

    lr_df = pd.DataFrame({"feature": features, "logreg_coef": coef})
    return uni_df.merge(lr_df, on="feature", how="left"), model_auc


def plot_outer_mass_correction_comparison(segments_df, sample_name, dp_filter, output_dir):
    if "outer_mass_raw" not in segments_df.columns or "outer_mass" not in segments_df.columns:
        return

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    axes[0].hist(segments_df["outer_mass_raw"].dropna(), bins=50, alpha=0.7, color="blue", edgecolor="black")
    axes[0].set_xlabel("Raw outer_mass", fontsize=12, fontweight="bold")
    axes[0].set_ylabel("Count (segments)", fontsize=12, fontweight="bold")
    axes[0].set_title("Before Purity Correction", fontsize=13, fontweight="bold")
    axes[0].grid(alpha=0.3); axes[0].set_xlim(0, 1)

    axes[1].hist(segments_df["outer_mass"].dropna(), bins=50, alpha=0.7, color="green", edgecolor="black")
    axes[1].set_xlabel("Corrected outer_mass", fontsize=12, fontweight="bold")
    axes[1].set_ylabel("Count (segments)", fontsize=12, fontweight="bold")
    axes[1].set_title("After Purity Correction", fontsize=13, fontweight="bold")
    axes[1].grid(alpha=0.3); axes[1].set_xlim(0, 1)

    axes[2].scatter(segments_df["outer_mass_raw"], segments_df["outer_mass"],
                    alpha=0.5, s=30, color="purple", edgecolor="black", linewidth=0.5)
    axes[2].plot([0, 1], [0, 1], "k--", linewidth=2, label="No change", alpha=0.5)
    mean_shift = np.abs(segments_df["outer_mass"] - segments_df["outer_mass_raw"]).mean()
    axes[2].text(0.05, 0.95, f"Mean shift: {mean_shift:.4f}", transform=axes[2].transAxes,
                 fontsize=10, va="top", bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.8))
    axes[2].set_xlabel("Raw outer_mass", fontsize=12, fontweight="bold")
    axes[2].set_ylabel("Corrected outer_mass", fontsize=12, fontweight="bold")
    axes[2].set_title("Raw vs Corrected", fontsize=13, fontweight="bold")
    axes[2].legend(); axes[2].grid(alpha=0.3)
    axes[2].set_xlim(0, 1); axes[2].set_ylim(0, 1)

    plt.suptitle(f"{sample_name} | outer_mass Purity Correction | DP≥{dp_filter}",
                 fontsize=14, fontweight="bold")
    plt.tight_layout()

    os.makedirs(output_dir, exist_ok=True)
    out_file = os.path.join(output_dir, f"outer_mass_correction_dp{dp_filter}.png")
    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"      ✓ Saved outer_mass comparison: {out_file}")


def plot_baf_correction_comparison(baf_df, sample_name, dp_filter, purity, output_dir):
    if "BAF_raw" not in baf_df.columns or "BAF_corrected" not in baf_df.columns:
        return

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    plot_df = baf_df.sample(10000, random_state=42) if len(baf_df) > 10000 else baf_df

    axes[0].hist(plot_df["BAF_raw"], bins=50, alpha=0.7, color="blue", edgecolor="black")
    axes[0].axvline(0.5, color="red", linestyle="--", linewidth=2, label="Het line")
    axes[0].set_xlabel("Raw BAF", fontsize=12, fontweight="bold")
    axes[0].set_ylabel("Count", fontsize=12, fontweight="bold")
    axes[0].set_title("Before Purity Correction", fontsize=13, fontweight="bold")
    axes[0].legend(); axes[0].grid(alpha=0.3); axes[0].set_xlim(0, 1)

    axes[1].hist(plot_df["BAF_corrected"], bins=50, alpha=0.7, color="green", edgecolor="black")
    axes[1].axvline(0.5, color="red", linestyle="--", linewidth=2, label="Het line")
    axes[1].set_xlabel("Corrected BAF", fontsize=12, fontweight="bold")
    axes[1].set_ylabel("Count", fontsize=12, fontweight="bold")
    axes[1].set_title(f"After Correction (purity={purity:.3f})", fontsize=13, fontweight="bold")
    axes[1].legend(); axes[1].grid(alpha=0.3); axes[1].set_xlim(0, 1)

    axes[2].scatter(plot_df["BAF_raw"], plot_df["BAF_corrected"],
                    alpha=0.3, s=5, color="purple", rasterized=True)
    axes[2].plot([0, 1], [0, 1], "k--", linewidth=2, label="No change", alpha=0.5)
    axes[2].axhline(0.5, color="red", linestyle=":", alpha=0.3)
    axes[2].axvline(0.5, color="red", linestyle=":", alpha=0.3)
    axes[2].set_xlabel("Raw BAF", fontsize=12, fontweight="bold")
    axes[2].set_ylabel("Corrected BAF", fontsize=12, fontweight="bold")
    axes[2].set_title("Raw vs Corrected", fontsize=13, fontweight="bold")
    axes[2].legend(); axes[2].grid(alpha=0.3)
    axes[2].set_xlim(0, 1); axes[2].set_ylim(0, 1)

    plt.suptitle(f"{sample_name} | BAF Purity Correction | DP≥{dp_filter}",
                 fontsize=14, fontweight="bold")
    plt.tight_layout()

    os.makedirs(output_dir, exist_ok=True)
    out_file = os.path.join(output_dir, f"baf_correction_dp{dp_filter}.png")
    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"      ✓ Saved correction plot: {out_file}")


def create_input_data_scarHRD(segments_df, sample):
    df = segments_df.copy()

    if len(df) == 0:
        print(f"⚠️  WARNING: No segments for sample {sample}")
        return pd.DataFrame(columns=[
            "SampleID", "Chromosome", "Start_position", "End_position",
            "total_cn", "A_cn", "B_cn", "ploidy",
        ])

    df["totalCN"] = pd.to_numeric(df["totalCN"], errors="coerce")
    df["Outer>mid_mass"] = df["outer_mass"] > df["mid_mass"]
    df["Outer<mid_mass"] = df["outer_mass"] < df["mid_mass"]
    df["Outer=mid_mass"] = df["outer_mass"] == df["mid_mass"]

    cols = ["outer_mass", "mid_mass"]
    df["max_metric"] = df[cols].replace(0, np.nan).idxmax(axis=1)

    is_even    = (df["totalCN"] % 2 == 0) & df["totalCN"].notna()
    is_odd     = (df["totalCN"] % 2 == 1) & df["totalCN"].notna()
    is_loh1    = (df["totalCN"] == 1)
    is_cn0     = (df["totalCN"] == 0)
    is_odd_gt1 = is_odd & (df["totalCN"] > 1)

    df["A_cn"] = np.nan
    df["B_cn"] = np.nan

    mask = is_even & df["max_metric"].eq("outer_mass") & (df["outer_mass"] != 0)
    df.loc[mask, "A_cn"] = df.loc[mask, "totalCN"]
    df.loc[mask, "B_cn"] = 0

    mask = is_even & df["max_metric"].eq("mid_mass") & (df["mid_mass"] != 0)
    df.loc[mask, "A_cn"] = df.loc[mask, "totalCN"] / 2
    df.loc[mask, "B_cn"] = df.loc[mask, "totalCN"] / 2

    df.loc[is_loh1, "A_cn"] = 1
    df.loc[is_loh1, "B_cn"] = 0
    df.loc[is_cn0,  "A_cn"] = 0
    df.loc[is_cn0,  "B_cn"] = 0

    mask = is_odd_gt1 & df["max_metric"].eq("outer_mass") & (df["outer_mass"] != 0)
    df.loc[mask, "A_cn"] = df.loc[mask, "totalCN"]
    df.loc[mask, "B_cn"] = 0

    mask = is_odd_gt1 & df["max_metric"].eq("mid_mass") & (df["mid_mass"] != 0)
    df.loc[mask, "A_cn"] = np.floor(df.loc[mask, "totalCN"] / 2)
    df.loc[mask, "B_cn"] = np.ceil(df.loc[mask, "totalCN"] / 2)

    mask = is_even & df["max_metric"].isna()
    df.loc[mask, "A_cn"] = df.loc[mask, "totalCN"] / 2
    df.loc[mask, "B_cn"] = df.loc[mask, "totalCN"] / 2

    mask = is_odd & df["max_metric"].isna()
    df.loc[mask, "A_cn"] = df.loc[mask, "totalCN"] - 1
    df.loc[mask, "B_cn"] = 1

    df["A_cn"]   = df["A_cn"].round().astype("Int64")
    df["B_cn"]   = df["B_cn"].round().astype("Int64")
    df["ploidy"] = df["totalCN"].mean()

    return pd.DataFrame({
        "SampleID":       sample,
        "Chromosome":     df["Chromosome"],
        "Start_position": df["start"],
        "End_position":   df["end"],
        "total_cn":       df["totalCN"],
        "A_cn":           df["A_cn"],
        "B_cn":           df["B_cn"],
        "ploidy":         df["ploidy"],
    })


def classify_cn_match(row):
    if pd.isna(row["nMajor"]) or pd.isna(row["nMinor"]):
        return "no_prediction"
    if pd.isna(row["nMajor_gold"]) or pd.isna(row["nMinor_gold"]):
        return "no_gold_standard"
    total_pred = row["nMajor"] + row["nMinor"]
    total_gold = row["nMajor_gold"] + row["nMinor_gold"]
    if row["nMajor"] == row["nMajor_gold"] and row["nMinor"] == row["nMinor_gold"]:
        return "allelic_match"
    if row["nMajor"] == row["nMinor_gold"] and row["nMinor"] == row["nMajor_gold"]:
        return "allelic_match"
    if total_pred == total_gold:
        return "totalCN_match"
    return "no_match"


def create_input_data_scarHRD_with_gold(segments_df):
    df = segments_df.copy()
    if "sample_name" in df.columns:
        sample = df["sample_name"].iloc[0] if len(df) > 0 else "SAMPLE"
    elif "sample" in df.columns:
        sample = df["sample"].iloc[0] if len(df) > 0 else "SAMPLE"
    else:
        sample = "SAMPLE"
    scarhrd_input = create_input_data_scarHRD(df, sample)
    return scarhrd_input, df


# ============================================================
# SCRIPT 2 HELPER FUNCTIONS
# ============================================================

def assert_required_columns(df, required_cols, df_name):
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"{df_name} is missing required columns: {missing}")


def load_manifest(manifest_file):
    if not os.path.exists(manifest_file):
        raise FileNotFoundError(
            f"Manifest file not found: {manifest_file}\nRun script 1 first."
        )
    manifest = pd.read_csv(manifest_file)
    if len(manifest) == 0:
        raise ValueError("Manifest file is empty")
    return manifest


def deduplicate_manifest(manifest_df):
    subset_cols = ["sample_name", "sample", "file_type", "baf_folder", "segments_file", "baf_hits_file"]
    existing = [c for c in subset_cols if c in manifest_df.columns]
    return manifest_df.drop_duplicates(subset=existing).reset_index(drop=True)


# ============================================================
# SCRIPT 3 HELPER FUNCTIONS
# ============================================================

def classify_final_state_from_acn_bcn(a_cn, b_cn):
    if pd.isna(a_cn) or pd.isna(b_cn):
        return "Unknown"
    if a_cn == 0 and b_cn == 0:
        return "Deletion"
    if b_cn == 0:
        return "LoH"
    if a_cn == b_cn:
        return "Balanced"
    return "AI"


def worst_tier(series):
    tier_rank = {
        "tier1": 1, "tier2": 2, "tier3": 3,
        "tier3_heuristic": 4, "tier3_ai_default": 5,
        "no_baf": 6, "unknown": 7,
    }
    inv_rank = {v: k for k, v in tier_rank.items()}
    vals = [x for x in series.dropna().astype(str).tolist() if x in tier_rank]
    if not vals:
        return "unknown"
    return inv_rank[max(tier_rank[x] for x in vals)]


def build_scarHRD_input_using_thresholds(segments_df, sample, thresholds):
    """
    Infer A_cn / B_cn from total_cn + outer_mass_raw using pre-computed thresholds.
    Returns (output_df, annotated_df) — both pre-merge.
    """
    df = segments_df.copy()

    if len(df) == 0:
        empty_cols = ["SampleID", "Chromosome", "Start_position", "End_position",
                      "total_cn", "A_cn", "B_cn", "ploidy"]
        return pd.DataFrame(columns=empty_cols), df

    if "total_cn_raw" not in df.columns:
        df["total_cn_raw"] = pd.to_numeric(df["total_cn"], errors="coerce")

    df["totalCN_raw"] = pd.to_numeric(df["total_cn_raw"], errors="coerce")
    df["totalCN"]     = np.floor(df["totalCN_raw"] + 0.5).astype("Int64")
    df["total_cn_used_for_decision"] = df["totalCN"]
    df["total_cn_rounding_rule"]     = "classical_rounding_floor(x+0.5)"

    # When purity correction is active (PURITY_FILE set), use the purity-corrected
    # outer_mass so that classification and calibration are on the same BAF scale
    # (thresholds were derived from outer_mass_purity_corrected, which equals the
    # pipeline's outer_mass column when apply_purity_correction was applied).
    # Without purity file, use outer_mass_raw to match the raw-calibrated thresholds.
    _purity_active = PURITY_FILE is not None
    if _purity_active and "outer_mass" in df.columns:
        score_col = "outer_mass"
        print("ℹ️  Purity correction active — using corrected 'outer_mass' for classification")
    elif "outer_mass_raw" in df.columns:
        score_col = "outer_mass_raw"
    elif "outer_mass_raw_x" in df.columns:
        score_col = "outer_mass_raw_x"
        print("ℹ️  Using 'outer_mass_raw_x' (merged CSV with _x/_y suffixes)")
    elif "outer_mass_raw_y" in df.columns:
        score_col = "outer_mass_raw_y"
        print("ℹ️  Using 'outer_mass_raw_y' (merged CSV with _x/_y suffixes)")
    elif "outer_mass" in df.columns:
        score_col = "outer_mass"
        print("⚠️  'outer_mass_raw' not found — falling back to corrected 'outer_mass'")
    elif "outer_mass_x" in df.columns:
        score_col = "outer_mass_x"
        print("⚠️  Using 'outer_mass_x' (merged CSV with _x/_y suffixes, corrected)")
    else:
        raise ValueError("Neither 'outer_mass_raw' nor 'outer_mass' found in input dataframe")

    df["score_column_used"] = score_col
    df["outer_mass_used"]   = df[score_col]

    df["A_cn"]              = np.nan
    df["B_cn"]              = np.nan
    df["confidence"]        = np.nan
    df["tier"]              = "unknown"
    df["threshold_applied"] = np.nan
    df["decision_rule"]     = "unknown"

    print(f"\n{'='*70}")
    print("THRESHOLD-BASED ALLELE-SPECIFIC CN CLASSIFICATION")
    print(f"Sample: {sample}")
    print(f"{'='*70}")

    for total_cn in sorted(df["totalCN"].dropna().unique()):
        total_cn   = int(total_cn)
        cn_mask    = df["totalCN"] == total_cn
        cn_segments = df[cn_mask]
        n_segments = cn_mask.sum()
        print(f"\nProcessing CN={total_cn}: {n_segments} segments")

        if total_cn == 0:
            df.loc[cn_mask, ["A_cn", "B_cn", "confidence", "tier"]] = [0, 0, 0.99, "tier1"]
            df.loc[cn_mask, "threshold_applied"] = np.nan
            df.loc[cn_mask, "decision_rule"] = "CN=0 -> homozygous deletion (0:0)"
            print("  → 0:0 (homozygous deletion)")

        elif total_cn == 1:
            df.loc[cn_mask, ["A_cn", "B_cn", "confidence", "tier"]] = [1, 0, 0.95, "tier1"]
            df.loc[cn_mask, "threshold_applied"] = np.nan
            df.loc[cn_mask, "decision_rule"] = "CN=1 -> mandatory LoH (1:0)"
            print("  → 1:0 (LoH)")

        elif total_cn == 2:
            threshold = thresholds.get("2", {}).get("loh_vs_balanced", 0.93)
            for idx in cn_segments.index:
                outer_mass = df.loc[idx, score_col]
                df.loc[idx, "threshold_applied"] = threshold
                if pd.isna(outer_mass):
                    df.loc[idx, ["A_cn", "B_cn", "confidence", "tier"]] = [1, 1, 0.50, "no_baf"]
                    df.loc[idx, "decision_rule"] = f"CN=2; no BAF -> default Balanced (1:1); threshold={threshold:.3f}"
                elif outer_mass > threshold:
                    df.loc[idx, ["A_cn", "B_cn", "confidence", "tier"]] = [2, 0, 0.90, "tier1"]
                    df.loc[idx, "decision_rule"] = f"CN=2; outer_mass>{threshold:.3f} -> LoH (2:0)"
                else:
                    df.loc[idx, ["A_cn", "B_cn", "confidence", "tier"]] = [1, 1, 0.88, "tier1"]
                    df.loc[idx, "decision_rule"] = f"CN=2; outer_mass<={threshold:.3f} -> Balanced (1:1)"
            loh = ((df.loc[cn_mask, "A_cn"] == 2) & (df.loc[cn_mask, "B_cn"] == 0)).sum()
            bal = ((df.loc[cn_mask, "A_cn"] == 1) & (df.loc[cn_mask, "B_cn"] == 1)).sum()
            print(f"  → LoH (2:0): {loh}, Balanced (1:1): {bal}")

        elif total_cn == 3:
            threshold = thresholds.get("3", {}).get("loh_vs_ai", 0.92)
            for idx in cn_segments.index:
                outer_mass = df.loc[idx, score_col]
                df.loc[idx, "threshold_applied"] = threshold
                if pd.isna(outer_mass):
                    df.loc[idx, ["A_cn", "B_cn", "confidence", "tier"]] = [2, 1, 0.50, "no_baf"]
                    df.loc[idx, "decision_rule"] = f"CN=3; no BAF -> default AI (2:1); threshold={threshold:.3f}"
                elif outer_mass > threshold:
                    df.loc[idx, ["A_cn", "B_cn", "confidence", "tier"]] = [3, 0, 0.75, "tier2"]
                    df.loc[idx, "decision_rule"] = f"CN=3; outer_mass>{threshold:.3f} -> LoH (3:0)"
                else:
                    df.loc[idx, ["A_cn", "B_cn", "confidence", "tier"]] = [2, 1, 0.70, "tier2"]
                    df.loc[idx, "decision_rule"] = f"CN=3; outer_mass<={threshold:.3f} -> AI (2:1)"
            loh = ((df.loc[cn_mask, "A_cn"] == 3) & (df.loc[cn_mask, "B_cn"] == 0)).sum()
            ai  = ((df.loc[cn_mask, "A_cn"] == 2) & (df.loc[cn_mask, "B_cn"] == 1)).sum()
            print(f"  → LoH (3:0): {loh}, AI (2:1): {ai}")

        elif total_cn == 4:
            loh_threshold    = thresholds.get("4", {}).get("loh_vs_balanced", 0.965)
            ai_bal_threshold = thresholds.get("4", {}).get("ai_vs_balanced",  0.85)
            for idx in cn_segments.index:
                outer_mass = df.loc[idx, score_col]
                if pd.isna(outer_mass):
                    df.loc[idx, ["A_cn", "B_cn", "confidence", "tier"]] = [2, 2, 0.50, "no_baf"]
                    df.loc[idx, "threshold_applied"] = np.nan
                    df.loc[idx, "decision_rule"] = (
                        f"CN=4; no BAF -> default Balanced (2:2); "
                        f"loh_threshold={loh_threshold:.3f}, ai_bal_threshold={ai_bal_threshold:.3f}"
                    )
                elif outer_mass > loh_threshold:
                    df.loc[idx, ["A_cn", "B_cn", "confidence", "tier"]] = [4, 0, 0.75, "tier2"]
                    df.loc[idx, "threshold_applied"] = loh_threshold
                    df.loc[idx, "decision_rule"] = f"CN=4; outer_mass>{loh_threshold:.3f} -> LoH (4:0)"
                elif outer_mass > ai_bal_threshold:
                    df.loc[idx, ["A_cn", "B_cn", "confidence", "tier"]] = [3, 1, 0.65, "tier2"]
                    df.loc[idx, "threshold_applied"] = ai_bal_threshold
                    df.loc[idx, "decision_rule"] = (
                        f"CN=4; outer_mass<={loh_threshold:.3f} and >{ai_bal_threshold:.3f} -> AI (3:1)"
                    )
                else:
                    df.loc[idx, ["A_cn", "B_cn", "confidence", "tier"]] = [2, 2, 0.72, "tier2"]
                    df.loc[idx, "threshold_applied"] = ai_bal_threshold
                    df.loc[idx, "decision_rule"] = f"CN=4; outer_mass<={ai_bal_threshold:.3f} -> Balanced (2:2)"
            loh = ((df.loc[cn_mask, "A_cn"] == 4) & (df.loc[cn_mask, "B_cn"] == 0)).sum()
            ai  = ((df.loc[cn_mask, "A_cn"] == 3) & (df.loc[cn_mask, "B_cn"] == 1)).sum()
            bal = ((df.loc[cn_mask, "A_cn"] == 2) & (df.loc[cn_mask, "B_cn"] == 2)).sum()
            print(f"  → LoH (4:0): {loh}, AI (3:1): {ai}, Balanced (2:2): {bal}")

        elif total_cn == 5:
            threshold = thresholds.get("5", {}).get("loh_vs_ai", 0.96)
            for idx in cn_segments.index:
                outer_mass = df.loc[idx, score_col]
                df.loc[idx, "threshold_applied"] = threshold
                if pd.isna(outer_mass):
                    df.loc[idx, ["A_cn", "B_cn", "confidence", "tier"]] = [3, 2, 0.50, "no_baf"]
                    df.loc[idx, "decision_rule"] = f"CN=5; no BAF -> default AI (3:2); threshold={threshold:.3f}"
                elif outer_mass > threshold:
                    df.loc[idx, ["A_cn", "B_cn", "confidence", "tier"]] = [5, 0, 0.70, "tier3"]
                    df.loc[idx, "decision_rule"] = f"CN=5; outer_mass>{threshold:.3f} -> LoH (5:0)"
                else:
                    df.loc[idx, ["A_cn", "B_cn", "confidence", "tier"]] = [3, 2, 0.60, "tier3"]
                    df.loc[idx, "decision_rule"] = f"CN=5; outer_mass<={threshold:.3f} -> AI (3:2)"
            loh = ((df.loc[cn_mask, "A_cn"] == 5) & (df.loc[cn_mask, "B_cn"] == 0)).sum()
            ai  = ((df.loc[cn_mask, "A_cn"] == 3) & (df.loc[cn_mask, "B_cn"] == 2)).sum()
            print(f"  → LoH (5:0): {loh}, AI (3:2): {ai}")

        elif total_cn == 6:
            threshold = thresholds.get("6", {}).get("ai_vs_balanced", 0.80)
            for idx in cn_segments.index:
                outer_mass = df.loc[idx, score_col]
                df.loc[idx, "threshold_applied"] = threshold
                if pd.isna(outer_mass):
                    df.loc[idx, ["A_cn", "B_cn", "confidence", "tier"]] = [4, 2, 0.50, "no_baf"]
                    df.loc[idx, "decision_rule"] = f"CN=6; no BAF -> default AI (4:2); threshold={threshold:.3f}"
                elif outer_mass > threshold:
                    df.loc[idx, ["A_cn", "B_cn", "confidence", "tier"]] = [4, 2, 0.55, "tier3"]
                    df.loc[idx, "decision_rule"] = f"CN=6; outer_mass>{threshold:.3f} -> AI (4:2)"
                else:
                    df.loc[idx, ["A_cn", "B_cn", "confidence", "tier"]] = [3, 3, 0.60, "tier3"]
                    df.loc[idx, "decision_rule"] = f"CN=6; outer_mass<={threshold:.3f} -> Balanced (3:3)"
            ai  = ((df.loc[cn_mask, "A_cn"] == 4) & (df.loc[cn_mask, "B_cn"] == 2)).sum()
            bal = ((df.loc[cn_mask, "A_cn"] == 3) & (df.loc[cn_mask, "B_cn"] == 3)).sum()
            print(f"  → AI (4:2): {ai}, Balanced (3:3): {bal}")

        elif total_cn == 7:
            heuristic_threshold = 0.95
            for idx in cn_segments.index:
                outer_mass = df.loc[idx, score_col]
                df.loc[idx, "threshold_applied"] = heuristic_threshold
                if pd.isna(outer_mass):
                    df.loc[idx, ["A_cn", "B_cn", "confidence", "tier"]] = [4, 3, 0.40, "no_baf"]
                    df.loc[idx, "decision_rule"] = f"CN=7; no BAF -> default AI (4:3); heuristic_threshold={heuristic_threshold:.3f}"
                elif outer_mass > heuristic_threshold:
                    df.loc[idx, ["A_cn", "B_cn", "confidence", "tier"]] = [7, 0, 0.60, "tier3_heuristic"]
                    df.loc[idx, "decision_rule"] = f"CN=7; outer_mass>{heuristic_threshold:.3f} -> heuristic LoH (7:0)"
                else:
                    df.loc[idx, ["A_cn", "B_cn", "confidence", "tier"]] = [4, 3, 0.50, "tier3_ai_default"]
                    df.loc[idx, "decision_rule"] = f"CN=7; outer_mass<={heuristic_threshold:.3f} -> default AI (4:3)"
            loh = ((df.loc[cn_mask, "A_cn"] == 7) & (df.loc[cn_mask, "B_cn"] == 0)).sum()
            ai  = ((df.loc[cn_mask, "A_cn"] == 4) & (df.loc[cn_mask, "B_cn"] == 3)).sum()
            print(f"  → LoH (7:0): {loh}, AI (4:3): {ai}")

        elif total_cn == 8:
            threshold = thresholds.get("8", {}).get("ai_vs_balanced", 0.81)
            for idx in cn_segments.index:
                outer_mass = df.loc[idx, score_col]
                df.loc[idx, "threshold_applied"] = threshold
                if pd.isna(outer_mass):
                    df.loc[idx, ["A_cn", "B_cn", "confidence", "tier"]] = [6, 2, 0.50, "no_baf"]
                    df.loc[idx, "decision_rule"] = f"CN=8; no BAF -> default AI (6:2); threshold={threshold:.3f}"
                elif outer_mass > threshold:
                    df.loc[idx, ["A_cn", "B_cn", "confidence", "tier"]] = [6, 2, 0.50, "tier3"]
                    df.loc[idx, "decision_rule"] = f"CN=8; outer_mass>{threshold:.3f} -> AI (6:2)"
                else:
                    df.loc[idx, ["A_cn", "B_cn", "confidence", "tier"]] = [4, 4, 0.55, "tier3"]
                    df.loc[idx, "decision_rule"] = f"CN=8; outer_mass<={threshold:.3f} -> Balanced (4:4)"
            ai  = ((df.loc[cn_mask, "A_cn"] == 6) & (df.loc[cn_mask, "B_cn"] == 2)).sum()
            bal = ((df.loc[cn_mask, "A_cn"] == 4) & (df.loc[cn_mask, "B_cn"] == 4)).sum()
            print(f"  → AI (6:2): {ai}, Balanced (4:4): {bal}")

        else:
            heuristic_threshold = 0.95
            if total_cn % 2 == 1:
                a_cn = (total_cn + 1) // 2
                b_cn = (total_cn - 1) // 2
            else:
                a_cn = total_cn // 2 + 1
                b_cn = total_cn // 2 - 1

            for idx in cn_segments.index:
                outer_mass = df.loc[idx, score_col]
                df.loc[idx, "threshold_applied"] = heuristic_threshold
                if pd.isna(outer_mass):
                    df.loc[idx, ["A_cn", "B_cn", "confidence", "tier"]] = [a_cn, b_cn, 0.40, "no_baf"]
                    df.loc[idx, "decision_rule"] = (
                        f"CN>8; no BAF -> default AI ({a_cn}:{b_cn}); heuristic_threshold={heuristic_threshold:.3f}"
                    )
                elif outer_mass > heuristic_threshold:
                    df.loc[idx, ["A_cn", "B_cn", "confidence", "tier"]] = [total_cn, 0, 0.60, "tier3_heuristic"]
                    df.loc[idx, "decision_rule"] = f"CN>8; outer_mass>{heuristic_threshold:.3f} -> heuristic LoH ({total_cn}:0)"
                else:
                    df.loc[idx, ["A_cn", "B_cn", "confidence", "tier"]] = [a_cn, b_cn, 0.50, "tier3_ai_default"]
                    df.loc[idx, "decision_rule"] = f"CN>8; outer_mass<={heuristic_threshold:.3f} -> default AI ({a_cn}:{b_cn})"

            loh     = ((df.loc[cn_mask, "A_cn"] == total_cn) & (df.loc[cn_mask, "B_cn"] == 0)).sum()
            ai      = cn_mask.sum() - loh
            ai_state = (f"{(total_cn+1)//2}:{(total_cn-1)//2}" if total_cn % 2 == 1
                        else f"{total_cn//2+1}:{total_cn//2-1}")
            print(f"  → LoH ({total_cn}:0): {loh}, AI ({ai_state}): {ai}")

    df["A_cn"]   = df["A_cn"].round().astype("Int64")
    df["B_cn"]   = df["B_cn"].round().astype("Int64")
    df["ploidy"] = df["totalCN"].mean()
    df["final_class"] = df.apply(
        lambda r: classify_final_state_from_acn_bcn(r["A_cn"], r["B_cn"]), axis=1
    )

    print(f"\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")
    print(f"Total segments processed: {len(df)}")
    if len(df) > 0:
        print(f"Ploidy: {df['ploidy'].iloc[0]:.2f}")
        print("\nAllelic state distribution:")
        for (a, b), count in df.groupby(["A_cn", "B_cn"]).size().sort_values(ascending=False).head(10).items():
            pct = 100 * count / len(df)
            state_type = classify_final_state_from_acn_bcn(a, b)
            print(f"  {a}:{b} ({state_type:8s}): {count:5d} ({pct:5.1f}%)")
        print("\nConfidence tier distribution:")
        for tier in ["tier1", "tier2", "tier3", "tier3_heuristic", "tier3_ai_default", "no_baf"]:
            tc = (df["tier"] == tier).sum()
            if tc > 0:
                print(f"  {tier:20s}: {tc:5d} ({100*tc/len(df):5.1f}%), "
                      f"avg conf: {df.loc[df['tier']==tier,'confidence'].mean():.2f}")
        print("\nOverall class distribution:")
        class_counts = df["final_class"].value_counts()
        for cls in ["Deletion", "LoH", "AI", "Balanced", "Unknown"]:
            if cls in class_counts:
                count = class_counts[cls]
                print(f"  {cls:9s}: {count:5d} ({100*count/len(df):5.1f}%)")

    df_renamed = df.rename(columns={
        "chromosome": "Chromosome",
        "start": "Start_position",
        "end": "End_position",
    }).copy()
    df_renamed["SampleID"] = sample

    output_df = df_renamed[[
        "SampleID", "Chromosome", "Start_position", "End_position",
        "total_cn", "A_cn", "B_cn", "ploidy",
    ]].copy()

    return output_df, df


def merge_adjacent_segments_by_allelic_state(annotated_df, sample, max_gap_bp=1):
    """Merge contiguous segments that share the same A_cn and B_cn."""
    df = annotated_df.copy()

    if len(df) == 0:
        empty = pd.DataFrame(columns=[
            "SampleID", "Chromosome", "Start_position", "End_position",
            "total_cn", "A_cn", "B_cn", "ploidy",
        ])
        return empty, df

    for col in ["chromosome", "start", "end", "A_cn", "B_cn"]:
        if col not in df.columns:
            raise ValueError(f"Missing required column for merging: {col}")

    df = df.sort_values(["chromosome", "start", "end"]).reset_index(drop=True)

    blocks = []
    current_rows = []

    def flush_block(rows):
        if not rows:
            return None
        block = pd.DataFrame(rows).copy()
        a_cn  = block["A_cn"].iloc[0]
        b_cn  = block["B_cn"].iloc[0]
        chrom = block["chromosome"].iloc[0]
        return {
            "sample":    block["sample"].iloc[0] if "sample" in block.columns else sample,
            "dp_filter": block["dp_filter"].iloc[0] if "dp_filter" in block.columns else np.nan,
            "chromosome": chrom,
            "start": block["start"].min(),
            "end":   block["end"].max(),
            "total_cn": int(a_cn + b_cn) if pd.notna(a_cn) and pd.notna(b_cn) else np.nan,
            "A_cn": a_cn,
            "B_cn": b_cn,
            "ploidy":     block["ploidy"].iloc[0] if "ploidy" in block.columns else np.nan,
            "confidence": block["confidence"].min() if "confidence" in block.columns else np.nan,
            "tier":       worst_tier(block["tier"]) if "tier" in block.columns else "unknown",
            "final_class": classify_final_state_from_acn_bcn(a_cn, b_cn),
            "n_original_segments": len(block),
            "merged_from_segment_ids": (
                ",".join(block["segment_id"].astype(str).tolist())
                if "segment_id" in block.columns else ""
            ),
            "score_column_used": block["score_column_used"].iloc[0] if "score_column_used" in block.columns else "",
            "outer_mass_used_mean":   block["outer_mass_used"].mean()   if "outer_mass_used" in block.columns else np.nan,
            "outer_mass_used_median": block["outer_mass_used"].median() if "outer_mass_used" in block.columns else np.nan,
            "inner_mass_mean": block["inner_mass"].mean() if "inner_mass" in block.columns else np.nan,
            "baf_mean_mean":   block["baf_mean"].mean()  if "baf_mean" in block.columns else np.nan,
            "baf_sd_mean":     block["baf_sd"].mean()    if "baf_sd" in block.columns else np.nan,
            "threshold_applied_min": block["threshold_applied"].min() if "threshold_applied" in block.columns else np.nan,
            "threshold_applied_max": block["threshold_applied"].max() if "threshold_applied" in block.columns else np.nan,
            "decision_rule_summary": (
                "; ".join(block["decision_rule"].dropna().unique().tolist())
                if "decision_rule" in block.columns else ""
            ),
            "total_cn_rounding_rule": (
                block["total_cn_rounding_rule"].iloc[0]
                if "total_cn_rounding_rule" in block.columns else ""
            ),
        }

    for _, row in df.iterrows():
        row_dict = row.to_dict()
        if not current_rows:
            current_rows.append(row_dict)
            continue

        prev = current_rows[-1]
        same_chrom = prev["chromosome"] == row_dict["chromosome"]
        same_acn   = prev["A_cn"] == row_dict["A_cn"]
        same_bcn   = prev["B_cn"] == row_dict["B_cn"]
        adjacent   = row_dict["start"] <= prev["end"] + max_gap_bp

        if same_chrom and same_acn and same_bcn and adjacent:
            current_rows.append(row_dict)
        else:
            merged_row = flush_block(current_rows)
            if merged_row:
                blocks.append(merged_row)
            current_rows = [row_dict]

    merged_row = flush_block(current_rows)
    if merged_row:
        blocks.append(merged_row)

    merged_annotated_df = pd.DataFrame(blocks)

    merged_scarhrd_df = merged_annotated_df.rename(columns={
        "chromosome": "Chromosome",
        "start": "Start_position",
        "end": "End_position",
    })[["sample", "Chromosome", "Start_position", "End_position",
        "total_cn", "A_cn", "B_cn", "ploidy"]].copy()
    merged_scarhrd_df.insert(0, "SampleID", sample)

    n_before = len(annotated_df)
    n_after  = len(merged_annotated_df)
    print(f"  Merged {n_before} → {n_after} segments "
          f"({n_before - n_after} segments collapsed, "
          f"reduction {100*(n_before-n_after)/max(n_before,1):.1f}%)")

    return merged_scarhrd_df, merged_annotated_df


# ============================================================
# STAGE 1: PER-SAMPLE BAF ANALYSIS
# ============================================================

def run_script1():
    unified_logger = UnifiedLogger(
        log_file=UNIFIED_LOG_FILE,
        sample_id=SAMPLE_ID_FILTER,
        verbose_console=False,
    )
    _real_stdout = sys.stdout
    sys.stdout = unified_logger

    try:
        print(f"{'='*70}")
        if PARALLEL_MODE:
            print(f"PARALLEL MODE - Processing sample: {SAMPLE_ID_FILTER}")
            print(f"Long name filter: {LONG_NAME_FILTER}")
            print(f"Folder override:  {FOLDER_PATH_OVERRIDE}")
        else:
            print("SEQUENTIAL MODE - Processing all samples")
        print(f"{'='*70}")
        print(f"Unified log: {UNIFIED_LOG_FILE}")
        print(f"Output dir:  {OUTPUT_DIR}")
        print(f"{'='*70}\n")

        if not PARALLEL_MODE:
            print(f"\n{'='*70}")
            print("SAVING SCRIPT BACKUP")
            print(f"{'='*70}")
            latest_file = os.path.join(SCRIPTS_BACKUP_DIR, "pipeline_latest.py")
            try:
                shutil.copy(__file__, latest_file)
                print(f"✓ Latest copy saved: {latest_file}")
            except Exception as e:
                print(f"⚠️  Could not save script: {e}")
            print(f"{'='*70}\n")

        purity_mapping, ploidy_mapping = load_purity_mappings(PURITY_FILE)

        all_segments_list = []
        manifest_rows     = []

        print("\n" + "=" * 70)
        print("SEGMENT-LEVEL ANALYSIS")
        print("=" * 70)

        for folder_idx, folder_path in enumerate(FOLDER_PATHS):
            folder_label = os.path.basename(folder_path)

            print(f'\n{"="*70}')
            print(f"FOLDER {folder_idx + 1}/{len(FOLDER_PATHS)}: {folder_label}")
            print(f"Path: {folder_path}")
            print(f'{"="*70}')

            if PARALLEL_MODE:
                samples_dict = {SAMPLE_ID_FILTER: LONG_NAME_FILTER}
            else:
                samples_dict = create_sample_name_dict_custom(folder_path)

            segments_master = load_segments_for_samples(samples_dict, RAW_SEGMENTS_DIR)

            if len(segments_master) == 0:
                print("  ⚠️  No segments loaded, skipping folder")
                continue

            for sample_id, long_name in samples_dict.items():
                actual_suffix = detect_baf_suffix(folder_path, long_name)
                for file_type in TYPES_OF_FILES:
                    print(f'\n{"="*70}')
                    print(f"Sample: {sample_id} | BAF suffix: {actual_suffix}")
                    print(f'{"="*70}')

                    # Tag outputs by the detected BAF file type so the manifest
                    # is meaningful when both TSO and ASCAT files are present.
                    if actual_suffix == ".tumour.tso500_sparse_downsample.txt":
                        file_label = "ascat_downsampled"
                    else:
                        file_label = "original_tso"
                    stem = build_output_stem(sample_id, file_label, folder_label)

                    prs = Presentation()
                    prs.slide_width  = Inches(11.69)
                    prs.slide_height = Inches(8.27)

                    segments = segments_master[segments_master["SampleID"] == sample_id].copy()
                    print(f"original segments length: {len(segments)}")
                    if len(segments) == 0:
                        print("  ⚠️  No segments for this sample, skipping")
                        continue

                    segments = segments.sort_values(["chromosome", "start", "end"]).reset_index(drop=True)
                    segments["segment_id"]     = np.arange(len(segments), dtype=int)
                    segments["segment_len_bp"] = segments["end"] - segments["start"] + 1

                    sample_baf_hit_list  = []
                    sample_segments_list = []
                    sample_purity_used   = np.nan
                    sample_purity_note   = "no_purity_file" if not purity_mapping else "no_purity_data"

                    for dp in DP_VALUES:
                        print(f"\n  Processing DP≥{dp}...")

                        baf_file = os.path.join(folder_path, f"{long_name}{actual_suffix}")
                        print(f"    Loading BAF from: {baf_file}")

                        if not os.path.exists(baf_file):
                            print("    ⚠️  BAF file not found, skipping")
                            continue

                        baf_raw = normalize_baf_df(pd.read_csv(baf_file, sep="\t", compression="infer"))
                        print(f"    loaded positions: {len(baf_raw):,}")

                        baf_filtered = baf_raw[
                            (baf_raw["DP"] >= dp) & baf_raw["_ok"]
                        ].copy()
                        print(f"    after DP≥{dp}: {len(baf_filtered):,}")

                        if len(baf_filtered) == 0:
                            print("    ⚠️  No data after DP filter")
                            continue

                        baf_filtered["chromosome"] = baf_filtered["CHROM"].astype(str)
                        baf_filtered["BAF_n"] = (
                            baf_filtered["Alt_Reads"] /
                            (baf_filtered["Ref_Reads"] + baf_filtered["Alt_Reads"])
                        )

                        baf_df = prepare_baf_origin(
                            baf_filtered,
                            baf_col="BAF_n",
                            het_chr_col="chromosome",
                            baf_missing="./.",
                            depth=dp,
                        )

                        n_positions = len(baf_df)
                        print(f"    BAF positions (raw): {n_positions:,}")

                        if n_positions == 0:
                            print("    ⚠️  No BAF data")
                            continue

                        if not purity_mapping:
                            print("    ⚠️  No purity data available - skipping correction")
                            baf_df["BAF_raw"]       = baf_df["BAF_n"]
                            baf_df["BAF_corrected"] = baf_df["BAF_n"]
                            baf_df["purity_used"]   = np.nan
                        else:
                            sample_purity = purity_mapping.get(long_name) or purity_mapping.get(sample_id)
                            if pd.isna(sample_purity) or sample_purity is None:
                                print(f"    ⚠️  WARNING: No purity found for {sample_id} - proceeding without correction")
                                baf_df["BAF_raw"]       = baf_df["BAF_n"]
                                baf_df["BAF_corrected"] = baf_df["BAF_n"]
                                baf_df["purity_used"]   = np.nan
                            else:
                                print(f"    ✓ Found purity: {sample_purity:.4f}")
                                sample_purity_used = sample_purity
                                sample_purity_note = "ok"
                                baf_df = apply_purity_correction(baf_df, sample_purity, baf_col="BAF_n")
                                comparison_dir = os.path.join(BAF_CORRECTION_PLOTS_DIR, stem)
                                plot_baf_correction_comparison(baf_df, sample_id, dp, sample_purity, comparison_dir)

                        baf_df["BAF_n"] = baf_df["BAF_corrected"]
                        print(f"    BAF positions (corrected): {len(baf_df):,}")

                        print("    Mapping BAF to segments...")
                        segments_with_baf, baf_hits = map_baf_to_segments_FINE(
                            segments, baf_df,
                            seg_chr="chromosome", seg_start="start", seg_end="end",
                            baf_chr="CHROM", baf_pos="POS",
                            baf_val_corrected="BAF_corrected", baf_val_raw="BAF_raw",
                            n_bins=N_BINS,
                        )

                        print("    Calculating segment-level BAF summary metrics...")
                        segment_metrics = calculate_segment_baf_summary_metrics(baf_hits)

                        segments_final = segments_with_baf.merge(
                            segment_metrics, on="segment_id", how="left"
                        )

                        comparison_dir = os.path.join(BAF_CORRECTION_PLOTS_DIR, stem)
                        plot_outer_mass_correction_comparison(segments_final, sample_id, dp, comparison_dir)

                        segments_final["sample_name"]          = sample_id
                        segments_final["sample"]               = long_name
                        segments_final["file_type"]            = file_label
                        segments_final["dp_filter"]            = dp
                        segments_final["n_baf_positions_total"] = n_positions
                        segments_final["baf_folder"]           = folder_label
                        segments_final["purity_used"]          = (
                            baf_df["purity_used"].iloc[0]
                            if "purity_used" in baf_df.columns and len(baf_df) > 0
                            else np.nan
                        )

                        if len(baf_hits) > 0:
                            baf_hits["sample_name"]          = sample_id
                            baf_hits["sample"]               = long_name
                            baf_hits["file_type"]            = file_label
                            baf_hits["dp_filter"]            = dp
                            baf_hits["n_baf_positions_total"] = n_positions
                            baf_hits["baf_folder"]           = folder_label
                            baf_hits["purity_used"]          = (
                                segments_final["purity_used"].iloc[0]
                                if len(segments_final) > 0 else np.nan
                            )

                        sample_segments_list.append(segments_final)
                        sample_baf_hit_list.append(baf_hits)
                        all_segments_list.append(segments_final)

                    # Ensure REFcounts / ALTcounts are present for Script 2
                    if "Ref_Reads" in baf_df.columns and "REFcounts" not in baf_df.columns:
                        baf_df["REFcounts"] = baf_df["Ref_Reads"]
                    if "Alt_Reads" in baf_df.columns and "ALTcounts" not in baf_df.columns:
                        baf_df["ALTcounts"] = baf_df["Alt_Reads"]

                    # Save BAF hits
                    baf_hits_output_file = os.path.join(
                        BAF_HITS_PER_SAMPLE_DIR, f"alldp_baf_hits_{stem}.csv"
                    )
                    non_empty_hits = [x for x in sample_baf_hit_list if len(x) > 0]
                    if non_empty_hits:
                        all_baf_hits_df = pd.concat(non_empty_hits, ignore_index=True)
                        all_baf_hits_df.to_csv(baf_hits_output_file, index=False)
                        print(f"  ✓ Saved BAF hits: {baf_hits_output_file} ({len(all_baf_hits_df):,} rows)")
                    else:
                        all_baf_hits_df = pd.DataFrame()
                        print(f"  ⚠️  No non-empty BAF hits for sample {sample_id}")

                    # Save segments and generate scarHRD inputs
                    segments_output_file = os.path.join(
                        SEGMENTS_PER_SAMPLE_DIR, f"segments_{stem}.csv"
                    )
                    if sample_segments_list:
                        sample_segments_df = pd.concat(sample_segments_list, ignore_index=True)
                        sample_segments_df.to_csv(segments_output_file, index=False)
                        print(f"  ✓ Saved segments: {segments_output_file} ({len(sample_segments_df):,} rows)")

                        print(f'\n{"="*70}')
                        print(f"GENERATING scarHRD INPUT - {sample_id} - {file_label}")
                        print(f'{"="*70}')

                        segments_auto = sample_segments_df.copy()
                        segments_auto["totalCN"]   = segments_auto["total_cn"]
                        segments_auto["mid_mass"]  = segments_auto["inner_mass"]
                        segments_auto["Chromosome"] = (
                            segments_auto["chromosome"].astype(str).str.strip()
                            .str.replace("^chr", "", regex=True)
                            .str.replace(r"\.0$", "", regex=True)
                        )
                        segments_auto = segments_auto[
                            segments_auto["Chromosome"].isin([str(i) for i in range(1, 23)])
                        ].drop_duplicates(keep="first")
                        print(f"  Autosomal segment rows: {len(segments_auto):,}")

                        for dp in sorted(segments_auto["dp_filter"].dropna().unique().tolist()):
                            print(f"\n  Processing DP≥{dp}...")
                            df_sample = segments_auto[segments_auto["dp_filter"] == dp].copy()
                            if len(df_sample) == 0:
                                print(f"    ⚠️  No segments for DP≥{dp}")
                                continue
                            try:
                                scarhrd_input = create_input_data_scarHRD(df_sample, sample_id)
                                if len(scarhrd_input) > 0:
                                    scarhrd_file = os.path.join(
                                        SCARHRD_INPUT_DIR,
                                        f"{stem}_dp{dp}_scarHRD_input.txt",
                                    )
                                    scarhrd_input.to_csv(scarhrd_file, sep="\t", index=False)
                                    print(f"    ✓ Saved: {scarhrd_file} ({len(scarhrd_input):,} rows)")
                                else:
                                    print(f"    ⚠️  No output for DP≥{dp}")
                            except Exception as e:
                                print(f"    ✗ Error generating scarHRD input for DP≥{dp}: {e}")

                        manifest_rows.append({
                            "sample_name":    sample_id,
                            "sample":         long_name,
                            "file_type":      file_label,
                            "baf_folder":     folder_label,
                            "segment_format": (
                                sample_segments_df["segment_format"].iloc[0]
                                if "segment_format" in sample_segments_df.columns and len(sample_segments_df) > 0
                                else "unknown"
                            ),
                            "segments_file":       segments_output_file,
                            "baf_hits_file":       baf_hits_output_file if len(all_baf_hits_df) > 0 else "",
                            "purity_used":         sample_purity_used,
                            "purity_note":         sample_purity_note,
                            "dp_filters":          ",".join(map(str, sorted(
                                sample_segments_df["dp_filter"].dropna().unique().tolist()
                            ))),
                            "scarhrd_output_prefix": os.path.join(SCARHRD_INPUT_DIR, stem),
                        })
                    else:
                        print(f"  ⚠️  No processed segments for sample {sample_id}")

        # Save global combined segment file
        if all_segments_list:
            all_segments_df = pd.concat(all_segments_list, ignore_index=True)
            global_output_file = os.path.join(OUTPUT_DIR, "segments_all_samples_combined.csv")
            all_segments_df.to_csv(global_output_file, index=False)
            print(f"\n✓ Saved combined segment dataset: {global_output_file} ({len(all_segments_df):,} rows)")
        else:
            print("\n⚠️  No segments collected - skipping global save")

        # Save manifest
        if manifest_rows:
            manifest_df = pd.DataFrame(manifest_rows)
            manifest_df.to_csv(MANIFEST_FILE, index=False)
            print(f"\n✓ Saved pipeline manifest: {MANIFEST_FILE} ({len(manifest_df):,} rows)")
        else:
            print("\n⚠️  No manifest rows collected - manifest not written")

        # Classify BAF hits by CN match (only in sequential mode and only when allelic CN is present)
        if not PARALLEL_MODE:
            print(f"\n{'='*70}")
            print("CHECKING FOR CLASSIFICATION CAPABILITY")
            print(f"{'='*70}")

            baf_hits_dir      = Path(BAF_HITS_PER_SAMPLE_DIR)
            segments_dir      = Path(SEGMENTS_PER_SAMPLE_DIR)
            classification_dir = Path(os.path.join(OUTPUT_DIR, "baf_hits_segments_classified"))
            classification_dir.mkdir(exist_ok=True, parents=True)

            file_list_baf = list(baf_hits_dir.glob("alldp_baf_hits_*.csv"))
            file_list_seg = list(segments_dir.glob("segments_*.csv"))

            print(f"\nFound {len(file_list_baf)} BAF hit files")
            print(f"Found {len(file_list_seg)} segment files")

            if not file_list_baf or not file_list_seg:
                print("⚠️  No files to process, skipping classification")
            else:
                test_seg = pd.read_csv(file_list_seg[0])
                has_allelic_cn = (
                    "nMajor" in test_seg.columns and
                    "nMinor" in test_seg.columns and
                    test_seg["nMajor"].notna().any()
                )

                if not has_allelic_cn:
                    print(f"\n⚠️  CLASSIFICATION SKIPPED")
                    print("Reason: Segment files do NOT contain usable nMajor/nMinor columns")
                else:
                    print(f"\n✓ ALLELIC CN DATA DETECTED - PROCEEDING WITH CLASSIFICATION")

                    seg_file_lookup = {
                        f.name.replace("segments_", "").replace(".csv", ""): f
                        for f in file_list_seg
                    }
                    processed = 0
                    skipped   = 0

                    for i, baf_file in enumerate(file_list_baf):
                        print(f"\n[{i+1}/{len(file_list_baf)}] Processing: {baf_file.name}")
                        core = baf_file.name.replace("alldp_baf_hits_", "").replace(".csv", "")

                        if core not in seg_file_lookup:
                            print("  ⚠️  No matching segment file")
                            skipped += 1
                            continue

                        try:
                            baf_temp = pd.read_csv(baf_file)
                            seg_temp = pd.read_csv(seg_file_lookup[core])

                            seg_temp["totalCN"]    = seg_temp["total_cn"]
                            seg_temp["mid_mass"]   = seg_temp["inner_mass"]
                            seg_temp["Chromosome"] = seg_temp["chromosome"]

                            scarhrd_input, seg_temp2 = create_input_data_scarHRD_with_gold(seg_temp)
                            seg_temp2["nMajor_gold"] = seg_temp2["nMajor"]
                            seg_temp2["nMinor_gold"] = seg_temp2["nMinor"]
                            seg_temp2["nMajor"]      = seg_temp2["A_cn"]
                            seg_temp2["nMinor"]      = seg_temp2["B_cn"]
                            seg_temp2["cn_classification"] = seg_temp2.apply(classify_cn_match, axis=1)

                            bigmerged = baf_temp.merge(
                                seg_temp2[[
                                    "sample", "segment_id", "dp_filter",
                                    "nMajor", "nMinor", "nMajor_gold", "nMinor_gold",
                                    "cn_classification", "outer_mass", "inner_mass",
                                ]],
                                on=["sample", "segment_id", "dp_filter"],
                                how="left",
                                suffixes=("", "_seg"),
                            )

                            output_file = classification_dir / f"{core}_classified.csv"
                            bigmerged.to_csv(output_file, index=False)
                            print(f"  ✓ Saved: {output_file.name}")
                            processed += 1

                        except Exception as e:
                            print(f"  ✗ Error: {e}")
                            traceback.print_exc()
                            skipped += 1

                        if (i + 1) % 10 == 0:
                            gc.collect()

                    print(f"\nCLASSIFICATION COMPLETE: {processed}/{len(file_list_baf)} processed, {skipped} skipped")

        print(f"\n{'='*70}")
        print("✓ STAGE 1 COMPLETE")
        print(f"{'='*70}")

    finally:
        sys.stdout = _real_stdout
        unified_logger.close()


# ============================================================
# STAGE 2: CONSOLIDATION
# ============================================================

def run_script2():
    print("╔" + "=" * 74 + "╗")
    print("║              BAF Analysis Workflow - CONSOLIDATION              ║")
    print("╚" + "=" * 74 + "╝")
    print(f"Manifest file:     {MANIFEST_FILE}")
    print(f"Output file:       {MERGED_OUTPUT_FILE}")
    print(f"Enhanced features: {'✓ ENABLED' if ENHANCED_FEATURES_AVAILABLE else '✗ DISABLED'}")
    print()

    manifest_df = load_manifest(MANIFEST_FILE)
    manifest_df = deduplicate_manifest(manifest_df)
    print(f"✓ Loaded manifest with {len(manifest_df):,} rows\n")

    assert_required_columns(manifest_df, [
        "sample_name", "sample", "file_type", "baf_folder", "segments_file", "baf_hits_file",
    ], "manifest_df")

    print("Loading canonical segment datasets...\n")
    all_merged_list  = []
    processing_summary = []

    for idx, row in manifest_df.iterrows():
        sample_name = row["sample_name"]
        sample      = row["sample"]
        file_type   = row["file_type"]
        baf_folder  = row["baf_folder"]
        seg_file    = row["segments_file"]
        baf_file    = row["baf_hits_file"] if isinstance(row["baf_hits_file"], str) else ""

        print(f"[{idx + 1}/{len(manifest_df)}] sample_name={sample_name} | file_type={file_type} | folder={baf_folder}")

        try:
            if not os.path.exists(seg_file):
                raise FileNotFoundError(f"Segment file not found: {seg_file}")

            seg_df = pd.read_csv(seg_file)
            print(f"  Loaded processed segments: {len(seg_df):,} rows")

            assert_required_columns(seg_df, [
                "sample", "sample_name", "segment_id", "dp_filter", "baf_folder",
                "chromosome", "start", "end", "total_cn",
                "outer_mass", "inner_mass",
                "baf_mean", "baf_median", "baf_sd", "baf_q25", "baf_q75", "baf_iqr",
            ], "processed segment file")

            merged = seg_df.copy()

            if ENHANCED_FEATURES_AVAILABLE:
                if baf_file and os.path.exists(baf_file):
                    print("  Loading BAF hits for enhanced features...")
                    baf_temp = pd.read_csv(baf_file)
                    print(f"  Loaded BAF hits: {len(baf_temp):,} rows")

                    if "REFcounts" not in baf_temp.columns:
                        baf_temp["REFcounts"] = baf_temp.get("Ref_Reads", baf_temp.get("REF", np.nan))
                    if "ALTcounts" not in baf_temp.columns:
                        baf_temp["ALTcounts"] = baf_temp.get("Alt_Reads", baf_temp.get("ALT", np.nan))
                    if "chromosome" not in baf_temp.columns and "CHROM" in baf_temp.columns:
                        baf_temp["chromosome"] = baf_temp["CHROM"]
                    if "position" not in baf_temp.columns and "POS" in baf_temp.columns:
                        baf_temp["position"] = baf_temp["POS"]

                    missing_needed = [c for c in ["REFcounts", "ALTcounts"] if c not in baf_temp.columns]
                    if missing_needed:
                        print(f"  ⚠️  Missing columns for enhanced features: {missing_needed} — skipping")
                    else:
                        print("  Adding enhanced features...")
                        merged = add_enhanced_features_to_merged(merged, baf_temp)
                else:
                    print("  ⚠️  BAF hits file missing; skipping add_enhanced_features_to_merged")

                print("  Adding length-normalized features...")
                merged = add_length_normalized_features(merged)
            else:
                print("  ⚠️  Enhanced features disabled")

            merged["manifest_sample_name"] = sample_name
            merged["manifest_sample"]      = sample
            merged["manifest_file_type"]   = file_type
            merged["manifest_baf_folder"]  = baf_folder

            all_merged_list.append(merged)
            processing_summary.append({
                "sample_name": sample_name, "sample": sample,
                "file_type": file_type, "baf_folder": baf_folder,
                "segment_rows": len(seg_df), "merged_rows": len(merged), "status": "ok",
            })
            print(f"  Final merged dataset: {len(merged):,} rows × {len(merged.columns)} columns")
            print("  ✓ Completed\n")

        except Exception as e:
            print(f"  ✗ Error processing manifest row {idx}: {e}")
            traceback.print_exc()
            processing_summary.append({
                "sample_name": sample_name, "sample": sample,
                "file_type": file_type, "baf_folder": baf_folder,
                "segment_rows": np.nan, "merged_rows": np.nan, "status": f"error: {e}",
            })
            print()

    if not all_merged_list:
        print("\n⚠️  No datasets were successfully processed")
        print("\n╔" + "=" * 74 + "╗")
        print("║                      PROCESSING COMPLETE                      ║")
        print("╚" + "=" * 74 + "╝")
        return

    print(f"\n{'=' * 70}")
    print("CONCATENATING ALL DATASETS")
    print(f"{'=' * 70}")

    merged_all = pd.concat(all_merged_list, ignore_index=True)
    print(f"Combined dataset: {len(merged_all):,} rows")
    print(f"Unique samples:   {merged_all['sample'].nunique() if 'sample' in merged_all.columns else 'N/A'}")
    print(f"Total columns:    {len(merged_all.columns)}")

    if "segment_format" in merged_all.columns:
        print("\nSegment format distribution:")
        for fmt, count in merged_all["segment_format"].value_counts(dropna=False).items():
            print(f"  {str(fmt):15s}: {count:>10,} ({count/len(merged_all)*100:>5.1f}%)")

    if "nMajor" in merged_all.columns and "nMinor" in merged_all.columns:
        nan_nmajor = merged_all["nMajor"].isna().sum()
        print(f"\nRows with NaN nMajor: {nan_nmajor:,} ({nan_nmajor/len(merged_all)*100:.1f}%)")

    if "dp_filter" in merged_all.columns:
        print("\nDP filter distribution:")
        for dp, count in merged_all["dp_filter"].value_counts(dropna=False).sort_index().items():
            print(f"  DP={dp}: {count:>10,} ({count/len(merged_all)*100:>5.1f}%)")

    merged_all.to_csv(MERGED_OUTPUT_FILE, index=False)
    print(f"\n✓ Saved combined dataset: {MERGED_OUTPUT_FILE}")

    pd.DataFrame(processing_summary).to_csv(SUMMARY_OUTPUT_FILE, index=False)
    print(f"✓ Saved processing summary: {SUMMARY_OUTPUT_FILE}")

    print(f"\n{'=' * 70}")
    print("COLUMN SUMMARY")
    print(f"{'=' * 70}")
    count_cols     = [c for c in merged_all.columns if c.startswith("count_") and "raw" not in c]
    ratio_cols     = [c for c in merged_all.columns if c.startswith("ratio_") and "raw" not in c]
    count_raw_cols = [c for c in merged_all.columns if c.startswith("count_raw_")]
    ratio_raw_cols = [c for c in merged_all.columns if c.startswith("ratio_raw_")]
    print(f"Corrected BAF count columns: {len(count_cols)}")
    print(f"Corrected BAF ratio columns: {len(ratio_cols)}")
    print(f"Raw BAF count columns:       {len(count_raw_cols)}")
    print(f"Raw BAF ratio columns:       {len(ratio_raw_cols)}")
    print(f"Total BAF bin columns:       {len(count_cols)+len(ratio_cols)+len(count_raw_cols)+len(ratio_raw_cols)}")

    if "outer_mass" in merged_all.columns:
        print(f"\n  outer_mass range: [{merged_all['outer_mass'].min():.3f}, {merged_all['outer_mass'].max():.3f}]")
    if "inner_mass" in merged_all.columns:
        print(f"  inner_mass range: [{merged_all['inner_mass'].min():.3f}, {merged_all['inner_mass'].max():.3f}]")
    if "baf_mean" in merged_all.columns:
        print(f"  baf_mean range:   [{merged_all['baf_mean'].min():.3f}, {merged_all['baf_mean'].max():.3f}]")

    print("\n╔" + "=" * 74 + "╗")
    print("║                      STAGE 2 COMPLETE                        ║")
    print("╚" + "=" * 74 + "╝")


# ============================================================
# STAGE 3: ALLELE-SPECIFIC CN INFERENCE → scarHRD INPUTS
# ============================================================

def run_script3():
    print(f"{'='*70}")
    print("BAF-TO-HRD PIPELINE — STAGE 3: ALLELE-SPECIFIC CN INFERENCE")
    print(f"{'='*70}")
    print(f"Input file:      {MERGED_OUTPUT_FILE}")
    print(f"Thresholds file: {THRESHOLDS_FILE}")
    print(f"Output dir:      {SCARHRD_THRESHOLDS_DIR}")
    print(f"DP filters:      {DP_VALUES}")
    print(f"Max gap (merge): {MAX_GAP_BP_FOR_MERGE} bp")
    print()

    df = pd.read_csv(MERGED_OUTPUT_FILE)
    df["total_cn_raw"] = pd.to_numeric(df["total_cn"], errors="coerce")
    df["total_cn"]     = np.floor(df["total_cn_raw"] + 0.5).astype("Int64")

    with open(THRESHOLDS_FILE, "r") as f:
        thresholds_raw = json.load(f)

    # Detect JSON format:
    #   Nested (new): {dp: {cn: {comparison_key: threshold}}}
    #   Flat   (old): {cn: {comparison_key: threshold}}
    # Use DP_VALUES membership as the discriminator: if any outer key is a known
    # DP filter value, the file is nested.  Checking the value structure is
    # unreliable when a DP key has an empty dict (no thresholds were fitted).
    _is_nested = any(
        k.lstrip('-').isdigit() and int(k) in DP_VALUES
        for k in thresholds_raw
    )
    if _is_nested:
        thresholds_by_dp = {int(k): v for k, v in thresholds_raw.items()}
        print(f"✓ Loaded per-DP thresholds for DP filters: {sorted(thresholds_by_dp.keys())}")
    else:
        # Flat format — replicate the same thresholds for every DP filter
        thresholds_by_dp = {dp: thresholds_raw for dp in DP_VALUES}
        print(f"⚠️  Loaded flat (single-DP) thresholds — applied to all DP filters: {DP_VALUES}")
        print("    Re-run calibrate_thresholds.py to generate per-DP thresholds.")

    # Build (sample, baf_folder) combination list.
    # The same sample UUID can appear under multiple baf_folder values (e.g. one entry
    # per downsampling backbone such as JBLAB231 and JBLAB17041).  Looping on sample
    # alone would mix both backbones into a single file, doubling the segment count.
    _has_folder = "baf_folder" in df.columns
    if _has_folder:
        combos = sorted(df[["sample", "baf_folder"]].dropna(subset=["sample"]).drop_duplicates().values.tolist())
    else:
        combos = [(s, None) for s in sorted(df["sample"].dropna().unique())]

    print(
        f"There are {len(combos)} sample×backbone combinations; "
        f"with {len(DP_VALUES)} DP filters × 2 merge states this gives "
        f"{len(combos) * len(DP_VALUES) * 2} total rows across {len(combos)} output files.\n"
    )

    def _folder_label(baf_folder):
        """Short backbone label extracted from baf_folder for filenames."""
        if not baf_folder:
            return ""
        import re
        m = re.search(r'JBLAB\w+', str(baf_folder))
        return f"_{m.group(0)}" if m else f"_{baf_folder.split('_')[-1]}"

    for sample, baf_folder in combos:
        label     = _folder_label(baf_folder)
        sample_id = f"{sample}{label}"   # e.g. uuid_JBLAB17041

        # Accumulate all DP × merge_state rows into a single file per sample×backbone.
        # Columns: standard scarHRD columns + dp_filter + merge_state.
        # SampleID encodes dp and merge_state so scarHRD output rows are identifiable.
        all_rows = []

        for dp in DP_VALUES:
            mask = (df["sample"] == sample) & (df["dp_filter"] == dp)
            if _has_folder and baf_folder is not None:
                mask &= (df["baf_folder"] == baf_folder)
            data_to_analyse = df[mask].copy()

            if len(data_to_analyse) == 0:
                print(f"⚠️  No data for {sample_id} at DP {dp}, skipping")
                continue

            if dp in thresholds_by_dp:
                thresholds = thresholds_by_dp[dp]
            else:
                nearest = min(thresholds_by_dp.keys(), key=lambda k: abs(k - dp))
                print(f"  ⚠️  No thresholds for DP={dp}, using nearest DP={nearest}")
                thresholds = thresholds_by_dp[nearest]

            sid_dp = f"{sample_id}_dp{dp}"

            output_df_premerge, annotated_df_premerge = build_scarHRD_input_using_thresholds(
                segments_df=data_to_analyse,
                sample=f"{sid_dp}_premerge",
                thresholds=thresholds,
            )
            output_df_premerge["dp_filter"]   = dp
            output_df_premerge["merge_state"] = "premerge"

            merged_scarhrd_df, _ = merge_adjacent_segments_by_allelic_state(
                annotated_df_premerge, sample=f"{sid_dp}_postmerge",
                max_gap_bp=MAX_GAP_BP_FOR_MERGE,
            )
            merged_scarhrd_df["dp_filter"]   = dp
            merged_scarhrd_df["merge_state"] = "postmerge"

            all_rows.extend([output_df_premerge, merged_scarhrd_df])

            scarhrd_cols = ["SampleID", "Chromosome", "Start_position",
                            "End_position", "total_cn", "A_cn", "B_cn", "ploidy"]
            for df_out, merge_state in [(output_df_premerge, "premerge"),
                                        (merged_scarhrd_df,  "postmerge")]:
                fname = f"{sample_id}_dp{dp}_{merge_state}_scarHRD_input.txt"
                df_out[scarhrd_cols].to_csv(
                    os.path.join(SCARHRD_THRESHOLDS_DIR, fname),
                    sep="\t", index=False,
                )

            print(f"✓ {sample_id} DP={dp}: "
                  f"pre={len(output_df_premerge)} segs → post={len(merged_scarhrd_df)} segs")

    print(f"\n{'='*70}")
    print("✓ STAGE 3 COMPLETE")
    print(f"{'='*70}")


# ============================================================
# MAIN
# ============================================================

def main():
    args = _parse_args()
    apply_config(args)

    for _d in [
        BASE_DIR, SEGMENTS_PER_SAMPLE_DIR, BAF_HITS_PER_SAMPLE_DIR,
        BAF_CORRECTION_PLOTS_DIR, SCARHRD_INPUT_DIR, SCARHRD_THRESHOLDS_DIR,
        SCRIPTS_BACKUP_DIR,
    ]:
        os.makedirs(_d, exist_ok=True)

    print(f"\n{'='*70}")
    print("BAF-TO-HRD PIPELINE — INTEGRATED")
    print(f"DATE:          {DATE}")
    print(f"OUTPUT:        {OUTPUT_DIR}")
    print(f"STAGES TO RUN: {sorted(STAGES_TO_RUN)}")
    print(f"{'='*70}\n")

    if 1 in STAGES_TO_RUN:
        print(f"\n{'='*70}")
        print("STAGE 1: PER-SAMPLE BAF ANALYSIS")
        print(f"{'='*70}\n")
        run_script1()

    if PARALLEL_MODE:
        print("\nParallel mode: stages 2 and 3 are skipped (run sequentially after all samples complete).")
        return

    if 2 in STAGES_TO_RUN:
        print(f"\n{'='*70}")
        print("STAGE 2: CONSOLIDATION")
        print(f"{'='*70}\n")
        run_script2()

    if 3 in STAGES_TO_RUN:
        print(f"\n{'='*70}")
        print("STAGE 3: ALLELE-SPECIFIC CN INFERENCE")
        print(f"{'='*70}\n")
        run_script3()

    print(f"\n{'='*70}")
    print("✓ PIPELINE COMPLETE")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
