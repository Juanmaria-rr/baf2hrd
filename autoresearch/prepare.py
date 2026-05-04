#!/usr/bin/env python3
"""
prepare.py — run once before experiments begin.

Builds val_data.csv: segments with pre-computed BAF bin ratios (20 bins)
and ASCAT ground-truth allele-specific copy numbers (nMajor / nMinor).

These are the segments on which every experiment is scored.

Requires:
  - Stage 2 merged CSV (merged_segments_baf_metrics_ALL_with20bins_enhancedFeatures.csv)
    produced by baf_to_hrd_pipeline.py.  Rows where nMajor / nMinor are non-null
    correspond to ASCAT-profiled samples and serve as the validation set.
  - (Optional) a separate ASCAT segment table if the merged CSV does not already
    carry nMajor / nMinor — set ASCAT_SEGMENTS_FILE below.

Output:
  val_data.csv   — one row per segment × DP filter, columns:
      sample, sample_name, segment_id, dp_filter, chromosome, start, end,
      total_cn,
      ratio_0.00-0.05 … ratio_0.95-1.00  (20 corrected BAF bin ratios)
      ratio_raw_0.00-0.05 … ratio_raw_0.95-1.00  (20 raw BAF bin ratios)
      outer_mass_raw, inner_mass_raw,    (uncorrected — same scale as TSO500 application)
      outer_mass, inner_mass,            (purity-corrected — kept for reference)
      baf_mean, baf_sd, baf_median, baf_iqr, n_positions,
      nMajor_true, nMinor_true,          (ASCAT ground truth)
      class_true                          (LoH / Balanced / AI / Deletion)
"""

import os
import sys
import numpy as np
import pandas as pd

# ============================================================
# PATHS — adjust if needed
# ============================================================

MERGED_CSV = (
    "/storage/scratch01/groups/co/cn_extra/alleleSpecific/tso_samples/20260406"
    "/merged_segments_baf_metrics_ALL_with20bins_enhancedFeatures.csv"
)
# NOTE: The 20260421 run contains only TSO500 samples (nMajor/nMinor are all NaN).
# The 20260406 run contains 1,162 ASCAT-downsampled samples with nMajor/nMinor
# populated for all rows — these are the validation segments for the autoresearch loop.
# Column names in 20260406 use _x/_y suffixes (legacy pandas merge artifact):
#   outer_mass_y = purity-corrected outer_mass (20-bin recomputed)
#   outer_mass_x = original 4-bin outer_mass from Stage 1

# Optional separate ASCAT table.  Must have columns:
#   sample (or SampleID), chromosome, start, end, nMajor, nMinor
# Set to None to rely entirely on the nMajor/nMinor columns in MERGED_CSV.
ASCAT_SEGMENTS_FILE = None

OUTPUT_FILE = os.path.join(os.path.dirname(__file__), "val_data.csv")

# Minimum BAF positions a segment must have to be included
MIN_BAF_POSITIONS = 5

# DP filter to use for validation.  Must match the dp_filter used in
# calibrate_thresholds.py (PARAMS['dp_filter'] = 8) so that the outer_mass
# values seen during threshold calibration match those used for evaluation.
DP_FILTER_FOR_VAL = 8

# Optional: path to a plain-text file (one sample ID per line) listing the
# samples to USE for validation.  When set, only those samples are included in
# val_data.csv — use split_samples.py to generate this file so that validation
# samples are never seen during threshold calibration.
# Set to None to use all ASCAT-profiled samples (original in-sample behaviour).
VALIDATION_SAMPLES_FILE = os.path.join(os.path.dirname(__file__), "validation_samples.txt")


# ============================================================
# HELPERS
# ============================================================

def classify_true_state(nMajor, nMinor):
    if pd.isna(nMajor) or pd.isna(nMinor):
        return "Unknown"
    a, b = int(nMajor), int(nMinor)
    if a == 0 and b == 0:
        return "Deletion"
    if b == 0:
        return "LoH"
    if a == b:
        return "Balanced"
    return "AI"


def load_merged(path):
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"Merged CSV not found: {path}\n"
            "Run baf_to_hrd_pipeline.py (stages 1 + 2) first."
        )
    print(f"Loading merged CSV: {path}")
    df = pd.read_csv(path, low_memory=False)
    print(f"  {len(df):,} rows, {len(df.columns)} columns")
    return df


def merge_with_ascat(df, ascat_file):
    print(f"Loading ASCAT segment file: {ascat_file}")
    ascat = pd.read_csv(ascat_file)

    rename = {"SampleID": "sample", "Chr": "chromosome", "chr": "chromosome",
              "startpos": "start", "endpos": "end"}
    ascat.rename(columns=rename, inplace=True)

    ascat["chromosome"] = (
        ascat["chromosome"].astype(str)
        .str.replace("^chr", "", regex=True)
        .str.replace(r"\.0$", "", regex=True)
        .str.strip()
    )

    for col in ["chromosome", "start", "end"]:
        if col not in ascat.columns:
            raise ValueError(f"ASCAT file missing column: {col}")

    ascat = ascat[["sample", "chromosome", "start", "end", "nMajor", "nMinor"]].copy()
    ascat.rename(columns={"nMajor": "nMajor_ascat", "nMinor": "nMinor_ascat"}, inplace=True)

    df["chromosome"] = df["chromosome"].astype(str).str.strip()
    ascat["chromosome"] = ascat["chromosome"].astype(str).str.strip()
    df["start"] = pd.to_numeric(df["start"], errors="coerce").astype("Int64")
    df["end"]   = pd.to_numeric(df["end"],   errors="coerce").astype("Int64")
    ascat["start"] = pd.to_numeric(ascat["start"], errors="coerce").astype("Int64")
    ascat["end"]   = pd.to_numeric(ascat["end"],   errors="coerce").astype("Int64")

    merged = df.merge(ascat, on=["sample", "chromosome", "start", "end"], how="left")

    # Prefer ascat columns, fall back to existing nMajor/nMinor
    for col in ["nMajor", "nMinor"]:
        ascat_col = f"{col}_ascat"
        if ascat_col in merged.columns:
            merged[col] = merged[ascat_col].combine_first(merged.get(col, pd.Series(dtype=float)))
            merged.drop(columns=[ascat_col], inplace=True)

    n_matched = merged["nMajor"].notna().sum()
    print(f"  Matched {n_matched:,} / {len(merged):,} segments with ASCAT ground truth")
    return merged


# ============================================================
# MAIN
# ============================================================

def main():
    print("=" * 60)
    print("PREPARING VALIDATION DATASET")
    print("=" * 60)

    df = load_merged(MERGED_CSV)

    if ASCAT_SEGMENTS_FILE is not None:
        df = merge_with_ascat(df, ASCAT_SEGMENTS_FILE)

    # Filter to the chosen DP filter
    if "dp_filter" in df.columns:
        df = df[df["dp_filter"] == DP_FILTER_FOR_VAL].copy()
        print(f"\nFiltered to DP={DP_FILTER_FOR_VAL}: {len(df):,} rows")

    # Keep only segments with ASCAT ground truth
    has_gt = df["nMajor"].notna() & df["nMinor"].notna()
    df_val = df[has_gt].copy()
    print(f"Segments with ASCAT ground truth: {len(df_val):,} / {len(df):,}")

    # Optional sample filter for honest train/val split
    if VALIDATION_SAMPLES_FILE is not None:
        if not os.path.exists(VALIDATION_SAMPLES_FILE):
            raise FileNotFoundError(f"VALIDATION_SAMPLES_FILE not found: {VALIDATION_SAMPLES_FILE}")
        with open(VALIDATION_SAMPLES_FILE) as _f:
            val_samples = set(l.strip() for l in _f if l.strip() and not l.startswith('#'))
        _sample_col = next((c for c in ('sample', 'SampleID', 'sample_name') if c in df_val.columns), None)
        if _sample_col is None:
            raise ValueError("Cannot find a sample ID column in the merged CSV.")
        _before = len(df_val)
        df_val = df_val[df_val[_sample_col].isin(val_samples)].copy()
        print(f"Filtered to {len(val_samples):,} validation samples: "
              f"{len(df_val):,} / {_before:,} segments retained")

    if len(df_val) == 0:
        print("\nERROR: No ground-truth segments found.")
        print("Make sure the merged CSV was built from ASCAT-format samples")
        print("(segment_format == 'ascat') or set ASCAT_SEGMENTS_FILE.")
        sys.exit(1)

    # Filter segments with enough BAF data
    if "n_baf_points" in df_val.columns:
        before = len(df_val)
        df_val = df_val[df_val["n_baf_points"] >= MIN_BAF_POSITIONS].copy()
        print(f"After requiring ≥{MIN_BAF_POSITIONS} BAF positions: {len(df_val):,} (dropped {before - len(df_val):,})")
    elif "n_positions" in df_val.columns:
        before = len(df_val)
        df_val = df_val[df_val["n_positions"] >= MIN_BAF_POSITIONS].copy()
        print(f"After requiring ≥{MIN_BAF_POSITIONS} BAF positions: {len(df_val):,} (dropped {before - len(df_val):,})")

    # Round ground-truth CN
    df_val["nMajor_true"] = pd.to_numeric(df_val["nMajor"], errors="coerce").round().astype("Int64")
    df_val["nMinor_true"] = pd.to_numeric(df_val["nMinor"], errors="coerce").round().astype("Int64")
    df_val["class_true"]  = df_val.apply(
        lambda r: classify_true_state(r["nMajor_true"], r["nMinor_true"]), axis=1
    )

    # ── Legacy 20260406 column suffix cleanup ────────────────────────────────
    # The 20260406 merged CSV was built by merging a segment file (with 4-bin outer_mass)
    # against freshly recomputed 20-bin metrics. Overlapping column names received _x/_y.
    #
    # Convention used here:
    #   outer_mass / inner_mass / outer_mass_raw / inner_mass_raw
    #     → keep _y  (20-bin recomputed in Stage 2; _x is the obsolete 4-bin version)
    #   ratio_* / count_* / baf_iqr and all other _x/_y pairs
    #     → keep _x  (Stage 1 pre-computed; both are equivalent 20-bin values)
    prefer_y = {"outer_mass", "inner_mass", "outer_mass_raw", "inner_mass_raw"}
    for col in list(df_val.columns):
        if col.endswith("_x"):
            base = col[:-2]
            y_col = base + "_y"
            if y_col in df_val.columns:
                if base in prefer_y:
                    # Keep _y (20-bin), drop _x (4-bin)
                    df_val = df_val.drop(columns=[col]).rename(columns={y_col: base})
                else:
                    # Keep _x, drop _y
                    df_val = df_val.drop(columns=[y_col]).rename(columns={col: base})
    print(f"  After suffix cleanup: {len(df_val.columns)} columns")

    # Identify 20 bin ratio columns AFTER suffix cleanup so names are final
    bin_cols_corr = sorted([c for c in df_val.columns if c.startswith("ratio_") and "raw" not in c])
    bin_cols_raw  = sorted([c for c in df_val.columns if c.startswith("ratio_raw_")])
    print(f"  BAF bin columns: {len(bin_cols_corr)} corrected + {len(bin_cols_raw)} raw")

    # Core columns to keep
    keep_cols = [
        "sample", "sample_name", "segment_id", "dp_filter",
        "chromosome", "start", "end", "total_cn",
        "outer_mass_raw", "inner_mass_raw",   # uncorrected — calibration/application scale
        "outer_mass", "inner_mass",            # purity-corrected — kept for reference
        "baf_mean", "baf_sd", "baf_median", "baf_iqr",
        "nMajor_true", "nMinor_true", "class_true",
    ]
    optional_cols = ["n_baf_points", "n_positions", "segment_format", "purity_used"]
    existing_core     = [c for c in keep_cols     if c in df_val.columns]
    existing_optional = [c for c in optional_cols if c in df_val.columns]

    final_cols = existing_core + existing_optional + bin_cols_corr + bin_cols_raw
    df_out = df_val[final_cols].copy()

    df_out.to_csv(OUTPUT_FILE, index=False)
    print(f"\n✓ Saved: {OUTPUT_FILE}")
    print(f"  Rows: {len(df_out):,}")
    print(f"  Unique samples: {df_out['sample'].nunique() if 'sample' in df_out.columns else 'N/A'}")

    print("\nGround-truth class distribution:")
    for cls, count in df_out["class_true"].value_counts().items():
        print(f"  {cls:9s}: {count:6,} ({100*count/len(df_out):5.1f}%)")

    print("\nTotal CN distribution:")
    df_out["total_cn_int"] = pd.to_numeric(df_out["total_cn"], errors="coerce").round().astype("Int64")
    for cn, count in df_out["total_cn_int"].value_counts().sort_index().items():
        print(f"  CN={cn}: {count:6,} ({100*count/len(df_out):5.1f}%)")

    print("\n" + "=" * 60)
    print("DONE. Run experiment.py to start experimenting.")
    print("=" * 60)


if __name__ == "__main__":
    main()
