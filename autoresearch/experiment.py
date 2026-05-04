#!/usr/bin/env python3
"""
experiment.py — THE file agents modify.

DO NOT modify: prepare.py, evaluate.py, baf_to_hrd_pipeline.py

──────────────────────────────────────────────────────────────
Goal
──────────────────────────────────────────────────────────────
Maximise val_allelic_acc — the fraction of segments where the predicted
(A_cn, B_cn) exactly matches the ASCAT ground truth.

Secondary metric: val_class_acc — fraction with the correct allelic state
class (LoH / Balanced / AI / Deletion), regardless of exact copy numbers.

──────────────────────────────────────────────────────────────
On this HPC cluster: submit via sbatch, do NOT run directly
──────────────────────────────────────────────────────────────
  sbatch autoresearch_baf_hrd/run_experiment.sbatch

Results are appended to results/log.jsonl after each run.

──────────────────────────────────────────────────────────────
What you can change
──────────────────────────────────────────────────────────────
  1. THRESHOLDS          — outer_mass_raw decision boundaries per total CN
  2. OUTER_LOWER / OUTER_UPPER  — the BAF interval that defines outer_mass_raw
                                   (recomputed from the 20 ratio_raw_* fine bins)
  3. HEURISTIC_THRESHOLD — threshold for CN=7 and CN>8
  4. classify_segment()  — the function itself; use any columns available
                           in val_data.csv (bin ratios, baf_mean, baf_sd, etc.)

──────────────────────────────────────────────────────────────
Key signal: outer_mass_raw vs outer_mass
──────────────────────────────────────────────────────────────
IMPORTANT: this experiment uses outer_mass_raw (uncorrected BAF tails), NOT
the purity-corrected outer_mass.  This matches the scale used by Stage 3 of
the production pipeline for TSO500 samples, which also has no purity model.

  outer_mass_raw = computed from ratio_raw_* bins (raw BAF, no correction)
  outer_mass     = computed from ratio_* bins (purity-corrected BAF, ASCAT only)

Both are available in val_data.csv; recompute_outer_mass() always uses raw bins.

──────────────────────────────────────────────────────────────
Available columns in val_data.csv
──────────────────────────────────────────────────────────────
  ratio_raw_0.00-0.05 … ratio_raw_0.95-1.00   raw BAF bin ratios  (20 bins)
  ratio_0.00-0.05 … ratio_0.95-1.00           corrected BAF bin ratios (for reference)
  outer_mass_raw, inner_mass_raw               pre-computed from raw bins (boundary 0.25/0.75)
  outer_mass, inner_mass                       purity-corrected versions (ASCAT scale, reference)
  baf_mean, baf_sd, baf_median, baf_iqr
  total_cn, nMajor_true, nMinor_true, class_true
"""

import os
import sys
import json
import time
import numpy as np
import pandas as pd
from datetime import datetime

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "pipeline"))
from baf_to_hrd_pipeline import classify_final_state_from_acn_bcn

VAL_DATA  = os.path.join(os.path.dirname(__file__), "val_data.csv")
RESULTS_LOG = os.path.join(os.path.dirname(__file__), "results", "log.jsonl")

# ============================================================
# PARAMETERS — modify these
# ============================================================

THRESHOLDS = {
    # Calibrated per-DP=8 Youden's-index thresholds from calibrate_thresholds.py
    # (re-run 2026-04-29 on both JBLAB231+JBLAB17041 backbones, 201,218 segments)
    "2": {"loh_vs_balanced": 0.8537},
    "3": {"loh_vs_ai":       0.8824},
    "4": {"loh_vs_balanced": 0.8261, "loh_vs_ai": 0.8571, "ai_vs_balanced": 0.8333},
    "5": {"loh_vs_ai":       0.8444},
    "6": {"loh_vs_balanced": 1.0000, "loh_vs_ai": 1.0000, "ai_vs_balanced": 0.8400},
    "8": {"ai_vs_balanced":  0.8276},
}

# Boundary for recomputing outer_mass from the 20 fine bins.
# outer_mass = sum of ratio_{bin} for bins fully below OUTER_LOWER
#              or fully above OUTER_UPPER.
# Pipeline default: 0.25 / 0.75
OUTER_LOWER = 0.25
OUTER_UPPER = 0.75

# Heuristic fallback for CN=7 and CN>8
HEURISTIC_THRESHOLD = 0.95

# Short description of what this experiment changes (shown in log)
EXPERIMENT_NOTES = "baseline — calibrated DP=8 thresholds"


# ============================================================
# outer_mass recomputation from fine bins — DO NOT MODIFY
# (Change OUTER_LOWER / OUTER_UPPER above instead)
# ============================================================

def recompute_outer_mass(df, lower=OUTER_LOWER, upper=OUTER_UPPER):
    """
    Recompute outer_mass_raw from the 20 raw BAF bin ratios using custom
    boundary values.  Bins are 0.05-wide: 0.00-0.05, 0.05-0.10, …
    A bin is "outer" if its upper edge <= lower or its lower edge >= upper.

    Uses ratio_raw_* columns (uncorrected BAF) to match the scale used by
    Stage 3 of the production pipeline for TSO500 samples.  Falls back to
    the pre-computed outer_mass_raw column when fine bins are absent.
    """
    bin_cols = sorted([c for c in df.columns if c.startswith("ratio_raw_")])
    if not bin_cols:
        # Fall back to the pre-computed column (same scale)
        return df["outer_mass_raw"].values

    outer_sum = pd.Series(np.zeros(len(df)), index=df.index)
    for col in bin_cols:
        # Parse boundaries from column name e.g. "ratio_raw_0.00-0.05"
        label = col.replace("ratio_raw_", "")
        try:
            lo, hi = map(float, label.split("-"))
        except ValueError:
            continue
        if hi <= lower or lo >= upper:
            vals = pd.to_numeric(df[col], errors="coerce").fillna(0)
            outer_sum = outer_sum + vals

    return outer_sum.values


# ============================================================
# CLASSIFICATION FUNCTION — agents rewrite this
# ============================================================

def classify_segment(total_cn, outer_mass):
    """
    Classify a single segment into (A_cn, B_cn).

    Parameters
    ----------
    total_cn   : int
    outer_mass : float or NaN  — recomputed with OUTER_LOWER / OUTER_UPPER

    Returns
    -------
    (A_cn, B_cn) : tuple of int
    """
    if pd.isna(total_cn):
        return (pd.NA, pd.NA)

    total_cn = int(total_cn)

    if total_cn == 0:
        return (0, 0)

    if total_cn == 1:
        return (1, 0)

    if total_cn == 2:
        t = THRESHOLDS["2"]["loh_vs_balanced"]
        if pd.isna(outer_mass):
            return (1, 1)
        return (2, 0) if outer_mass > t else (1, 1)

    if total_cn == 3:
        t = THRESHOLDS["3"]["loh_vs_ai"]
        if pd.isna(outer_mass):
            return (2, 1)
        return (3, 0) if outer_mass > t else (2, 1)

    if total_cn == 4:
        t_loh = THRESHOLDS["4"]["loh_vs_balanced"]
        t_ai  = THRESHOLDS["4"]["ai_vs_balanced"]
        if pd.isna(outer_mass):
            return (2, 2)
        if outer_mass > t_loh:
            return (4, 0)
        if outer_mass > t_ai:
            return (3, 1)
        return (2, 2)

    if total_cn == 5:
        t = THRESHOLDS["5"]["loh_vs_ai"]
        if pd.isna(outer_mass):
            return (3, 2)
        return (5, 0) if outer_mass > t else (3, 2)

    if total_cn == 6:
        t = THRESHOLDS["6"]["ai_vs_balanced"]
        if pd.isna(outer_mass):
            return (4, 2)
        return (4, 2) if outer_mass > t else (3, 3)

    if total_cn == 7:
        t = HEURISTIC_THRESHOLD
        if pd.isna(outer_mass):
            return (4, 3)
        return (7, 0) if outer_mass > t else (4, 3)

    if total_cn == 8:
        t = THRESHOLDS["8"]["ai_vs_balanced"]
        if pd.isna(outer_mass):
            return (6, 2)
        return (6, 2) if outer_mass > t else (4, 4)

    # CN > 8: heuristic
    if pd.isna(outer_mass):
        a = (total_cn + 1) // 2 if total_cn % 2 == 1 else total_cn // 2 + 1
        b = (total_cn - 1) // 2 if total_cn % 2 == 1 else total_cn // 2 - 1
        return (a, b)
    if outer_mass > HEURISTIC_THRESHOLD:
        return (total_cn, 0)
    a = (total_cn + 1) // 2 if total_cn % 2 == 1 else total_cn // 2 + 1
    b = (total_cn - 1) // 2 if total_cn % 2 == 1 else total_cn // 2 - 1
    return (a, b)


# ============================================================
# EVALUATION — DO NOT MODIFY BELOW THIS LINE
# ============================================================

def run_evaluation():
    if not os.path.exists(VAL_DATA):
        raise FileNotFoundError(
            f"val_data.csv not found: {VAL_DATA}\nRun prepare.py first."
        )

    t0 = time.time()
    df = pd.read_csv(VAL_DATA, low_memory=False)
    print(f"Loaded val_data.csv: {len(df):,} segments")

    # Recompute outer_mass with current OUTER_LOWER / OUTER_UPPER
    df["outer_mass_exp"] = recompute_outer_mass(df, lower=OUTER_LOWER, upper=OUTER_UPPER)

    # Clamp total_cn to integer
    df["total_cn_int"] = (
        pd.to_numeric(df["total_cn"], errors="coerce")
        .apply(lambda x: int(np.floor(x + 0.5)) if pd.notna(x) else pd.NA)
    )

    # Run classification
    predictions = df.apply(
        lambda r: classify_segment(r["total_cn_int"], r["outer_mass_exp"]), axis=1
    )
    df["A_cn_pred"] = [p[0] for p in predictions]
    df["B_cn_pred"] = [p[1] for p in predictions]

    df["A_cn_true"] = pd.to_numeric(df["nMajor_true"], errors="coerce")
    df["B_cn_true"] = pd.to_numeric(df["nMinor_true"], errors="coerce")

    # Allelic match: exact (A, B) or flipped (B, A)
    exact   = (df["A_cn_pred"] == df["A_cn_true"]) & (df["B_cn_pred"] == df["B_cn_true"])
    flipped = (df["A_cn_pred"] == df["B_cn_true"]) & (df["B_cn_pred"] == df["A_cn_true"])
    df["allelic_match"] = exact | flipped

    # Class match: LoH / Balanced / AI / Deletion
    df["class_pred"] = df.apply(
        lambda r: classify_final_state_from_acn_bcn(r["A_cn_pred"], r["B_cn_pred"]), axis=1
    )
    df["class_match"] = df["class_pred"] == df["class_true"]

    n_total           = len(df)
    val_allelic_acc   = df["allelic_match"].mean()
    val_class_acc     = df["class_match"].mean()
    elapsed           = time.time() - t0

    print(f"\n{'='*50}")
    print(f"val_allelic_acc : {val_allelic_acc:.4f}  ({df['allelic_match'].sum():,} / {n_total:,})")
    print(f"val_class_acc   : {val_class_acc:.4f}  ({df['class_match'].sum():,} / {n_total:,})")
    print(f"{'='*50}")

    # Per-CN breakdown
    print("\nPer-CN allelic accuracy:")
    per_cn = {}
    for cn in sorted(df["total_cn_int"].dropna().unique()):
        mask   = df["total_cn_int"] == cn
        n      = mask.sum()
        acc    = df.loc[mask, "allelic_match"].mean()
        per_cn[int(cn)] = round(float(acc), 4)
        print(f"  CN={int(cn):2d}  n={n:5,}  acc={acc:.3f}")

    # Class breakdown
    print("\nPer-class accuracy:")
    per_class = {}
    for cls in ["Deletion", "LoH", "AI", "Balanced"]:
        mask = df["class_true"] == cls
        if mask.sum() == 0:
            continue
        acc = df.loc[mask, "allelic_match"].mean()
        per_class[cls] = round(float(acc), 4)
        print(f"  {cls:9s}  n={mask.sum():5,}  allelic_acc={acc:.3f}")

    # Outer_mass boundary used
    outer_boundary_info = {
        "OUTER_LOWER": OUTER_LOWER,
        "OUTER_UPPER": OUTER_UPPER,
    }

    # Append result to log
    log_entry = {
        "timestamp":      datetime.now().isoformat(),
        "val_allelic_acc": round(float(val_allelic_acc), 6),
        "val_class_acc":   round(float(val_class_acc),   6),
        "n_segments":      n_total,
        "elapsed_s":       round(elapsed, 2),
        "thresholds":      THRESHOLDS,
        "outer_boundary":  outer_boundary_info,
        "heuristic_threshold": HEURISTIC_THRESHOLD,
        "per_cn_allelic_acc":  per_cn,
        "per_class_allelic_acc": per_class,
        "notes":           EXPERIMENT_NOTES,
    }

    os.makedirs(os.path.dirname(RESULTS_LOG), exist_ok=True)
    with open(RESULTS_LOG, "a") as f:
        f.write(json.dumps(log_entry) + "\n")

    print(f"\n✓ Result appended to {RESULTS_LOG}")
    print(f"  val_allelic_acc = {val_allelic_acc:.4f}   (target: higher is better)")

    return val_allelic_acc, val_class_acc


if __name__ == "__main__":
    run_evaluation()
