#!/usr/bin/env python3
"""
compute_purity_corrected_outer_mass.py

Analytical purity correction for BAF outer_mass calibration.

Takes the merged segments calibration CSV (which already contains per-segment
raw BAF histogram bins and sample purity) and adds a new column:

    outer_mass_purity_corrected

Computed by applying the pipeline purity correction formula
    baf_corr = 0.5 ± (baf - 0.5) / purity
to each raw histogram bin center, then summing the ratios that map to outer
bins (BAF ≤ 0.25 or ≥ 0.75) after correction.

Context: The per-sample ASCAT CN segment files that would be needed to re-run
Stage 1 from scratch are no longer available individually. The calibration CSV
already contains all 20 raw histogram bins per segment, so the purity-corrected
outer_mass can be derived analytically without re-running Stage 1–2.
"""

import argparse
import os

import numpy as np
import pandas as pd

# ──────────────────────────────────────────────────────────────────────────────
# Bin structure — must match baf_to_hrd_pipeline.py
# ──────────────────────────────────────────────────────────────────────────────
N_BINS = 20
BIN_EDGES = np.linspace(0, 1, N_BINS + 1)
BIN_LABELS = [f"{BIN_EDGES[i]:.2f}-{BIN_EDGES[i+1]:.2f}" for i in range(N_BINS)]
BIN_CENTERS = (BIN_EDGES[:-1] + BIN_EDGES[1:]) / 2   # 0.025, 0.075, …, 0.975
BIN_WIDTH = float(BIN_EDGES[1] - BIN_EDGES[0])

# Outer bins: upper edge ≤ 0.25  OR  lower edge ≥ 0.75
# Indices 0-4 (0.00-0.25) and 15-19 (0.75-1.00)
OUTER_MASK = np.array(
    [BIN_EDGES[i + 1] <= 0.25 + 1e-9 or BIN_EDGES[i] >= 0.75 - 1e-9
     for i in range(N_BINS)],
    dtype=bool,
)

# Raw ratio columns present in the calibration CSV (JBLAB17041 source)
RAW_RATIO_COLS = [f"ratio_raw_{lbl}_x" for lbl in BIN_LABELS]


def compute_corrected_outer_mass(
    raw_ratios: np.ndarray, purities: np.ndarray
) -> np.ndarray:
    """
    Vectorised purity-corrected outer_mass from raw histogram bins.

    Parameters
    ----------
    raw_ratios : (N, 20) float64 — per-segment raw BAF histogram ratios
    purities   : (N,)   float64 — per-sample tumour purity in (0, 1]

    Returns
    -------
    outer_mass_corr : (N,) float64
    """
    # Invalid purity → no correction (treat as pure, purity=1)
    p = np.where((purities > 0) & (purities <= 1), purities, 1.0)

    c  = BIN_CENTERS[np.newaxis, :]   # (1, 20)
    p2 = p[:, np.newaxis]             # (N,  1)

    # Shift BAF values toward 0 and 1
    c_corr = np.where(
        c >= 0.5,
        0.5 + (c - 0.5) / p2,
        0.5 - (0.5 - c) / p2,
    )
    # Clip to [0, 1) so floor mapping stays within bounds
    c_corr = np.clip(c_corr, 0.0, 1.0 - 1e-12)   # (N, 20)

    # Index of the corrected bin each original bin maps to
    corr_idx = np.clip(
        np.floor(c_corr / BIN_WIDTH).astype(np.int32), 0, N_BINS - 1
    )   # (N, 20)

    # Sum raw ratios whose corrected bin falls in the outer region
    return (raw_ratios * OUTER_MASK[corr_idx]).sum(axis=1)   # (N,)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input", required=True,
        help="Merged calibration CSV with ratio_raw_*_x columns and purity column",
    )
    parser.add_argument(
        "--output-dir", required=True,
        help="Directory where the enriched CSV will be written",
    )
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    print("=" * 70)
    print("PURITY-CORRECTED outer_mass — analytical computation")
    print("=" * 70)
    print(f"Input  : {args.input}")
    print(f"Outdir : {args.output_dir}")

    print("\nLoading CSV…")
    df = pd.read_csv(args.input)
    print(f"  {len(df):,} rows  ×  {len(df.columns)} columns")

    # Validate required columns
    missing = [c for c in RAW_RATIO_COLS if c not in df.columns]
    if missing:
        raise ValueError(
            f"Missing raw ratio columns (first 3): {missing[:3]}\n"
            "Expected columns like 'ratio_raw_0.00-0.05_x' in the calibration CSV."
        )
    if "purity" not in df.columns:
        raise ValueError("'purity' column not found in the input CSV.")

    raw_ratios = df[RAW_RATIO_COLS].fillna(0.0).to_numpy(dtype=np.float64)
    purities   = df["purity"].to_numpy(dtype=np.float64)

    n_invalid = int(((purities <= 0) | (purities > 1) | np.isnan(purities)).sum())
    print(f"\n  Purity range  : [{np.nanmin(purities):.3f}, {np.nanmax(purities):.3f}]")
    print(f"  Invalid purity: {n_invalid:,} rows (will use raw outer_mass = no correction)")

    print("\nComputing outer_mass_purity_corrected…")
    df["outer_mass_purity_corrected"] = compute_corrected_outer_mass(raw_ratios, purities)

    # Sanity check against the raw value
    if "outer_mass_raw_x" in df.columns:
        delta = df["outer_mass_purity_corrected"] - df["outer_mass_raw_x"]
        print(f"  Delta vs outer_mass_raw_x — mean: {delta.mean():.4f}, "
              f"abs mean: {delta.abs().mean():.4f}")
        print(f"  Fraction increased: {(delta > 0).mean():.1%}")

    out_path = os.path.join(
        args.output_dir,
        "merged_segments_baf_metrics_ALL_purity_corrected.csv",
    )
    print(f"\nWriting {out_path}…")
    df.to_csv(out_path, index=False)
    print(f"  {len(df):,} rows written.")
    print("\nDone.")


if __name__ == "__main__":
    main()
