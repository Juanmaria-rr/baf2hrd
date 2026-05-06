#!/usr/bin/env python3
"""
split_segments.py — Extract per-sample segment CSVs from the combined segtable.

Reads the combined TSO500 segment file and writes one CSV per sample to an
output directory named YYYYMMDD_segments under tso_samples_root.

Usage:
    python split_segments.py --samples JBLAB001 JBLAB002 ...
    python split_segments.py --all
    python split_segments.py --samples JBLAB001 --force

Output files:
    <outdir>/<SAMPLE>_100kb_abs_segtable_free_purity_filtered_unique.csv

Exit codes:
    0 — all requested samples written (or already exist with --no-force)
    1 — one or more requested samples not found in the combined file
"""

import argparse
import os
import sys
from datetime import date

import pandas as pd

COMBINED_FILE = (
    "/storage/scratch01/groups/co/cn_extra/alleleSpecific/tso_samples/"
    "unique_copyright_100kb_abs_segtable_free_purity_filtered_unique.csv"
)
TSO_SAMPLES_ROOT = "/storage/scratch01/groups/co/cn_extra/alleleSpecific/tso_samples"
SUFFIX = "_100kb_abs_segtable_free_purity_filtered_unique.csv"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("--samples", nargs="+", metavar="SAMPLE", help="Sample IDs to extract.")
    group.add_argument("--all", action="store_true", help="Extract all samples in the combined file.")
    p.add_argument("--combined", default=COMBINED_FILE, metavar="FILE",
                   help="Path to combined segment CSV (default: %(default)s).")
    p.add_argument("--outdir", default=None, metavar="DIR",
                   help="Output directory. Default: <tso_samples_root>/YYYYMMDD_segments/")
    p.add_argument("--force", action="store_true",
                   help="Overwrite existing per-sample CSVs.")
    return p.parse_args()


def main():
    args = parse_args()

    outdir = args.outdir or os.path.join(TSO_SAMPLES_ROOT, f"{date.today().strftime('%Y%m%d')}_segments")
    os.makedirs(outdir, exist_ok=True)

    print(f"[split_segments] Combined file : {args.combined}")
    print(f"[split_segments] Output dir    : {outdir}")

    df = pd.read_csv(args.combined)
    if "sample" not in df.columns:
        sys.exit(f"[ERROR] Column 'sample' not found in {args.combined}. Columns: {df.columns.tolist()}")

    available = set(df["sample"].unique())
    requested = available if args.all else set(args.samples)

    missing = requested - available
    if missing:
        print(f"[ERROR] {len(missing)} sample(s) not found in combined file:", file=sys.stderr)
        for s in sorted(missing):
            print(f"  - {s}", file=sys.stderr)

    found = requested & available
    print(f"[split_segments] Samples requested: {len(requested)}  |  found: {len(found)}  |  missing: {len(missing)}")

    n_written = 0
    n_skipped = 0
    for sample in sorted(found):
        out_path = os.path.join(outdir, f"{sample}{SUFFIX}")
        if os.path.exists(out_path) and not args.force:
            print(f"  [SKIP]  {sample}  (already exists — use --force to overwrite)")
            n_skipped += 1
            continue
        sample_df = df[df["sample"] == sample].copy()
        sample_df.to_csv(out_path, index=False)
        print(f"  [OK]    {sample}  ({len(sample_df)} segments → {out_path})")
        n_written += 1

    print(f"[split_segments] Done. Written: {n_written}  Skipped: {n_skipped}  Missing: {len(missing)}")
    print(f"[split_segments] raw_segments_dir={outdir}")

    if missing:
        sys.exit(1)


if __name__ == "__main__":
    main()
