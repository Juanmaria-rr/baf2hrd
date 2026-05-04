#!/usr/bin/env python3
"""
consolidate_results_parallelization.py
Consolidate results from parallel jobs and generate scarHRD inputs
"""

import pandas as pd
import glob
import os
import numpy as np

DATE = '20260316'
output_dir = f"/storage/scratch01/groups/co/cn_extra/alleleSpecific/tso_samples/{DATE}"

print("="*70)
print("CONSOLIDATING PARALLEL RESULTS")
print("="*70)

# Find all segment files
segment_files = sorted(glob.glob(f"{output_dir}/segments_per_sample/segments_*.csv"))
print(f"Found {len(segment_files)} segment files")

if len(segment_files) == 0:
    print("ERROR: No segment files found!")
    exit(1)

# Load all
all_segments = []
for i, file in enumerate(segment_files, 1):
    if i % 100 == 0:
        print(f"  Loading {i}/{len(segment_files)}...")
    df = pd.read_csv(file)
    all_segments.append(df)

all_segments_df = pd.concat(all_segments, ignore_index=True)

# Save consolidated
consolidated_file = f"{output_dir}/all_segments_consolidated.csv"
all_segments_df.to_csv(consolidated_file, index=False)

print(f"\n✓ Consolidated:")
print(f"  Total segments: {len(all_segments_df):,}")
print(f"  Unique samples: {all_segments_df['sample_name'].nunique()}")
print(f"  Saved: {consolidated_file}")

# Now run scarHRD generation...
# (Copy the scarHRD section from your original script)

print("\n✓ CONSOLIDATION COMPLETE!")