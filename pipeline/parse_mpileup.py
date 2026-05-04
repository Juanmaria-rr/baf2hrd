#!/usr/bin/env python3
"""
parse_mpileup.py — Convert raw bcftools-mpileup TSV to pipeline-ready _qc.tsv

Stage 0b in the BAF-to-HRD pipeline:
  Input:  {sample}_allele_counts_snp6_wrapper_qualityMetrics_q30Pass.tsv
          (no header; columns produced by bcftools query:
           CHROM  POS  REF  ALT  DP  AD  ADF  ADR  SP  MQBZ  BQBZ)

  Output: {sample}_qc.tsv
          (header: CHROM  POS  REF  ALT  DP  Ref_Reads  Alt_Reads  BAF
                   GT_inferred  problematic  problematic_reason)

The output schema is identical to what baf_to_hrd_pipeline.py expects for
TSO500 samples.  normalize_baf_df() will recognise the 'problematic_reason'
column and apply the _ok filter automatically.

Logic extracted from phasingTSO.ipynb cell 24 and cell 17.

Usage (standalone):
    python parse_mpileup.py <input_tsv> [output_qc_tsv]

    If output_qc_tsv is omitted the output is written next to the input:
        {dir}/{sample}_qc.tsv
    where {sample} is derived by stripping the long suffix from the filename.

Usage (from orchestrate_pipeline.sh):
    python parse_mpileup.py \\
        /path/to/SAMPLE_allele_counts_snp6_wrapper_qualityMetrics_q30Pass.tsv \\
        /path/to/out/SAMPLE_qc.tsv
"""

import sys
import numpy as np
import pandas as pd
from pathlib import Path

RAW_SUFFIX = "_allele_counts_snp6_wrapper_qualityMetrics_q30Pass.tsv"
RAW_COLS   = ["CHROM", "POS", "REF", "ALT", "DP", "AD", "ADF", "ADR", "SP", "MQBZ", "BQBZ"]


# ── helpers ──────────────────────────────────────────────────────────────────

def _parse_ad(x):
    """Return (ref_count, alt_count) from bcftools AD string e.g. '12,3' or '12,3,0'."""
    parts = str(x).split(",")
    ref = int(parts[0]) if len(parts) > 0 and parts[0].isdigit() else 0
    alt = int(parts[1]) if len(parts) > 1 and parts[1].isdigit() else 0
    return ref, alt


def _parse_allele_count(field, allele_idx=1):
    """Return the count for a given allele index from a comma-separated AD-style field."""
    parts = str(field).split(",")
    if len(parts) > allele_idx and parts[allele_idx].isdigit():
        return int(parts[allele_idx])
    return 0


def _infer_gt(baf):
    if pd.isna(baf):
        return "./."
    if baf < 0.01:
        return "0/0"
    if baf > 0.99:
        return "1/1"
    return "0/1"


def _normalise_alt(alt_str):
    """Replace reference-only placeholder <*> with N; keep first allele for multi-allelics."""
    s = str(alt_str) if pd.notna(alt_str) else "N"
    if s == "<*>":
        return "N"
    # keep only first allele (e.g. "G,<*>" → "G")
    return s.split(",")[0]


# ── main conversion ───────────────────────────────────────────────────────────

def convert(input_path: Path, output_path: Path) -> None:
    print(f"[parse_mpileup] Reading  : {input_path}")

    df = pd.read_csv(
        input_path,
        sep="\t",
        names=RAW_COLS,
        dtype={"CHROM": "string", "REF": "string", "ALT": "string"},
        na_values=["NA"],       # do NOT treat '.' as NaN — it is a valid ALT value
    )
    print(f"[parse_mpileup] Loaded   : {len(df):,} raw positions")

    # ── 1. Parse AD → Ref_Reads / Alt_Reads ──────────────────────────────────
    ad_parsed       = df["AD"].apply(_parse_ad)
    df["Ref_Reads"] = [x[0] for x in ad_parsed]
    df["Alt_Reads"] = [x[1] for x in ad_parsed]

    # DP: use the explicit DP field; fall back to Ref+Alt sum
    df["DP"] = (
        pd.to_numeric(df["DP"], errors="coerce")
        .fillna(df["Ref_Reads"] + df["Alt_Reads"])
        .astype(int)
    )

    # ── 2. Normalise ALT ─────────────────────────────────────────────────────
    df["ALT"] = df["ALT"].apply(_normalise_alt)

    # ── 3. BAF ───────────────────────────────────────────────────────────────
    total    = df["Ref_Reads"] + df["Alt_Reads"]
    df["BAF"] = np.where(total > 0, df["Alt_Reads"] / total, 0.0)

    # ── 4. GT_inferred ───────────────────────────────────────────────────────
    df["GT_inferred"] = df["BAF"].apply(_infer_gt)

    # ── 5. Strand-balance QC → problematic_reason ────────────────────────────
    # Flag ALT positions that have zero forward OR zero reverse strand support
    # (classic strand-bias artefact).  REF-only sites (Alt_Reads == 0) are clean.
    df["_alt_fwd"] = df["ADF"].apply(lambda x: _parse_allele_count(x, 1))
    df["_alt_rev"] = df["ADR"].apply(lambda x: _parse_allele_count(x, 1))

    def _problematic_reason(row):
        if row["Alt_Reads"] > 0:
            if row["_alt_fwd"] == 0 or row["_alt_rev"] == 0:
                return "strand_bias"
        return None

    df["problematic_reason"] = df.apply(_problematic_reason, axis=1)
    df["problematic"]        = df["problematic_reason"].notna().astype(int)

    # ── 6. Select output columns ─────────────────────────────────────────────
    out_cols = [
        "CHROM", "POS", "REF", "ALT", "DP",
        "Ref_Reads", "Alt_Reads", "BAF",
        "GT_inferred", "problematic", "problematic_reason",
    ]
    out = df[out_cols].copy()

    n_problematic = df["problematic"].sum()
    print(f"[parse_mpileup] Positions: {len(out):,}  "
          f"(flagged strand-bias: {n_problematic:,})")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output_path, sep="\t", index=False)
    print(f"[parse_mpileup] Written  : {output_path}")


# ── CLI ───────────────────────────────────────────────────────────────────────

def _derive_output(input_path: Path) -> Path:
    """Derive {dir}/{sample}_qc.tsv from the long input filename."""
    name = input_path.name
    if name.endswith(RAW_SUFFIX):
        sample = name[: -len(RAW_SUFFIX)]
    else:
        # fallback: strip .tsv and append _qc.tsv
        sample = name.removesuffix(".tsv")
    return input_path.parent / f"{sample}_qc.tsv"


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    inp = Path(sys.argv[1])
    out = Path(sys.argv[2]) if len(sys.argv) >= 3 else _derive_output(inp)

    if not inp.exists():
        print(f"[ERROR] Input file not found: {inp}", file=sys.stderr)
        sys.exit(1)

    convert(inp, out)
    print("[parse_mpileup] Done.")
