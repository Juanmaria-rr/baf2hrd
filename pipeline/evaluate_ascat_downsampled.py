#!/usr/bin/env python3
"""
evaluate_ascat_downsampled.py

Correlates pipeline-computed HRD scores (from ASCAT-downsampled samples run
through the threshold classifier) against gold-standard ASCAT scarHRD scores.

Produces a FacetGrid scatter plot (correlation_ascat_downsampled.png) with one
panel per metric × DP filter × backbone, each annotated with n / Pearson r / p.

──────────────────────────────────────────────────────────────
SampleID formats (verified 2026-04-29 against actual CSV files)
──────────────────────────────────────────────────────────────

DOWNSAMPLED (Stage 4 on ASCAT-downsampled samples):
  {uuid1}_vs_{uuid2}_{backbone_suffix}_{dp_tag}
  e.g.  002484bc-72fc-11e8-b373-b4da36d49785_vs_275b947c-c09f-11e7-949e-ae70460a7042_17041_dp8

  Parsed with rsplit('_', 2):
    str[0] = {uuid1}_vs_{uuid2}          ← join key
    str[1] = {backbone_suffix}           ← 17041 or 231; backbone detected via InputFile
    str[2] = dp{N}  → strip 'dp' → N    ← dp_filter (0, 2, 4, 8)

GOLD STANDARD (ASCAT run on matched tumor-normal):
  {uuid1}_vs_{uuid2}_gs_scarHRD_input.txt
  e.g.  002484bc-72fc-11e8-b373-b4da36d49785_vs_275b947c-c09f-11e7-949e-ae70460a7042_gs_scarHRD_input.txt

  Join key = SampleID.replace('_gs_scarHRD_input.txt', '')
           = {uuid1}_vs_{uuid2}

Both join keys are the bare UUID pair, enabling a clean inner join on sample_name.
Gold standard has no dp filter column — each sample is one row.

──────────────────────────────────────────────────────────────
Inputs
──────────────────────────────────────────────────────────────
  DOWNSAMPLED_CSV   -- ALL_scarHRD_scores.csv from Stage 4 on ASCAT-downsampled samples
  GOLD_STANDARD_CSV -- ALL_scarHRD_scores.csv from the ASCAT gold-standard run

Output
──────────────────────────────────────────────────────────────
  correlation_ascat_downsampled.png  (saved to current working directory)
"""

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

# ============================================================
# PATHS
# ============================================================
DOWNSAMPLED_CSV = (
    "/storage/scratch01/groups/co/cn_extra/alleleSpecific/tso_samples"
    "/20260406/outputs_scarHRD_thresholds/ALL_scarHRD_scores.csv"
)
GOLD_STANDARD_CSV = (
    "/storage/scratch01/groups/co/cn_extra/alleleSpecific/tso_samples"
    "/20260310/goldStandard_outputs_scarHRD/ALL_scarHRD_scores.csv"
)

# ============================================================
# LOAD DOWNSAMPLED RESULTS
# ============================================================
downsampled = pd.read_csv(DOWNSAMPLED_CSV, sep=",")
df = downsampled.copy()

_sid_parts        = df['SampleID'].str.rsplit('_', 2)
# str[0] = {uuid1}_vs_{uuid2}  (join key)
# str[1] = backbone suffix (17041 / 231 etc.) — backbone detected via InputFile below
# str[2] = dp{N}  →  strip 'dp' prefix to get numeric depth filter
df['sample_name'] = _sid_parts.str[0]
df['dp_filter']   = _sid_parts.str[2].str.replace('dp', '', regex=False).astype(str)
df['backbone'] = np.where(
    df['InputFile'].str.contains('JBLAB231', na=False),
    'JBLAB231',
    np.where(
        df['InputFile'].str.contains('JBLAB17041', na=False),
        'JBLAB17041',
        'Unknown'
    )
)
df['LoH']       = df['HRD']
df['HRD_Total'] = df['HRD-sum']

downsamp = df.copy()

# ============================================================
# LOAD GOLD STANDARD
# ============================================================
goldStandard = pd.read_csv(GOLD_STANDARD_CSV, sep=",")
df = goldStandard.copy()


# Bug 3 fix: strip the scarHRD filename suffix to recover the UUID pair join key
# SampleID format: {uuid1}_vs_{uuid2}_gs_scarHRD_input.txt
df['sample_name'] = df['SampleID'].str.replace('_gs_scarHRD_input\\.txt$', '', regex=True)
df['backbone'] = np.where(
    df['InputFile'].str.contains('JBLAB231', na=False),
    'JBLAB231',
    np.where(
        df['InputFile'].str.contains('JBLAB17041', na=False),
        'JBLAB17041',
        'Unknown'
    )
)
# Bug 5 fix: gold standard has no dp filter — dp column removed
df['LoH']       = df['HRD']
df['HRD_Total'] = df['HRD-sum']

goldStand = df.copy()
goldStand.columns = ['gs_' + c for c in goldStand.columns]

# ============================================================
# MERGE
# ============================================================
merged = downsamp.merge(
    goldStand,
    left_on='sample_name', right_on='gs_sample_name',
    how='inner'
)
# Bug 4 fix: dp_filter already correctly extracted from downsampled SampleID above;
# the old double-strip (split('dp')[1] on an already-stripped string) produced NaN.
print(f"Merged rows: {len(merged)}")

# ============================================================
# PLOT
# ============================================================
metric_pairs = [
    ("LoH",          "gs_LoH"),
    ("Telomeric AI", "gs_Telomeric AI"),
    ("LST",          "gs_LST"),
    ("HRD_Total",    "gs_HRD_Total"),
]

df_plot = merged.drop_duplicates(subset=["SampleID", "dp_filter", "backbone"]).copy()

all_metric_cols = [col for pair in metric_pairs for col in pair]
df_plot[["dp_filter"] + all_metric_cols] = (
    df_plot[["dp_filter"] + all_metric_cols].apply(pd.to_numeric, errors="coerce")
)
df_plot = df_plot.dropna(subset=all_metric_cols + ["dp_filter"])

long_data = []
for calc_metric, gs_metric in metric_pairs:
    temp = df_plot[["SampleID", "dp_filter", "backbone", calc_metric, gs_metric]].copy()
    temp["Metric"]       = calc_metric.replace("_", " ")
    temp["Calculated"]   = temp[calc_metric]
    temp["GoldStandard"] = temp[gs_metric]
    long_data.append(temp[["SampleID", "dp_filter", "backbone", "Metric", "Calculated", "GoldStandard"]])

long = pd.concat(long_data, ignore_index=True).dropna(subset=["Calculated", "GoldStandard"])

dp_order       = sorted(long["dp_filter"].unique())
metric_order   = [p[0].replace("_", " ") for p in metric_pairs]
backbone_order = sorted(long["backbone"].unique())

long['row_var'] = 'DP≥' + long['dp_filter'].astype(str) + ' | ' + long['backbone']
row_order = [f"DP≥{dp} | {bb}" for dp in dp_order for bb in backbone_order]


def scatter_reg_corr(data, **kws):
    ax = plt.gca()
    x, y = data["GoldStandard"], data["Calculated"]
    n = len(data)

    sns.scatterplot(x=x, y=y, s=30, alpha=0.8, edgecolor="black", linewidth=0.4, ax=ax)

    if len(x) >= 3 and np.std(x) > 0 and np.std(y) > 0:
        sns.regplot(x=x, y=y, scatter=False, ci=None,
                    line_kws={"color": "crimson", "linewidth": 2}, ax=ax)
        r, p = pearsonr(x, y)
        label = f"n = {n}\nr = {r:.2f}\np = {p:.1e}"
    else:
        label = f"n = {n}\nr = NA\np = NA"

    ax.text(0.05, 0.92, label, transform=ax.transAxes, fontsize=9,
            verticalalignment="top",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                      alpha=0.9, edgecolor="lightgray"))


g = sns.FacetGrid(
    long,
    row="row_var", col="Metric",
    row_order=row_order, col_order=metric_order,
    margin_titles=True, sharex=False, sharey=False,
    height=4.0, aspect=1.3,
)
g.map_dataframe(scatter_reg_corr)
g.set_axis_labels("Gold Standard\n(ASCAT CN)", "Calculated")
g.set_titles(col_template="{col_name}", row_template="{row_name}")

sns.despine(trim=True)
plt.tight_layout()
plt.savefig("correlation_ascat_downsampled.png", dpi=150, bbox_inches="tight")
plt.show()
print("Saved: correlation_ascat_downsampled.png")
