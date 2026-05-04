#!/usr/bin/env python3
"""
evaluate_tso_vs_myriad.py

Correlates pipeline-computed HRD scores (original TSO500 samples) against
clinical gold-standard GIS scores from the MYRIAD NGS panel.

SampleID format from launchScarHRD_2_sbatch.sbatch:
  {sample_id}_dp{N}_{postmerge|premerge}
  e.g.  JBLAB341_dp8_postmerge

Parsing:
  sample_name = SampleID.split('_', n=4)[0]                      → JBLAB341
  dp_filter   = SampleID.split('_', n=4)[1].replace('dp', '')    → 8

Inputs:
  TSO_CSV      -- ALL_scarHRD_scores.csv from Stage 4 on TSO500 samples
  MYRIAD_CSV   -- GIS scores CSV with columns: sample, GIS
"""

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

# ============================================================
# PATHS
# ============================================================
TSO_CSV = (
    "/storage/scratch01/groups/co/cn_extra/alleleSpecific/tso_samples"
    "/20260421/outputs_scarHRD_thresholds/ALL_scarHRD_scores.csv"
)
MYRIAD_CSV = (
    "/storage/scratch01/groups/co/cn_extra/alleleSpecific/tso_samples"
    "/MYRIAD_GIS_SCORES/gis_scores.csv"
)

# ============================================================
# FONT SIZE CONFIGURATION
# ============================================================
TITLE_SIZE   = 20
AXIS_LABEL   = 14
TICK_SIZE    = 18
ANNOT_SIZE   = 12
MARGIN_TITLE = 18

plt.rcParams.update({
    'font.size':       TICK_SIZE,
    'axes.titlesize':  TITLE_SIZE,
    'axes.labelsize':  AXIS_LABEL,
    'xtick.labelsize': TICK_SIZE,
    'ytick.labelsize': TICK_SIZE,
})

# ============================================================
# LOAD TSO RESULTS
# ============================================================
tso_samples = pd.read_csv(TSO_CSV, sep=",")
df = tso_samples.copy()

df['sample_name'] = df['SampleID'].str.split('_', n=4).str[0]
df['dp_filter']   = df['SampleID'].str.split('_', n=4).str[1].str.replace("dp", "")
df['LoH']         = df['HRD']
df['HRD_Total']   = df['HRD-sum']

tso_hrd = df.copy()

# ============================================================
# LOAD MYRIAD GOLD STANDARD
# ============================================================
gis_scores = pd.read_csv(MYRIAD_CSV, sep=",")
df_scores = gis_scores.copy()
df_scores['sample_name'] = df_scores['sample'].str.split('-', n=1).str[0]

hrd_ref = df_scores[['sample_name', 'GIS']].copy()
hrd_ref.columns = ['gs_' + c for c in hrd_ref.columns]

# ============================================================
# MERGE
# ============================================================
merged = tso_hrd.merge(
    hrd_ref,
    left_on='sample_name', right_on='gs_sample_name',
    how='inner'
)
print(f"Merged rows: {len(merged)}")

# ============================================================
# PLOT
# ============================================================
metric_pairs = [
    ("LoH",          "gs_GIS"),
    ("Telomeric AI", "gs_GIS"),
    ("LST",          "gs_GIS"),
    ("HRD_Total",    "gs_GIS"),
]

df_plot = merged.drop_duplicates(subset=["sample_name", "dp_filter"]).copy()

all_metric_cols = list(dict.fromkeys([col for pair in metric_pairs for col in pair]))
cols_to_numeric = ["dp_filter"] + all_metric_cols
df_plot[cols_to_numeric] = df_plot[cols_to_numeric].apply(pd.to_numeric, errors="coerce")
df_plot = df_plot.dropna(subset=all_metric_cols + ["dp_filter"])

long_data = []
for calc_metric, gs_metric in metric_pairs:
    temp = df_plot[["SampleID", "dp_filter", calc_metric, gs_metric]].copy()
    temp["Metric"]       = calc_metric.replace("_", " ")
    temp["Calculated"]   = temp[calc_metric]
    temp["GoldStandard"] = temp[gs_metric]
    long_data.append(temp[["SampleID", "dp_filter", "Metric", "Calculated", "GoldStandard"]])

long = pd.concat(long_data, ignore_index=True).dropna(subset=["Calculated", "GoldStandard"])

dp_order     = sorted(long["dp_filter"].unique())
metric_order = [p[0].replace("_", " ") for p in metric_pairs]

long['row_var'] = 'DP≥' + long['dp_filter'].astype(str) + ' | '


def scatter_reg_corr(data, **kws):
    ax = plt.gca()
    x, y = data["GoldStandard"], data["Calculated"]
    metric = data["Metric"].iloc[0]
    n = len(data)

    sns.scatterplot(x=x, y=y, s=40, alpha=0.8,
                    edgecolor="black", linewidth=0.4, ax=ax)

    if len(x) >= 3 and np.std(x) > 0 and np.std(y) > 0:
        sns.regplot(x=x, y=y, scatter=False, ci=None,
                    line_kws={"color": "crimson", "linewidth": 2}, ax=ax)
        r, p = pearsonr(x, y)
        label = f"n = {n}\nr = {r:.2f}\np = {p:.1e}"
    else:
        label = f"n = {n}\nr = NA\np = NA"

    ax.text(0.5, 0.92, label,
            transform=ax.transAxes, fontsize=ANNOT_SIZE,
            verticalalignment="top",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                      alpha=0.9, edgecolor="lightgray"))

    ax.tick_params(axis='both', labelsize=TICK_SIZE)

    if metric == "HRD Total":
        ax.axvline(42, color="black", linestyle="--", linewidth=1.5)
        ax.axhline(42, color="black", linestyle="--", linewidth=1.5)
        ax.text(42, ax.get_ylim()[1] * 0.95, "42", ha="right", va="top",    fontsize=TICK_SIZE)
        ax.text(ax.get_xlim()[1] * 0.95, 42, "42", ha="right", va="bottom", fontsize=TICK_SIZE)


g = sns.FacetGrid(
    long,
    row="row_var", col="Metric",
    col_order=metric_order,
    margin_titles=True, sharex=False, sharey=False,
    height=4.5, aspect=1.3,
)
g.map_dataframe(scatter_reg_corr)
g.set_axis_labels("Gold Standard\nTSO-500-HRD", "Calculated", fontsize=AXIS_LABEL)
g.set_titles(col_template="{col_name}", row_template="{row_name}", size=MARGIN_TITLE)

sns.despine(trim=True)
plt.tight_layout()
plt.savefig("correlation_tso_vs_myriad.png", dpi=150, bbox_inches="tight")
plt.show()
print("Saved: correlation_tso_vs_myriad.png")
