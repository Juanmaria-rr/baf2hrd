import matplotlib
matplotlib.use('Agg')

import argparse
import base64
import glob
import io
import json
import math
import os
import sys
import warnings
from datetime import datetime

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import numpy as np
import pandas as pd
import yaml


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description="Generate QC report for BAF-to-HRD pipeline")
    p.add_argument("--config", help="YAML config file (pipeline_config_YYYYMMDD.yaml)")
    p.add_argument("--base-dir", dest="base_dir", help="Direct override for base_dir")
    return p.parse_args()


def load_config(path):
    with open(path) as f:
        return yaml.safe_load(f)


def resolve_base_dir(args, cfg):
    if args.base_dir:
        return args.base_dir
    if cfg.get("run", {}).get("base_dir"):
        return cfg["run"]["base_dir"]
    root = cfg["paths"]["tso_samples_root"]
    date = cfg["run"]["date"]
    return os.path.join(root, date)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def shorten_name(name, maxlen=20):
    name = name.split("/")[-1]
    for prefix in ("CNIO_", "TSO500_", "Sample_", "sample_"):
        if name.startswith(prefix):
            name = name[len(prefix):]
    if len(name) > maxlen:
        name = name[:maxlen - 3] + "..."
    return name


def fig_to_b64(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=120, bbox_inches="tight")
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("ascii")


def round_cn(x):
    return int(math.floor(x + 0.5))


# ---------------------------------------------------------------------------
# Data loaders
# ---------------------------------------------------------------------------

def load_merged_csv(base_dir):
    path = os.path.join(
        base_dir,
        "merged_segments_baf_metrics_ALL_with20bins_enhancedFeatures.csv",
    )
    if not os.path.exists(path):
        warnings.warn(f"Merged CSV not found: {path}")
        return None
    df = pd.read_csv(path)
    df["total_cn_int"] = df["total_cn"].apply(round_cn)
    df["sample_label"] = df["sample_name"].apply(shorten_name)
    return df


def load_thresholds(thresholds_file):
    if not thresholds_file or not os.path.exists(thresholds_file):
        warnings.warn(f"Thresholds file not found: {thresholds_file}")
        return None
    with open(thresholds_file) as f:
        raw = json.load(f)
    # Nested: {dp: {cn: {comparison: thr}}} — second-level values are dicts of floats.
    # Flat:   {cn: {comparison: thr}}       — second-level values are floats.
    first_val = next(iter(raw.values()))
    if not isinstance(first_val, dict):
        return {"flat": raw}
    inner_val = next(iter(first_val.values()))
    if isinstance(inner_val, dict):
        return {"nested": raw}
    return {"flat": raw}


def get_thresholds_for_dp(thresholds, dp):
    """Return {cn_str: {comparison: value}} for a given dp (int or str)."""
    if thresholds is None:
        return {}
    if "flat" in thresholds:
        return thresholds["flat"]
    nested = thresholds["nested"]
    key = str(dp)
    return nested.get(key, {})


def load_scarHRD_inputs(base_dir):
    """Return dict keyed by (sample, dp, merge_state) → DataFrame."""
    pattern = os.path.join(base_dir, "scarHRD_inputs_thresholds", "*_scarHRD_input.txt")
    files = glob.glob(pattern)
    result = {}
    for fp in files:
        fname = os.path.basename(fp)
        # {sample}_dp{dp}_{merge_state}_scarHRD_input.txt
        try:
            stem = fname.replace("_scarHRD_input.txt", "")
            parts = stem.rsplit("_", 1)  # split off merge_state
            merge_state = parts[1]
            rest = parts[0]
            dp_idx = rest.rfind("_dp")
            sample = rest[:dp_idx]
            dp = int(rest[dp_idx + 3:])
        except Exception:
            warnings.warn(f"Cannot parse filename: {fname}")
            continue
        try:
            df = pd.read_csv(fp, sep="\t")
            result[(sample, dp, merge_state)] = df
        except Exception as e:
            warnings.warn(f"Cannot read {fp}: {e}")
    return result


# ---------------------------------------------------------------------------
# Plot 1: BAF coverage (n_positions per segment per sample)
# ---------------------------------------------------------------------------

def plot_baf_coverage(df, dp_values):
    if df is None:
        return None

    n_dp = len(dp_values)
    n_samples = df["sample_label"].nunique()
    fig, axes = plt.subplots(1, n_dp, figsize=(max(4, n_samples * 0.35) * n_dp, 6), sharey=False)
    if n_dp == 1:
        axes = [axes]

    for ax, dp in zip(axes, dp_values):
        sub = df[df["dp_filter"] == dp].copy()
        if sub.empty:
            ax.set_title(f"DP={dp}\n(no data)")
            continue

        order = (
            sub.groupby("sample_label")["n_positions"]
            .median()
            .sort_values(ascending=False)
            .index.tolist()
        )
        grouped = [sub[sub["sample_label"] == s]["n_positions"].dropna().values for s in order]

        ax.boxplot(grouped, tick_labels=order, patch_artist=True,
                   boxprops=dict(facecolor="steelblue", alpha=0.6),
                   medianprops=dict(color="navy", linewidth=2),
                   flierprops=dict(marker=".", markersize=3, alpha=0.4))
        ax.axhline(10, color="red", linestyle="--", linewidth=1.2, label="n=10 (sparse)")
        ax.set_yscale("log")
        ax.set_title(f"DP filter = {dp}")
        ax.set_xlabel("Sample")
        ax.set_ylabel("n_positions (log)" if dp == dp_values[0] else "")
        ax.tick_params(axis="x", rotation=45)
        ax.legend(fontsize=7)

    fig.suptitle("BAF coverage per segment (n_positions)", fontsize=13, fontweight="bold")
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Plot 2: outer_mass_raw by total_cn with threshold lines
# ---------------------------------------------------------------------------

COMPARISON_COLORS = {
    "loh_vs_balanced": "dodgerblue",
    "loh_vs_ai": "darkorange",
    "ai_vs_balanced": "forestgreen",
}


def plot_outer_mass(df, dp_values, thresholds):
    if df is None:
        return None

    n_dp = len(dp_values)
    n_samples = df["sample_label"].nunique()
    fig, axes = plt.subplots(1, n_dp, figsize=(max(5, n_samples * 0.2) * n_dp, 6), sharey=True)
    if n_dp == 1:
        axes = [axes]

    cn_range = list(range(2, 9))

    for ax, dp in zip(axes, dp_values):
        sub = df[(df["dp_filter"] == dp) & (df["total_cn_int"].between(2, 8))].copy()
        thr = get_thresholds_for_dp(thresholds, dp)

        positions = []
        data_by_cn = []
        for i, cn in enumerate(cn_range):
            vals = sub[sub["total_cn_int"] == cn]["outer_mass_raw"].dropna().values
            if len(vals) > 0:
                positions.append(i + 1)
                data_by_cn.append(vals)

        if data_by_cn:
            parts = ax.violinplot(data_by_cn, positions=positions, showmedians=True, widths=0.7)
            for pc, pos in zip(parts["bodies"], positions):
                cn = cn_range[pos - 1]
                color = plt.cm.tab10(cn % 10)
                pc.set_facecolor(color)
                pc.set_alpha(0.6)

        # Threshold lines
        drawn_labels = set()
        for i, cn in enumerate(cn_range):
            cn_str = str(cn)
            if cn_str not in thr:
                continue
            for comp, val in thr[cn_str].items():
                color = COMPARISON_COLORS.get(comp, "grey")
                label = comp if comp not in drawn_labels else None
                ax.hlines(val, i + 0.6, i + 1.4, colors=color, linestyles="--",
                          linewidth=1.5, label=label)
                ax.text(i + 1.45, val, f"{val:.3f}", fontsize=6, va="center", color=color)
                drawn_labels.add(comp)

        ax.set_xticks(range(1, len(cn_range) + 1))
        ax.set_xticklabels(cn_range)
        ax.set_xlabel("total_cn")
        ax.set_ylabel("outer_mass_raw" if dp == dp_values[0] else "")
        ax.set_ylim(0, 1.05)
        ax.set_title(f"DP filter = {dp}")
        if drawn_labels:
            ax.legend(fontsize=7, loc="lower right")

    fig.suptitle(
        "outer_mass_raw distribution by total_cn with calibrated thresholds",
        fontsize=12, fontweight="bold",
    )
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Plot 3: Threshold margin heatmap
# ---------------------------------------------------------------------------

def compute_margin(row, thr_for_dp):
    cn = row["total_cn_int"]
    om = row["outer_mass_raw"]
    if pd.isna(om) or cn < 2 or cn > 8:
        return np.nan

    cn_str = str(cn)
    cn_thr = thr_for_dp.get(cn_str, {})

    if cn == 2:
        keys = ["loh_vs_balanced"]
    elif cn == 3:
        keys = ["loh_vs_ai"]
    elif cn == 4:
        keys = ["loh_vs_balanced", "ai_vs_balanced"]
    elif cn == 5:
        keys = ["loh_vs_ai"]
    elif cn == 6:
        keys = ["ai_vs_balanced"]
    else:  # 7, 8
        return abs(om - 0.95)

    margins = [abs(om - cn_thr[k]) for k in keys if k in cn_thr]
    return min(margins) if margins else np.nan


def plot_margin_heatmap(df, dp_values, thresholds):
    if df is None:
        return None

    n_dp = len(dp_values)
    samples_all = df["sample_label"].unique().tolist()
    n_samples = len(samples_all)
    cn_range = list(range(2, 9))

    fig, axes = plt.subplots(1, n_dp, figsize=(max(6, n_samples * 0.45) * n_dp, 5))
    if n_dp == 1:
        axes = [axes]

    cmap = plt.cm.RdYlGn
    norm = mcolors.Normalize(vmin=0, vmax=0.3)

    for ax, dp in zip(axes, dp_values):
        sub = df[df["dp_filter"] == dp].copy()
        thr = get_thresholds_for_dp(thresholds, dp)

        sub["margin"] = sub.apply(lambda r: compute_margin(r, thr), axis=1)

        samples = sorted(
            samples_all,
            key=lambda s: sub[sub["sample_label"] == s]["margin"].median()
        )

        grid = np.full((len(cn_range), len(samples)), np.nan)
        counts = np.zeros((len(cn_range), len(samples)), dtype=int)

        for si, s in enumerate(samples):
            for ci, cn in enumerate(cn_range):
                cell = sub[(sub["sample_label"] == s) & (sub["total_cn_int"] == cn)]["margin"].dropna()
                if len(cell) > 0:
                    grid[ci, si] = cell.median()
                    counts[ci, si] = len(cell)

        masked = np.ma.masked_invalid(grid)
        im = ax.imshow(masked, aspect="auto", cmap=cmap, norm=norm,
                       interpolation="nearest", origin="upper")

        for ci in range(len(cn_range)):
            for si in range(len(samples)):
                n = counts[ci, si]
                if n > 0:
                    ax.text(si, ci, str(n), ha="center", va="center",
                            fontsize=6, color="black")

        ax.set_xticks(range(len(samples)))
        ax.set_xticklabels(samples, rotation=45, ha="right", fontsize=7)
        ax.set_yticks(range(len(cn_range)))
        ax.set_yticklabels([f"CN={cn}" for cn in cn_range], fontsize=8)
        ax.set_title(f"DP filter = {dp}")

    plt.colorbar(im, ax=axes[-1], label="Median margin (0=uncertain, 0.3+=clear)")
    fig.suptitle(
        "Median threshold margin (classification uncertainty)",
        fontsize=12, fontweight="bold",
    )
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Plot 4: Allelic state composition per sample (DP=8 premerge)
# ---------------------------------------------------------------------------

def classify_allelic_state(row):
    a, b, total = row["A_cn"], row["B_cn"], row["total_cn"]
    if total == 0:
        return "Deletion"
    if min(a, b) == 0:
        return "LoH"
    if a == b:
        return "Balanced"
    return "AI"


def plot_allelic_composition(scarHRD_data, dp_values):
    # Prefer DP=8, else highest available
    target_dp = None
    for dp in sorted(dp_values, reverse=True):
        if any(k[1] == dp and k[2] == "premerge" for k in scarHRD_data):
            target_dp = dp
            break

    if target_dp is None:
        warnings.warn("No premerge scarHRD input files found for any DP — skipping Plot 4")
        return None

    keys = [k for k in scarHRD_data if k[1] == target_dp and k[2] == "premerge"]
    if not keys:
        warnings.warn(f"No premerge files for DP={target_dp}")
        return None

    frames = []
    for k in keys:
        df = scarHRD_data[k].copy()
        df["_sample"] = k[0]
        frames.append(df)
    combined = pd.concat(frames, ignore_index=True)

    combined["seg_len"] = combined["End_position"] - combined["Start_position"]
    combined["state"] = combined.apply(classify_allelic_state, axis=1)
    combined["sample_label"] = combined["_sample"].apply(shorten_name)

    state_order = ["LoH", "Balanced", "AI", "Deletion"]
    colors = {"LoH": "crimson", "Balanced": "steelblue", "AI": "darkorange", "Deletion": "mediumpurple"}

    fractions = (
        combined.groupby(["sample_label", "state"])["seg_len"]
        .sum()
        .unstack(fill_value=0)
        .reindex(columns=state_order, fill_value=0)
    )
    fractions = fractions.div(fractions.sum(axis=1), axis=0)
    fractions = fractions.sort_values("LoH", ascending=True)

    n_samples = len(fractions)
    fig, ax = plt.subplots(figsize=(8, max(4, n_samples * 0.35)))

    left = np.zeros(n_samples)
    for state in state_order:
        vals = fractions[state].values
        ax.barh(range(n_samples), vals, left=left, color=colors[state], label=state, height=0.7)
        left += vals

    ax.set_yticks(range(n_samples))
    ax.set_yticklabels(fractions.index.tolist(), fontsize=8)
    ax.set_xlabel("Fraction of segment length")
    ax.set_xlim(0, 1)
    ax.legend(loc="lower right", fontsize=9)
    ax.set_title(
        f"Allelic state composition by segment length (DP={target_dp}, premerge)",
        fontsize=12, fontweight="bold",
    )
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Plot 5: Pre vs post-merge segment count
# ---------------------------------------------------------------------------

def plot_segment_counts(scarHRD_data, dp_values):
    if not scarHRD_data:
        return None

    records = []
    for (sample, dp, merge_state), df in scarHRD_data.items():
        records.append({
            "sample": sample,
            "sample_label": shorten_name(sample),
            "dp": dp,
            "merge_state": merge_state,
            "n_segs": len(df),
        })
    rec_df = pd.DataFrame(records)

    pre = rec_df[rec_df["merge_state"] == "premerge"][["sample_label", "dp", "n_segs"]].rename(columns={"n_segs": "pre"})
    post = rec_df[rec_df["merge_state"] == "postmerge"][["sample_label", "dp", "n_segs"]].rename(columns={"n_segs": "post"})

    if pre.empty or post.empty:
        warnings.warn("Pre or post-merge files missing — skipping Plot 5")
        return None

    merged = pre.merge(post, on=["sample_label", "dp"], how="inner")
    if merged.empty:
        warnings.warn("No matched pre/post-merge pairs — skipping Plot 5")
        return None

    fig, ax = plt.subplots(figsize=(7, 6))

    dp_colors = {dp: plt.cm.Set1(i / max(len(dp_values) - 1, 1)) for i, dp in enumerate(sorted(dp_values))}

    for dp_val in sorted(dp_values):
        sub = merged[merged["dp"] == dp_val]
        if sub.empty:
            continue
        sc = ax.scatter(
            sub["pre"], sub["post"],
            c=[dp_colors[dp_val]] * len(sub),
            label=f"DP={dp_val}",
            alpha=0.8, s=60, edgecolors="white", linewidth=0.5,
        )
        for _, row in sub.iterrows():
            ratio = row["post"] / row["pre"] if row["pre"] > 0 else 1.0
            ax.annotate(
                f"{ratio:.2f}",
                (row["pre"], row["post"]),
                textcoords="offset points", xytext=(4, 2),
                fontsize=6, alpha=0.7,
            )

    max_val = max(merged["pre"].max(), merged["post"].max()) * 1.05
    ax.plot([0, max_val], [0, max_val], "k--", linewidth=1, alpha=0.5, label="y=x")
    ax.set_xlabel("Segment count (premerge)")
    ax.set_ylabel("Segment count (postmerge)")
    ax.set_xlim(0, max_val)
    ax.set_ylim(0, max_val)
    ax.legend(fontsize=9)
    ax.set_title("Segment count: pre vs post-merge", fontsize=12, fontweight="bold")
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# HTML assembly
# ---------------------------------------------------------------------------

PLOT_DESCRIPTIONS = {
    "baf_coverage": (
        "BAF coverage per segment",
        "Look for samples with many segments near or below n=10 (red line) — "
        "those may lack enough heterozygous positions for reliable outer_mass estimation. "
        "The log scale makes low-coverage outliers visible."
    ),
    "outer_mass": (
        "outer_mass_raw distribution by total_cn",
        "Violin bodies show the spread of outer_mass values within each CN class. "
        "Dashed lines mark calibrated thresholds: blue=loh_vs_balanced, orange=loh_vs_ai, green=ai_vs_balanced. "
        "Overlapping distributions near threshold lines indicate ambiguous classification regions."
    ),
    "margin_heatmap": (
        "Threshold margin heatmap",
        "Red cells indicate segments where outer_mass falls close to (or on) a classification boundary — "
        "small perturbations would flip the allelic state call. "
        "Numbers in each cell show how many segments contribute to the median. "
        "Empty cells mean no segments of that CN were observed for that sample."
    ),
    "allelic_composition": (
        "Allelic state composition (DP=8, premerge)",
        "Fraction of total segment length in each allelic state, weighted by genomic span. "
        "High LoH fraction (red) is expected for HRD-positive tumours. "
        "Large Deletion bars may indicate poor segmentation or low-purity samples."
    ),
    "segment_counts": (
        "Pre vs post-merge segment counts",
        "Points below the diagonal indicate merging reduced segment count. "
        "Labels show compression ratio (post/pre). "
        "Extreme compression (ratio <0.5) may indicate over-aggressive merging or broad genomic regions with uniform state."
    ),
}


def build_html(plots, base_dir, date_str):
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    plot_keys = ["baf_coverage", "outer_mass", "margin_heatmap", "allelic_composition", "segment_counts"]

    sections = []
    for key in plot_keys:
        title, caption = PLOT_DESCRIPTIONS[key]
        if key not in plots or plots[key] is None:
            sections.append(f"""
            <div class="card skipped">
                <h2>{title}</h2>
                <p class="skip-msg">Plot skipped — required data not found.</p>
            </div>
            """)
            continue
        b64 = fig_to_b64(plots[key])
        sections.append(f"""
        <div class="card">
            <h2>{title}</h2>
            <figure>
                <img src="data:image/png;base64,{b64}" alt="{title}" />
                <figcaption>{caption}</figcaption>
            </figure>
        </div>
        """)

    sections_html = "\n".join(sections)
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8" />
<meta name="viewport" content="width=device-width, initial-scale=1.0" />
<title>QC Report — BAF-to-HRD Pipeline</title>
<style>
  * {{ box-sizing: border-box; margin: 0; padding: 0; }}
  body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
          background: #f5f7fa; color: #2d3748; }}
  header {{ background: #1a202c; color: #e2e8f0; padding: 24px 40px; }}
  header h1 {{ font-size: 1.6rem; font-weight: 700; }}
  header p {{ margin-top: 6px; font-size: 0.9rem; color: #a0aec0; }}
  main {{ max-width: 1600px; margin: 32px auto; padding: 0 32px; }}
  .card {{ background: #fff; border-radius: 8px; box-shadow: 0 1px 4px rgba(0,0,0,.12);
           margin-bottom: 32px; padding: 24px; }}
  .card h2 {{ font-size: 1.1rem; font-weight: 600; margin-bottom: 16px;
              border-bottom: 2px solid #edf2f7; padding-bottom: 8px; color: #2d3748; }}
  .card img {{ max-width: 100%; height: auto; display: block; margin: 0 auto; }}
  figcaption {{ margin-top: 12px; font-size: 0.82rem; color: #718096;
                line-height: 1.55; max-width: 900px; }}
  .skipped {{ opacity: 0.6; }}
  .skip-msg {{ color: #e53e3e; font-size: 0.9rem; margin-top: 8px; }}
  footer {{ text-align: center; padding: 24px; font-size: 0.78rem; color: #a0aec0; }}
</style>
</head>
<body>
<header>
  <h1>QC Report — BAF-to-HRD Pipeline</h1>
  <p>Date: {date_str} &nbsp;|&nbsp; base_dir: {base_dir}</p>
</header>
<main>
{sections_html}
</main>
<footer>Generated: {timestamp}</footer>
</body>
</html>
"""


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    if args.config:
        cfg = load_config(args.config)
    else:
        cfg = {"run": {"date": "", "base_dir": None}, "paths": {"tso_samples_root": ""}, "processing": {"dp_values": [0, 2, 4, 8]}, "input": {"thresholds_file": None}}

    base_dir = resolve_base_dir(args, cfg)
    dp_values = cfg.get("processing", {}).get("dp_values", [0, 2, 4, 8])
    thresholds_file = cfg.get("input", {}).get("thresholds_file")
    date_str = cfg.get("run", {}).get("date", "")

    print(f"base_dir:        {base_dir}")
    print(f"dp_values:       {dp_values}")
    print(f"thresholds_file: {thresholds_file}")

    df = load_merged_csv(base_dir)
    thresholds = load_thresholds(thresholds_file)
    scarHRD_data = load_scarHRD_inputs(base_dir)
    print(f"scarHRD files loaded: {len(scarHRD_data)}")

    plots = {}

    print("Generating Plot 1: BAF coverage...")
    try:
        plots["baf_coverage"] = plot_baf_coverage(df, dp_values)
    except Exception as e:
        warnings.warn(f"Plot 1 failed: {e}")
        plots["baf_coverage"] = None

    print("Generating Plot 2: outer_mass distribution...")
    try:
        plots["outer_mass"] = plot_outer_mass(df, dp_values, thresholds)
    except Exception as e:
        warnings.warn(f"Plot 2 failed: {e}")
        plots["outer_mass"] = None

    print("Generating Plot 3: margin heatmap...")
    try:
        plots["margin_heatmap"] = plot_margin_heatmap(df, dp_values, thresholds)
    except Exception as e:
        warnings.warn(f"Plot 3 failed: {e}")
        plots["margin_heatmap"] = None

    print("Generating Plot 4: allelic composition...")
    try:
        plots["allelic_composition"] = plot_allelic_composition(scarHRD_data, dp_values)
    except Exception as e:
        warnings.warn(f"Plot 4 failed: {e}")
        plots["allelic_composition"] = None

    print("Generating Plot 5: segment counts...")
    try:
        plots["segment_counts"] = plot_segment_counts(scarHRD_data, dp_values)
    except Exception as e:
        warnings.warn(f"Plot 5 failed: {e}")
        plots["segment_counts"] = None

    html = build_html(plots, base_dir, date_str)
    out_path = os.path.join(base_dir, "qc_report.html")
    os.makedirs(base_dir, exist_ok=True)
    with open(out_path, "w") as f:
        f.write(html)
    print(f"QC report written to: {out_path}")

    for fig in plots.values():
        if fig is not None:
            plt.close(fig)


if __name__ == "__main__":
    main()
