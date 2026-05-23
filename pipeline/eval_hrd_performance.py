import matplotlib
matplotlib.use('Agg')

import argparse
import base64
import glob
import io
import os
import re
import sys
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import pearsonr
from sklearn.metrics import confusion_matrix, cohen_kappa_score

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description="HRD performance evaluation across all runs")
    p.add_argument("--tso-samples-root", required=True,
                   help="Root dir containing YYYYMMDD run subdirs (e.g. .../tso_samples)")
    p.add_argument("--gis-scores", required=True,
                   help="CSV with columns: sample, GIS")
    p.add_argument("--current-run", dest="current_run", default=None,
                   help="YYYYMMDD of the run to highlight in per-run plots")
    p.add_argument("--output", default=None,
                   help="Output HTML path (default: <tso-samples-root>/hrd_performance_report.html)")
    p.add_argument("--min-run", default="20260511",
                   help="Exclude runs earlier than this date (YYYYMMDD, default: 20260511)")
    p.add_argument("--exclude-samples", nargs="*", default=[],
                   help="Sample base IDs to exclude from analysis (e.g. JBLAB17037 JBLAB296)")
    return p.parse_args()


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def discover_runs(root, min_run="20260511"):
    base_runs, purity_runs = {}, {}

    for f in sorted(glob.glob(os.path.join(root, "20??????", "outputs_scarHRD", "ALL_scarHRD_scores.csv"))):
        run_date = os.path.basename(os.path.dirname(os.path.dirname(f)))
        if run_date >= min_run:
            base_runs[run_date] = f

    for f in sorted(glob.glob(os.path.join(root, "20??????_purity", "outputs_scarHRD", "ALL_scarHRD_scores.csv"))):
        dir_name  = os.path.basename(os.path.dirname(os.path.dirname(f)))
        date_part = dir_name.replace('_purity', '')
        if date_part >= min_run:
            purity_runs[dir_name] = f

    return base_runs, purity_runs


def load_scores(run_files):
    frames = []
    for run_date, path in run_files.items():
        df = pd.read_csv(path)
        df['run_date'] = run_date
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


def parse_sample_ids(df):
    sid = df['SampleID'].str.strip('"')

    # sample_name_hrd: first _-delimited token → e.g. JBLAB17058-HRD
    df['sample_name_hrd'] = sid.str.split('_').str[0]

    # sample_name: strip known -HRD suffix, then take base ID before first -
    df['sample_name'] = (
        df['sample_name_hrd']
        .str.replace(r'-HRD$', '', regex=True)
        .str.split('-').str[0]
    )

    # HRD_panel: whether the first token ends in -HRD
    df['HRD_panel'] = df['sample_name_hrd'].str.contains(r'-HRD$', regex=True)
    df['HRD_panel'] = df['HRD_panel'].map({True: 'HRD', False: 'non-HRD_panel'})

    # dp_filter: regex, position-independent
    df['dp_filter'] = sid.str.extract(r'_dp(\d+)_', expand=False).astype(float)

    # merge_stage: premerge / postmerge
    df['merge_stage'] = sid.str.extract(r'_((?:pre|post)merge)$', expand=False)

    df['LoH'] = df['HRD']
    df['HRD_Total'] = df['HRD-sum']

    return df


def load_gis(path):
    gis = pd.read_csv(path)
    gis['sample_name'] = gis['sample'].str.split('-').str[0]
    return gis[['sample_name', 'GIS']].rename(columns={'sample_name': 'gs_sample_name', 'GIS': 'gs_GIS'})


# ---------------------------------------------------------------------------
# Parsing QC audit
# ---------------------------------------------------------------------------

def build_parsing_audit(df, gis):
    audit = df[['SampleID', 'sample_name_hrd', 'sample_name', 'HRD_panel',
                'dp_filter', 'merge_stage', 'run_date']].drop_duplicates()

    # Flag unexpected patterns
    audit['dp_ok']    = audit['dp_filter'].notna()
    audit['stage_ok'] = audit['merge_stage'].notna()

    gis_names = set(gis['gs_sample_name'])
    audit['gis_match'] = audit['sample_name'].isin(gis_names)

    return audit


def parsing_audit_html(audit, gis):
    gis_names = set(gis['gs_sample_name'])
    pipeline_names = set(audit['sample_name'].unique())

    n_total    = audit['sample_name_hrd'].nunique()
    n_matched  = audit[audit['gis_match']]['sample_name_hrd'].nunique()
    n_no_dp    = (~audit['dp_ok']).sum()
    n_no_stage = (~audit['stage_ok']).sum()

    unmatched_pipeline = sorted(pipeline_names - gis_names)
    unmatched_gis      = sorted(gis_names - pipeline_names)

    match_rate = n_matched / n_total if n_total else 0
    warn_class = 'warning' if match_rate < 0.7 else 'ok'

    rows = audit.sort_values(['run_date', 'sample_name_hrd']).to_dict('records')
    table_rows = ''.join(
        f"<tr class='{'ok' if r['gis_match'] else 'fail'}'>"
        f"<td>{r['SampleID']}</td>"
        f"<td>{r['sample_name_hrd']}</td>"
        f"<td>{r['sample_name']}</td>"
        f"<td>{r['HRD_panel']}</td>"
        f"<td class='{'ok' if r['dp_ok'] else 'fail'}'>{r['dp_filter'] if r['dp_ok'] else '⚠ not found'}</td>"
        f"<td class='{'ok' if r['stage_ok'] else 'fail'}'>{r['merge_stage'] if r['stage_ok'] else '⚠ not found'}</td>"
        f"<td class='{('ok' if r['gis_match'] else 'fail')}'>{('✓' if r['gis_match'] else '✗')}</td>"
        f"<td>{r['run_date']}</td>"
        f"</tr>"
        for r in rows
    )

    unmatched_pipeline_html = (
        '<ul>' + ''.join(f'<li>{s}</li>' for s in unmatched_pipeline) + '</ul>'
        if unmatched_pipeline else '<em>none</em>'
    )
    unmatched_gis_html = (
        '<ul>' + ''.join(f'<li>{s}</li>' for s in unmatched_gis) + '</ul>'
        if unmatched_gis else '<em>none</em>'
    )

    return f"""
<section id="parsing-audit">
  <h2>Parsing audit</h2>

  <div class="summary-box {warn_class}">
    <strong>GIS match rate:</strong> {n_matched}/{n_total} unique samples
    ({match_rate*100:.0f}%)
    {'<span class="badge warning">⚠ below 70%</span>' if match_rate < 0.7 else ''}
    {'<span class="badge fail">⚠ dp_filter not parsed in '+str(n_no_dp)+' rows</span>' if n_no_dp else ''}
    {'<span class="badge fail">⚠ merge_stage not parsed in '+str(n_no_stage)+' rows</span>' if n_no_stage else ''}
  </div>

  <details>
    <summary>Pipeline samples without GIS match ({len(unmatched_pipeline)})</summary>
    {unmatched_pipeline_html}
  </details>
  <details>
    <summary>GIS samples without pipeline data ({len(unmatched_gis)})</summary>
    {unmatched_gis_html}
  </details>

  <details>
    <summary>Full parsing table ({len(rows)} rows)</summary>
    <div class="table-wrap">
      <table>
        <thead><tr>
          <th>SampleID</th><th>sample_name_hrd</th><th>sample_name</th>
          <th>HRD_panel</th><th>dp_filter</th><th>merge_stage</th>
          <th>GIS match</th><th>run_date</th>
        </tr></thead>
        <tbody>{table_rows}</tbody>
      </table>
    </div>
  </details>
</section>
"""


# ---------------------------------------------------------------------------
# Scatter plots
# ---------------------------------------------------------------------------

METRIC_PAIRS = [
    ("LoH",         "gs_GIS"),
    ("Telomeric AI", "gs_GIS"),
    ("LST",          "gs_GIS"),
    ("HRD_Total",    "gs_GIS"),
]

TITLE_SIZE   = 14
AXIS_LABEL   = 10
TICK_SIZE    = 9
ANNOT_SIZE   = 8
MARGIN_TITLE = 10

plt.rcParams.update({
    'font.size':      TICK_SIZE,
    'axes.titlesize': TITLE_SIZE,
    'axes.labelsize': AXIS_LABEL,
    'xtick.labelsize': TICK_SIZE,
    'ytick.labelsize': TICK_SIZE,
})


def _scatter_ax(ax, x, y, metric_label, n_samples):
    sns.scatterplot(x=x, y=y, s=30, alpha=0.8,
                    edgecolor='black', linewidth=0.3, ax=ax)
    if len(x) >= 3 and np.std(x) > 0 and np.std(y) > 0:
        sns.regplot(x=x, y=y, scatter=False, ci=None,
                    line_kws={'color': 'crimson', 'linewidth': 1.5}, ax=ax)
        r, p = pearsonr(x, y)
        label = f"n={n_samples}\nr={r:.2f}\np={p:.1e}"
    else:
        label = f"n={n_samples}\nr=NA"

    ax.text(0.97, 0.97, label, transform=ax.transAxes,
            fontsize=ANNOT_SIZE, va='top', ha='right',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      alpha=0.9, edgecolor='lightgray'))

    if 'HRD Total' in metric_label or metric_label == 'HRD_Total':
        ax.axvline(42, color='black', linestyle='--', linewidth=1)
        ax.axhline(42, color='black', linestyle='--', linewidth=1)


def make_scatter_figure(data, title):
    """
    data: DataFrame with columns Calculated, GoldStandard, Metric, dp_filter, merge_stage
    Returns base64-encoded PNG.
    """
    metric_order = [p[0].replace('_', ' ') for p in METRIC_PAIRS]
    dp_values    = sorted(data['dp_filter'].dropna().unique())
    stages       = sorted(data['merge_stage'].dropna().unique())

    n_rows = len(dp_values) * len(stages)
    n_cols = len(metric_order)

    if n_rows == 0 or n_cols == 0:
        fig, ax = plt.subplots(1, 1, figsize=(6, 2))
        ax.text(0.5, 0.5, 'No data (dp_filter / merge_stage not parsed)',
                ha='center', va='center', transform=ax.transAxes, fontsize=10)
        ax.axis('off')
        fig.suptitle(title, fontsize=TITLE_SIZE)
        return _fig_to_b64(fig)

    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(4 * n_cols, 3.5 * n_rows),
                             squeeze=False)

    row = 0
    for dp in dp_values:
        for stage in stages:
            for col, metric in enumerate(metric_order):
                ax = axes[row][col]
                sub = data[
                    (data['dp_filter'] == dp) &
                    (data['merge_stage'] == stage) &
                    (data['Metric'] == metric)
                ].dropna(subset=['Calculated', 'GoldStandard'])

                if sub.empty:
                    ax.axis('off')
                    col_label = metric if col == 0 else ''
                    ax.set_title(f"DP≥{int(dp)} {stage}\n{metric}", fontsize=MARGIN_TITLE)
                    row_label_done = True
                    continue

                _scatter_ax(ax, sub['GoldStandard'], sub['Calculated'],
                            metric, len(sub))
                ax.set_xlabel('Gold Standard (GIS)', fontsize=AXIS_LABEL)
                ax.set_ylabel(metric, fontsize=AXIS_LABEL)
                if col == 0:
                    ax.set_title(f"DP≥{int(dp)} | {stage}\n{metric}",
                                 fontsize=MARGIN_TITLE)
                else:
                    ax.set_title(metric, fontsize=MARGIN_TITLE)
            row += 1

    fig.suptitle(title, fontsize=TITLE_SIZE, y=1.01)
    sns.despine()
    plt.tight_layout()
    return _fig_to_b64(fig)


def build_long(df):
    frames = []
    for calc_col, gs_col in METRIC_PAIRS:
        if calc_col not in df.columns:
            continue
        tmp = df[['sample_name_hrd', 'dp_filter', 'merge_stage',
                  'HRD_panel', 'run_date', calc_col, gs_col]].copy()
        tmp['Metric']       = calc_col.replace('_', ' ')
        tmp['Calculated']   = pd.to_numeric(tmp[calc_col], errors='coerce')
        tmp['GoldStandard'] = pd.to_numeric(tmp[gs_col],   errors='coerce')
        frames.append(tmp[['sample_name_hrd', 'dp_filter', 'merge_stage',
                            'HRD_panel', 'run_date', 'Metric',
                            'Calculated', 'GoldStandard']])
    return pd.concat(frames, ignore_index=True)


# ---------------------------------------------------------------------------
# Confusion matrices
# ---------------------------------------------------------------------------

def make_confusion_figure(df, title):
    df = df.copy()
    df['HRD_class'] = df['HRD_Total'].ge(42)
    df['GIS_class']  = df['gs_GIS'].ge(42)

    dp_values = sorted(df['dp_filter'].dropna().unique().astype(int))
    stages    = sorted(df['merge_stage'].dropna().unique())

    n_cols = len(dp_values)
    n_rows = len(stages)

    if n_rows == 0 or n_cols == 0:
        fig, ax = plt.subplots(1, 1, figsize=(6, 2))
        ax.text(0.5, 0.5, 'No data (dp_filter / merge_stage not parsed)',
                ha='center', va='center', transform=ax.transAxes, fontsize=10)
        ax.axis('off')
        fig.suptitle(title, fontsize=TITLE_SIZE)
        return _fig_to_b64(fig)

    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(4 * n_cols, 4 * n_rows),
                             squeeze=False)

    for ri, stage in enumerate(stages):
        for ci, dp in enumerate(dp_values):
            ax = axes[ri][ci]
            sub = df[
                (df['dp_filter'] == dp) &
                (df['merge_stage'] == stage)
            ].drop_duplicates(subset=['sample_name_hrd'])

            if len(sub) < 2 or sub['GIS_class'].nunique() < 2:
                ax.set_title(f"DP≥{dp} | {stage}\n(insufficient data)", fontsize=9)
                ax.axis('off')
                continue

            y_true = sub['GIS_class'].astype(int)
            y_pred = sub['HRD_class'].astype(int)

            cm = confusion_matrix(y_true, y_pred)
            tn, fp, fn, tp = cm.ravel()
            sens  = tp / (tp + fn) if (tp + fn) > 0 else float('nan')
            spec  = tn / (tn + fp) if (tn + fp) > 0 else float('nan')
            kappa = cohen_kappa_score(y_true, y_pred)

            # Color by correctness (diagonal=1, off-diagonal=0), annotate with counts
            correct_mask = np.eye(cm.shape[0], dtype=float)
            sns.heatmap(correct_mask, annot=cm, fmt='d', cmap='RdYlGn',
                        vmin=0, vmax=1,
                        xticklabels=['neg', 'pos'],
                        yticklabels=['neg', 'pos'],
                        linewidths=0.5, ax=ax, cbar=False,
                        annot_kws={'size': 11})
            ax.set_xlabel('Gold standard (GIS)', fontsize=8)
            ax.set_ylabel('Predicted (HRD_Total)', fontsize=8)
            ax.set_title(
                f"DP≥{dp} | {stage}  (n={len(sub)})\n"
                f"Sens={sens:.2f}  Spec={spec:.2f}  κ={kappa:.2f}",
                fontsize=9
            )

    fig.suptitle(title, fontsize=TITLE_SIZE, y=1.01)
    plt.tight_layout()
    return _fig_to_b64(fig)


# ---------------------------------------------------------------------------
# HTML helpers
# ---------------------------------------------------------------------------

def _fig_to_b64(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=130, bbox_inches='tight')
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode()


def img_tag(b64):
    return f'<img src="data:image/png;base64,{b64}" style="max-width:100%;"/>'


CSS = """
body { font-family: sans-serif; max-width: 1400px; margin: auto; padding: 1em; }
h1   { border-bottom: 2px solid #333; padding-bottom: .3em; }
h2   { border-bottom: 1px solid #aaa; margin-top: 2em; }
h3   { color: #444; }
section { margin-bottom: 3em; }
.summary-box { padding: .6em 1em; border-radius: 4px; margin: .8em 0;
               border: 1px solid #ccc; background: #f9f9f9; }
.summary-box.warning { background: #fff3cd; border-color: #ffc107; }
.badge { display: inline-block; padding: .2em .5em; border-radius: 3px;
         font-size: .85em; margin-left: .4em; }
.badge.warning { background: #ffc107; color: #333; }
.badge.fail    { background: #dc3545; color: #fff; }
.table-wrap { overflow-x: auto; max-height: 400px; overflow-y: auto; }
table { border-collapse: collapse; font-size: .82em; width: 100%; }
th, td { border: 1px solid #ddd; padding: .3em .6em; white-space: nowrap; }
th { background: #f0f0f0; position: sticky; top: 0; }
tr.ok  td { background: #f8fff8; }
tr.fail td { background: #fff5f5; }
td.ok  { color: #198754; font-weight: bold; }
td.fail { color: #dc3545; font-weight: bold; }
details { margin: .5em 0; }
summary { cursor: pointer; color: #0077cc; }
"""


def render_html(sections, generated_at):
    body = '\n'.join(sections)
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8"/>
  <title>HRD Performance Report</title>
  <style>{CSS}</style>
</head>
<body>
  <h1>HRD Performance Report</h1>
  <p>Generated: {generated_at}</p>
  {body}
</body>
</html>"""


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def _build_analysis_sections(merged, run_files, section_id, section_label, current_run):
    sections = []

    # ── Combined ──────────────────────────────────────────────────────────
    long_all = build_long(merged)
    sections.append(f'<section id="{section_id}-combined"><h2>{section_label} — combined</h2>')
    for hrd_grp, grp in long_all.groupby('HRD_panel'):
        b64 = make_scatter_figure(grp, f"Scatter — {hrd_grp} (all runs)")
        sections.append(f'<h3>{hrd_grp}</h3>' + img_tag(b64))
    for hrd_grp, grp in merged.groupby('HRD_panel'):
        b64 = make_confusion_figure(grp, f"Confusion — {hrd_grp} (all runs)")
        sections.append(f'<h3>Confusion matrix — {hrd_grp}</h3>' + img_tag(b64))
    sections.append('</section>')

    # ── Per-run ───────────────────────────────────────────────────────────
    sections.append(f'<section id="{section_id}-per-run"><h2>{section_label} — per run</h2>')
    for run_date in sorted(run_files.keys()):
        highlight = ' ★ current' if run_date == current_run else ''
        sections.append(f'<h3>Run {run_date}{highlight}</h3>')
        run_merged = merged[merged['run_date'] == run_date]
        if run_merged.empty:
            sections.append('<p><em>No GIS matches for this run.</em></p>')
            continue
        long_run = build_long(run_merged)
        for hrd_grp, grp in long_run.groupby('HRD_panel'):
            b64 = make_scatter_figure(grp, f"Scatter — {hrd_grp} | run {run_date}")
            sections.append(f'<h4>{hrd_grp}</h4>' + img_tag(b64))
        for hrd_grp, grp in run_merged.groupby('HRD_panel'):
            b64 = make_confusion_figure(grp, f"Confusion — {hrd_grp} | run {run_date}")
            sections.append(f'<h4>Confusion — {hrd_grp}</h4>' + img_tag(b64))
    sections.append('</section>')

    return sections


def main():
    args = parse_args()

    out_path = args.output or os.path.join(args.tso_samples_root,
                                           'hrd_performance_report.html')

    # ── Load data ─────────────────────────────────────────────────────────
    base_runs, purity_runs = discover_runs(args.tso_samples_root, min_run=args.min_run)
    if not base_runs:
        sys.exit(f"No runs found under {args.tso_samples_root} (min_run={args.min_run})")

    print(f"Min run filter  : {args.min_run}")
    print(f"Base runs found : {sorted(base_runs.keys())}")
    print(f"Purity runs     : {sorted(purity_runs.keys())}")
    if args.exclude_samples:
        print(f"Excluded samples: {sorted(args.exclude_samples)}")

    all_runs = {**base_runs, **purity_runs}
    raw = load_scores(all_runs)
    df  = parse_sample_ids(raw)
    if args.exclude_samples:
        before = len(df)
        df = df[~df['sample_name'].isin(args.exclude_samples)]
        print(f"Rows after exclusion: {len(df):,} (removed {before - len(df):,})")
    gis = load_gis(args.gis_scores)

    # ── Parsing audit (all runs) ──────────────────────────────────────────
    audit     = build_parsing_audit(df, gis)
    audit_html = parsing_audit_html(audit, gis)

    # ── Merge with GIS ────────────────────────────────────────────────────
    merged = df.merge(gis, left_on='sample_name', right_on='gs_sample_name', how='inner')
    merged = merged.drop_duplicates(subset=['sample_name_hrd', 'dp_filter', 'merge_stage', 'run_date'])

    merged_base   = merged[merged['run_date'].isin(base_runs)]
    merged_purity = merged[merged['run_date'].isin(purity_runs)]

    sections = [audit_html]

    # ── Base runs ─────────────────────────────────────────────────────────
    sections += _build_analysis_sections(
        merged_base, base_runs,
        section_id='base', section_label='Base runs',
        current_run=args.current_run,
    )

    # ── Purity-corrected runs ─────────────────────────────────────────────
    if purity_runs and not merged_purity.empty:
        sections += _build_analysis_sections(
            merged_purity, purity_runs,
            section_id='purity', section_label='Purity-corrected runs',
            current_run=args.current_run + '_purity' if args.current_run else None,
        )

    # ── Write HTML ────────────────────────────────────────────────────────
    html = render_html(sections, datetime.now().strftime('%Y-%m-%d %H:%M'))
    with open(out_path, 'w') as f:
        f.write(html)

    print(f"Report written to: {out_path}")


if __name__ == '__main__':
    main()
