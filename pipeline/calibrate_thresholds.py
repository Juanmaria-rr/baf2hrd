#!/usr/bin/env python3
"""
calibrate_thresholds.py — Threshold calibration for allele-specific CN classification.

Uses ASCAT-downsampled samples (which carry ground-truth nMajor/nMinor) to find the
optimal outer_mass decision boundaries that separate allelic states for each total CN
value.  The resulting thresholds_lookup.json is consumed by Stage 3 of
baf_to_hrd_pipeline.py.

Origin: section "How to decide for thresholds in order to implement a rule-base
framework for inputting to scarHRD" in lastVersion_hrdPredictor.ipynb.

Usage:
    python calibrate_thresholds.py

Outputs (all under PARAMS['outdir']):
    threshold_*.png                 — distribution + ROC per comparison
    confusion_*.png                 — confusion matrix per comparison
    threshold_summary_complete.csv  — one row per comparison
    thresholds_lookup.json          — machine-readable lookup consumed by Stage 3
    separability_summary_all_cn.png — Cohen's d / AUC overview across all CN values
"""

import argparse
import os
import json

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.metrics import roc_curve, auc, accuracy_score, confusion_matrix

# ============================================================
# CONFIGURATION
# ============================================================

PARAMS = {
    # Stage 2 merged CSV produced by baf_to_hrd_pipeline.py (must contain ASCAT rows
    # with nMajor/nMinor ground truth to calibrate thresholds from).
    'input_file': (
        '/storage/scratch01/groups/co/cn_extra/alleleSpecific/tso_samples/20260406'
        '/merged_segments_baf_metrics_ALL_with20bins_enhancedFeatures.csv'
    ),
    # DP filters to calibrate — one independent threshold set per filter.
    # The same DP filter must be used when running Stage 3 on TSO samples so that
    # outer_mass values are computed on the same BAF position subset as calibration.
    'dp_filters': [0, 2, 4, 8],
    # Minimum segments required per allelic state to include it in a comparison
    'min_samples': 50,
    # BAF feature used as the decision variable.
    # We calibrate on outer_mass_raw (uncorrected BAF) rather than the
    # purity-corrected outer_mass because TSO500 samples have no purity estimate —
    # calibrating on a corrected quantity and applying to an uncorrected one
    # compares different scales.  outer_mass_raw is computable for both cohorts.
    #
    # Auto-detect logic (None = auto):
    #   new pipeline  → 'outer_mass_raw'
    #   legacy 20260406 CSV (pandas merge suffix artifact) → 'outer_mass_raw_y'
    'feature': None,
    # Root output directory — per-DP results go into {outdir}/dp{N}/ subdirectories.
    # The combined thresholds_lookup.json is saved directly under outdir.
    'outdir': (
        '/storage/scratch01/groups/co/cn_extra/alleleSpecific/tso_samples/20260406'
        '/thresholds_analysis'
    ),
    'random_state': 42,
    # Optional: path to a plain-text file (one sample ID per line) listing the
    # samples to USE for calibration.  When set, all other samples are excluded.
    # Use split_samples.py to generate this file for an honest train/val split.
    # Set to None to use all ASCAT-profiled samples (original behaviour).
    'calibration_samples_file': None,
}

# CLI overrides — allow non-interactive runs without editing PARAMS
_p = argparse.ArgumentParser(add_help=False)
_p.add_argument('--input',   default=None, help='Path to merged segments CSV')
_p.add_argument('--feature', default=None, help='Column to use as decision variable')
_p.add_argument('--outdir',  default=None, help='Output directory (overrides PARAMS)')
_cli, _ = _p.parse_known_args()
if _cli.input:   PARAMS['input_file'] = _cli.input
if _cli.feature: PARAMS['feature']    = _cli.feature
if _cli.outdir:  PARAMS['outdir']     = _cli.outdir

# ============================================================
# LOAD DATA
# ============================================================

print("\n" + "=" * 70)
print("LOADING DATA")
print("=" * 70)

df_full = pd.read_csv(PARAMS['input_file'])
print(f"✓ Loaded: {len(df_full):,} segments, {len(df_full.columns)} columns")

# Optional sample filter for honest train/val split
if PARAMS.get('calibration_samples_file'):
    _cal_file = PARAMS['calibration_samples_file']
    if not os.path.exists(_cal_file):
        raise FileNotFoundError(f"calibration_samples_file not found: {_cal_file}")
    with open(_cal_file) as _f:
        _cal_samples = set(l.strip() for l in _f if l.strip() and not l.startswith('#'))
    _sample_col = next((c for c in ('sample', 'SampleID', 'sample_name') if c in df_full.columns), None)
    if _sample_col is None:
        raise ValueError("Cannot find a sample ID column in the merged CSV.")
    _before = len(df_full)
    df_full = df_full[df_full[_sample_col].isin(_cal_samples)].copy()
    print(f"✓ Filtered to {len(_cal_samples):,} calibration samples: "
          f"{len(df_full):,} / {_before:,} segments retained")

# Auto-detect the raw outer_mass column.
# We use the *uncorrected* (raw) outer_mass so that calibration and TSO500
# application are on the same BAF scale (TSO500 samples have no purity estimate).
# New pipeline produces 'outer_mass_raw'; legacy 20260406 CSV has 'outer_mass_raw_y'
# (pandas merge suffix artifact from the Stage 2 v2 script).
if PARAMS['feature'] is None:
    if 'outer_mass_raw' in df_full.columns:
        PARAMS['feature'] = 'outer_mass_raw'
    elif 'outer_mass_raw_y' in df_full.columns:
        PARAMS['feature'] = 'outer_mass_raw_y'
        print("⚠️  Using 'outer_mass_raw_y' (legacy 20260406 merge artifact). "
              "New pipeline outputs will use 'outer_mass_raw'.")
    else:
        raise ValueError(
            "Cannot find raw outer_mass column. Expected 'outer_mass_raw' (new pipeline) "
            "or 'outer_mass_raw_y' (legacy 20260406 CSV). "
            f"Available columns: {[c for c in df_full.columns if 'mass' in c]}"
        )
print(f"✓ Feature column: {PARAMS['feature']}")
print(f"✓ DP filters to calibrate: {PARAMS['dp_filters']}")

os.makedirs(PARAMS['outdir'], exist_ok=True)

# ============================================================
# SEGMENT COUNT SUMMARY BY CN
# ============================================================

def summarize_cn_distribution(df_filtered):
    """
    Print per-CN / per-allelic-state segment counts and return a dict of states.

    Returns
    -------
    cn_states_dict : dict
        {total_cn: {'loh': [...], 'ai': [...], 'balanced': [...], 'total_segments': int}}
        Each inner list contains ((nMajor, nMinor), count) tuples.
    """
    print("\n" + "=" * 70)
    print("SEGMENT DISTRIBUTION BY TOTAL CN AND ALLELIC STATE")
    print("=" * 70)

    total_segments = len(df_filtered)
    cn_states_dict = {}

    for total_cn in sorted(df_filtered['total_cn'].unique()):
        cn_data = df_filtered[df_filtered['total_cn'] == total_cn]
        cn_pct = 100 * len(cn_data) / total_segments

        print(f"\n{'=' * 70}")
        print(f"Total CN = {total_cn}: {len(cn_data):,} segments ({cn_pct:.1f}% of all)")
        print(f"{'=' * 70}")

        state_counts = cn_data.groupby(['nMajor', 'nMinor']).size().sort_values(ascending=False)

        loh_states, ai_states, balanced_states = [], [], []

        for (nmaj, nmin), count in state_counts.items():
            pct = 100 * count / len(cn_data)
            if nmaj == 0 or nmin == 0:
                class_type = 'LoH'
                loh_states.append(((nmaj, nmin), count))
            elif nmaj == nmin:
                class_type = 'Balanced'
                balanced_states.append(((nmaj, nmin), count))
            else:
                class_type = 'AI'
                ai_states.append(((nmaj, nmin), count))
            print(f"  {nmaj}+{nmin} ({class_type:8s}): {count:6,} ({pct:5.1f}%)")

        cn_states_dict[total_cn] = {
            'loh': loh_states,
            'ai': ai_states,
            'balanced': balanced_states,
            'total_segments': len(cn_data),
        }

        loh_count = sum(c for _, c in loh_states)
        ai_count  = sum(c for _, c in ai_states)
        bal_count = sum(c for _, c in balanced_states)
        print(f"\n  Class Summary:")
        print(f"    LoH:      {loh_count:6,} ({100 * loh_count / len(cn_data):5.1f}%)")
        print(f"    AI:       {ai_count:6,} ({100 * ai_count  / len(cn_data):5.1f}%)")
        print(f"    Balanced: {bal_count:6,} ({100 * bal_count / len(cn_data):5.1f}%)")

    print("\n" + "=" * 70)
    return cn_states_dict

# ============================================================
# SINGLE-COMPARISON ANALYSIS
# ============================================================

def analyze_cn_state(df_filtered, total_cn, state_pairs,
                     feature,
                     comparison_name=None,
                     outdir=None):
    """
    Fit and evaluate an outer_mass threshold for one binary comparison.

    Parameters
    ----------
    df_filtered   : DataFrame   — filtered merged CSV (one dp_filter, all ground-truth rows)
    total_cn      : int         — total copy number to analyse
    state_pairs   : list        — exactly two (nMajor, nMinor) tuples to compare
    feature       : str         — column used as decision variable
    comparison_name : str       — label used in output filenames
    outdir        : str         — directory for saved figures

    Returns
    -------
    dict with threshold, AUC, Cohen's d, accuracy, etc.; or None if insufficient data.
    """
    if comparison_name is None:
        comparison_name = f"cn{total_cn}"
    if outdir is None:
        outdir = PARAMS['outdir']

    print("\n" + "=" * 70)
    print(f"ANALYSIS: Total CN = {total_cn}")
    print(f"Comparison: {state_pairs[0]} vs {state_pairs[1]}")
    print("=" * 70)

    cn_data = df_filtered[df_filtered['total_cn'] == total_cn].copy()
    if len(cn_data) == 0:
        print(f"No segments with total CN = {total_cn}")
        return None

    # --- separate states ---
    def _classify(nmaj, nmin):
        if nmaj == 0 or nmin == 0:
            return 'LoH'
        if nmaj == nmin:
            return 'Balanced'
        return 'AI'

    state_data   = {}
    state_labels = {}
    for nmaj, nmin in state_pairs:
        mask = (cn_data['nMajor'] == nmaj) & (cn_data['nMinor'] == nmin)
        state_data[(nmaj, nmin)]   = cn_data[mask].copy()
        state_labels[(nmaj, nmin)] = _classify(nmaj, nmin)

    print(f"\nState distribution:")
    for (nmaj, nmin), data in state_data.items():
        pct = 100 * len(data) / len(cn_data)
        print(f"  {nmaj}+{nmin} ({state_labels[(nmaj, nmin)]:8s}): {len(data):6,} ({pct:5.1f}%)")

    if len(state_pairs) != 2:
        print(f"⚠️  Need exactly 2 states (got {len(state_pairs)})")
        return None

    # drop NaN
    clean_data = {}
    for state, data in state_data.items():
        clean = data.dropna(subset=[feature])
        removed = len(data) - len(clean)
        if removed:
            print(f"  Removed {removed} NaN rows from {state}")
        clean_data[state] = clean

    state1, state2 = state_pairs
    data1, data2   = clean_data[state1], clean_data[state2]
    label1, label2 = state_labels[state1], state_labels[state2]

    if len(data1) == 0 or len(data2) == 0:
        print("❌ No valid data after removing NaN")
        return None

    # --- distribution statistics ---
    def _desc(s):
        return {
            'mean': s.mean(), 'median': s.median(), 'std': s.std(),
            'min': s.min(), 'max': s.max(),
            'q1': s.quantile(0.25), 'q3': s.quantile(0.75),
        }

    stats1, stats2 = _desc(data1[feature]), _desc(data2[feature])

    print(f"\n{state1[0]}+{state1[1]} ({label1}):")
    for k, v in stats1.items():
        print(f"  {k:6s}: {v:.4f}")
    print(f"\n{state2[0]}+{state2[1]} ({label2}):")
    for k, v in stats2.items():
        print(f"  {k:6s}: {v:.4f}")

    # --- assign positive / negative class (LoH > AI > Balanced) ---
    priority = {'LoH': 0, 'AI': 1, 'Balanced': 2}
    if priority[label1] <= priority[label2]:
        pos_state, neg_state = state1, state2
        pos_data,  neg_data  = data1, data2
        pos_label, neg_label = label1, label2
    else:
        pos_state, neg_state = state2, state1
        pos_data,  neg_data  = data2, data1
        pos_label, neg_label = label2, label1

    y_true = np.concatenate([np.ones(len(pos_data)), np.zeros(len(neg_data))])
    feature_values = np.concatenate([pos_data[feature].values, neg_data[feature].values])

    # remove any remaining non-finite values
    valid = np.isfinite(feature_values)
    if not valid.all():
        print(f"⚠️  Removed {(~valid).sum()} NaN/Inf values")
        y_true, feature_values = y_true[valid], feature_values[valid]

    # --- ROC / Youden ---
    fpr, tpr, roc_thresholds = roc_curve(y_true, feature_values)
    roc_auc = auc(fpr, tpr)
    youden_idx = np.argmax(tpr - fpr)
    optimal_threshold = roc_thresholds[youden_idx]

    print(f"\nAUC-ROC: {roc_auc:.4f}")
    print(f"Optimal threshold (Youden): {optimal_threshold:.4f}")
    print(f"  Sensitivity ({pos_label}): {tpr[youden_idx]:.4f}")
    print(f"  Specificity ({neg_label}): {1 - fpr[youden_idx]:.4f}")

    # --- Cohen's d ---
    pooled_std = np.sqrt((stats1['std'] ** 2 + stats2['std'] ** 2) / 2)
    cohens_d   = abs(stats1['mean'] - stats2['mean']) / pooled_std

    print(f"\nCohen's d: {cohens_d:.4f}", end="")
    if cohens_d > 0.8:
        print(" → LARGE effect (easily separable)")
        separability = "HIGH"
    elif cohens_d > 0.5:
        print(" → MEDIUM effect (moderately separable)")
        separability = "MEDIUM"
    else:
        print(" → SMALL effect (difficult to separate)")
        separability = "LOW"

    mw_stat, mw_p = stats.mannwhitneyu(data1[feature], data2[feature], alternative='two-sided')
    print(f"Mann-Whitney p-value: {mw_p:.2e}")

    # --- accuracy at optimal threshold ---
    y_pred_opt = (feature_values > optimal_threshold).astype(int)
    acc_opt    = accuracy_score(y_true, y_pred_opt)

    # ============================================================
    # VISUALISATIONS
    # ============================================================

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # 1. Distributions
    ax = axes[0, 0]
    ax.hist(data2[feature], bins=40, alpha=0.6,
            label=f'{state2[0]}+{state2[1]} {label2}',
            density=True, color='blue', edgecolor='black', linewidth=0.5)
    ax.hist(data1[feature], bins=40, alpha=0.6,
            label=f'{state1[0]}+{state1[1]} {label1}',
            density=True, color='red', edgecolor='black', linewidth=0.5)
    ax.axvline(optimal_threshold, color='black', linestyle='--', linewidth=2.5,
               label=f'Optimal={optimal_threshold:.3f}')
    ax.axvline(stats1['mean'], color='red',  linestyle=':', linewidth=2, alpha=0.7)
    ax.axvline(stats2['mean'], color='blue', linestyle=':', linewidth=2, alpha=0.7)
    ax.set_xlabel(feature, fontsize=12)
    ax.set_ylabel('Density', fontsize=12)
    ax.set_title(f'CN={total_cn}: {state1} vs {state2}\n{feature} Distribution',
                 fontsize=13, fontweight='bold')
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    textstr = f"Cohen's d = {cohens_d:.3f}\nAUC = {roc_auc:.3f}\nSeparability: {separability}"
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=11,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    # 2. ROC
    ax = axes[0, 1]
    ax.plot(fpr, tpr, linewidth=2.5, label=f'AUC={roc_auc:.4f}', color='darkblue')
    ax.plot([0, 1], [0, 1], 'k--', linewidth=1.5, label='Random')
    ax.scatter(fpr[youden_idx], tpr[youden_idx], s=150, c='red', marker='o',
               edgecolors='black', linewidth=2, zorder=5,
               label=f'Optimal={optimal_threshold:.3f}')
    ax.set_xlabel('False Positive Rate', fontsize=12)
    ax.set_ylabel('True Positive Rate', fontsize=12)
    ax.set_title(f'ROC Curve: {pos_label} Detection', fontsize=13, fontweight='bold')
    ax.legend(fontsize=9, loc='lower right')
    ax.grid(alpha=0.3)
    ax.set_xlim([-0.02, 1.02])
    ax.set_ylim([-0.02, 1.02])

    # 3. Violin
    ax = axes[1, 0]
    parts = ax.violinplot([data2[feature], data1[feature]], positions=[1, 2],
                          widths=0.6, showmeans=True, showmedians=True)
    for pc, col in zip(parts['bodies'], ['blue', 'red']):
        pc.set_facecolor(col)
        pc.set_alpha(0.6)
        pc.set_edgecolor('black')
        pc.set_linewidth(1.5)
    ax.axhline(optimal_threshold, color='black', linestyle='--', linewidth=2.5,
               label=f'Threshold={optimal_threshold:.3f}')
    ax.set_xticks([1, 2])
    ax.set_xticklabels(
        [f'{state2[0]}+{state2[1]}\n{label2}', f'{state1[0]}+{state1[1]}\n{label1}'],
        fontsize=10)
    ax.set_ylabel(feature, fontsize=12)
    ax.set_title('Distribution Comparison', fontsize=13, fontweight='bold')
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3, axis='y')
    for pos, n in zip([1, 2], [len(data2), len(data1)]):
        ax.text(pos, -0.05, f'n={n:,}', ha='center', fontsize=9,
                transform=ax.get_xaxis_transform())

    # 4. Performance vs threshold
    ax = axes[1, 1]
    thresh_range = np.linspace(
        max(0.70, feature_values.min()),
        min(0.98, feature_values.max()),
        100)
    accs, precs, recs = [], [], []
    for t in thresh_range:
        yp = (feature_values > t).astype(int)
        accs.append(accuracy_score(y_true, yp))
        tp = ((yp == 1) & (y_true == 1)).sum()
        fp = ((yp == 1) & (y_true == 0)).sum()
        fn = ((yp == 0) & (y_true == 1)).sum()
        precs.append(tp / (tp + fp) if (tp + fp) > 0 else 0)
        recs.append(tp / (tp + fn) if (tp + fn) > 0 else 0)
    ax.plot(thresh_range, accs,  linewidth=2.5, label='Accuracy',           color='green')
    ax.plot(thresh_range, precs, linewidth=2.5, label=f'{pos_label} Precision', color='blue')
    ax.plot(thresh_range, recs,  linewidth=2.5, label=f'{pos_label} Recall',    color='red')
    ax.axvline(optimal_threshold, color='black', linestyle='--', linewidth=2.5,
               alpha=0.7, label=f'Optimal={optimal_threshold:.3f}')
    ax.set_xlabel('Threshold Value', fontsize=12)
    ax.set_ylabel('Performance Metric', fontsize=12)
    ax.set_title('Performance vs Threshold', fontsize=13, fontweight='bold')
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    ax.set_ylim([0, 1.05])

    plt.tight_layout()
    out_png = os.path.join(
        outdir,
        f'threshold_{comparison_name}_{state1[0]}{state1[1]}_vs_{state2[0]}{state2[1]}.png')
    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {out_png}")
    plt.close()

    # confusion matrix
    fig2, ax2 = plt.subplots(figsize=(8, 6))
    cm = confusion_matrix(y_true, y_pred_opt)
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', ax=ax2,
                xticklabels=[neg_label, pos_label],
                yticklabels=[neg_label, pos_label],
                cbar_kws={'label': 'Count'},
                annot_kws={'fontsize': 14, 'fontweight': 'bold'})
    ax2.set_xlabel('Predicted', fontsize=13)
    ax2.set_ylabel('True', fontsize=13)
    ax2.set_title(
        f'CN={total_cn}: {state1} vs {state2}\n'
        f'Threshold={optimal_threshold:.3f}, Acc={acc_opt * 100:.1f}%',
        fontsize=14, fontweight='bold')
    for i in range(2):
        for j in range(2):
            row_total = cm[i, :].sum()
            pct = 100 * cm[i, j] / row_total if row_total > 0 else 0
            ax2.text(j + 0.5, i + 0.7, f'({pct:.1f}%)',
                     ha='center', va='center', fontsize=11, color='darkred')
    plt.tight_layout()
    out_cm = os.path.join(
        outdir,
        f'confusion_{comparison_name}_{state1[0]}{state1[1]}_vs_{state2[0]}{state2[1]}.png')
    plt.savefig(out_cm, dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {out_cm}")
    plt.close()

    return {
        'total_cn':           total_cn,
        'comparison':         f'{state1} vs {state2}',
        'state1':             f'{state1[0]}+{state1[1]}',
        'state2':             f'{state2[0]}+{state2[1]}',
        'label1':             label1,
        'label2':             label2,
        'feature':            feature,
        'optimal_threshold':  optimal_threshold,
        'auc':                roc_auc,
        'cohens_d':           cohens_d,
        'separability':       separability,
        'accuracy':           acc_opt,
        'sensitivity':        tpr[youden_idx],
        'specificity':        1 - fpr[youden_idx],
        'n_state1':           len(data1),
        'n_state2':           len(data2),
        'n_total':            len(feature_values),
    }

# ============================================================
# DETERMINE COMPARISONS TO RUN
# ============================================================

def determine_comparisons(cn_states_dict, min_samples=50):
    """
    For each CN, build the list of binary comparisons that have enough data.

    Priority order: LoH vs Balanced → LoH vs AI → AI vs Balanced.
    A state is included only if it has >= min_samples segments.

    Returns list of (total_cn, state_pairs, comparison_name).
    """
    comparisons = []

    for total_cn, info in cn_states_dict.items():
        if info['total_segments'] < min_samples * 2:
            print(f"\n⚠️  CN={total_cn}: only {info['total_segments']} segments, skipping")
            continue

        def _best(lst):
            candidates = [(s, c) for s, c in lst if c >= min_samples]
            if not candidates:
                return None
            return max(candidates, key=lambda x: x[1])[0]

        loh = _best(info['loh'])
        ai  = _best(info['ai'])
        bal = _best(info['balanced'])

        if loh and bal:
            comparisons.append((total_cn, [loh, bal], f'cn{total_cn}_loh_vs_bal'))
        if loh and ai:
            comparisons.append((total_cn, [loh, ai],  f'cn{total_cn}_loh_vs_ai'))
        if ai and bal:
            comparisons.append((total_cn, [ai, bal],  f'cn{total_cn}_ai_vs_bal'))

        # fallback: any 2 states from DIFFERENT classes when none of the above triggered.
        # Same-class comparisons (e.g. AI_vs_AI: 5:2 vs 4:3) are excluded because
        # outer_mass cannot distinguish different AI sub-states.
        if not any(c[0] == total_cn for c in comparisons):
            def _class(state):
                nmaj, nmin = state
                if nmin == 0:
                    return 'loh'
                if nmaj == nmin:
                    return 'balanced'
                return 'ai'

            all_states = [
                (s, c, _class(s))
                for lst in [info['loh'], info['ai'], info['balanced']]
                for s, c in lst if c >= min_samples
            ]
            all_states.sort(key=lambda x: x[1], reverse=True)
            chosen = None
            for i in range(len(all_states)):
                for j in range(i + 1, len(all_states)):
                    if all_states[i][2] != all_states[j][2]:
                        chosen = (all_states[i][0], all_states[j][0])
                        break
                if chosen:
                    break
            if chosen:
                s1, s2 = chosen
                comparisons.append((
                    total_cn, [s1, s2],
                    f'cn{total_cn}_{s1[0]}{s1[1]}_vs_{s2[0]}{s2[1]}'
                ))
            else:
                print(f"\n⚠️  CN={total_cn}: only same-class states available — no valid comparison")

    return comparisons

# ============================================================
# MAIN
# ============================================================

# ============================================================
# PER-DP CALIBRATION FUNCTION
# ============================================================

def run_calibration_for_dp(df_full, dp, feature, outdir_dp, min_samples):
    """
    Run full threshold calibration for one DP filter.

    Selects all rows at dp_filter == dp that carry ASCAT ground-truth
    (nMajor and nMinor non-null), regardless of which baf_folder they
    came from.  Using all available backbones maximises segment counts
    and keeps calibration data at exactly the same BAF-depth cut as
    will be applied during classification.

    Returns
    -------
    threshold_lookup : dict  {cn: {comparison_key: threshold}}
    summary_df       : pd.DataFrame  one row per comparison (or None if no results)
    """
    print("\n" + "█" * 70)
    print(f"  DP FILTER = {dp}")
    print("█" * 70)

    df_filtered = df_full[
        (df_full['dp_filter'] == dp) &
        df_full['nMajor'].notna() &
        df_full['nMinor'].notna()
    ].copy()

    n_folders = df_filtered['baf_folder'].nunique() if 'baf_folder' in df_filtered.columns else '?'
    print(f"✓ Segments: {len(df_filtered):,}  (from {n_folders} folder(s) with ground truth)")
    if len(df_filtered) == 0:
        print(f"  ⚠️  No ground-truth data for dp={dp} — skipping")
        return {}, None

    os.makedirs(outdir_dp, exist_ok=True)

    cn_states_dict    = summarize_cn_distribution(df_filtered)
    comparisons_to_run = determine_comparisons(cn_states_dict, min_samples=min_samples)

    print(f"\nWill run {len(comparisons_to_run)} comparisons:")
    for total_cn, state_pairs, name in comparisons_to_run:
        print(f"  CN={total_cn}: {state_pairs[0]} vs {state_pairs[1]} ({name})")

    results = []
    for total_cn, state_pairs, comparison_name in comparisons_to_run:
        print("\n" + "-" * 70)
        print(f"CN={total_cn} — {comparison_name}")
        print("-" * 70)
        result = analyze_cn_state(
            df_filtered,
            total_cn=total_cn,
            state_pairs=state_pairs,
            feature=feature,
            comparison_name=comparison_name,
            outdir=outdir_dp,
        )
        if result:
            results.append(result)

    if not results:
        print(f"  ⚠️  No comparisons passed for dp={dp}")
        return {}, None

    col_order = [
        'total_cn', 'comparison', 'state1', 'state2', 'label1', 'label2',
        'n_state1', 'n_state2', 'feature', 'optimal_threshold', 'accuracy',
        'sensitivity', 'specificity', 'auc', 'cohens_d', 'separability',
    ]
    summary_df = pd.DataFrame(results)[col_order]
    summary_df['dp_filter'] = dp

    summary_file = os.path.join(outdir_dp, 'threshold_summary_complete.csv')
    summary_df.to_csv(summary_file, index=False)
    print(f"\n✓ Summary saved: {summary_file}")

    # Build threshold lookup for this DP
    threshold_lookup = {}
    for _, row in summary_df.iterrows():
        cn   = int(row['total_cn'])
        comp = str(row['comparison']).lower()
        l1   = row['label1'].lower()
        l2   = row['label2'].lower()

        if 'loh_vs_bal' in comp or ('loh' in (l1, l2) and 'balanced' in (l1, l2)):
            key = 'loh_vs_balanced'
        elif 'loh_vs_ai' in comp or ('loh' in (l1, l2) and 'ai' in (l1, l2)):
            key = 'loh_vs_ai'
        elif 'ai_vs_bal' in comp or ('ai' in (l1, l2) and 'balanced' in (l1, l2)):
            key = 'ai_vs_balanced'
        else:
            key = f"{row['label1']}_vs_{row['label2']}"

        threshold_lookup.setdefault(cn, {})[key] = float(row['optimal_threshold'])

    # Per-DP separability plot
    labels_plot = [
        f"CN={r['total_cn']}\n{r['state1']} vs {r['state2']}" for r in results
    ]
    x     = np.arange(len(labels_plot))
    width = 0.35
    fig, ax = plt.subplots(figsize=(14, 7))
    bars1 = ax.bar(x - width / 2, [r['cohens_d'] for r in results], width,
                   label="Cohen's d", color='steelblue', edgecolor='black', linewidth=1.5)
    bars2 = ax.bar(x + width / 2, [r['auc']      for r in results], width,
                   label='AUC-ROC',   color='coral',     edgecolor='black', linewidth=1.5)
    ax.axhline(0.8, color='green',  linestyle='--', linewidth=2, alpha=0.7,
               label='Large effect (d>0.8)')
    ax.axhline(0.5, color='orange', linestyle='--', linewidth=2, alpha=0.7,
               label='Medium effect (d>0.5)')
    ax.set_xlabel('Comparison', fontsize=13)
    ax.set_ylabel('Value', fontsize=13)
    ax.set_title(f'DP≥{dp} — Separability: Cohen\'s d & AUC-ROC',
                 fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(labels_plot, rotation=45, ha='right', fontsize=9)
    ax.legend(fontsize=10)
    ax.grid(alpha=0.3, axis='y')
    ax.set_ylim([0, 1.1])
    for bars in [bars1, bars2]:
        for bar in bars:
            h = bar.get_height()
            if h > 0.05:
                ax.text(bar.get_x() + bar.get_width() / 2., h, f'{h:.2f}',
                        ha='center', va='bottom', fontsize=8, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(outdir_dp, 'separability_summary_all_cn.png'),
                dpi=300, bbox_inches='tight')
    plt.close()

    return threshold_lookup, summary_df


# ============================================================
# MAIN — loop over all DP filters
# ============================================================

print("\n" + "=" * 70)
print("COMPREHENSIVE THRESHOLD CALIBRATION — ALL DP FILTERS")
print("=" * 70)

all_thresholds = {}   # {dp: {cn: {comparison_key: threshold}}}
all_summaries  = []

for dp in PARAMS['dp_filters']:
    outdir_dp = os.path.join(PARAMS['outdir'], f'dp{dp}')
    thresholds_dp, summary_dp = run_calibration_for_dp(
        df_full,
        dp          = dp,
        feature     = PARAMS['feature'],
        outdir_dp   = outdir_dp,
        min_samples = PARAMS['min_samples'],
    )
    if thresholds_dp:
        all_thresholds[dp] = thresholds_dp
    if summary_dp is not None:
        all_summaries.append(summary_dp)

# ============================================================
# COMBINED OUTPUTS
# ============================================================

print("\n" + "=" * 70)
print("COMBINED OUTPUTS")
print("=" * 70)

# Nested JSON: {dp: {cn: {comparison_key: threshold}}}
threshold_file = os.path.join(PARAMS['outdir'], 'thresholds_lookup.json')
with open(threshold_file, 'w') as f:
    json.dump(all_thresholds, f, indent=2)
print(f"\n✓ Nested thresholds saved: {threshold_file}")
print("  Structure: {dp_filter → {cn → {comparison_key → threshold}}}")

print("\nThreshold values by DP and CN:")
print("-" * 70)
for dp in sorted(all_thresholds):
    print(f"\n  DP≥{dp}:")
    for cn in sorted(all_thresholds[dp]):
        for k, v in all_thresholds[dp][cn].items():
            print(f"    CN={cn}  {k:25s}: {v:.4f}")

# Combined summary CSV
if all_summaries:
    combined_summary = pd.concat(all_summaries, ignore_index=True)
    combined_file = os.path.join(PARAMS['outdir'], 'threshold_summary_all_dp.csv')
    combined_summary.to_csv(combined_file, index=False)
    print(f"\n✓ Combined summary saved: {combined_file}")

    # Cross-DP separability overview
    fig, axes = plt.subplots(1, len(PARAMS['dp_filters']),
                             figsize=(6 * len(PARAMS['dp_filters']), 6),
                             sharey=True)
    if len(PARAMS['dp_filters']) == 1:
        axes = [axes]
    for ax, dp in zip(axes, PARAMS['dp_filters']):
        dp_rows = combined_summary[combined_summary['dp_filter'] == dp]
        if len(dp_rows) == 0:
            ax.set_title(f'DP≥{dp}\n(no data)')
            continue
        labels = [f"CN={r['total_cn']}\n{r['state1']}v{r['state2']}"
                  for _, r in dp_rows.iterrows()]
        x = np.arange(len(labels))
        w = 0.35
        ax.bar(x - w/2, dp_rows['cohens_d'].values, w,
               label="Cohen's d", color='steelblue', edgecolor='black', linewidth=1)
        ax.bar(x + w/2, dp_rows['auc'].values, w,
               label='AUC-ROC',   color='coral',     edgecolor='black', linewidth=1)
        ax.axhline(0.8, color='green',  linestyle='--', linewidth=1.5, alpha=0.7)
        ax.axhline(0.5, color='orange', linestyle='--', linewidth=1.5, alpha=0.7)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=7)
        ax.set_title(f'DP≥{dp}', fontsize=13, fontweight='bold')
        ax.set_ylim([0, 1.1])
        ax.grid(alpha=0.3, axis='y')
        if ax == axes[0]:
            ax.set_ylabel('Value')
            ax.legend(fontsize=8)
    plt.suptitle('Separability Across DP Filters', fontsize=15, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(PARAMS['outdir'], 'separability_cross_dp.png'),
                dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✓ Cross-DP plot saved: {os.path.join(PARAMS['outdir'], 'separability_cross_dp.png')}")

print("\n" + "=" * 70)
print("CALIBRATION COMPLETE")
print("=" * 70)

print("\n" + "=" * 70)
print("ALL ANALYSES COMPLETE")
print("=" * 70)
print(f"\nGenerated files in: {PARAMS['outdir']}/")
print(f"  threshold_*.png                 — distribution + ROC per comparison")
print(f"  confusion_*.png                 — confusion matrix per comparison")
print(f"  threshold_summary_complete.csv  — one row per comparison")
print(f"  thresholds_lookup.json          — consumed by baf_to_hrd_pipeline.py Stage 3")
print(f"  separability_summary_all_cn.png — overview across all CN values")
