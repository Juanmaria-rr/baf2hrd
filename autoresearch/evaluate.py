#!/usr/bin/env python3
"""
evaluate.py — read-only utility for inspecting experiment results.

DO NOT MODIFY — agents should not change this file.

Usage:
    python evaluate.py              # print leaderboard of all runs
    python evaluate.py --best       # print only the best run
    python evaluate.py --last N     # print the last N runs
    python evaluate.py --plot       # save a PNG of val_allelic_acc over time
"""

import os
import sys
import json
import argparse
import numpy as np
import pandas as pd

RESULTS_LOG = os.path.join(os.path.dirname(__file__), "results", "log.jsonl")


def load_results():
    if not os.path.exists(RESULTS_LOG):
        print(f"No results yet. Run experiment.py first.")
        return pd.DataFrame()

    rows = []
    with open(RESULTS_LOG) as f:
        for line in f:
            line = line.strip()
            if line:
                try:
                    rows.append(json.loads(line))
                except json.JSONDecodeError:
                    continue

    if not rows:
        return pd.DataFrame()

    df = pd.DataFrame(rows)
    df["rank"] = df["val_allelic_acc"].rank(ascending=False, method="min").astype(int)
    df = df.sort_values("val_allelic_acc", ascending=False).reset_index(drop=True)
    return df


def print_leaderboard(df, n=None):
    if df.empty:
        return
    show = df.head(n) if n else df
    print(f"\n{'Rank':>4}  {'val_allelic_acc':>15}  {'val_class_acc':>13}  "
          f"{'n_seg':>7}  {'elapsed_s':>9}  {'timestamp':>19}  notes")
    print("-" * 110)
    for _, row in show.iterrows():
        print(
            f"{int(row.get('rank', 0)):>4}  "
            f"{row['val_allelic_acc']:>15.4f}  "
            f"{row.get('val_class_acc', float('nan')):>13.4f}  "
            f"{int(row.get('n_segments', 0)):>7,}  "
            f"{row.get('elapsed_s', float('nan')):>9.1f}  "
            f"{str(row.get('timestamp',''))[:19]:>19}  "
            f"{str(row.get('notes',''))[:50]}"
        )
    print()


def print_best(df):
    if df.empty:
        return
    best = df.iloc[0]
    print("\n=== BEST RUN ===")
    print(f"  val_allelic_acc : {best['val_allelic_acc']:.4f}")
    print(f"  val_class_acc   : {best.get('val_class_acc', float('nan')):.4f}")
    print(f"  timestamp       : {best.get('timestamp','')}")
    print(f"  notes           : {best.get('notes','')}")
    print("\n  thresholds:")
    for cn, vals in best.get("thresholds", {}).items():
        print(f"    CN={cn}: {vals}")
    ob = best.get("outer_boundary", {})
    print(f"\n  outer_mass boundary: [{ob.get('OUTER_LOWER','?')}, {ob.get('OUTER_UPPER','?')}]")
    print(f"  heuristic_threshold: {best.get('heuristic_threshold','?')}")

    print("\n  per-CN allelic accuracy:")
    for cn, acc in sorted(best.get("per_cn_allelic_acc", {}).items(), key=lambda x: int(x[0])):
        print(f"    CN={cn}: {acc:.3f}")


def plot_progress(df):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available — skipping plot")
        return

    if df.empty:
        return

    df_time = df.sort_values("timestamp").reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(df_time.index + 1, df_time["val_allelic_acc"], marker="o", linewidth=1.5)
    ax.axhline(df_time["val_allelic_acc"].max(), color="red", linestyle="--",
               alpha=0.5, label=f"best={df_time['val_allelic_acc'].max():.4f}")
    ax.set_xlabel("Experiment #")
    ax.set_ylabel("val_allelic_acc")
    ax.set_title("HRD Pipeline — Allelic Accuracy over Experiments")
    ax.legend()
    ax.grid(alpha=0.3)
    plt.tight_layout()

    out = os.path.join(os.path.dirname(__file__), "results", "progress.png")
    plt.savefig(out, dpi=150)
    plt.close()
    print(f"✓ Saved plot: {out}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--best",      action="store_true")
    parser.add_argument("--last",      type=int, default=None)
    parser.add_argument("--plot",      action="store_true")
    args = parser.parse_args()

    df = load_results()
    if df.empty:
        return

    if args.best:
        print_best(df)
    elif args.plot:
        plot_progress(df)
    else:
        n = args.last if args.last else None
        if args.last:
            df_show = df.sort_values("timestamp").tail(args.last).copy()
            df_show["rank"] = df_show["val_allelic_acc"].rank(ascending=False, method="min").astype(int)
        else:
            df_show = df
        print_leaderboard(df_show, n=None)
        print(f"Total experiments: {len(df)}")
        print(f"Best val_allelic_acc: {df['val_allelic_acc'].max():.4f}")
        print(f"Latest val_allelic_acc: {df.sort_values('timestamp').iloc[-1]['val_allelic_acc']:.4f}")


if __name__ == "__main__":
    main()
