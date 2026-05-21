#!/usr/bin/env python3
"""
log_processed_samples.py

Scan all pipeline run directories from 2026-05-11 onwards and generate
PROCESSED_SAMPLES.tsv in the repo root.

One row per (run_date, sample_id) at the representative dp
(dp4 postmerge preferred, falling back to dp2 / dp0 / dp8).

Run manually at any time:
    conda run -n shapeit4 python pipeline/log_processed_samples.py
"""

import argparse
import csv
import json
import subprocess
from pathlib import Path

PIPELINE_REPO = Path(__file__).resolve().parent.parent

def _parse_args():
    p = argparse.ArgumentParser(description="Generate PROCESSED_SAMPLES.tsv from pipeline run dirs.")
    p.add_argument(
        "--tso-root",
        type=Path,
        default=PIPELINE_REPO.parent / "tso_samples",
        help="Root directory containing all pipeline run subdirectories (default: ../tso_samples relative to repo root)",
    )
    return p.parse_args()

_ARGS = _parse_args()
TSO_SAMPLES_ROOT = _ARGS.tso_root
OUTPUT_TSV = PIPELINE_REPO / "PROCESSED_SAMPLES.tsv"
MIN_RUN_DATE = "20260511"
DP_PRIORITY = [4, 2, 0, 8]


def _git_commit_at_date(date_str: str) -> str:
    date_end = f"{date_str[:4]}-{date_str[4:6]}-{date_str[6:8]}T23:59:59"
    r = subprocess.run(
        ["git", "-C", str(PIPELINE_REPO), "log", "--format=%H", f"--before={date_end}", "-1"],
        capture_output=True, text=True,
    )
    return r.stdout.strip() or "unknown"


def _find_config(date_str: str) -> str:
    p = PIPELINE_REPO / "configs" / f"pipeline_config_{date_str}.yaml"
    return str(p) if p.exists() else "unknown"


def _load_run_info(run_date: str) -> tuple:
    """Return (git_commit, config_file) for the given run date.

    Reads configs/run_info_YYYYMMDD.json if present (written by orchestrator
    for runs from this pipeline); otherwise reconstructs from git history.
    """
    info_path = PIPELINE_REPO / "configs" / f"run_info_{run_date}.json"
    if info_path.exists():
        with open(info_path) as f:
            info = json.load(f)
        return info.get("git_commit", "unknown"), info.get("config_file", "unknown")
    return _git_commit_at_date(run_date), _find_config(run_date)


def _parse_scores(scores_path: Path) -> dict:
    """Return {sample_name: score_dict} at the best available dp (postmerge).

    sample_name preserves the -HRD suffix when present so base and HRD BAMs
    appear as separate rows.
    """
    by_sample: dict[str, dict] = {}

    with open(scores_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            sid = row["SampleID"].strip('"')
            if "_results_dp" not in sid:
                continue
            sample_name, rest = sid.split("_results_dp", 1)
            dp_str, merge_type = rest.split("_", 1)
            if merge_type != "postmerge":
                continue
            dp = int(dp_str)
            by_sample.setdefault(sample_name, {})[dp] = {
                "hrd":     row["HRD"].strip('"'),
                "tai":     row["Telomeric AI"].strip('"'),
                "lst":     row["LST"].strip('"'),
                "hrd_sum": row["HRD-sum"].strip('"'),
                "dp_used": dp,
            }

    result = {}
    for sample_name, dp_scores in by_sample.items():
        for dp in DP_PRIORITY:
            if dp in dp_scores:
                result[sample_name] = dp_scores[dp]
                break
    return result


def main():
    rows = []

    for run_dir in sorted(TSO_SAMPLES_ROOT.iterdir()):
        if not run_dir.is_dir():
            continue
        run_date = run_dir.name
        if not (run_date.isdigit() and len(run_date) == 8):
            continue
        if run_date < MIN_RUN_DATE:
            continue

        scores_path = run_dir / "outputs_scarHRD" / "ALL_scarHRD_scores.csv"
        if not scores_path.exists():
            continue

        git_commit, config_file = _load_run_info(run_date)
        scores = _parse_scores(scores_path)

        for sample_name, score in sorted(scores.items()):
            rows.append({
                "run_date":    run_date,
                "sample_id":   sample_name,
                "git_commit":  git_commit,
                "config_file": config_file,
                "output_dir":  str(run_dir),
                "dp_used":     score["dp_used"],
                "hrd":         score["hrd"],
                "tai":         score["tai"],
                "lst":         score["lst"],
                "hrd_sum":     score["hrd_sum"],
            })

    fields = [
        "run_date", "sample_id", "git_commit", "config_file",
        "output_dir", "dp_used", "hrd", "tai", "lst", "hrd_sum",
    ]
    with open(OUTPUT_TSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"[log_processed_samples] {len(rows)} rows → {OUTPUT_TSV}")


if __name__ == "__main__":
    main()
