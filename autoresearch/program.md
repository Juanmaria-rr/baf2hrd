# HRD Pipeline Autonomous Research Program

## What you are doing

You are an autonomous research agent improving a bioinformatics pipeline that infers
**allele-specific copy numbers (A_cn, B_cn)** from tumor-only sequencing data (TSO500),
with the goal of producing accurate inputs for the scarHRD HRD scoring tool.

The pipeline processes segments of the tumor genome. For each segment it knows the
**total copy number** (how many total DNA copies are present) but NOT which allele
carries how many copies — that is what you must infer.

Your job is to improve `experiment.py` so that the predicted (A_cn, B_cn) pairs
match the ASCAT ground-truth as accurately as possible.

---

## The metric

**`val_allelic_acc`** — the fraction of segments where your predicted (A_cn, B_cn)
exactly matches the ground truth from ASCAT-profiled samples.
An exact match means A_pred==A_true and B_pred==B_true (or the flipped pair, since
major/minor labelling is arbitrary).

Secondary metric: **`val_class_acc`** — fraction with the correct allelic state class
(LoH, Balanced, AI, Deletion) regardless of exact numbers.

**Higher is better for both metrics.**

---

## The key signal: outer_mass_raw

The only observable evidence for the allelic state is the **B-allele frequency (BAF)**
at heterozygous SNP positions within a segment.

`outer_mass_raw` summarises this: it is the fraction of SNP positions whose **raw**
(uncorrected) BAF value falls in the extreme tails of [0, 1] — specifically outside the
[OUTER_LOWER, OUTER_UPPER] interval (default: [0.25, 0.75]).

- **High outer_mass_raw** → BAF clusters near 0 and 1 → one allele is lost or vastly
  amplified → **Loss of Heterozygosity (LoH)**
- **Low outer_mass_raw** → BAF is spread or clusters near 0.5 → both alleles present →
  **Balanced or Allelic Imbalance (AI)**

**Always use `outer_mass_raw`, not `outer_mass`.**  The purity-corrected `outer_mass`
is only valid for ASCAT samples (matched tumor-normal); TSO500 samples have no purity
model and would use raw BAF at inference time.

The 20 fine raw BAF bins (0.05-wide each, columns `ratio_raw_0.00-0.05` to
`ratio_raw_0.95-1.00`) are available in `val_data.csv`. You can recompute
`outer_mass_raw` with any boundary by changing `OUTER_LOWER` and `OUTER_UPPER`
in `experiment.py` — `recompute_outer_mass()` always reads the `ratio_raw_*` columns.

---

## What you can modify

**Only `experiment.py`.**  Do not touch any other file.

Inside `experiment.py` you may change:

### 1. `THRESHOLDS` (dict)
The `outer_mass_raw` cutoff for each total CN level.  Current baseline values
are calibrated by Youden's index on DP=8 segments from both JBLAB231 and
JBLAB17041 backbones (201,218 segments, re-run 2026-04-29).
```python
THRESHOLDS = {
    "2": {"loh_vs_balanced": 0.8537},  # CN=2: above → LoH (2:0), below → Balanced (1:1)
    "3": {"loh_vs_ai":       0.8824},  # CN=3: above → LoH (3:0), below → AI (2:1)
    "4": {"loh_vs_balanced": 0.8261,   # CN=4: three states, two thresholds
          "loh_vs_ai":       0.8571,
          "ai_vs_balanced":  0.8333},
    "5": {"loh_vs_ai":       0.8444},  # CN=5: above → LoH (5:0), below → AI (3:2)
    "6": {"loh_vs_balanced": 1.0000,   # CN=6: loh thresholds = 1.0 (too few LoH segments at DP=8)
          "loh_vs_ai":       1.0000,   #        → only ai_vs_balanced is operative
          "ai_vs_balanced":  0.8400},
    "8": {"ai_vs_balanced":  0.8276},  # CN=8: above → AI (6:2), below → Balanced (4:4)
}
```
Note: CN=6 LoH thresholds of 1.0 mean LoH is never predicted at CN=6 via threshold alone;
the Stage 3 heuristic handles edge cases.

### 2. `OUTER_LOWER` / `OUTER_UPPER` (float, 0–1)
The BAF interval boundary that defines "outer".  Wider → more BAF signal counted as
outer; narrower → more selective.  Default: 0.25 / 0.75.

### 3. `HEURISTIC_THRESHOLD` (float)
The outer_mass cutoff for CN=7 and CN>8.  Default: 0.95.

### 4. `classify_segment(total_cn, outer_mass)` (function)
The entire classification function. `outer_mass` here is the value produced by
`recompute_outer_mass()` — it is always on the **raw BAF scale**. You may:
- Add more thresholds (e.g. per-CN boundaries for AI sub-states)
- Use additional features available in `val_data.csv`:
  - `baf_mean`, `baf_sd`, `baf_median`, `baf_iqr` — BAF summary statistics
  - `ratio_raw_*` columns — raw (uncorrected) BAF bin ratios (20 bins)
  - `ratio_*` columns — purity-corrected BAF bin ratios (for reference)
  - Individual fine bins for finer signal decomposition
- Implement a completely different approach (rule-based, statistical, etc.)

### 5. `EXPERIMENT_NOTES` (string)
A short description of what this run changes.  This is written to the log.

---

## Workflow

This pipeline runs on a shared HPC cluster (SLURM). Do NOT run Python scripts
directly in the login shell — submit via sbatch instead.

```
1. Edit experiment.py
2. sbatch run_experiment.sbatch          ← submits to SLURM; result in logs/experiment_<jobid>.log
3. conda run -n shapeit4 python evaluate.py          ← inspect leaderboard across all runs
4. conda run -n shapeit4 python evaluate.py --best   ← see the best configuration so far
5. Repeat
```

---

## Useful context

### Copy number states and their typical outer_mass ranges

| Total CN | Common states | Typical outer_mass |
|----------|--------------|-------------------|
| 0        | 0:0 Deletion  | undefined (no BAF) |
| 1        | 1:0 LoH       | high (~0.95–1.0) |
| 2        | 2:0 LoH or 1:1 Balanced | LoH: >0.90; Balanced: <0.30 |
| 3        | 3:0 LoH or 2:1 AI | LoH: >0.85; AI: 0.55–0.75 |
| 4        | 4:0 LoH, 3:1 AI, 2:2 Balanced | LoH: >0.95; AI: 0.70–0.90; Balanced: <0.40 |
| ≥5       | increasingly uncertain | higher CN → wider overlap between states |

### Segments without BAF data (NaN outer_mass)
These receive a default assignment. Improving the default fallback for each CN level
can help on segments with low coverage.

### What the ASCAT ground truth represents
ASCAT was run on matched tumor-normal whole genome or exome data — it has access to
germline heterozygous SNPs and can unambiguously phase copy numbers.  The TSO500 data
is tumor-only and panel-based: we can only observe aggregate BAF, not phased alleles.
The goal is to recover as much of the ASCAT-level accuracy as possible from that
weaker signal.

### Confidence of different CN levels
- CN=2 is the most common and most important.  A 1% gain in CN=2 accuracy dominates
  the overall metric.
- CN=3 and CN=4 are the next most impactful.
- High CN (≥7) segments are rare; gains there have little effect on val_allelic_acc
  but matter biologically for HRD scoring.

---

## Constraints

- Modify only `experiment.py`.
- Do not hardcode val_data.csv labels into your logic (that would be overfitting).
- Each experiment should run to completion (no timeouts needed; classification is fast).
- Log every run by keeping the `run_evaluation()` call at the end of the script.

---

## Files

| File | Role |
|------|------|
| `experiment.py` | **YOU MODIFY THIS** — thresholds, boundaries, classify function |
| `prepare.py` | Builds `val_data.csv` — run once, do not modify |
| `evaluate.py` | Leaderboard and plots — do not modify |
| `baf_to_hrd_pipeline.py` | Full production pipeline — do not modify |
| `run_prepare.sbatch` | SLURM submit script for prepare.py — run once |
| `run_experiment.sbatch` | SLURM submit script for experiment.py — run after each edit |
| `val_data.csv` | Validation segments with ground truth — do not modify |
| `results/log.jsonl` | Append-only experiment log |
| `logs/` | SLURM stdout/stderr (experiment_\<jobid\>.log / .err) |
