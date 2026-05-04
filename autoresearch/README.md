# autoresearch_baf_hrd

Autonomous research setup for improving the allele-specific copy-number classification
step of the BAF-to-HRD pipeline, inspired by
[karpathy/autoresearch](https://github.com/karpathy/autoresearch).

---

## The problem

TSO500 is a tumor-only targeted panel. It produces copy-number segments with a
**total CN** per segment but no allele breakdown — there is no matched germline to
distinguish which allele carries extra or missing copies.

scarHRD (which computes HRD scores for clinical use) requires **allele-specific CN**:
how many copies of each allele (A_cn, B_cn) are present per segment.

This repo optimises the threshold-based classifier that infers (A_cn, B_cn) from the
shape of the B-allele frequency (BAF) distribution within each segment.

---

## Key signal: `outer_mass_raw`

For each segment, the pipeline counts the fraction of heterozygous SNP positions whose
BAF falls outside a central interval — in the tails near 0 and 1.  This fraction is
the primary indicator of Loss of Heterozygosity:

```
outer_mass high  →  BAF piles up at 0 and 1  →  one allele gone  →  LoH
outer_mass low   →  BAF spreads or sits at 0.5  →  both alleles  →  Balanced / AI
```

**Two variants exist — this loop uses `outer_mass_raw`:**

| Column | BAF scale | When to use |
|--------|-----------|-------------|
| `outer_mass_raw` | Uncorrected raw BAF | Production pipeline (TSO500, no purity model) — **use this** |
| `outer_mass` | Purity-corrected BAF | ASCAT samples only; not available for TSO500 |

`experiment.py` always uses `outer_mass_raw` via `recompute_outer_mass()`, which reads
the `ratio_raw_*` fine bins. The 20 raw fine bins (5% intervals) allow agents to
recompute with any boundary they choose.

---

## Metric

**`val_allelic_acc`** — fraction of segments where predicted (A_cn, B_cn) exactly
matches the ASCAT ground truth (or its mirror image).

Ground truth comes from samples that were also profiled with ASCAT on matched
tumor-normal data, giving unambiguous allele-specific CN per segment.

---

## Files

| File | Purpose | Modify? |
|------|---------|---------|
| `experiment.py` | Thresholds, boundaries, classify function | **YES — agents only modify this** |
| `prepare.py` | Builds `val_data.csv` from pipeline output | No |
| `evaluate.py` | Leaderboard and progress plot | No |
| `baf_to_hrd_pipeline.py` | Full production pipeline (source) | No |
| `program.md` | Agent instructions (read this first) | No |
| `run_prepare.sbatch` | SLURM job for prepare.py — submit once | No |
| `run_experiment.sbatch` | SLURM job for experiment.py — submit after each edit | No |
| `val_data.csv` | Validation set (created by prepare.py) | No |
| `results/log.jsonl` | Append-only experiment log | No |
| `logs/` | SLURM stdout/stderr for each job | No |

---

## Quickstart

> **HPC note:** this cluster uses SLURM. Do NOT run Python scripts directly in
> the login shell — use `sbatch` for compute jobs.

### 1. Run the full pipeline to generate Stage 2 output

The merged CSV for the ASCAT-downsampled samples (20260406 run) already exists:
```
/storage/scratch01/groups/co/cn_extra/alleleSpecific/tso_samples/20260406/
  merged_segments_baf_metrics_ALL_with20bins_enhancedFeatures.csv
```
No action needed unless you re-run Stage 1/2.

### 2. Build the validation set (once)

```bash
sbatch autoresearch_baf_hrd/run_prepare.sbatch
```

This creates `autoresearch_baf_hrd/val_data.csv` with segments from ASCAT-profiled
samples (nMajor/nMinor present) at DP=8, including 20 corrected + 20 raw BAF bin
columns and ground-truth class labels.

### 3. Run the baseline experiment

```bash
sbatch autoresearch_baf_hrd/run_experiment.sbatch
```

Check the log: `autoresearch_baf_hrd/logs/experiment_<jobid>.log`
```
val_allelic_acc : 0.XXXX  (N / M segments)
val_class_acc   : 0.XXXX
✓ Result appended to results/log.jsonl
```

### 4. Inspect results

```bash
conda run -n shapeit4 python autoresearch_baf_hrd/evaluate.py            # full leaderboard
conda run -n shapeit4 python autoresearch_baf_hrd/evaluate.py --best     # best run details
conda run -n shapeit4 python autoresearch_baf_hrd/evaluate.py --last 5   # last 5 runs
conda run -n shapeit4 python autoresearch_baf_hrd/evaluate.py --plot     # save results/progress.png
```

### 5. Iterate

```
1. Edit experiment.py  (THRESHOLDS, classify_segment, OUTER_LOWER/UPPER, etc.)
2. sbatch autoresearch_baf_hrd/run_experiment.sbatch
3. Check logs/ and results/log.jsonl
4. Repeat
```

### 6. Ask an agent to improve it

Point an agent (Claude Code, or any LLM-based agent) at `program.md`:

```
Read program.md and then autonomously improve experiment.py to maximise
val_allelic_acc.  Submit via sbatch after each change and check the result.
Only modify experiment.py.
```

---

## What agents can try

The following directions are all fair game inside `experiment.py`:

- **Threshold tuning** — adjust the `outer_mass` cutoff for each CN level
- **Boundary tuning** — change `OUTER_LOWER` / `OUTER_UPPER` (e.g. 0.2/0.8 instead
  of 0.25/0.75) to redefine what counts as "outer"
- **Multi-feature rules** — use `baf_sd`, `baf_mean`, `baf_iqr` in addition to
  `outer_mass` (all available in `val_data.csv`)
- **Per-CN logic** — different decision rules for different CN levels
- **Raw vs corrected BAF** — the `ratio_raw_*` columns reflect uncorrected BAF;
  mixing them with corrected values may help for low-purity samples
- **Fine-bin features** — aggregate individual bins differently (e.g. inner 5 bins
  vs outer 5 bins vs specific peaks)

---

## Relationship to the full pipeline

`baf_to_hrd_pipeline.py` (in the parent directory) is the production code.
The classification logic optimised here corresponds exactly to Stage 3
(`run_script3` / `build_scarHRD_input_using_thresholds`).

Once a better configuration is found in `experiment.py`, the improved thresholds
and `classify_segment` logic can be ported back into `baf_to_hrd_pipeline.py`
(specifically into `build_scarHRD_input_using_thresholds` and the `THRESHOLDS`
parameter block).

---

## Dependencies

Same as `baf_to_hrd_pipeline.py`: `numpy`, `pandas`, `matplotlib`.
No GPU required — all experiments run on CPU in seconds.
