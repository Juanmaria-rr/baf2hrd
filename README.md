# baf2hrd

> **Work in progress** — pipeline under active development.

End-to-end pipeline for inferring allele-specific copy numbers (A_cn, B_cn)
from TSO500 tumor-only sequencing data, producing inputs for
[scarHRD](https://github.com/sztup/scarHRD) HRD scoring.

## Overview

```
BAM files
   └── Pre-processing: split combined segtable → per-sample CSVs
   └── Stage 0: mpileup → _qc.tsv
   └── Stage 1: per-sample BAF analysis → segments + manifest
   └── Stage 2: consolidation → merged dataset
   └── Stage 3: allele-specific CN inference → scarHRD inputs
   │             ↑ optimised by autoresearch/
   └── Stage 4: scarHRD → HRD scores
   └── QC report → {base_dir}/qc_report.html
   └── HRD performance report → {tso_samples_root}/hrd_performance_report.html
```

## Structure

```
pipeline/     Python source (main pipeline + utilities)
  split_segments.py        ← pre-processing: combined segtable → per-sample CSVs
  qc_report.py             ← post-pipeline QC report (HTML autocontenido)
  log_processed_samples.py ← genera PROCESSED_SAMPLES.tsv con trazabilidad de runs
  eval_hrd_performance.py  ← evaluación de performance vs GIS gold standard (HTML)
slurm/        SLURM submission scripts
  orchestrate_pipeline.sh  ← entry point; encadena Stage 0 → 4 → QC report
  run_index_bam.sbatch     ← indexado de BAM (lanzado automáticamente si .bai ausente/vacío)
  run_qc_report.sbatch     ← lanzamiento ad hoc del QC report
  run_eval_hrd.sbatch      ← lanzamiento ad hoc del informe de evaluación vs GIS
configs/      Per-run YAML configuration (auto-generated per run)
  run_info_YYYYMMDD.json   ← trazabilidad: git_commit + config_file por run (auto-generated)
autoresearch/ Autonomous loop for improving Stage 3 classification
PROCESSED_SAMPLES.tsv      ← log de muestras procesadas (generado por log_processed_samples.py)
```

## Usage

### Full run from BAMs (recommended)

```bash
bash slurm/orchestrate_pipeline.sh \
  --bam-dir /path/to/rg_bam_YYYYMMDD_batchNN \
  --tsv-dir /path/to/YYYYMMDD_batchNN_mpileup/results
```

When purity data is available, pass `--purity-file` to launch the base and
purity-corrected runs **in parallel** from the same Stage 0 dependency:

```bash
bash slurm/orchestrate_pipeline.sh \
  --bam-dir  /path/to/all_bam_tso \
  --tsv-dir  /path/to/YYYYMMDD_mpileup_results \
  --segments-dir /path/to/YYYYMMDD_segments \
  --purity-file /path/to/purity_ploidy_pipeline.tsv
```

This submits two independent pipeline chains simultaneously:
- `YYYYMMDD/` — base run (raw BAF, raw thresholds)
- `YYYYMMDD_purity/` — purity-corrected BAF + purity-calibrated thresholds

The orchestrator will:
1. Run `split_segments.py` to generate per-sample segment CSVs from the combined segtable
2. Auto-generate `configs/pipeline_config_YYYYMMDD.yaml` (and `_purity.yaml` if `--purity-file` set)
3. Write `configs/run_info_YYYYMMDD.json` (git commit hash + config path) for traceability
4. **Pre-flight**: validate `.bai` index for every BAM — if absent or empty, submit `run_index_bam.sbatch` and chain mpileup after it
5. Submit Stage 0 (mpileup → `_qc.tsv`) per sample; samples that already have a `_qc.tsv` are skipped
6. Submit Stages 1–3 (base, and purity in parallel if `--purity-file` set) dependent on Stage 0
7. Submit Stage 4 (scarHRD) for each pipeline chain
8. Submit QC report → `{base_dir}/qc_report.html`
9. Submit HRD eval report once all scarHRD jobs complete → `tso_samples/hrd_performance_report.html`

The `pipeline_manifest.csv` in each run directory includes a `purity_note` column:
`ok` / `no_purity_data` / `no_purity_file` — records why correction was or wasn't applied per sample.

If `--bam-dir` / `--tsv-dir` are omitted, defaults to the original `rg_bam` and
`20260225_GoldStandard_mpileup_stats_q30/results` directories.

### Processing new BAMs from local storage (batches)

BAMs are ~7–12 GB each. Transfer and process in batches to keep peak storage
bounded. Once Stage 0 produces the `_qc.tsv` for a batch, the BAMs can be
deleted — Stages 1–3 only need the TSV files.

```bash
# 1. Transfer batch
rsync -avzP --progress /local/path/batch01/*.mdup.rg.bam \
  user@cluster:/storage/.../tso_samples/rg_bam_20260505_batch01/

# 2. Run pipeline for this batch
bash slurm/orchestrate_pipeline.sh \
  --bam-dir /storage/.../tso_samples/rg_bam_20260505_batch01 \
  --tsv-dir /storage/.../tso_samples/20260505_batch01_mpileup/results

# 3. After Stage 0 completes, free the BAM storage
rm /storage/.../tso_samples/rg_bam_20260505_batch01/*.bam

# 4. Repeat for batch02, batch03, ...
```

Batch directories follow the convention `rg_bam_YYYYMMDD_batchNN` (BAMs) and
`YYYYMMDD_batchNN_mpileup/results` (TSV output).

### QC report (ad hoc)

El orchestrator lanza el QC report automáticamente al final de cada run. Para regenerarlo
sobre un run ya procesado:

```bash
sbatch --export=ALL,CONFIG=configs/pipeline_config_YYYYMMDD.yaml \
  slurm/run_qc_report.sbatch
```

O directamente sin SLURM:

```bash
conda run -n shapeit4 python pipeline/qc_report.py \
  --config configs/pipeline_config_YYYYMMDD.yaml
```

El HTML resultante se escribe en `{base_dir}/qc_report.html` e incluye:

| Plot | Qué muestra |
|------|-------------|
| BAF coverage | `n_positions` por segmento y muestra — detecta cobertura esparsa |
| outer_mass por CN | Distribución de la feature de decisión con líneas de threshold calibradas |
| Heatmap de margen | Qué segmentos caen cerca del límite de clasificación (inciertos) |
| Composición alélica | Fracción LoH / Balanced / AI / Deletion por muestra (DP=8) |
| Pre/post-merge | Compresión de segmentos tras el merge por adyacencia |

### HRD performance report

Generates a self-contained HTML report comparing pipeline HRD scores against
Myriad GIS gold standard across all production runs (≥ 20260511). Includes
scatter plots (LoH / TAI / LST / HRD_Total vs GIS), confusion matrices
(threshold 42), and a SampleID parsing audit.

```bash
sbatch --export=ALL,REPO_DIR="$(pwd)",CURRENT_RUN=20260520 \
  slurm/run_eval_hrd.sbatch
```

Outputs:
- `{tso_samples_root}/hrd_performance_report.html` — all samples
- `{tso_samples_root}/hrd_performance_report_filtered.html` — blacklisted samples excluded (generated automatically if `configs/blacklist.txt` is non-empty)

Key options (passed via `--export` or editable in the sbatch):

| Variable | Default | Description |
|----------|---------|-------------|
| `CURRENT_RUN` | _(none)_ | Highlights a specific run with ★ in per-run plots |
| `MIN_RUN` | `20260511` | Excludes runs before this date (pre-pipeline runs) |
| `BLACKLIST_FILE` | `configs/blacklist.txt` | Path to sample blacklist for filtered report |

**Sample blacklist** — edit `configs/blacklist.txt` to add or remove samples
from the filtered report. One base ID per line; lines starting with `#` are
ignored. The filtered report is regenerated automatically on every eval run.

### Purity correction

When a `purity_file` is provided in the config (tab-separated, columns:
`sample`, `purity`, `ploidy`), the pipeline applies per-sample BAF correction
in Stage 1 before computing segment-level metrics:

```
BAF_corrected = 0.5 ± (BAF_raw - 0.5) / purity
```

Both the raw and corrected BAF values are retained in outputs (`BAF_raw`,
`BAF_corrected`, `purity_used`). Samples without a purity entry proceed
without correction. The `-HRD` BAM variant is automatically mapped to the
same purity as the base sample.

Set `purity_file` in the run config:

```yaml
input:
  purity_file: "/path/to/purity_ploidy.tsv"  # columns: sample, purity, ploidy
```

### Processing log

`PROCESSED_SAMPLES.tsv` in the repo root tracks every sample processed since
2026-05-11 (first production run). Columns: `run_date`, `sample_id`,
`git_commit`, `config_file`, `output_dir`, `dp_used`, `hrd`, `tai`, `lst`,
`hrd_sum`. Score is dp4 postmerge (fallback dp2/dp0/dp8).

Regenerate at any time:

```bash
conda run -n shapeit4 python pipeline/log_processed_samples.py [--tso-root /path/to/tso_samples]
```

`--tso-root` defaults to `../tso_samples` relative to the repo root.

Git commit is read from `configs/run_info_YYYYMMDD.json` (written automatically
by the orchestrator) or reconstructed from `git log` for retroactive runs.

### Manual invocation (single stage or debugging)

```bash
# Stages 1–3 only (TSV files already exist)
python pipeline/baf_to_hrd_pipeline.py --config configs/pipeline_config_YYYYMMDD.yaml

# Single stage
python pipeline/baf_to_hrd_pipeline.py --config configs/pipeline_config_YYYYMMDD.yaml --stages 3

# Split segments only (generate per-sample CSVs without running the pipeline)
python pipeline/split_segments.py --samples JBLAB001 JBLAB002
python pipeline/split_segments.py --all
python pipeline/split_segments.py --all --outdir /custom/output/dir
```

## Segment files

Per-sample CN segment CSVs are extracted from the combined segtable:
```
/storage/.../tso_samples/unique_copyright_100kb_abs_segtable_free_purity_filtered_unique.csv
```

`split_segments.py` reads this file and writes one CSV per sample to a
date-stamped directory (`tso_samples/YYYYMMDD_segments/`). The orchestrator
runs this automatically at the start of each batch. Use `--force` to overwrite
existing CSVs.

Required format (one CSV per sample, naming convention below):
```
{SAMPLE}_100kb_abs_segtable_free_purity_filtered_unique.csv

chromosome,start,end,segVal,sample
1,900001,211100000,2.02,JBLAB17012
```

## Status

| Stage | Status |
|-------|--------|
| Pre-processing — segment split | functional |
| 0 — BAM → TSV | functional |
| 1 — BAF analysis | functional |
| 2 — Consolidation | functional |
| 3 — Allele-specific CN | functional · being optimised |
| 4 — scarHRD scoring | functional |
| QC report | functional · auto-lanzado por orchestrator |
| HRD performance report | functional · ad hoc via run_eval_hrd.sbatch |
| autoresearch loop | active |

Current `val_allelic_acc` baseline: **0.922** (ASCAT ground truth, DP≥8).
