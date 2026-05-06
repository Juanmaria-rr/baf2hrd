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
                 ↑ optimised by autoresearch/
```

## Structure

```
pipeline/     Python source (main pipeline + utilities)
  split_segments.py        ← pre-processing: combined segtable → per-sample CSVs
slurm/        SLURM submission scripts
configs/      Per-run YAML configuration (auto-generated per run)
autoresearch/ Autonomous loop for improving Stage 3 classification
```

## Usage

### Full run from BAMs (recommended)

```bash
bash slurm/orchestrate_pipeline.sh \
  --bam-dir /path/to/rg_bam_YYYYMMDD_batchNN \
  --tsv-dir /path/to/YYYYMMDD_batchNN_mpileup/results
```

The orchestrator will:
1. Run `split_segments.py` to generate per-sample segment CSVs from the combined segtable
2. Auto-generate `configs/pipeline_config_YYYYMMDD.yaml` with the correct paths
3. Submit Stage 0 (mpileup → `_qc.tsv`) per sample via SLURM
4. Submit Stages 1–3 as a dependent job once all Stage 0 jobs complete

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
| autoresearch loop | active |

Current `val_allelic_acc` baseline: **0.922** (ASCAT ground truth, DP≥8).
