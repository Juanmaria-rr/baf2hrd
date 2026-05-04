# baf2hrd

> **Work in progress** — pipeline under active development.

End-to-end pipeline for inferring allele-specific copy numbers (A_cn, B_cn)
from TSO500 tumor-only sequencing data, producing inputs for
[scarHRD](https://github.com/sztup/scarHRD) HRD scoring.

## Overview

```
BAM files
   └── Stage 0: mpileup → _qc.tsv
   └── Stage 1: per-sample BAF analysis → segments + manifest
   └── Stage 2: consolidation → merged dataset
   └── Stage 3: allele-specific CN inference → scarHRD inputs
                 ↑ optimised by autoresearch/
```

## Structure

```
pipeline/     Python source (main pipeline + utilities)
slurm/        SLURM submission scripts
configs/      Per-run YAML configuration
autoresearch/ Autonomous loop for improving Stage 3 classification
```

## Usage

```bash
# Full pipeline (Stages 1–3)
python pipeline/baf_to_hrd_pipeline.py --config configs/pipeline_config_YYYYMMDD.yaml

# Single stage
python pipeline/baf_to_hrd_pipeline.py --config configs/pipeline_config_YYYYMMDD.yaml --stages 3

# From BAMs (orchestrates Stage 0 + pipeline via SLURM)
bash slurm/orchestrate_pipeline.sh
```

Copy and fill `configs/pipeline_config_template.yaml` for each new run.

## Status

| Stage | Status |
|-------|--------|
| 0 — BAM → TSV | functional |
| 1 — BAF analysis | functional |
| 2 — Consolidation | functional |
| 3 — Allele-specific CN | functional · being optimised |
| autoresearch loop | active |

Current `val_allelic_acc` baseline: **0.922** (ASCAT ground truth, DP≥8).
