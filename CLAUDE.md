# CLAUDE.md — baf_hrd_pipeline

## Proyecto

Pipeline end-to-end para inferir alelo-específicos copy numbers (A_cn, B_cn)
a partir de datos TSO500 tumor-only, con el objetivo de calcular HRD scores
mediante scarHRD.

## Estructura

```
baf_hrd_pipeline/
├── pipeline/               # Código Python de producción
│   ├── baf_to_hrd_pipeline.py   ← pipeline principal (Stage 1–3)
│   ├── calibrate_thresholds.py  ← calibración de thresholds por DP
│   ├── parse_mpileup.py         ← Stage 0b: BAM → _qc.tsv
│   ├── split_samples.py         ← split 80/20 calibración/validación
│   ├── functions.py             ← utilidades compartidas
│   ├── correct_REFALT_*.py      ← preprocesado ASCAT downsampled
│   └── evaluate_*.py            ← evaluación de resultados
├── slurm/                  # Scripts de envío a SLURM
│   ├── orchestrate_pipeline.sh  ← entry point principal (BAM → HRD)
│   ├── run_pipeline.sbatch      ← Stages 1–3 completos
│   └── run_*.sbatch             ← stages individuales
├── configs/                # Configuración por run (YAML)
│   ├── pipeline_config_template.yaml
│   └── pipeline_config_YYYYMMDD.yaml
└── autoresearch/           # Loop autónomo de mejora de Stage 3
    ├── experiment.py       ← ÚNICO fichero que modifican los agentes
    ├── prepare.py
    ├── evaluate.py
    └── README.md / program.md
```

## Flujo completo

```
BAMs → slurm/orchestrate_pipeline.sh
         └── Stage 0: mpileup + parse_mpileup → _qc.tsv
         └── Stage 1: BAF analysis → segments + manifest
         └── Stage 2: consolidation → merged CSV
         └── Stage 3: allele-specific CN → scarHRD inputs
                       ↑ optimizado por autoresearch/
```

## Reglas para agentes

- **Sólo modificar `autoresearch/experiment.py`** en el loop de autoresearch.
- No modificar `pipeline/` sin instrucción explícita.
- Toda ejecución en este clúster HPC va por `sbatch`, nunca directamente en login shell.
- Configuración de runs en `configs/pipeline_config_YYYYMMDD.yaml`.

## Métrica objetivo

`val_allelic_acc` — fracción de segmentos donde (A_cn, B_cn) predicho
coincide exactamente con el ground truth ASCAT. Baseline actual: **0.9221**.
