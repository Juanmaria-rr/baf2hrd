#!/bin/bash
#SBATCH --job-name=autoresearch_experiment
#SBATCH --output=logs/experiment_%j.log
#SBATCH --error=logs/experiment_%j.err
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=0:30:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jmroldan@cnio.es

# ============================================================
# run_experiment.sbatch
#
# Run experiment.py — re-submit after every change to experiment.py.
# Results are appended to autoresearch_baf_hrd/results/log.jsonl.
#
# Requires val_data.csv to already exist (run run_prepare.sbatch first).
#
# Workflow:
#   1. Edit experiment.py (THRESHOLDS, OUTER_LOWER/UPPER, classify_segment, etc.)
#   2. sbatch autoresearch_baf_hrd/run_experiment.sbatch
#   3. Check logs/experiment_<jobid>.log for val_allelic_acc
#   4. Review results/log.jsonl for the full history
#   5. Repeat
#
# Submit:
#   sbatch autoresearch_baf_hrd/run_experiment.sbatch
# ============================================================

set -euo pipefail

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
mkdir -p "${SCRIPT_DIR}/logs"
mkdir -p "${SCRIPT_DIR}/results"

VAL_DATA="${SCRIPT_DIR}/val_data.csv"
if [[ ! -f "${VAL_DATA}" ]]; then
    echo "ERROR: val_data.csv not found at ${VAL_DATA}"
    echo "Run run_prepare.sbatch first."
    exit 1
fi

echo "=== run_experiment.sbatch  job=${SLURM_JOB_ID}  $(date) ==="
echo "val_data.csv: $(wc -l < "${VAL_DATA}") lines"

conda run -n shapeit4 python "${SCRIPT_DIR}/experiment.py"

echo "=== DONE  $(date) ==="
echo ""
echo "--- Last result ---"
_LAST_JSON=$(tail -1 "${SCRIPT_DIR}/results/log.jsonl")
python3 -c "
import json, sys
d = json.loads(sys.argv[1])
print(f'  val_allelic_acc : {d[\"val_allelic_acc\"]:.4f}')
print(f'  val_class_acc   : {d[\"val_class_acc\"]:.4f}')
print(f'  notes           : {d[\"notes\"]}')
print(f'  timestamp       : {d[\"timestamp\"]}')
" "${_LAST_JSON}"
