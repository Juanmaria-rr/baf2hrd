#!/bin/bash
#SBATCH --job-name=recalibrate
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00

# ============================================================
# run_recalibrate.sbatch
#
# Re-runs calibrate_thresholds.py restricted to the 80%
# calibration samples produced by run_split_samples.sbatch.
# Writes honest thresholds to:
#   .../20260406/thresholds_analysis_honest/
#
# Submit after run_split_samples.sbatch completes:
#   sbatch --dependency=afterok:<split_jobid> run_recalibrate.sbatch
#
# After this job, update THRESHOLDS in experiment.py with the
# DP=8 values from thresholds_analysis_honest/thresholds_lookup.json
# and re-submit autoresearch/run_experiment.sbatch.
# ============================================================

set -euo pipefail

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
mkdir -p "${REPO_DIR}/logs"

echo "============================================================"
echo "Threshold re-calibration (calibration set only)"
echo "DATE:    $(date)"
echo "JOB ID:  ${SLURM_JOB_ID:-NA}"
echo "============================================================"

# Patch calibrate_thresholds.py at runtime via a small wrapper so we
# don't permanently modify the source file.  The wrapper overrides:
#   - PARAMS['calibration_samples_file']  → calibration_samples.txt
#   - PARAMS['outdir']                    → thresholds_analysis_honest/
python3 - <<'PYEOF'
import importlib.util, sys, os

# Point to source
spec = importlib.util.spec_from_file_location(
    "calibrate_thresholds",
    os.path.join(os.path.dirname(os.path.abspath(__file__)) if '__file__' in dir() else os.getcwd(),
                 "calibrate_thresholds.py")
)

# We cannot import as a module (top-level code runs on import), so we
# patch sys.argv and run it as a script via exec after patching PARAMS.
# The patching happens by injecting a sitecustomize-style override.

script_dir = os.environ.get('SLURM_SUBMIT_DIR', os.getcwd())
repo_dir = os.path.dirname(script_dir)

cal_file  = os.path.join(script_dir, "autoresearch", "calibration_samples.txt")
out_dir   = ("/storage/scratch01/groups/co/cn_extra/alleleSpecific/tso_samples/20260406"
             "/thresholds_analysis_honest")

if not os.path.exists(cal_file):
    print(f"ERROR: calibration_samples.txt not found: {cal_file}")
    print("Run run_split_samples.sbatch first.")
    sys.exit(1)

print(f"Calibration samples file : {cal_file}")
print(f"Output directory         : {out_dir}")

# Read and patch the script source
src_file = os.path.join(repo_dir, "pipeline", "calibrate_thresholds.py")
with open(src_file) as f:
    src = f.read()

# Inject overrides immediately after PARAMS definition (before any code runs)
patch = f"""
PARAMS['calibration_samples_file'] = {repr(cal_file)}
PARAMS['outdir'] = {repr(out_dir)}
"""
# Insert after the closing brace of PARAMS dict
insert_marker = "    'calibration_samples_file': None,\n}"
patched_src   = src.replace(insert_marker, insert_marker + "\n" + patch, 1)
if patched_src == src:
    print("ERROR: PARAMS patch marker not found — calibrate_thresholds.py may have changed")
    sys.exit(1)

exec(compile(patched_src, src_file, 'exec'), {'__file__': src_file, '__name__': '__main__'})
PYEOF

echo "============================================================"
echo "Recalibration complete: $(date)"
echo "Output: /storage/scratch01/groups/co/cn_extra/alleleSpecific/tso_samples/20260406/thresholds_analysis_honest/"
echo "============================================================"
