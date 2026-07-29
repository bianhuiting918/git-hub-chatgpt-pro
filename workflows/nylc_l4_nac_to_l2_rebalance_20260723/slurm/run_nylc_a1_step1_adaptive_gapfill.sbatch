#!/usr/bin/env bash
#SBATCH -J nylc_a1adgap
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -c 1
#SBATCH --mem-per-cpu=2500M
#SBATCH --time=24:00:00
#SBATCH --array=0-7
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-nylc_a1adgap-%A_%a.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-nylc_a1adgap-%A_%a.err

set -euo pipefail
: "${A1_ADAPTIVE_GAPFILL_CODE_ROOT:?set immutable code root}"
: "${A1_ADAPTIVE_GAPFILL_GITHUB_COMMIT:?set immutable GitHub commit}"
CODE_ROOT="$A1_ADAPTIVE_GAPFILL_CODE_ROOT"
EXPECTED_COMMIT="$A1_ADAPTIVE_GAPFILL_GITHUB_COMMIT"
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
export PYTHONPYCACHEPREFIX="${SLURM_TMPDIR:-/tmp}/nylc_a1_adaptive_pycache_${SLURM_ARRAY_JOB_ID:-manual}_${SLURM_ARRAY_TASK_ID:-manual}"
cd "$CODE_ROOT"
test "$(cat GITHUB_COMMIT)" = "$EXPECTED_COMMIT"
test "$(git rev-parse HEAD)" = "$EXPECTED_COMMIT"
sha256sum -c SNAPSHOT_SHA256.tsv
"$PY" tests/test_nylc_a1_step1_adaptive_gapfill.py
"$PY" -m py_compile scripts/prepare_audit_nylc_a1_step1_adaptive_gapfill.py
"$PY" scripts/prepare_audit_nylc_a1_step1_adaptive_gapfill.py --mode describe
bash -n scripts/run_nylc_a1_step1_adaptive_gapfill.sh
export PYTHONDONTWRITEBYTECODE=1
exec bash scripts/run_nylc_a1_step1_adaptive_gapfill.sh
