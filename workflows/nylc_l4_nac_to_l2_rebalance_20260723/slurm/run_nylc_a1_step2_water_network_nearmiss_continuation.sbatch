#!/usr/bin/env bash
#SBATCH -J nylc_a2nmcont
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -c 1
#SBATCH --mem-per-cpu=2500M
#SBATCH -t 18:00:00
#SBATCH --array=0-7
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-nylc_a2nmcont-%A_%a.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-nylc_a2nmcont-%A_%a.err

set -euo pipefail
: "${A1_STEP2_NEARMISS_CODE_ROOT:?set immutable code root}"
: "${A1_STEP2_NEARMISS_GITHUB_COMMIT:?set immutable GitHub commit}"
CODE_ROOT="$A1_STEP2_NEARMISS_CODE_ROOT"
EXPECTED_COMMIT="$A1_STEP2_NEARMISS_GITHUB_COMMIT"
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
cd "$CODE_ROOT"
test "$(tr -d '\r\n' < GITHUB_COMMIT)" = "$EXPECTED_COMMIT"
test -s SNAPSHOT_SHA256.tsv
sha256sum -c SNAPSHOT_SHA256.tsv
export PYTHONDONTWRITEBYTECODE=1
export PYTHONPYCACHEPREFIX="${SLURM_TMPDIR:-/tmp}/nylc_a1_step2_nearmiss_pycache_${SLURM_ARRAY_JOB_ID:-manual}_${SLURM_ARRAY_TASK_ID:-manual}"
"$PY" tests/test_nylc_a1_launch_environment_contract.py
bash -n scripts/run_nylc_a1_step2_water_network_nearmiss_continuation.sh
exec bash scripts/run_nylc_a1_step2_water_network_nearmiss_continuation.sh
