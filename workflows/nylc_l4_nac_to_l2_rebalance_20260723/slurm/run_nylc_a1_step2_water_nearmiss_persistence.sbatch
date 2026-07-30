#!/usr/bin/env bash
#SBATCH -J nylc_a2wpersist
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -c 1
#SBATCH --mem-per-cpu=2500M
#SBATCH -t 12:00:00
#SBATCH --array=0-8
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-nylc_a2wpersist-%A_%a.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-nylc_a2wpersist-%A_%a.err

set -euo pipefail
: "${A1_STEP2_WATER_PERSIST_CODE_ROOT:?set immutable code root}"
: "${A1_STEP2_WATER_PERSIST_GITHUB_COMMIT:?set immutable GitHub commit}"
CODE_ROOT="$A1_STEP2_WATER_PERSIST_CODE_ROOT"
EXPECTED_COMMIT="$A1_STEP2_WATER_PERSIST_GITHUB_COMMIT"
cd "$CODE_ROOT"
test "$(cat GITHUB_COMMIT)" = "$EXPECTED_COMMIT"
sha256sum -c SNAPSHOT_SHA256.tsv
export PYTHONDONTWRITEBYTECODE=1
export PYTHONPYCACHEPREFIX="${SLURM_TMPDIR:-/tmp}/nylc_a1_step2_water_persist_pycache_${SLURM_ARRAY_JOB_ID:-manual}_${SLURM_ARRAY_TASK_ID:-manual}"
exec bash scripts/run_nylc_a1_step2_water_nearmiss_persistence.sh
