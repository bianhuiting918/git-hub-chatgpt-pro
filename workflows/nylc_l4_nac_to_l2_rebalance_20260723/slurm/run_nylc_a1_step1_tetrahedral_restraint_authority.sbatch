#!/usr/bin/env bash
#SBATCH -J nylc_tauth
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -c 1
#SBATCH --mem-per-cpu=2500M
#SBATCH --time=02:00:00
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-nylc_tauth-%j.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-nylc_tauth-%j.err

set -euo pipefail
: "$A1_TETRA_AUTH_CODE_ROOT"
: "$A1_TETRA_AUTH_GITHUB_COMMIT"
CODE_ROOT="$A1_TETRA_AUTH_CODE_ROOT"
EXPECTED_COMMIT="$A1_TETRA_AUTH_GITHUB_COMMIT"
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
export PYTHONDONTWRITEBYTECODE=1
set +u
CACHE_BASE="$SLURM_TMPDIR"
set -u
if [[ -z "$CACHE_BASE" ]]; then CACHE_BASE=/tmp; fi
export PYTHONPYCACHEPREFIX="$CACHE_BASE/nylc_pycache_"$SLURM_JOB_ID
cd "$CODE_ROOT"
test "$(tr -d '\r\n' < GITHUB_COMMIT)" = "$EXPECTED_COMMIT"
sha256sum -c SNAPSHOT_SHA256.tsv
"$PY" scripts/prepare_audit_nylc_a1_step1_tetrahedral_bridge.py --mode describe
bash -n scripts/run_nylc_a1_step1_tetrahedral_restraint_authority.sh
exec bash scripts/run_nylc_a1_step1_tetrahedral_restraint_authority.sh
