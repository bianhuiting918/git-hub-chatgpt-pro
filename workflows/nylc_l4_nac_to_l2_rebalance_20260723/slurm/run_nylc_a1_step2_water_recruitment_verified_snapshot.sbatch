#!/usr/bin/env bash
#SBATCH -J nylc_a2qmfix
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -c 1
#SBATCH --mem-per-cpu=2500M
#SBATCH -t 08:00:00
#SBATCH --array=0-3
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-nylc_a2qmfix-%A_%a.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-nylc_a2qmfix-%A_%a.err

set -euo pipefail
: "${A1_STEP2_WATER_RECRUIT_CODE_ROOT:?set immutable code root}"
: "${A1_STEP2_WATER_RECRUIT_GITHUB_COMMIT:?set immutable GitHub commit}"
CODE_ROOT="$A1_STEP2_WATER_RECRUIT_CODE_ROOT"
EXPECTED_COMMIT="$A1_STEP2_WATER_RECRUIT_GITHUB_COMMIT"
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
cd "$CODE_ROOT"
test "$(cat GITHUB_COMMIT)" = "$EXPECTED_COMMIT"
test -s SNAPSHOT_SHA256.tsv
sha256sum -c SNAPSHOT_SHA256.tsv
"$PY" tests/test_nylc_a1_step2_water_recruitment.py
"$PY" scripts/prepare_audit_nylc_a1_step2_water_recruitment.py --mode describe
bash -n scripts/run_nylc_a1_step2_water_recruitment.sh
export PYTHONDONTWRITEBYTECODE=1
exec bash scripts/run_nylc_a1_step2_water_recruitment.sh
