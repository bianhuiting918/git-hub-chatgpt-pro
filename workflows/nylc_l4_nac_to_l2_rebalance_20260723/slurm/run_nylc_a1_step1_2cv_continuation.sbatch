#!/usr/bin/env bash
#SBATCH -J nylc_a12cvc
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -c 1
#SBATCH --mem-per-cpu=2500M
#SBATCH --time=36:00:00
#SBATCH --array=0-7
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-nylc_a12cvc-%A_%a.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-nylc_a12cvc-%A_%a.err

set -euo pipefail
: "${A1_2CV_CONT_CODE_ROOT:?set immutable code root}"
: "${A1_2CV_CONT_GITHUB_COMMIT:?set immutable GitHub commit}"
CODE_ROOT="$A1_2CV_CONT_CODE_ROOT"
EXPECTED_COMMIT="$A1_2CV_CONT_GITHUB_COMMIT"
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
export PYTHONDONTWRITEBYTECODE=1
cd "$CODE_ROOT"
test "$(cat GITHUB_COMMIT)" = "$EXPECTED_COMMIT"
sha256sum -c SNAPSHOT_SHA256.tsv
"$PY" scripts/prepare_audit_nylc_a1_step1_2cv_continuation.py --mode describe
bash -n scripts/run_nylc_a1_step1_2cv_continuation.sh
exec bash scripts/run_nylc_a1_step1_2cv_continuation.sh
