#!/bin/bash
#SBATCH -J l4tol2_build5
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 01:00:00
#SBATCH -o prepare_candidates_%j.out
#SBATCH -e prepare_candidates_%j.err

set -euo pipefail
TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
REPO="$TASK_ROOT/repo"
FLOW="$REPO/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
PYTHON=/work/home/acshdt1dks/anaconda3/bin/python3.9
cd "$REPO"
started="$(date -Is)"
printf '%s\t%s\t%s\t%s\n' "$started" prepare_candidates START "job_id=${SLURM_JOB_ID:-NA}" >> "$TASK_ROOT/run_history.tsv"
"$PYTHON" "$FLOW/scripts/prepare_candidates.py" \
  --manifest "$FLOW/manifests/candidates.json" \
  --task-root "$TASK_ROOT" \
  > "$TASK_ROOT/audit/prepare_candidates.stdout.json"
rc=$?
printf '%s\t%s\t%s\t%s\n' "$(date -Is)" prepare_candidates "rc=$rc" "job_id=${SLURM_JOB_ID:-NA}" >> "$TASK_ROOT/run_history.tsv"
exit "$rc"
