#!/bin/bash
#SBATCH -J l4tol2_charge3
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=3G
#SBATCH -t 01:00:00
#SBATCH -o charge_correction_%j.out
#SBATCH -e charge_correction_%j.err

set -euo pipefail
TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
REPO="$TASK_ROOT/repo"
FLOW="$REPO/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
PYTHON=/work/home/acshdt1dks/anaconda3/bin/python3.9
cd "$REPO"
printf '%s\t%s\t%s\t%s\n' "$(date -Is)" charge_correction_rebuild START "job_id=${SLURM_JOB_ID:-NA};candidate_count=3" >> "$TASK_ROOT/run_history.tsv"
"$PYTHON" "$FLOW/scripts/prepare_candidates.py" \
  --manifest "$FLOW/manifests/candidates.json" \
  --task-root "$TASK_ROOT" \
  --candidate-id nylc_c18_11854ps \
  --candidate-id nylc_c23_29684ps \
  --candidate-id nyl50_c18_70792ps \
  --summary-name prepare_charge_correction_summary.tsv \
  > "$TASK_ROOT/audit/prepare_charge_correction.stdout.json"
rc=$?
printf '%s\t%s\t%s\t%s\n' "$(date -Is)" charge_correction_rebuild "rc=$rc" "job_id=${SLURM_JOB_ID:-NA}" >> "$TASK_ROOT/run_history.tsv"
exit "$rc"
