#!/bin/bash
#SBATCH -J nylc_t267_select
#SBATCH -p xahdnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem-per-cpu=3500M
#SBATCH -t 00:20:00
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/true_thr267_select_%j.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/true_thr267_select_%j.err

set -euo pipefail
TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
FLOW="$TASK_ROOT/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
PYTHON=/work/home/acshdt1dks/anaconda3/bin/python3.9
source_root="$TASK_ROOT/candidates/nylc_C18_trueT267_recapture/recapture_release2"
output_root="$TASK_ROOT/candidates/nylc_C18_trueT267_recapture/selected_low_energy_nac/release2_${SLURM_JOB_ID}"
mkdir -p "$output_root" "$TASK_ROOT/logs"

record() {
  local state=$1 detail=$2 stamp
  stamp="$(date -Is)"
  (
    flock -x 9
    printf '%s\t%s\t%s\t%s\n' "$stamp" true_thr267_low_energy_nac_selection "$state" "$detail" >&9
  ) 9>>"$TASK_ROOT/run_history.tsv"
  (
    flock -x 8
    "$PYTHON" - "$stamp" "$state" "$detail" "${SLURM_JOB_ID}" >&8 <<'PY'
import json,sys
stamp,state,detail,job=sys.argv[1:]
print(json.dumps({
  "timestamp":stamp,
  "event":"true_thr267_low_energy_nac_selection",
  "state":state,
  "detail":detail,
  "candidate_id":"nylc_C18_trueT267_recapture",
  "slurm_job_id":job,
},sort_keys=True))
PY
  ) 8>>"$TASK_ROOT/run_history.jsonl"
}

failure=1
trap 'rc=$?; if [[ $failure -ne 0 ]]; then record FAIL "rc=$rc;source=release2"; fi; exit $rc' EXIT
record START "source=release2;rule=lowest_potential_among_true_NAC"
"$FLOW/scripts/select_true_thr267_low_energy_nac.sh" "$source_root" "$output_root"
record PASS_TECHNICAL "output=$output_root/selected_lowest_potential_nac.json"
failure=0
trap - EXIT
