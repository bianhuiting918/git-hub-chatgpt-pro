#!/bin/bash
#SBATCH -J nylc_t267_gs
#SBATCH -p xahdnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --gres=dcu:1
#SBATCH --mem-per-cpu=3500M
#SBATCH -t 00:20:00
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/true_thr267_gs_%j.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/true_thr267_gs_%j.err

set -euo pipefail
TASK=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
FLOW="$TASK/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
PY=/work/home/acshdt1dks/anaconda3/bin/python3.9
job="${SLURM_JOB_ID:?}"
record() {
  local state=$1 detail=$2 stamp
  stamp="$(date -Is)"
  (
    flock -x 9
    printf '%s\ttrue_thr267_unrestrained_gs_selection\t%s\tjob=%s;%s\n' "$stamp" "$state" "$job" "$detail" >&9
  ) 9>>"$TASK/run_history.tsv"
  (
    flock -x 8
    "$PY" - "$stamp" "$state" "$detail" "$job" >&8 <<'PY'
import json,sys
stamp,state,detail,job=sys.argv[1:]
print(json.dumps({"timestamp":stamp,"event":"true_thr267_unrestrained_gs_selection","state":state,"detail":detail,"slurm_job_id":job,"candidate_id":"nylc_C18_trueT267_freeGS"},sort_keys=True))
PY
  ) 8>>"$TASK/run_history.jsonl"
}
failure=1
trap 'rc=$?; if [[ $failure -ne 0 ]]; then record FAIL "rc=$rc"; fi; exit $rc' EXIT
record START "source=61705307_seed26703_time8ps"
"$FLOW/scripts/select_true_thr267_unrestrained_gs.sh"
record PASS "manifest=$TASK/candidates/nylc_C18_trueT267_freeGS/source_manifest.json"
failure=0
trap - EXIT
