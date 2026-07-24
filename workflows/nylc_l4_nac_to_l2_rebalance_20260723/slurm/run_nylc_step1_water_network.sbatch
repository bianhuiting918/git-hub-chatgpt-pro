#!/bin/bash
#SBATCH -J nylc_s1_waternet
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem-per-cpu=1500M
#SBATCH -t 00:15:00
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/step1_waternet_%j.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/step1_waternet_%j.err

set -euo pipefail
TASK=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
FLOW="$TASK/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
GRO="$TASK/candidates/nylc_C18_trueT267_freeGS/qmmm/step1_dftb3_preflight_post61710861_job61712692/selected_free_l2_nac.gro"
job="${SLURM_JOB_ID:?}"
OUT="$TASK/candidates/nylc_C18_trueT267_freeGS/qmmm/step1_water_network_job_${job}"
PY=/work/home/acshdt1dks/anaconda3/bin/python3.9

record() {
  local state=$1 detail=$2 stamp
  stamp="$(date -Is)"
  printf '%s\tstep1_water_network\t%s\tcandidate=nylc_C18_trueT267_freeGS;job=%s;%s\n' "$stamp" "$state" "$job" "$detail" >> "$TASK/run_history.tsv"
  "$PY" - "$stamp" "$state" "$detail" "$job" >> "$TASK/run_history.jsonl" <<'PY'
import json,sys
stamp,state,detail,job=sys.argv[1:]
print(json.dumps({"timestamp":stamp,"event":"step1_water_network","state":state,"detail":detail,"candidate_id":"nylc_C18_trueT267_freeGS","slurm_job_id":job},sort_keys=True))
PY
}

on_error() {
  local rc=$?
  trap - ERR
  set +e
  record FAIL_TECHNICAL "exit_code=$rc;out=$OUT"
  exit "$rc"
}
trap on_error ERR

record START "source=434ps_lowest_potential_strict_NAC;dftb3_preflight=61712692"
mkdir "$OUT"
"$PY" "$FLOW/scripts/analyze_nylc_step1_water_network.py" --gro "$GRO" \
  --output "$OUT/step1_water_network_audit.json" > "$OUT/analyze.stdout" 2> "$OUT/analyze.stderr"
status=$("$PY" -c 'import json,sys;print(json.load(open(sys.argv[1]))["status"])' "$OUT/step1_water_network_audit.json")
if [[ "$status" == PASS_STRUCTURAL_GS_HYPOTHESIS ]]; then
  record PASS_STRUCTURAL "status=$status;audit=$OUT/step1_water_network_audit.json;not_TS=true"
else
  record NOT_EVALUATED "status=$status;audit=$OUT/step1_water_network_audit.json"
fi
