#!/bin/bash
#SBATCH -J nylc_s1_prodqm_audit
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem-per-cpu=1000M
#SBATCH -t 00:15:00
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/step1_prodqm_reaudit_%j.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/step1_prodqm_reaudit_%j.err

set -euo pipefail
TASK=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
FLOW="$TASK/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
PARENT="$TASK/candidates/nylc_C18_trueT267_freeGS/qmmm/step1_production_qm_smoke_job_61716670"
job="${SLURM_JOB_ID:?}"
OUT="$TASK/candidates/nylc_C18_trueT267_freeGS/qmmm/step1_production_qm_smoke_reaudit_job_${job}"
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python

record() {
  local state=$1 detail=$2 stamp
  stamp="$(date -Is)"
  printf '%s\tstep1_production_qm_smoke_reaudit\t%s\tcandidate=nylc_C18_trueT267_freeGS;job=%s;%s\n' "$stamp" "$state" "$job" "$detail" >> "$TASK/run_history.tsv"
  "$PY" - "$stamp" "$state" "$detail" "$job" >> "$TASK/run_history.jsonl" <<'PY'
import json,sys
stamp,state,detail,job=sys.argv[1:]
print(json.dumps({"timestamp":stamp,"event":"step1_production_qm_smoke_reaudit","state":state,"detail":detail,"candidate_id":"nylc_C18_trueT267_freeGS","slurm_job_id":job},sort_keys=True))
PY
}

record START "parent_job=61716670;reason=Amber_Run_spaced_done_parser"
mkdir "$OUT"
set +e
"$PY" "$FLOW/scripts/audit_nylc_step1_production_qm_smoke.py" \
  --region-audit "$PARENT/input/production_qm_region_audit.json" \
  --core "$PARENT/core_one_step.out" --network "$PARENT/network_one_step.out" \
  --output "$OUT/production_qm_smoke_result.json"
rc=$?
set -e
status=$("$PY" -c 'import json,sys;print(json.load(open(sys.argv[1]))["status"])' "$OUT/production_qm_smoke_result.json")
if [[ "$rc" -ne 0 || "$status" != PASS_PRODUCTION_QM_REGION_SMOKE ]]; then
  record FAIL_TECHNICAL "status=$status;audit=$OUT/production_qm_smoke_result.json"
  exit 31
fi
cp "$OUT/production_qm_smoke_result.json" "$OUT/PASS.json"
record PASS_TECHNICAL "status=$status;audit=$OUT/PASS.json;parent_job=61716670;interpretation=not_TS_RC_PMF"
