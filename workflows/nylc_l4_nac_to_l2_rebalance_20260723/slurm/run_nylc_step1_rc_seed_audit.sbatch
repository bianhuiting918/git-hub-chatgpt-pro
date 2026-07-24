#!/bin/bash
#SBATCH -J nylc_s1_rcseed
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem-per-cpu=2500M
#SBATCH -t 00:20:00
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/step1_rcseed_%j.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/step1_rcseed_%j.err

set -euo pipefail
TASK=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
FLOW="$TASK/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
job="${SLURM_JOB_ID:?}"
OUT="$TASK/candidates/nylc_C18_trueT267_freeGS/qmmm/step1_rc_seed_audit_job_${job}"
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python

record() {
  local state=$1 detail=$2 stamp
  stamp="$(date -Is)"
  printf '%s\tstep1_rc_seed_audit\t%s\tcandidate=nylc_C18_trueT267_freeGS;job=%s;%s\n' "$stamp" "$state" "$job" "$detail" >> "$TASK/run_history.tsv"
  "$PY" - "$stamp" "$state" "$detail" "$job" >> "$TASK/run_history.jsonl" <<'PY'
import json,sys
stamp,state,detail,job=sys.argv[1:]
print(json.dumps({"timestamp":stamp,"event":"step1_rc_seed_audit","state":state,"detail":detail,"candidate_id":"nylc_C18_trueT267_freeGS","slurm_job_id":job},sort_keys=True))
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

record START "sources=61712692,61718715;scope=geometry_and_candidate_basins"
mkdir "$OUT"
"$PY" "$FLOW/scripts/analyze_nylc_step1_rc_seed.py" --task-root "$TASK" \
  --output "$OUT/step1_rc_seed_audit.json" > "$OUT/analyze.stdout" 2> "$OUT/analyze.stderr"
status=$("$PY" -c 'import json,sys;print(json.load(open(sys.argv[1]))["status"])' "$OUT/step1_rc_seed_audit.json")
if [[ "$status" != PASS_STEP1_RC_SEED_AUDIT ]]; then
  record FAIL_TECHNICAL "status=$status;audit=$OUT/step1_rc_seed_audit.json"
  exit 31
fi
cp "$OUT/step1_rc_seed_audit.json" "$OUT/PASS.json"
record PASS_TECHNICAL "status=$status;audit=$OUT/PASS.json;interpretation=not_validated_RC_TS_PMF"
