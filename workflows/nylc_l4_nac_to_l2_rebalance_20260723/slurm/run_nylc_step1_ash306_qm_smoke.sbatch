#!/bin/bash
#SBATCH -J nylc_s1_ash306
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem-per-cpu=2500M
#SBATCH -t 03:00:00
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/step1_ash306_%j.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/step1_ash306_%j.err

set -euo pipefail
TASK=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
FLOW="$TASK/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
SRC="$TASK/candidates/nylc_C18_trueT267_freeGS/qmmm/ash306_full_system_preflight_job_61718715"
job="${SLURM_JOB_ID:?}"
OUT="$TASK/candidates/nylc_C18_trueT267_freeGS/qmmm/step1_ash306_qm_smoke_job_${job}"
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
AMBER_RUNTIME=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/nylc_gyaq_pa66_l2_nac_qmmm_20260723/amber_runtime

record() {
  local state=$1 detail=$2 stamp
  stamp="$(date -Is)"
  printf '%s\tstep1_ash306_qm_smoke\t%s\tcandidate=nylc_C18_trueT267_freeGS;job=%s;%s\n' "$stamp" "$state" "$job" "$detail" >> "$TASK/run_history.tsv"
  "$PY" - "$stamp" "$state" "$detail" "$job" >> "$TASK/run_history.jsonl" <<'PY'
import json,sys
stamp,state,detail,job=sys.argv[1:]
print(json.dumps({"timestamp":stamp,"event":"step1_ash306_qm_smoke","state":state,"detail":detail,"candidate_id":"nylc_C18_trueT267_freeGS","slurm_job_id":job},sort_keys=True))
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

record START "source_job=61718715;protonation=Asp306H_Asp308-;scope=one_step_numerical_smoke"
test -s "$SRC/amber/system.prmtop"
test -s "$SRC/amber/system.rst7"
mkdir "$OUT"
"$PY" "$FLOW/scripts/prepare_nylc_step1_ash306_qm_smoke.py" \
  --source-dir "$SRC/amber" --output-dir "$OUT/input" \
  > "$OUT/prepare.stdout" 2> "$OUT/prepare.stderr"

module purge >/dev/null 2>&1 || true
module load amber/2018-hpcx-gcc-7.3.1
export AMBERHOME="$AMBER_RUNTIME"
test -s "$AMBERHOME/dat/slko/3ob-3-1/C-C.skf"
sander -O -i "$OUT/input/core_one_step.in" -o "$OUT/core_one_step.out" \
  -p "$SRC/amber/system.prmtop" -c "$SRC/amber/system.rst7" \
  -r "$OUT/core_one_step.rst7" -inf "$OUT/core.mdinfo"
sander -O -i "$OUT/input/network_one_step.in" -o "$OUT/network_one_step.out" \
  -p "$SRC/amber/system.prmtop" -c "$SRC/amber/system.rst7" \
  -r "$OUT/network_one_step.rst7" -inf "$OUT/network.mdinfo"
"$PY" "$FLOW/scripts/audit_nylc_step1_ash306_qm_smoke.py" \
  --region-audit "$OUT/input/ash306_qm_region_audit.json" \
  --core "$OUT/core_one_step.out" --network "$OUT/network_one_step.out" \
  --output "$OUT/ash306_qm_smoke_result.json"
status=$("$PY" -c 'import json,sys;print(json.load(open(sys.argv[1]))["status"])' "$OUT/ash306_qm_smoke_result.json")
if [[ "$status" != PASS_ASH306_QM_REGION_SMOKE ]]; then
  record FAIL_TECHNICAL "status=$status;audit=$OUT/ash306_qm_smoke_result.json"
  exit 31
fi
cp "$OUT/ash306_qm_smoke_result.json" "$OUT/PASS.json"
record PASS_TECHNICAL "status=$status;audit=$OUT/PASS.json;interpretation=not_TS_RC_PMF"
