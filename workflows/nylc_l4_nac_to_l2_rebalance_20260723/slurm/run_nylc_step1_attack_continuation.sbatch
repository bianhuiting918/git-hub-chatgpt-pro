#!/bin/bash
#SBATCH -J nylc_s1_atkcont
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem-per-cpu=2500M
#SBATCH -t 04:00:00
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/step1_attack_cont_%j.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/step1_attack_cont_%j.err

set -euo pipefail
TASK=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
FLOW="$TASK/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
BASE="$TASK/candidates/nylc_C18_trueT267_freeGS/qmmm/step1_dftb3_preflight_post61710861_job61712692"
SOURCE="$TASK/candidates/nylc_C18_trueT267_freeGS/qmmm/step1_constrained_bracket_OD1_job_61723169/run/w01.rst7"
SOURCE_SHA=10d9b662cdd0bb2e945b34b798e4688373e43e1cca5737b8c86f7658c1975f38
job="${SLURM_JOB_ID:?}"
OUT="$TASK/candidates/nylc_C18_trueT267_freeGS/qmmm/step1_attack_continuation_job_$job"
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
AMBER_RUNTIME=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/nylc_gyaq_pa66_l2_nac_qmmm_20260723/amber_runtime

record() {
  local state=$1 detail=$2 stamp
  stamp="$(date -Is)"
  printf '%s\tstep1_attack_continuation\t%s\tcandidate=nylc_C18_trueT267_freeGS;job=%s;%s\n' "$stamp" "$state" "$job" "$detail" >> "$TASK/run_history.tsv"
  "$PY" - "$stamp" "$state" "$detail" "$job" >> "$TASK/run_history.jsonl" <<'PY'
import json,sys
stamp,state,detail,job=sys.argv[1:]
print(json.dumps({"timestamp":stamp,"event":"step1_attack_continuation","state":state,"detail":detail,"candidate_id":"nylc_C18_trueT267_freeGS","slurm_job_id":job},sort_keys=True))
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

record START "source_job=61723169;source_sha=$SOURCE_SHA;scope=attack_only_three_round_max"
test -s "$SOURCE"
actual="$(sha256sum "$SOURCE" | awk '{print $1}')"
if [[ "$actual" != "$SOURCE_SHA" ]]; then
  echo "source restart SHA mismatch" >&2
  exit 21
fi
mkdir "$OUT"
mkdir "$OUT/run"
"$PY" "$FLOW/scripts/prepare_nylc_step1_constrained_bracket.py" --source-dir "$BASE" --acceptor OD1 --output-dir "$OUT/input" > "$OUT/prepare.stdout" 2> "$OUT/prepare.stderr"

module purge >/dev/null 2>&1 || true
module load amber/2018-hpcx-gcc-7.3.1
export AMBERHOME="$AMBER_RUNTIME"
test -s "$AMBERHOME/dat/slko/3ob-3-1/C-C.skf"
current="$SOURCE"
passed=0
last_gate=
cd "$OUT/run"
for round in r01 r02 r03; do
  test -s "../input/w01.RST"
  sander -O -i "../input/w01.in" -o "${round}.out" -p "$BASE/system.prmtop" -c "$current" -ref "$current" -r "${round}.rst7" -inf "${round}.mdinfo"
  if ! grep -q "FINAL RESULTS" "${round}.out" || ! grep -q "5.  TIMINGS" "${round}.out"; then
    echo "Amber completion marker missing in $round" >&2
    exit 31
  fi
  if grep -Eq "SCC is not converged|SANDER BOMB|FATAL|NaN|forrtl|segmentation" "${round}.out"; then
    echo "Amber hard error in $round" >&2
    exit 32
  fi
  current="$OUT/run/${round}.rst7"
  last_gate="$OUT/run/${round}_gate.json"
  "$PY" "$FLOW/scripts/gate_nylc_step1_attack_prerequisite.py" --prmtop "$BASE/system.prmtop" --rst7 "$current" --output "$last_gate"
  gate_status=$("$PY" -c 'import json,sys;print(json.load(open(sys.argv[1]))["status"])' "$last_gate")
  if [[ "$gate_status" == PASS_ATTACK_PREREQUISITE ]]; then
    cp "$current" "$OUT/attack_seed.rst7"
    cp "$last_gate" "$OUT/PASS.json"
    record PASS_SCIENTIFIC "status=$gate_status;round=$round;audit=$OUT/PASS.json;interpretation=attack_seed_not_TS_PMF"
    passed=1
    break
  fi
done
cd "$OUT"
if [[ "$passed" -ne 1 ]]; then
  cp "$last_gate" "$OUT/SCIENTIFIC_FAIL.json"
  record FAIL_SCIENTIFIC "status=FAIL_SCIENTIFIC_ATTACK_PREREQUISITE;rounds=3;audit=$OUT/SCIENTIFIC_FAIL.json;proton_windows=NOT_RUN"
fi
exit 0
