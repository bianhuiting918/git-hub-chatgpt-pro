#!/bin/bash
#SBATCH -J nylc_s1_proton
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem-per-cpu=2500M
#SBATCH -t 08:00:00
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/step1_proton_%j.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/step1_proton_%j.err

set -euo pipefail
TASK=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
FLOW="$TASK/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
BASE="$TASK/candidates/nylc_C18_trueT267_freeGS/qmmm/step1_dftb3_preflight_post61710861_job61712692"
ATTACK="$TASK/candidates/nylc_C18_trueT267_freeGS/qmmm/step1_attack_continuation_job_61725801"
SEED="$ATTACK/attack_seed.rst7"
SEED_SHA=d5e607dd722be4fd6e5d33016afb4975d18b200a78d1fceedf5bfb9fed1b87bf
ACCEPTOR="${ACCEPTOR:?set ACCEPTOR=OD1_or_OD2}"
if [[ "$ACCEPTOR" != OD1 && "$ACCEPTOR" != OD2 ]]; then echo "invalid ACCEPTOR=$ACCEPTOR" >&2; exit 2; fi
job="${SLURM_JOB_ID:?}"
OUT="$TASK/candidates/nylc_C18_trueT267_freeGS/qmmm/step1_proton_bracket_${ACCEPTOR}_job_${job}"
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
AMBER_RUNTIME=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/nylc_gyaq_pa66_l2_nac_qmmm_20260723/amber_runtime

record() {
  local state=$1 detail=$2 stamp
  stamp="$(date -Is)"
  printf '%s\tstep1_proton_bracket\t%s\tcandidate=nylc_C18_trueT267_freeGS;job=%s;acceptor=%s;%s\n' "$stamp" "$state" "$job" "$ACCEPTOR" "$detail" >> "$TASK/run_history.tsv"
  "$PY" - "$stamp" "$state" "$detail" "$job" "$ACCEPTOR" >> "$TASK/run_history.jsonl" <<'PY'
import json,sys
stamp,state,detail,job,acceptor=sys.argv[1:]
print(json.dumps({"timestamp":stamp,"event":"step1_proton_bracket","state":state,"detail":detail,"candidate_id":"nylc_C18_trueT267_freeGS","slurm_job_id":job,"acceptor_hypothesis":acceptor},sort_keys=True))
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

record START "attack_job=61725801;attack_seed_sha=$SEED_SHA;scope=staged_proton_bracket_not_TS_PMF"
actual="$(sha256sum "$SEED" | awk '{print $1}')"
if [[ "$actual" != "$SEED_SHA" ]]; then echo "attack seed SHA mismatch" >&2; exit 21; fi
mkdir "$OUT"
mkdir "$OUT/run"
"$PY" "$FLOW/scripts/prepare_nylc_step1_proton_bracket.py" --source-manifest "$ATTACK/input/manifest.json" --attack-pass "$ATTACK/PASS.json" --attack-seed "$SEED" --acceptor "$ACCEPTOR" --output-dir "$OUT/input" > "$OUT/prepare.stdout" 2> "$OUT/prepare.stderr"

module purge >/dev/null 2>&1 || true
module load amber/2018-hpcx-gcc-7.3.1
export AMBERHOME="$AMBER_RUNTIME"
test -s "$AMBERHOME/dat/slko/3ob-3-1/C-C.skf"
current="$SEED"
cd "$OUT/run"
for stage in p00 p01 p02 p03; do
  sander -O -i "../input/${stage}.in" -o "${stage}.out" -p "$BASE/system.prmtop" -c "$current" -ref "$current" -r "${stage}.rst7" -inf "${stage}.mdinfo"
  if ! grep -q "FINAL RESULTS" "${stage}.out" || ! grep -q "5.  TIMINGS" "${stage}.out"; then
    echo "Amber completion marker missing in $stage" >&2
    exit 31
  fi
  if grep -Eq "SCC is not converged|SANDER BOMB|FATAL|NaN|forrtl|segmentation" "${stage}.out"; then
    echo "Amber hard error in $stage" >&2
    exit 32
  fi
  current="$OUT/run/${stage}.rst7"
  audit="$OUT/run/${stage}_audit.json"
  "$PY" "$FLOW/scripts/audit_nylc_step1_proton_stage.py" --manifest "$OUT/input/manifest.json" --prmtop "$BASE/system.prmtop" --rst7 "$current" --stage "$stage" --output "$audit"
  status=$("$PY" -c 'import json,sys;print(json.load(open(sys.argv[1]))["status"])' "$audit")
  expected=PASS_PROTON_BRACKET_STAGE
  if [[ "$stage" == p03 ]]; then expected=PASS_SCIENTIFIC_TETRAHEDRAL_SEED_REACHED; fi
  if [[ "$status" != "$expected" ]]; then
    cp "$audit" "$OUT/SCIENTIFIC_FAIL.json"
    record FAIL_SCIENTIFIC "stage=$stage;status=$status;audit=$OUT/SCIENTIFIC_FAIL.json"
    trap - ERR
    exit 0
  fi
done
cd "$OUT"
cp "$current" "$OUT/tetrahedral_seed.rst7"
cp "$OUT/run/p03_audit.json" "$OUT/PASS.json"
record PASS_SCIENTIFIC "status=PASS_SCIENTIFIC_TETRAHEDRAL_SEED_REACHED;audit=$OUT/PASS.json;interpretation=constrained_seed_not_TS_committor_PMF"
