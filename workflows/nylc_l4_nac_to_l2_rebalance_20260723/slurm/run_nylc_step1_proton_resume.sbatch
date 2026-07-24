#!/bin/bash
#SBATCH -J nylc_s1_presume
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem-per-cpu=2500M
#SBATCH -t 08:00:00
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/step1_proton_resume_%j.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/step1_proton_resume_%j.err

set -euo pipefail
TASK=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
FLOW="$TASK/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
BASE="$TASK/candidates/nylc_C18_trueT267_freeGS/qmmm/step1_dftb3_preflight_post61710861_job61712692"
ACCEPTOR="${ACCEPTOR:?set ACCEPTOR=OD1_or_OD2}"
SOURCE_JOB="${SOURCE_JOB:?set SOURCE_JOB=failed_p00_job_id}"
if [[ "$ACCEPTOR" != OD1 && "$ACCEPTOR" != OD2 ]]; then echo "invalid ACCEPTOR=$ACCEPTOR" >&2; exit 2; fi
if [[ ! "$SOURCE_JOB" =~ ^[0-9]+$ ]]; then echo "invalid SOURCE_JOB=$SOURCE_JOB" >&2; exit 2; fi
SOURCE="$TASK/candidates/nylc_C18_trueT267_freeGS/qmmm/step1_proton_bracket_${ACCEPTOR}_job_${SOURCE_JOB}"
SOURCE_FAIL="$SOURCE/SCIENTIFIC_FAIL.json"
SOURCE_RST="$SOURCE/run/p00.rst7"
job="${SLURM_JOB_ID:?}"
OUT="$TASK/candidates/nylc_C18_trueT267_freeGS/qmmm/step1_proton_resume_${ACCEPTOR}_job_${job}"
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
AMBER_RUNTIME=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/nylc_gyaq_pa66_l2_nac_qmmm_20260723/amber_runtime

record() {
  local state=$1 detail=$2 stamp
  stamp="$(date -Is)"
  printf '%s\tstep1_proton_resume\t%s\tcandidate=nylc_C18_trueT267_freeGS;job=%s;acceptor=%s;source_job=%s;%s\n' "$stamp" "$state" "$job" "$ACCEPTOR" "$SOURCE_JOB" "$detail" >> "$TASK/run_history.tsv"
  "$PY" - "$stamp" "$state" "$detail" "$job" "$ACCEPTOR" "$SOURCE_JOB" >> "$TASK/run_history.jsonl" <<'PY'
import json,sys
stamp,state,detail,job,acceptor,source_job=sys.argv[1:]
print(json.dumps({"timestamp":stamp,"event":"step1_proton_resume","state":state,"detail":detail,"candidate_id":"nylc_C18_trueT267_freeGS","slurm_job_id":job,"acceptor_hypothesis":acceptor,"source_job_id":source_job},sort_keys=True))
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

test -s "$SOURCE_FAIL"
test -s "$SOURCE_RST"
"$PY" - "$SOURCE_FAIL" "$ACCEPTOR" <<'PY'
import json,sys
data=json.load(open(sys.argv[1]))
if data.get("status") != "FAIL_SCIENTIFIC_PROTON_BRACKET_STAGE":
    raise SystemExit("source is not a scientific p00 failure")
if data.get("stage") != "p00":
    raise SystemExit("source failure is not p00")
if data.get("acceptor_hypothesis") != sys.argv[2]:
    raise SystemExit("source acceptor mismatch")
gates=data.get("gates",{})
if not gates.get("chemistry_valid") or not gates.get("proton_attached") or gates.get("stage_progress"):
    raise SystemExit("source does not meet safe continuation gates")
PY
source_sha="$(sha256sum "$SOURCE_RST" | awk '{print $1}')"
record START "source_endpoint_sha=$source_sha;scope=same_hamiltonian_p00_continuation_then_existing_schedule_not_TS_PMF"
mkdir "$OUT"
mkdir "$OUT/input" "$OUT/run"
cp -a "$SOURCE/input/." "$OUT/input/"
"$PY" - "$OUT/resume_manifest.json" "$SOURCE_JOB" "$ACCEPTOR" "$source_sha" <<'PY'
import json,sys
path,source_job,acceptor,source_sha=sys.argv[1:]
data={"schema_version":1,"status":"READY_FOR_P00_CONTINUATION","candidate_id":"nylc_C18_trueT267_freeGS","source_job_id":source_job,"acceptor_hypothesis":acceptor,"source_endpoint_sha256":source_sha,"continuation_rounds":["r01","r02","r03"],"interpretation":"Same-Hamiltonian minimization continuation; not a TS, committor, PMF, or barrier."}
open(path,"w").write(json.dumps(data,indent=2,sort_keys=True)+"\n")
PY

module purge >/dev/null 2>&1 || true
module load amber/2018-hpcx-gcc-7.3.1
export AMBERHOME="$AMBER_RUNTIME"
test -s "$AMBERHOME/dat/slko/3ob-3-1/C-C.skf"
check_amber() {
  local outfile=$1
  grep -q "FINAL RESULTS" "$outfile"
  grep -q "5.  TIMINGS" "$outfile"
  if grep -Eq "SCC is not converged|SANDER BOMB|FATAL|NaN|forrtl|segmentation" "$outfile"; then
    echo "Amber hard error in $outfile" >&2
    return 1
  fi
}

current="$SOURCE_RST"
cd "$OUT/run"
p00_pass=false
last_audit=
for round in r01 r02 r03; do
  prefix="p00_${round}"
  sander -O -i "../input/p00.in" -o "${prefix}.out" -p "$BASE/system.prmtop" -c "$current" -ref "$current" -r "${prefix}.rst7" -inf "${prefix}.mdinfo"
  check_amber "${prefix}.out"
  current="$OUT/run/${prefix}.rst7"
  last_audit="$OUT/run/${prefix}_audit.json"
  "$PY" "$FLOW/scripts/audit_nylc_step1_proton_stage.py" --manifest "$OUT/input/manifest.json" --prmtop "$BASE/system.prmtop" --rst7 "$current" --stage p00 --output "$last_audit"
  status=$("$PY" -c 'import json,sys;print(json.load(open(sys.argv[1]))["status"])' "$last_audit")
  if [[ "$status" == PASS_PROTON_BRACKET_STAGE ]]; then
    p00_pass=true
    break
  fi
done
if [[ "$p00_pass" != true ]]; then
  cp "$last_audit" "$OUT/SCIENTIFIC_FAIL.json"
  record FAIL_SCIENTIFIC "stage=p00_continuation;status=FAIL_SCIENTIFIC_PROTON_BRACKET_STAGE;audit=$OUT/SCIENTIFIC_FAIL.json"
  trap - ERR
  exit 0
fi

for stage in p01 p02 p03; do
  sander -O -i "../input/${stage}.in" -o "${stage}.out" -p "$BASE/system.prmtop" -c "$current" -ref "$current" -r "${stage}.rst7" -inf "${stage}.mdinfo"
  check_amber "${stage}.out"
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
