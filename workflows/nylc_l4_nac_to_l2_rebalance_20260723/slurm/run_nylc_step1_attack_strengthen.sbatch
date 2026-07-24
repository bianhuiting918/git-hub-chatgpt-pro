#!/bin/bash
#SBATCH -J nylc_s1_astr
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem-per-cpu=2500M
#SBATCH -t 08:00:00
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/step1_attack_strengthen_%j.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/step1_attack_strengthen_%j.err

set -euo pipefail
TASK=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
FLOW="$TASK/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
BASE="$TASK/candidates/nylc_C18_trueT267_freeGS/qmmm/step1_dftb3_preflight_post61710861_job61712692"
ACCEPTOR="${ACCEPTOR:?set ACCEPTOR=OD1_or_OD2}"
SOURCE_JOB="${SOURCE_JOB:?set SOURCE_JOB=failed_resume_job_id}"
if [[ "$ACCEPTOR" != OD1 && "$ACCEPTOR" != OD2 ]]; then echo "invalid ACCEPTOR=$ACCEPTOR" >&2; exit 2; fi
if [[ ! "$SOURCE_JOB" =~ ^[0-9]+$ ]]; then echo "invalid SOURCE_JOB=$SOURCE_JOB" >&2; exit 2; fi
SOURCE="$TASK/candidates/nylc_C18_trueT267_freeGS/qmmm/step1_proton_resume_${ACCEPTOR}_job_${SOURCE_JOB}"
SOURCE_FAIL="$SOURCE/SCIENTIFIC_FAIL.json"
SOURCE_RST="$SOURCE/run/p00_r03.rst7"
job="${SLURM_JOB_ID:?}"
OUT="$TASK/candidates/nylc_C18_trueT267_freeGS/qmmm/step1_attack_strengthen_${ACCEPTOR}_job_${job}"
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
AMBER_RUNTIME=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/nylc_gyaq_pa66_l2_nac_qmmm_20260723/amber_runtime

record() {
  local state=$1 detail=$2 stamp
  stamp="$(date -Is)"
  printf '%s\tstep1_attack_strengthen\t%s\tcandidate=nylc_C18_trueT267_freeGS;job=%s;acceptor=%s;source_job=%s;%s\n' "$stamp" "$state" "$job" "$ACCEPTOR" "$SOURCE_JOB" "$detail" >> "$TASK/run_history.tsv"
  "$PY" - "$stamp" "$state" "$detail" "$job" "$ACCEPTOR" "$SOURCE_JOB" >> "$TASK/run_history.jsonl" <<'PY'
import json,sys
stamp,state,detail,job,acceptor,source_job=sys.argv[1:]
print(json.dumps({"timestamp":stamp,"event":"step1_attack_strengthen","state":state,"detail":detail,"candidate_id":"nylc_C18_trueT267_freeGS","slurm_job_id":job,"acceptor_hypothesis":acceptor,"source_job_id":source_job},sort_keys=True))
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
source_sha="$(sha256sum "$SOURCE_RST" | awk '{print $1}')"
record START "source_endpoint_sha=$source_sha;change=attack_force_50_to_100_only;scope=constrained_seed_not_TS_PMF"
mkdir "$OUT"
mkdir "$OUT/run"
"$PY" "$FLOW/scripts/prepare_nylc_step1_attack_strengthen.py" --source-manifest "$SOURCE/input/manifest.json" --source-fail "$SOURCE_FAIL" --acceptor "$ACCEPTOR" --output-dir "$OUT/input" > "$OUT/prepare.stdout" 2> "$OUT/prepare.stderr"
for name in p01 p02 p03; do
  cp "$SOURCE/input/${name}.in" "$OUT/input/${name}.in"
  cp "$SOURCE/input/${name}.RST" "$OUT/input/${name}.RST"
done
cp "$SOURCE/input/manifest.json" "$OUT/input/proton_manifest.json"

"$PY" "$FLOW/scripts/audit_nylc_step1_proton_stage.py" --manifest "$OUT/input/manifest.json" --prmtop "$BASE/system.prmtop" --rst7 "$SOURCE_RST" --stage p00 --output "$OUT/source_reaudit.json"
"$PY" - "$OUT/source_reaudit.json" <<'PY'
import json,sys
data=json.load(open(sys.argv[1]))
g=data["gates"]
if not g["proton_attached"] or not g["attack_angle_valid"] or not g["chemistry_valid"]:
    raise SystemExit("source failed safe geometry re-audit")
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

cd "$OUT/run"
sander -O -i "../input/p00.in" -o "p00_strong.out" -p "$BASE/system.prmtop" -c "$SOURCE_RST" -ref "$SOURCE_RST" -r "p00_strong.rst7" -inf "p00_strong.mdinfo"
check_amber "p00_strong.out"
current="$OUT/run/p00_strong.rst7"
audit="$OUT/run/p00_strong_audit.json"
"$PY" "$FLOW/scripts/audit_nylc_step1_proton_stage.py" --manifest "$OUT/input/manifest.json" --prmtop "$BASE/system.prmtop" --rst7 "$current" --stage p00 --output "$audit"
status=$("$PY" -c 'import json,sys;print(json.load(open(sys.argv[1]))["status"])' "$audit")
if [[ "$status" != PASS_PROTON_BRACKET_STAGE ]]; then
  cp "$audit" "$OUT/SCIENTIFIC_FAIL.json"
  record FAIL_SCIENTIFIC "stage=p00_strong;status=$status;audit=$OUT/SCIENTIFIC_FAIL.json"
  trap - ERR
  exit 0
fi

for stage in p01 p02 p03; do
  sander -O -i "../input/${stage}.in" -o "${stage}.out" -p "$BASE/system.prmtop" -c "$current" -ref "$current" -r "${stage}.rst7" -inf "${stage}.mdinfo"
  check_amber "${stage}.out"
  current="$OUT/run/${stage}.rst7"
  audit="$OUT/run/${stage}_audit.json"
  "$PY" "$FLOW/scripts/audit_nylc_step1_proton_stage.py" --manifest "$OUT/input/proton_manifest.json" --prmtop "$BASE/system.prmtop" --rst7 "$current" --stage "$stage" --output "$audit"
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
