#!/bin/bash
#SBATCH -J nylc_s1_bracket
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem-per-cpu=2500M
#SBATCH -t 12:00:00
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/step1_bracket_%j.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/step1_bracket_%j.err

set -euo pipefail
TASK=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
FLOW="$TASK/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
SRC="$TASK/candidates/nylc_C18_trueT267_freeGS/qmmm/step1_dftb3_preflight_post61710861_job61712692"
ACCEPTOR="${ACCEPTOR:?set ACCEPTOR=OD1_or_OD2}"
if [[ "$ACCEPTOR" != OD1 && "$ACCEPTOR" != OD2 ]]; then echo "invalid ACCEPTOR=$ACCEPTOR" >&2; exit 2; fi
job="${SLURM_JOB_ID:?}"
OUT="$TASK/candidates/nylc_C18_trueT267_freeGS/qmmm/step1_constrained_bracket_${ACCEPTOR}_job_${job}"
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
AMBER_RUNTIME=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/nylc_gyaq_pa66_l2_nac_qmmm_20260723/amber_runtime

record() {
  local state=$1 detail=$2 stamp
  stamp="$(date -Is)"
  printf '%s\tstep1_constrained_bracket\t%s\tcandidate=nylc_C18_trueT267_freeGS;job=%s;acceptor=%s;%s\n' "$stamp" "$state" "$job" "$ACCEPTOR" "$detail" >> "$TASK/run_history.tsv"
  "$PY" - "$stamp" "$state" "$detail" "$job" "$ACCEPTOR" >> "$TASK/run_history.jsonl" <<'PY'
import json,sys
stamp,state,detail,job,acceptor=sys.argv[1:]
print(json.dumps({"timestamp":stamp,"event":"step1_constrained_bracket","state":state,"detail":detail,"candidate_id":"nylc_C18_trueT267_freeGS","slurm_job_id":job,"acceptor_hypothesis":acceptor},sort_keys=True))
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

record START "source_job=61712692;scope=constrained_tetrahedral_seed_not_TS_PMF"
test -s "$SRC/system.prmtop"
test -s "$SRC/selected_free_l2_nac.rst7"
mkdir "$OUT"
mkdir "$OUT/run"
"$PY" "$FLOW/scripts/prepare_nylc_step1_constrained_bracket.py" --source-dir "$SRC" --acceptor "$ACCEPTOR" --output-dir "$OUT/input" > "$OUT/prepare.stdout" 2> "$OUT/prepare.stderr"

module purge >/dev/null 2>&1 || true
module load amber/2018-hpcx-gcc-7.3.1
export AMBERHOME="$AMBER_RUNTIME"
test -s "$AMBERHOME/dat/slko/3ob-3-1/C-C.skf"
current="$SRC/selected_free_l2_nac.rst7"
cd "$OUT/run"
for window in w00 w01 w02 w03; do
  sander -O -i "../input/${window}.in" -o "${window}.out" -p "$SRC/system.prmtop" -c "$current" -r "${window}.rst7" -inf "${window}.mdinfo"
  current="$OUT/run/${window}.rst7"
done
cd "$OUT"
"$PY" "$FLOW/scripts/audit_nylc_step1_constrained_bracket.py" --prmtop "$SRC/system.prmtop" --input-dir "$OUT/input" --run-dir "$OUT/run" --output "$OUT/bracket_result.json"
technical=$("$PY" -c 'import json,sys;print(json.load(open(sys.argv[1]))["technical_status"])' "$OUT/bracket_result.json")
scientific=$("$PY" -c 'import json,sys;print(json.load(open(sys.argv[1]))["scientific_status"])' "$OUT/bracket_result.json")
if [[ "$technical" != PASS_TECHNICAL_STEP1_BRACKET ]]; then
  record FAIL_TECHNICAL "technical=$technical;scientific=$scientific;audit=$OUT/bracket_result.json"
  exit 31
fi
cp "$OUT/bracket_result.json" "$OUT/PASS.json"
record PASS_TECHNICAL "technical=$technical;scientific=$scientific;audit=$OUT/PASS.json;interpretation=not_TS_committor_PMF"
