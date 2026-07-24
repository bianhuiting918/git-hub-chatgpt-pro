#!/bin/bash
#SBATCH -J nylc_ash306_full
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem-per-cpu=2500M
#SBATCH -t 01:00:00
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/ash306_full_%j.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/ash306_full_%j.err

set -euo pipefail
TASK=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
FLOW="$TASK/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
SOURCE="$TASK/candidates/nylc_C18_trueT267_freeGS/qmmm/step1_dftb3_preflight_post61710861_job61712692/selected_free_l2_nac.gro"
SOURCE_BUILD="$TASK/candidates/nylc_C18_trueT267_freeGS/build"
PROBE="$TASK/candidates/nylc_C18_trueT267_freeGS/qmmm/ash306_chain_probe_job_61717760"
job="${SLURM_JOB_ID:?}"
OUT="$TASK/candidates/nylc_C18_trueT267_freeGS/qmmm/ash306_full_system_preflight_job_${job}"
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
GMX=/public/software/apps/Gromacs-DCU2/2022.1/mpi/bin/gmx_mpi

record() {
  local state=$1 detail=$2 stamp
  stamp="$(date -Is)"
  printf '%s\tash306_full_system_preflight\t%s\tcandidate=nylc_C18_trueT267_freeGS;job=%s;%s\n' "$stamp" "$state" "$job" "$detail" >> "$TASK/run_history.tsv"
  "$PY" - "$stamp" "$state" "$detail" "$job" >> "$TASK/run_history.jsonl" <<'PY'
import json,sys
stamp,state,detail,job=sys.argv[1:]
print(json.dumps({"timestamp":stamp,"event":"ash306_full_system_preflight","state":state,"detail":detail,"candidate_id":"nylc_C18_trueT267_freeGS","slurm_job_id":job},sort_keys=True))
PY
}

on_error() {
  local rc=$?
  trap - ERR
  set +e
  if [[ -d "$OUT" && ! -s "$OUT/FAIL.json" ]]; then
    "$PY" - "$rc" "$OUT" > "$OUT/FAIL.json" <<'PY'
import json,sys
print(json.dumps({"schema_version":1,"candidate_id":"nylc_C18_trueT267_freeGS","status":"FAIL_TECHNICAL_ASH306_FULL_SYSTEM_PREFLIGHT","exit_code":int(sys.argv[1]),"output_dir":sys.argv[2]},indent=2,sort_keys=True))
PY
  fi
  record FAIL_TECHNICAL "exit_code=$rc;out=$OUT"
  exit "$rc"
}
trap on_error ERR

record START "source=434ps_GS_hypothesis;probe_job=61717760;scope=build_grompp_parmed_only"
mkdir "$OUT"
"$PY" "$FLOW/scripts/build_nylc_ash306_full_system.py" \
  --source-gro "$SOURCE" --source-build-dir "$SOURCE_BUILD" --probe-dir "$PROBE" \
  --output-dir "$OUT/build" > "$OUT/build.stdout" 2> "$OUT/build.stderr"

set +u
source /work/home/acshdt1dks/opt/gromacs-fastest/env.sh
set -u
export GMXDATA=/public/software/apps/Gromacs-DCU2/2022.1/mpi/share/gromacs
export GROMACS_TOPDIR="$GMXDATA/top"
export GMX_MAXBACKUP=-1
(
  cd "$OUT/build"
  "$GMX" grompp -f grompp_preflight.mdp -c system.gro -p topol.top \
    -o preflight.tpr -po preflight.processed.mdp -maxwarn 0 > "$OUT/grompp.log" 2>&1
)
"$PY" "$FLOW/scripts/audit_nylc_ash306_full_system_preflight.py" \
  --build-dir "$OUT/build" --probe-pass "$PROBE/PASS.json" \
  --grompp-log "$OUT/grompp.log" --amber-output-dir "$OUT/amber" \
  --output "$OUT/preflight_audit.json" > "$OUT/audit.stdout" 2> "$OUT/audit.stderr"
status=$("$PY" -c 'import json,sys;print(json.load(open(sys.argv[1]))["status"])' "$OUT/preflight_audit.json")
if [[ "$status" != PASS_ASH306_FULL_SYSTEM_PREFLIGHT ]]; then
  cp "$OUT/preflight_audit.json" "$OUT/FAIL.json"
  record FAIL_SCIENTIFIC_GATE "status=$status;audit=$OUT/FAIL.json"
  exit 31
fi
cp "$OUT/preflight_audit.json" "$OUT/PASS.json"
record PASS_TECHNICAL "status=$status;audit=$OUT/PASS.json;interpretation=not_TS_RC_PMF"
