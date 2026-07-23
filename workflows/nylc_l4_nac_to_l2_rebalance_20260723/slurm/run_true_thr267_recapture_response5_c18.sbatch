#!/bin/bash
#SBATCH -J nylc_t267r5
#SBATCH -p xahdnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --gres=dcu:1
#SBATCH --mem-per-cpu=2000M
#SBATCH -t 01:00:00
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/true_thr267_response5_%j.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/true_thr267_response5_%j.err

set -euo pipefail
# SUPERSEDED_WRONG_NUCLEOPHILE_IDENTITY results are forbidden inputs.
TASK=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
FLOW="$TASK/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
PY=/work/home/acshdt1dks/anaconda3/bin/python3.9
GMX=/public/software/apps/Gromacs-DCU2/2022.1/mpi/bin/gmx_mpi
candidate=nylc_C18_trueT267_recapture
work="$TASK/candidates/$candidate/recapture_response5"
job="${SLURM_JOB_ID:?}"
task=0
mkdir -p "$TASK/logs"
source /work/home/acshdt1dks/opt/gromacs-fastest/env.sh
export GMX_MAXBACKUP=-1

record() {
  local state=$1 detail=$2 stamp
  stamp="$(date -Is)"
  (
    flock -x 9
    printf '%s\ttrue_thr267_recapture_response5\t%s\tcandidate=%s;job=%s;task=%s;%s\n' \
      "$stamp" "$state" "$candidate" "$job" "$task" "$detail" >&9
  ) 9>>"$TASK/run_history.tsv"
  (
    flock -x 8
    "$PY" - "$stamp" "$state" "$candidate" "$job" "$task" "$detail" >&8 <<'PY'
import json,sys
stamp,state,candidate,job,task,detail=sys.argv[1:]
print(json.dumps({
    "timestamp":stamp,
    "event":"true_thr267_recapture_response5",
    "state":state,
    "candidate":candidate,
    "slurm_array_job_id":job,
    "slurm_array_task_id":task,
    "detail":detail,
},sort_keys=True))
PY
  ) 8>>"$TASK/run_history.jsonl"
}

preflight_status="$("$PY" -c "import json;print(json.load(open('$work/preflight.json'))['status'])")"
if [[ "$preflight_status" != PASS_GROMPP_TRUE_THR267 ]]; then
  record BLOCKED "preflight_status=$preflight_status"
  exit 42
fi
for path in "$work/response5.log" "$work/response5.gro" "$work/response5.xtc"; do
  if [[ -e "$path" ]]; then
    record BLOCKED "refuse_overwrite=$path"
    exit 43
  fi
done

record START "preflight=PASS_GROMPP_TRUE_THR267;restrained=yes;duration_ps=100;protocol=response5"
set +e
mpirun -np 1 "$GMX" mdrun \
  -s "$work/response5.tpr" -deffnm "$work/response5" \
  -ntomp "${SLURM_CPUS_PER_TASK:-8}" -nb gpu -pme gpu -bonded gpu -pin on \
  >"$work/mdrun.stdout" 2>"$work/mdrun.stderr"
rc=$?
set -e
if [[ $rc -ne 0 ]]; then
  record TECHNICAL_FAIL "mdrun_rc=$rc"
  exit "$rc"
fi

set +e
"$PY" "$FLOW/scripts/audit_true_thr267_recapture_pilot.py" \
  --candidate "$candidate" --stage-dir recapture_response5 --stem response5 --min-response 0.20 --reference-gro "$TASK/candidates/$candidate/recapture_response3/response3.gro" >"$work/audit.stdout" 2>"$work/audit.stderr"
audit_rc=$?
set -e
audit_status="$("$PY" -c "import json;print(json.load(open('$work/pilot_audit.json'))['status'])" 2>/dev/null || echo MISSING_AUDIT)"
if [[ $audit_rc -ne 0 ]]; then
  record FAIL_GATE "audit_rc=$audit_rc;status=$audit_status"
  exit "$audit_rc"
fi
record PASS_TECHNICAL "status=$audit_status;scientific=NOT_EVALUATED_RESTRAINED"
