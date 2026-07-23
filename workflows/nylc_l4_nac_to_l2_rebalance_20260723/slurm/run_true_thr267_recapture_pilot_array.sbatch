#!/bin/bash
#SBATCH -J nylc_t267cap
#SBATCH -p xahdnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --gres=dcu:1
#SBATCH --mem-per-cpu=2000M
#SBATCH -t 01:00:00
#SBATCH --array=0-1%2
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/true_thr267_recapture_%A_%a.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/true_thr267_recapture_%A_%a.err

set -euo pipefail
# SUPERSEDED_WRONG_NUCLEOPHILE_IDENTITY results are forbidden inputs.
TASK=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
FLOW="$TASK/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
PY=/work/home/acshdt1dks/anaconda3/bin/python3.9
GMX=/public/software/apps/Gromacs-DCU2/2022.1/mpi/bin/gmx_mpi
candidates=(nylc_C18_trueT267_recapture nylc_C23_trueT267_recapture)
candidate="${candidates[${SLURM_ARRAY_TASK_ID:?}]}"
work="$TASK/candidates/$candidate/recapture_pilot"
job="${SLURM_ARRAY_JOB_ID:?}"
task="${SLURM_ARRAY_TASK_ID:?}"
mkdir -p "$TASK/logs"
source /work/home/acshdt1dks/opt/gromacs-fastest/env.sh
export GMX_MAXBACKUP=-1

record() {
  local state=$1 detail=$2 stamp
  stamp="$(date -Is)"
  (
    flock -x 9
    printf '%s\ttrue_thr267_recapture_pilot\t%s\tcandidate=%s;job=%s;task=%s;%s\n' \
      "$stamp" "$state" "$candidate" "$job" "$task" "$detail" >&9
  ) 9>>"$TASK/run_history.tsv"
  (
    flock -x 8
    "$PY" - "$stamp" "$state" "$candidate" "$job" "$task" "$detail" >&8 <<'PY'
import json,sys
stamp,state,candidate,job,task,detail=sys.argv[1:]
print(json.dumps({
    "timestamp":stamp,
    "event":"true_thr267_recapture_pilot",
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
for path in "$work/pilot.log" "$work/pilot.gro" "$work/pilot.xtc"; do
  if [[ -e "$path" ]]; then
    record BLOCKED "refuse_overwrite=$path"
    exit 43
  fi
done

record START "preflight=PASS_GROMPP_TRUE_THR267;restrained=yes;duration_ps=100"
set +e
mpirun -np 1 "$GMX" mdrun \
  -s "$work/pilot.tpr" -deffnm "$work/pilot" \
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
  --candidate "$candidate" >"$work/audit.stdout" 2>"$work/audit.stderr"
audit_rc=$?
set -e
audit_status="$("$PY" -c "import json;print(json.load(open('$work/pilot_audit.json'))['status'])" 2>/dev/null || echo MISSING_AUDIT)"
if [[ $audit_rc -ne 0 ]]; then
  record FAIL_GATE "audit_rc=$audit_rc;status=$audit_status"
  exit "$audit_rc"
fi
record PASS_TECHNICAL "status=$audit_status;scientific=NOT_EVALUATED_RESTRAINED"
