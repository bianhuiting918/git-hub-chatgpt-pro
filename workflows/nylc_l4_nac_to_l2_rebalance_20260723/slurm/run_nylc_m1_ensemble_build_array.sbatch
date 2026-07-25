#!/bin/bash
#SBATCH -J nylc_m1_build12
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem-per-cpu=2500M
#SBATCH -t 01:00:00
#SBATCH --array=0-11%4
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-%x-%A_%a.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-%x-%A_%a.err

set -euo pipefail

TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
FLOW="$TASK_ROOT/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
AUTHORITY="$FLOW/manifests/nylc_C18_m1_real_nac_ensemble.authority.json"
SELECTION="$TASK_ROOT/ensemble/selection_job_61813011/ensemble_selection.json"
LEGACY_BUILD="$TASK_ROOT/candidates/nylc_C18_trueT267_freeGS/build"
PYTHON=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
GMX=/public/software/apps/Gromacs-DCU2/2022.1/mpi/bin/gmx_mpi
JOB_ID="${SLURM_ARRAY_JOB_ID:?SLURM_ARRAY_JOB_ID is required}"
ARRAY_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"
EVENT=nylc_m1_real_nac_candidate_build

source /work/home/acshdt1dks/opt/gromacs-fastest/env.sh
export GMX_MAXBACKUP=-1
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"

candidate_id="$("$PYTHON" - "$SELECTION" "$ARRAY_ID" <<'PY'
import json
import sys

selection_path, array_id = sys.argv[1:]
selected = json.load(open(selection_path))["selected"]
index = int(array_id)
if index < 0 or index >= len(selected):
    raise SystemExit("array index is outside the frozen selection")
print(selected[index]["candidate_id"])
PY
)"
CANDIDATE_DIR="$TASK_ROOT/ensemble/candidates/$candidate_id/build_job_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
COMMAND="python prepare_nylc_m1_nac_ensemble.py --selection selection_job_61813011 --candidate-index $ARRAY_ID"

append_history() {
    local state="$1"
    local detail="$2"
    local timestamp
    timestamp="$(date '+%Y-%m-%dT%H:%M:%S%z')"
    (
        flock -x 9
        printf '%s\t%s\t%s\t%s\t%s\n' "$timestamp" "$EVENT" "$COMMAND" "$state" "$detail" >>"$TASK_ROOT/run_history.tsv"
        "$PYTHON" - "$timestamp" "$state" "$detail" "$JOB_ID" "$ARRAY_ID" "$candidate_id" "$CANDIDATE_DIR" >>"$TASK_ROOT/run_history.jsonl" <<'PY'
import json
import sys

timestamp, state, detail, job_id, array_id, candidate_id, candidate_dir = sys.argv[1:]
print(json.dumps({
    "timestamp": timestamp,
    "event": "nylc_m1_real_nac_candidate_build",
    "state": state,
    "candidate_id": candidate_id,
    "microstate": "M1_NalphaH2_Asp306H",
    "slurm_job_id": job_id,
    "slurm_array_task_id": array_id,
    "candidate_dir": candidate_dir,
    "detail": detail,
}, sort_keys=True))
PY
    ) 9>>"$TASK_ROOT/.run_history.lock"
}

if [[ -e "$CANDIDATE_DIR" ]]; then
    append_history FAIL_TECHNICAL_BUILD "reason=output_exists;refuse_overwrite=yes"
    exit 2
fi
append_history START "selection_job=61813011;candidate_index=$ARRAY_ID;source_is_real_unrestrained_NAC=yes"

set +e
"$PYTHON" "$FLOW/scripts/prepare_nylc_m1_nac_ensemble.py" --authority "$AUTHORITY" --selection "$SELECTION" --candidate-index "$ARRAY_ID" --candidate-dir "$CANDIDATE_DIR" --legacy-build-root "$LEGACY_BUILD" --gmx "$GMX" --em-mdp "$FLOW/mdp/em.mdp"
rc=$?
set -e

if [[ "$rc" -ne 0 ]]; then
    append_history FAIL_TECHNICAL_BUILD "exit_code=$rc;outputs_preserved=yes"
    exit "$rc"
fi
if [[ ! -s "$CANDIDATE_DIR/PASS.json" || ! -s "$CANDIDATE_DIR/build_audit.json" ]]; then
    append_history FAIL_TECHNICAL_BUILD "exit_code=3;reason=missing_pass_or_audit;outputs_preserved=yes"
    exit 3
fi

append_history PASS_TECHNICAL "exit_code=0;grompp_maxwarn=0;next_gate=double_precision_flexible_water_EM"
