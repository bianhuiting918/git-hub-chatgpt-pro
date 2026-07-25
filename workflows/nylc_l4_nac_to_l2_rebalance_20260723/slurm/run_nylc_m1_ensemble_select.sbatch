#!/bin/bash
#SBATCH -J nylc_m1_nac_select
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem-per-cpu=2500M
#SBATCH -t 02:00:00
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-%x-%j.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-%x-%j.err

set -euo pipefail

TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
FLOW="$TASK_ROOT/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
AUTHORITY="$FLOW/manifests/nylc_C18_m1_real_nac_ensemble.authority.json"
PYTHON=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
GMX=/public/software/apps/Gromacs-DCU2/2022.1/mpi/bin/gmx_mpi
JOB_ID="${SLURM_JOB_ID:?SLURM_JOB_ID is required}"
OUTPUT_ROOT="$TASK_ROOT/ensemble/selection_job_${JOB_ID}"
EVENT=nylc_m1_real_nac_ensemble_select
COMMAND="srun -n 1 -c 8 python select_nylc_m1_nac_ensemble.py --authority AUTHORITY --output-root OUTPUT_ROOT"

source /work/home/acshdt1dks/opt/gromacs-fastest/env.sh
export GMX_MAXBACKUP=-1
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"

append_history() {
    local state="$1"
    local detail="$2"
    local timestamp
    timestamp="$(date '+%Y-%m-%dT%H:%M:%S%z')"
    (
        flock -x 9
        printf '%s\t%s\t%s\t%s\t%s\n'             "$timestamp" "$EVENT" "$COMMAND" "$state" "$detail"             >>"$TASK_ROOT/run_history.tsv"
        "$PYTHON" - "$timestamp" "$state" "$detail" "$JOB_ID" "$OUTPUT_ROOT"             >>"$TASK_ROOT/run_history.jsonl" <<'PY'
import json
import sys

timestamp, state, detail, job_id, output_root = sys.argv[1:]
print(json.dumps({
    "timestamp": timestamp,
    "event": "nylc_m1_real_nac_ensemble_select",
    "state": state,
    "candidate_id": "nylc_C18_trueT267_freeGS",
    "microstate": "M1_NalphaH2_Asp306H",
    "slurm_job_id": job_id,
    "output_root": output_root,
    "detail": detail,
}, sort_keys=True))
PY
    ) 9>>"$TASK_ROOT/.run_history.lock"
}

if [[ -e "$OUTPUT_ROOT" ]]; then
    append_history FAIL "reason=output_exists;refuse_overwrite=yes"
    exit 2
fi
mkdir -p "$TASK_ROOT/ensemble"
append_history START "source_job=61780674;selection_audit_job=61801089;expected_events=72;expected_selected=12"

set +e
srun -n 1 -c "${SLURM_CPUS_PER_TASK:-8}" "$PYTHON"     "$FLOW/scripts/select_nylc_m1_nac_ensemble.py"     --authority "$AUTHORITY"     --output-root "$OUTPUT_ROOT"     --gmx "$GMX"
rc=$?
set -e

if [[ "$rc" -ne 0 ]]; then
    append_history FAIL "exit_code=$rc;reason=selection_script_failed;outputs_preserved=yes"
    exit "$rc"
fi
if [[ ! -s "$OUTPUT_ROOT/ensemble_selection.json" ]]; then
    append_history FAIL "exit_code=3;reason=missing_ensemble_selection_json;outputs_preserved=yes"
    exit 3
fi

append_history PASS "exit_code=0;selected=12;manifest=$OUTPUT_ROOT/ensemble_selection.json;scientific_status=SOURCE_SELECTION_ONLY_NOT_MD_VALIDATED"
