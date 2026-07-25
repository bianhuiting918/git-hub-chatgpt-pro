#!/bin/bash
#SBATCH -J nylc_m1_em12
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem-per-cpu=2500M
#SBATCH -t 02:00:00
#SBATCH --array=0-11%4
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-%x-%A_%a.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-%x-%A_%a.err

set -euo pipefail

TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
FLOW="$TASK_ROOT/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
SELECTION="$TASK_ROOT/ensemble/selection_job_61813011/ensemble_selection.json"
PYTHON=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
GMX=/work/home/acshdt1dks/opt/gmx-cp2k/gromacs-2023.1-cp2k-2024.1/bin/gmx_mpi_d
JOB_ID="${SLURM_ARRAY_JOB_ID:?SLURM_ARRAY_JOB_ID is required}"
ARRAY_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"
EVENT=nylc_m1_real_nac_em

source /work/home/acshdt1dks/opt/gmx-cp2k/env.sh
export GMX_MAXBACKUP=-1
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"

candidate_id="$("$PYTHON" - "$SELECTION" "$ARRAY_ID" <<'PY'
import json
import sys
selected = json.load(open(sys.argv[1]))["selected"]
print(selected[int(sys.argv[2])]["candidate_id"])
PY
)"
CANDIDATE_ROOT="$TASK_ROOT/ensemble/candidates/$candidate_id"
BUILD_ROOT="$CANDIDATE_ROOT/build_job_61813799_$ARRAY_ID"
BUILD="$BUILD_ROOT/build"
OUT="$CANDIDATE_ROOT/em_job_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
COMMAND="double_precision_flexible_water_EM candidate_index=$ARRAY_ID"

append_history() {
    local state="$1"
    local detail="$2"
    local timestamp
    timestamp="$(date '+%Y-%m-%dT%H:%M:%S%z')"
    (
        flock -x 9
        printf '%s\t%s\t%s\t%s\t%s\n' "$timestamp" "$EVENT" "$COMMAND" "$state" "$detail" >>"$TASK_ROOT/run_history.tsv"
        "$PYTHON" - "$timestamp" "$state" "$detail" "$JOB_ID" "$ARRAY_ID" "$candidate_id" "$OUT" >>"$TASK_ROOT/run_history.jsonl" <<'PY'
import json
import sys
timestamp, state, detail, job_id, array_id, candidate_id, output_root = sys.argv[1:]
print(json.dumps({
    "timestamp": timestamp,
    "event": "nylc_m1_real_nac_em",
    "state": state,
    "candidate_id": candidate_id,
    "microstate": "M1_NalphaH2_Asp306H",
    "slurm_job_id": job_id,
    "slurm_array_task_id": array_id,
    "output_root": output_root,
    "detail": detail,
}, sort_keys=True))
PY
    ) 9>>"$TASK_ROOT/.run_history.lock"
}

if [[ -e "$OUT" ]]; then
    append_history FAIL_TECHNICAL_EM "reason=output_exists;refuse_overwrite=yes"
    exit 2
fi
mkdir -p "$OUT"

if [[ ! -s "$BUILD_ROOT/PASS.json" ]]; then
    "$PYTHON" - "$OUT/NOT_EVALUATED.json" "$candidate_id" "$JOB_ID" "$ARRAY_ID" <<'PY'
import json
import sys
path, candidate_id, job_id, array_id = sys.argv[1:]
open(path, "w").write(json.dumps({
    "schema_version": 1,
    "status": "NOT_EVALUATED_BUILD_FAIL",
    "candidate_id": candidate_id,
    "slurm_job_id": job_id,
    "slurm_array_task_id": array_id,
    "reason": "candidate did not pass build_job_61813799",
}, indent=2, sort_keys=True) + "\n")
PY
    append_history NOT_EVALUATED_BUILD_FAIL "source_build_job=61813799;em_not_run=yes"
    exit 0
fi

append_history START "source_build_job=61813799;precision=double;integrator=cg;flexible_water=yes;scientific_status=NOT_EVALUATED"
(
    cd "$BUILD"
    "$GMX" grompp -f "$FLOW/mdp/em_cg_flexible_m1.mdp" -c system_M1.gro -r system_M1.gro -p topol.top -o "$OUT/run.tpr" -po "$OUT/run.processed.mdp" -maxwarn 0 >"$OUT/grompp.stdout" 2>"$OUT/grompp.stderr"
)
mpirun -np 1 "$GMX" mdrun -s "$OUT/run.tpr" -deffnm "$OUT/run" -ntomp "${SLURM_CPUS_PER_TASK:-8}" -nb cpu -pme cpu -bonded cpu -pin on >"$OUT/mdrun.stdout" 2>"$OUT/mdrun.stderr"

test -s "$OUT/run.gro"
fmax="$(awk '/Maximum force/{value=$4} END{print value}' "$OUT/run.log")"
test -n "$fmax"
hard=0
if grep -Eqi 'FATAL ERROR|LINCS WARNING|cannot be settled|SETTLE.*(warning|error|failed)|(^|[^[:alnum:]_])NaN([^[:alnum:]_]|$)' "$OUT/run.log" "$OUT/mdrun.stderr"; then
    hard=1
fi
status=FAIL_TECHNICAL_EM
if [[ "$hard" -eq 0 ]] && awk -v value="$fmax" 'BEGIN{exit !(value<=500.0)}'; then
    status=PASS_TECHNICAL_EM
fi

"$PYTHON" - "$OUT" "$candidate_id" "$JOB_ID" "$ARRAY_ID" "$status" "$fmax" "$hard" <<'PY'
import hashlib
import json
import pathlib
import sys
root, candidate_id, job_id, array_id, status, fmax, hard = sys.argv[1:]
root = pathlib.Path(root)
payload = {
    "schema_version": 1,
    "status": status,
    "candidate_id": candidate_id,
    "microstate": "M1_NalphaH2_Asp306H",
    "slurm_job_id": job_id,
    "slurm_array_task_id": array_id,
    "source_build_job_id": "61813799",
    "stage": "m1_em_cg_flexible_double",
    "integrator": "cg",
    "flexible_water": True,
    "precision": "double",
    "gromacs_version": "2023.1",
    "protein_and_l2_heavy_atom_restraints": True,
    "fmax_kj_mol_nm": float(fmax),
    "threshold_kj_mol_nm": 500.0,
    "hard_warning_count": int(hard),
    "scientific_nac_evidence": False,
    "run_gro_sha256": hashlib.sha256((root / "run.gro").read_bytes()).hexdigest(),
}
name = "PASS.json" if status == "PASS_TECHNICAL_EM" else "FAIL.json"
(root / name).write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
PY

if [[ "$status" != PASS_TECHNICAL_EM ]]; then
    append_history FAIL_TECHNICAL_EM "fmax=$fmax;hard_warning_count=$hard;outputs_preserved=yes"
    exit 31
fi
append_history PASS_TECHNICAL "fmax=$fmax;hard_warning_count=$hard;scientific_status=NOT_EVALUATED"
