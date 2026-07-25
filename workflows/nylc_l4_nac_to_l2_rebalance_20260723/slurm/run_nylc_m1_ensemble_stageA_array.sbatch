#!/bin/bash
#SBATCH -J nylc_m1_stageA
#SBATCH -p xahdnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --gres=dcu:1
#SBATCH --mem-per-cpu=2500M
#SBATCH -t 06:00:00
#SBATCH --array=0-35%8
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-%x-%A_%a.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-%x-%A_%a.err

set -euo pipefail

TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
FLOW="$TASK_ROOT/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
SELECTION="$TASK_ROOT/ensemble/selection_job_61813011/ensemble_selection.json"
PYTHON=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
GMX=/public/software/apps/Gromacs-DCU2/2022.1/mpi/bin/gmx_mpi
EM_ARRAY_JOB_ID="${EM_ARRAY_JOB_ID:?submit with --export=ALL,EM_ARRAY_JOB_ID=<completed EM array>}"
JOB_ID="${SLURM_ARRAY_JOB_ID:?SLURM_ARRAY_JOB_ID is required}"
TASK_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"
candidate_index=$((SLURM_ARRAY_TASK_ID / 3))
seed_index=$((SLURM_ARRAY_TASK_ID % 3))
seeds=(26711 26723 26737)
seed="${seeds[$seed_index]}"
EVENT=nylc_m1_real_nac_stageA

set +u
source /work/home/acshdt1dks/opt/gromacs-fastest/env.sh
set -u
export GMX_MAXBACKUP=-1
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"

candidate_id="$("$PYTHON" - "$SELECTION" "$candidate_index" <<'PY'
import json
import sys
selected = json.load(open(sys.argv[1]))["selected"]
print(selected[int(sys.argv[2])]["candidate_id"])
PY
)"
CANDIDATE_ROOT="$TASK_ROOT/ensemble/candidates/$candidate_id"
BUILD_ROOT="$CANDIDATE_ROOT/build_job_61813799_$candidate_index"
BUILD="$BUILD_ROOT/build"
EM_ROOT="$CANDIDATE_ROOT/em_job_${EM_ARRAY_JOB_ID}_$candidate_index"
RUN_ROOT="$CANDIDATE_ROOT/replicas/seed$seed/stageA_job_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
COMMAND="M1 release plus 100ps free StageA candidate_index=$candidate_index seed=$seed"

append_history() {
    local state="$1"
    local detail="$2"
    local timestamp
    timestamp="$(date '+%Y-%m-%dT%H:%M:%S%z')"
    (
        flock -x 9
        printf '%s\t%s\t%s\t%s\t%s\n' "$timestamp" "$EVENT" "$COMMAND" "$state" "$detail" >>"$TASK_ROOT/run_history.tsv"
        "$PYTHON" - "$timestamp" "$state" "$detail" "$JOB_ID" "$TASK_ID" "$candidate_id" "$seed" "$RUN_ROOT" >>"$TASK_ROOT/run_history.jsonl" <<'PY'
import json
import sys
timestamp, state, detail, job_id, task_id, candidate_id, seed, run_root = sys.argv[1:]
print(json.dumps({
    "timestamp": timestamp,
    "event": "nylc_m1_real_nac_stageA",
    "state": state,
    "candidate_id": candidate_id,
    "microstate": "M1_NalphaH2_Asp306H",
    "velocity_seed": int(seed),
    "slurm_job_id": job_id,
    "slurm_array_task_id": task_id,
    "run_root": run_root,
    "detail": detail,
}, sort_keys=True))
PY
    ) 9>>"$TASK_ROOT/.run_history.lock"
}

scan_log() {
    "$PYTHON" - "$1" "$2" <<'PY'
import json
import re
import sys
text = open(sys.argv[1], errors="replace").read()
audit = {
    "fatal": len(re.findall(r"(?i)fatal(?:\s+error)?", text)),
    "lincs_warning": len(re.findall(r"(?i)lincs\s+warning", text)),
    "settle_problem": len(re.findall(r"(?i)(?:cannot\s+be\s+settled|settle[^\n]*(?:warning|error|failed))", text)),
    "nan": len(re.findall(r"(?i)\bnan\b", text)),
    "finished_mdrun": bool(re.search(r"Finished mdrun", text)),
}
open(sys.argv[2], "w").write(json.dumps(audit, indent=2, sort_keys=True) + "\n")
raise SystemExit(
    0
    if audit["finished_mdrun"]
    and not any(audit[key] for key in ("fatal", "lincs_warning", "settle_problem", "nan"))
    else 1
)
PY
}

if [[ -e "$RUN_ROOT" ]]; then
    append_history FAIL_TECHNICAL_STAGEA "reason=output_exists;refuse_overwrite=yes"
    exit 2
fi
mkdir -p "$RUN_ROOT"

if [[ ! -s "$EM_ROOT/PASS.json" ]]; then
    "$PYTHON" - "$RUN_ROOT/NOT_EVALUATED.json" "$candidate_id" "$seed" "$EM_ARRAY_JOB_ID" <<'PY'
import json
import sys
path, candidate_id, seed, em_job = sys.argv[1:]
open(path, "w").write(json.dumps({
    "schema_version": 1,
    "status": "NOT_EVALUATED_EM_FAIL",
    "candidate_id": candidate_id,
    "velocity_seed": int(seed),
    "parent_em_array_job_id": em_job,
    "reason": "candidate lacks PASS_TECHNICAL_EM",
}, indent=2, sort_keys=True) + "\n")
PY
    append_history NOT_EVALUATED_EM_FAIL "parent_em_array=$EM_ARRAY_JOB_ID;stageA_not_run=yes"
    exit 0
fi

"$PYTHON" "$FLOW/scripts/render_nylc_m1_replica_mdp.py" --template "$FLOW/mdp/nvt50_m1.mdp" --output "$RUN_ROOT/nvt50.rendered.mdp" --set "gen-seed=$seed"
append_history START "parent_em_array=$EM_ARRAY_JOB_ID;seed=$seed;stages=restrained_release_then_free100ps;scientific_status=NOT_EVALUATED"

run_stage() {
    local name="$1"
    local mdp="$2"
    local coord="$3"
    local checkpoint="$4"
    local restrained="$5"
    local stage="$RUN_ROOT/$name"
    mkdir -p "$stage"
    local args=(grompp -f "$mdp" -c "$coord" -p "$BUILD/topol.top" -o "$stage/run.tpr" -po "$stage/run.processed.mdp" -maxwarn 0)
    if [[ "$restrained" == yes ]]; then
        args+=(-r "$BUILD/system_M1.gro")
    fi
    if [[ -n "$checkpoint" ]]; then
        args+=(-t "$checkpoint")
    fi
    (
        cd "$BUILD"
        "$GMX" "${args[@]}" >"$stage/grompp.stdout" 2>"$stage/grompp.stderr"
    )
    mpirun -np 1 "$GMX" mdrun -s "$stage/run.tpr" -deffnm "$stage/run" -ntomp "${SLURM_CPUS_PER_TASK:-8}" -nb gpu -pme gpu -bonded gpu -pin on >"$stage/mdrun.stdout" 2>"$stage/mdrun.stderr"
    test -s "$stage/run.gro"
    test -s "$stage/run.cpt"
    scan_log "$stage/run.log" "$stage/numerical_audit.json"
    "$PYTHON" - "$stage/PASS.json" "$candidate_id" "$seed" "$name" "$restrained" <<'PY'
import json
import sys
path, candidate_id, seed, stage, restrained = sys.argv[1:]
open(path, "w").write(json.dumps({
    "schema_version": 1,
    "status": "PASS_TECHNICAL_STAGE",
    "candidate_id": candidate_id,
    "velocity_seed": int(seed),
    "stage": stage,
    "restrained": restrained == "yes",
    "scientific_nac_evidence": stage == "npt300free_stageA" and restrained == "no",
}, indent=2, sort_keys=True) + "\n")
PY
}

failed=1
current=initialization
trap 'rc=$?; if [[ $failed -ne 0 ]]; then append_history FAIL_TECHNICAL_STAGEA "stage=$current;exit_code=$rc;outputs_preserved=yes"; fi; exit $rc' EXIT

current=nvt50
run_stage nvt50 "$RUN_ROOT/nvt50.rendered.mdp" "$EM_ROOT/run.gro" "" yes
current=nvt150
run_stage nvt150 "$FLOW/mdp/nvt150_m1.mdp" "$RUN_ROOT/nvt50/run.gro" "$RUN_ROOT/nvt50/run.cpt" yes
current=nvt300
run_stage nvt300 "$FLOW/mdp/nvt300_m1.mdp" "$RUN_ROOT/nvt150/run.gro" "$RUN_ROOT/nvt150/run.cpt" yes
current=npt300r
run_stage npt300r "$FLOW/mdp/npt300r_m1.mdp" "$RUN_ROOT/nvt300/run.gro" "$RUN_ROOT/nvt300/run.cpt" yes
current=npt300rel
run_stage npt300rel "$FLOW/mdp/npt300rel_m1.mdp" "$RUN_ROOT/npt300r/run.gro" "$RUN_ROOT/npt300r/run.cpt" yes

current=free_tpr_contract
FREE="$RUN_ROOT/npt300free_stageA"
mkdir -p "$FREE"
(
    cd "$BUILD"
    "$GMX" grompp -f "$FLOW/mdp/npt300free_m1_stageA.mdp" -c "$RUN_ROOT/npt300rel/run.gro" -t "$RUN_ROOT/npt300rel/run.cpt" -p "$BUILD/topol.top" -o "$FREE/run.tpr" -po "$FREE/run.processed.mdp" -maxwarn 0 >"$FREE/grompp.stdout" 2>"$FREE/grompp.stderr"
)
"$GMX" dump -s "$FREE/run.tpr" >"$FREE/tpr.dump" 2>"$FREE/tpr_dump.stderr"
"$PYTHON" - "$FREE/tpr.dump" "$FREE/run.processed.mdp" "$FREE/free_tpr_contract.json" <<'PY'
import json
import re
import sys

dump_path, mdp_path, output_path = sys.argv[1:]
dump = open(dump_path, errors="replace").read()
mdp = open(mdp_path, errors="replace").read()
position_counts = [int(value) for value in re.findall(r"#posres_xA\s*=\s*(\d+)", dump)]
position_restraints = sum(position_counts)
distance_restraints = len(re.findall(r"(?i)\bDISRES\b", dump))
define_present = any(
    line.split(";", 1)[0].strip().lower().startswith("define")
    for line in mdp.splitlines()
)
payload = {
    "position_restraints": 0,
    "distance_restraints": 0,
    "scientific_window": "fully_unrestrained_NPT_100ps",
    "observed_position_restraints": position_restraints,
    "observed_distance_restraint_tokens": distance_restraints,
    "processed_mdp_define_present": define_present,
}
open(output_path, "w").write(json.dumps(payload, indent=2, sort_keys=True) + "\n")
if position_restraints != 0 or distance_restraints != 0 or define_present:
    raise SystemExit("free Stage A TPR contains restraints")
PY

current=npt300free_stageA
mpirun -np 1 "$GMX" mdrun -s "$FREE/run.tpr" -deffnm "$FREE/run" -ntomp "${SLURM_CPUS_PER_TASK:-8}" -nb gpu -pme gpu -bonded gpu -pin on >"$FREE/mdrun.stdout" 2>"$FREE/mdrun.stderr"
test -s "$FREE/run.gro"
test -s "$FREE/run.cpt"
scan_log "$FREE/run.log" "$FREE/numerical_audit.json"
"$PYTHON" - "$RUN_ROOT/STAGE_A_COMPLETE.json" "$candidate_id" "$seed" "$JOB_ID" "$TASK_ID" "$EM_ARRAY_JOB_ID" <<'PY'
import hashlib
import json
import pathlib
import sys
path, candidate_id, seed, job_id, task_id, em_job = sys.argv[1:]
root = pathlib.Path(path).parent
free = root / "npt300free_stageA"
payload = {
    "schema_version": 1,
    "status": "PASS_TECHNICAL_STAGE_A",
    "scientific_status": "NOT_EVALUATED_PENDING_STAGE_A_AUDIT",
    "candidate_id": candidate_id,
    "velocity_seed": int(seed),
    "slurm_job_id": job_id,
    "slurm_array_task_id": task_id,
    "parent_em_array_job_id": em_job,
    "free_window_ps": [0.0, 100.0],
    "fully_unrestrained": True,
    "free_tpr_contract": str(free / "free_tpr_contract.json"),
    "final_gro_sha256": hashlib.sha256((free / "run.gro").read_bytes()).hexdigest(),
}
open(path, "w").write(json.dumps(payload, indent=2, sort_keys=True) + "\n")
PY

append_history PASS_TECHNICAL "seed=$seed;fully_unrestrained_100ps_complete=yes;scientific_status=NOT_EVALUATED_PENDING_STAGE_A_AUDIT"
failed=0
trap - EXIT
