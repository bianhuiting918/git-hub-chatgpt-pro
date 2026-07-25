#!/bin/bash
#SBATCH -J nylc_m1_stageA_recover
#SBATCH -p xahdnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --gres=dcu:1
#SBATCH --mem-per-cpu=2500M
#SBATCH -t 06:00:00
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-%x-%j.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-%x-%j.err

set -euo pipefail

TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
FLOW="$TASK_ROOT/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
SELECTION="$TASK_ROOT/ensemble/selection_job_61813011/ensemble_selection.json"
PYTHON=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
GMX=/public/software/apps/Gromacs-DCU2/2022.1/mpi/bin/gmx_mpi
EM_ARRAY_JOB_ID="${EM_ARRAY_JOB_ID:?EM_ARRAY_JOB_ID is required}"
FAILED_STAGEA_ARRAY_JOB_ID="${FAILED_STAGEA_ARRAY_JOB_ID:?FAILED_STAGEA_ARRAY_JOB_ID is required}"
RECOVERY_TASK_ID="${RECOVERY_TASK_ID:?RECOVERY_TASK_ID is required}"
JOB_ID="${SLURM_JOB_ID:?SLURM_JOB_ID is required}"
case "$RECOVERY_TASK_ID" in
    *[!0-9]*|'') echo "RECOVERY_TASK_ID must be an integer" >&2; exit 2 ;;
esac
if [[ "$RECOVERY_TASK_ID" -lt 0 || "$RECOVERY_TASK_ID" -gt 35 ]]; then
    echo "RECOVERY_TASK_ID must be within 0-35" >&2
    exit 2
fi
candidate_index=$((RECOVERY_TASK_ID / 3))
seed_index=$((RECOVERY_TASK_ID % 3))
seeds=(26711 26723 26737)
seed="${seeds[$seed_index]}"
candidate_id="$("$PYTHON" - "$SELECTION" "$candidate_index" <<'PY'
import json
import sys
print(json.load(open(sys.argv[1]))["selected"][int(sys.argv[2])]["candidate_id"])
PY
)"
CANDIDATE_ROOT="$TASK_ROOT/ensemble/candidates/$candidate_id"
BUILD="$CANDIDATE_ROOT/build_job_61813799_$candidate_index/build"
RUN_ROOT="$CANDIDATE_ROOT/replicas/seed$seed/stageA_resume_job_${FAILED_STAGEA_ARRAY_JOB_ID}_${RECOVERY_TASK_ID}"
EVENT=nylc_m1_real_nac_stageA_checkpoint_recovery
COMMAND="recover completed nvt50 checkpoint and finish Stage A candidate_index=$candidate_index seed=$seed"

set +u
source /work/home/acshdt1dks/opt/gromacs-fastest/env.sh
set -u
export GMX_MAXBACKUP=-1
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"

append_history() {
    local state="$1"
    local detail="$2"
    local timestamp
    timestamp="$(date '+%Y-%m-%dT%H:%M:%S%z')"
    (
        flock -x 9
        printf '%s\t%s\t%s\t%s\t%s\n' "$timestamp" "$EVENT" "$COMMAND" "$state" "$detail" >>"$TASK_ROOT/run_history.tsv"
        "$PYTHON" - "$timestamp" "$state" "$detail" "$JOB_ID" "$FAILED_STAGEA_ARRAY_JOB_ID" "$RECOVERY_TASK_ID" "$candidate_id" "$seed" "$RUN_ROOT" >>"$TASK_ROOT/run_history.jsonl" <<'PY'
import json
import sys
timestamp, state, detail, job_id, failed_job, task_id, candidate, seed, root = sys.argv[1:]
print(json.dumps({
    "timestamp": timestamp,
    "event": "nylc_m1_real_nac_stageA_checkpoint_recovery",
    "state": state,
    "detail": detail,
    "slurm_job_id": job_id,
    "failed_stageA_array_job_id": failed_job,
    "recovery_task_id": int(task_id),
    "candidate_id": candidate,
    "velocity_seed": int(seed),
    "run_root": root,
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
    0 if audit["finished_mdrun"]
    and not any(audit[key] for key in ("fatal", "lincs_warning", "settle_problem", "nan"))
    else 1
)
PY
}

if [[ ! -d "$RUN_ROOT" ]]; then
    append_history FAIL_TECHNICAL_RECOVERY "reason=failed_run_root_missing"
    exit 2
fi
[[ ! -e "$RUN_ROOT/STAGE_A_COMPLETE.json" ]]
if [[ -e "$RUN_ROOT/nvt50/RECOVERY_AUDIT.json" ]]; then
    append_history FAIL_TECHNICAL_RECOVERY "reason=recovery_already_attempted;refuse_overwrite=yes"
    exit 2
fi
for stage in nvt150 nvt300 npt300r npt300rel npt300free_stageA; do
    if [[ -e "$RUN_ROOT/$stage" ]]; then
        append_history FAIL_TECHNICAL_RECOVERY "reason=refuse_existing_later_stage;stage=$stage"
        exit 2
    fi
done
for name in run.gro run.cpt run.edr run.xtc run.log run.tpr; do
    test -s "$RUN_ROOT/nvt50/$name"
done
grep -q "Finished mdrun" "$RUN_ROOT/nvt50/run.log"
scan_log "$RUN_ROOT/nvt50/run.log" "$RUN_ROOT/nvt50/recovery_numerical_audit.json"
"$PYTHON" - "$RUN_ROOT/nvt50/RECOVERY_AUDIT.json" "$RUN_ROOT/nvt50/run.gro" "$RUN_ROOT/nvt50/run.cpt" "$JOB_ID" "$FAILED_STAGEA_ARRAY_JOB_ID" "$RECOVERY_TASK_ID" <<'PY'
import hashlib
import json
import pathlib
import sys
output, gro, cpt, job_id, failed_job, task_id = sys.argv[1:]
def digest(path):
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()
payload = {
    "schema_version": 1,
    "status": "PASS_AUDITED_CHECKPOINT_RECOVERY_SOURCE",
    "recovery_job_id": job_id,
    "failed_stageA_array_job_id": failed_job,
    "failed_stageA_task_id": int(task_id),
    "source_exit_code": 132,
    "source_failure_point": "post-Finished-mdrun HIP cleanup SIGILL",
    "source_gro_sha256": digest(gro),
    "source_cpt_sha256": digest(cpt),
    "scientific_evidence": False,
}
open(output, "w").write(json.dumps(payload, indent=2, sort_keys=True) + "\n")
PY
append_history START "validated_nvt50_checkpoint=yes;source_exit_code=132;source_outputs_preserved=yes"

run_stage() {
    local name="$1"
    local mdp="$2"
    local coord="$3"
    local checkpoint="$4"
    local stage="$RUN_ROOT/$name"
    mkdir "$stage"
    local args=(grompp -f "$mdp" -c "$coord" -p "$BUILD/topol.top" -o "$stage/run.tpr" -po "$stage/run.processed.mdp" -maxwarn 0 -r "$BUILD/system_M1.gro")
    if [[ -n "$checkpoint" ]]; then
        args+=(-t "$checkpoint")
    fi
    (
        cd "$BUILD"
        "$GMX" "${args[@]}" >"$stage/grompp.stdout" 2>"$stage/grompp.stderr"
    )
    mpirun -np 1 "$GMX" mdrun -s "$stage/run.tpr" -deffnm "$stage/run"         -ntomp "${SLURM_CPUS_PER_TASK:-8}" -nb gpu -pme gpu -bonded gpu -pin on         >"$stage/mdrun.stdout" 2>"$stage/mdrun.stderr"
    test -s "$stage/run.gro"
    test -s "$stage/run.cpt"
    scan_log "$stage/run.log" "$stage/numerical_audit.json"
    "$PYTHON" - "$stage/PASS.json" "$candidate_id" "$seed" "$name" "$JOB_ID" <<'PY'
import json
import sys
path, candidate, seed, stage, job_id = sys.argv[1:]
open(path, "w").write(json.dumps({
    "schema_version": 1,
    "status": "PASS_TECHNICAL_RECOVERY_STAGE",
    "candidate_id": candidate,
    "velocity_seed": int(seed),
    "stage": stage,
    "recovery_job_id": job_id,
    "restrained": True,
    "scientific_nac_evidence": False,
}, indent=2, sort_keys=True) + "\n")
PY
}

failed=1
current=nvt150
trap 'rc=$?; if [[ $failed -ne 0 ]]; then append_history FAIL_TECHNICAL_RECOVERY "stage=$current;exit_code=$rc;outputs_preserved=yes"; fi; exit $rc' EXIT

run_stage nvt150 "$FLOW/mdp/nvt150_m1.mdp" "$RUN_ROOT/nvt50/run.gro" "$RUN_ROOT/nvt50/run.cpt"
current=nvt300
run_stage nvt300 "$FLOW/mdp/nvt300_m1.mdp" "$RUN_ROOT/nvt150/run.gro" "$RUN_ROOT/nvt150/run.cpt"
current=npt300r
run_stage npt300r "$FLOW/mdp/npt300r_m1.mdp" "$RUN_ROOT/nvt300/run.gro" "$RUN_ROOT/nvt300/run.cpt"
current=npt300rel
run_stage npt300rel "$FLOW/mdp/npt300rel_m1.mdp" "$RUN_ROOT/npt300r/run.gro" "$RUN_ROOT/npt300r/run.cpt"

current=free_tpr_contract
FREE="$RUN_ROOT/npt300free_stageA"
mkdir "$FREE"
(
    cd "$BUILD"
    "$GMX" grompp -f "$FLOW/mdp/npt300free_m1_stageA.mdp"         -c "$RUN_ROOT/npt300rel/run.gro" -t "$RUN_ROOT/npt300rel/run.cpt"         -p "$BUILD/topol.top" -o "$FREE/run.tpr" -po "$FREE/run.processed.mdp"         -maxwarn 0 >"$FREE/grompp.stdout" 2>"$FREE/grompp.stderr"
)
"$GMX" dump -s "$FREE/run.tpr" >"$FREE/tpr.dump" 2>"$FREE/tpr_dump.stderr"
"$PYTHON" - "$FREE/tpr.dump" "$FREE/run.processed.mdp" "$FREE/free_tpr_contract.json" <<'PY'
import json
import re
import sys
dump_path, mdp_path, output_path = sys.argv[1:]
dump = open(dump_path, errors="replace").read()
mdp = open(mdp_path, errors="replace").read()
position_restraints = sum(int(value) for value in re.findall(r"#posres_xA\s*=\s*(\d+)", dump))
distance_restraints = len(re.findall(r"(?i)\bDISRES\b", dump))
define_values = []
for line in mdp.splitlines():
    body = line.split(";", 1)[0].strip()
    if "=" not in body:
        continue
    key, value = body.split("=", 1)
    if key.strip().lower() == "define" and value.strip():
        define_values.append(value.strip())
payload = {
    "position_restraints": position_restraints,
    "distance_restraints": distance_restraints,
    "processed_mdp_nonempty_define_values": define_values,
    "scientific_window": "fully_unrestrained_NPT_100ps",
}
open(output_path, "w").write(json.dumps(payload, indent=2, sort_keys=True) + "\n")
if position_restraints != 0 or distance_restraints != 0 or define_values:
    raise SystemExit("free Stage A recovery TPR contains restraints")
PY

current=npt300free_stageA
mpirun -np 1 "$GMX" mdrun -s "$FREE/run.tpr" -deffnm "$FREE/run"     -ntomp "${SLURM_CPUS_PER_TASK:-8}" -nb gpu -pme gpu -bonded gpu -pin on     >"$FREE/mdrun.stdout" 2>"$FREE/mdrun.stderr"
test -s "$FREE/run.gro"
test -s "$FREE/run.cpt"
scan_log "$FREE/run.log" "$FREE/numerical_audit.json"
"$PYTHON" - "$RUN_ROOT/STAGE_A_COMPLETE.json" "$candidate_id" "$seed" "$JOB_ID" "$FAILED_STAGEA_ARRAY_JOB_ID" "$RECOVERY_TASK_ID" <<'PY'
import hashlib
import json
import pathlib
import sys
path, candidate, seed, job_id, failed_job, task_id = sys.argv[1:]
root = pathlib.Path(path).parent
free = root / "npt300free_stageA"
payload = {
    "schema_version": 1,
    "status": "PASS_TECHNICAL_STAGE_A_RECOVERED",
    "scientific_status": "NOT_EVALUATED_PENDING_STAGE_A_AUDIT",
    "candidate_id": candidate,
    "velocity_seed": int(seed),
    "slurm_job_id": job_id,
    "recovered_failed_stageA_array_job_id": failed_job,
    "slurm_array_task_id": task_id,
    "free_window_ps": [0.0, 100.0],
    "fully_unrestrained": True,
    "free_tpr_contract": str(free / "free_tpr_contract.json"),
    "final_gro_sha256": hashlib.sha256((free / "run.gro").read_bytes()).hexdigest(),
}
open(path, "w").write(json.dumps(payload, indent=2, sort_keys=True) + "\n")
PY
append_history PASS_TECHNICAL "recovered_failed_stageA_array=$FAILED_STAGEA_ARRAY_JOB_ID;fully_unrestrained_100ps_complete=yes;scientific_status=NOT_EVALUATED_PENDING_STAGE_A_AUDIT"
failed=0
trap - EXIT
