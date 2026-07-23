#!/bin/bash
#SBATCH -J nylc_s1_pilot
#SBATCH -p xahdnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --gres=dcu:1
#SBATCH --mem-per-cpu=3500M
#SBATCH -t 02:00:00
#SBATCH --array=0-1%2
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/nylc_step1_pilot_%A_%a.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/nylc_step1_pilot_%A_%a.err

set -euo pipefail

TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
REPO="$TASK_ROOT/repo"
FLOW="$REPO/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
PYTHON=/work/home/acshdt1dks/anaconda3/bin/python3.9
GMX=/public/software/apps/Gromacs-DCU2/2022.1/mpi/bin/gmx_mpi
candidates=(nylc_c18_11854ps nylc_c23_29684ps)

task_index="${SLURM_ARRAY_TASK_ID:?array task id is required}"
candidate="${candidates[$task_index]}"
candidate_root="$TASK_ROOT/candidates/$candidate"
build="$candidate_root/build"
job_id="${SLURM_ARRAY_JOB_ID:?array job id is required}"
run_root="$candidate_root/runs_step1_gs/pilot_${job_id}_${task_index}"

mkdir -p "$TASK_ROOT/logs" "$run_root"
source /work/home/acshdt1dks/opt/gromacs-fastest/env.sh
export GMX_MAXBACKUP=-1

record() {
  local event=$1 state=$2 detail=$3
  local stamp
  stamp="$(date -Is)"
  (
    flock -x 9
    printf '%s\t%s\t%s\t%s\n' "$stamp" "$event" "$state" "$detail" >&9
  ) 9>>"$TASK_ROOT/run_history.tsv"
  (
    flock -x 8
    "$PYTHON" - "$stamp" "$event" "$state" "$detail" "$candidate" "$job_id" "$task_index" >&8 <<'PY'
import json
import sys
stamp, event, state, detail, candidate, job, task = sys.argv[1:]
print(json.dumps({
    "timestamp": stamp,
    "event": event,
    "state": state,
    "detail": detail,
    "candidate_id": candidate,
    "slurm_array_job_id": job,
    "slurm_array_task_id": task,
    "workflow": "nylc_step1_gs_pilot",
}, sort_keys=True))
PY
  ) 8>>"$TASK_ROOT/run_history.jsonl"
}

scan_log() {
  local log=$1 output=$2
  "$PYTHON" - "$log" "$output" <<'PY'
import json
import re
import sys
log_path, output_path = sys.argv[1:]
text = open(log_path, errors="replace").read()
patterns = {
    "fatal": r"(?i)fatal(?:\s+error)?",
    "lincs_warning": r"(?i)lincs\s+warning",
    "settle_problem": r"(?i)(?:cannot\s+be\s+settled|settle[^\n]*(?:warning|error|failed))",
    "nan": r"(?i)\bnan\b",
}
result = {name: len(re.findall(pattern, text)) for name, pattern in patterns.items()}
result["finished_mdrun"] = bool(re.search(r"Finished mdrun", text))
json.dump(result, open(output_path, "w"), indent=2)
open(output_path, "a").write("\n")
raise SystemExit(0 if result["finished_mdrun"] and not any(
    result[name] for name in ("fatal", "lincs_warning", "settle_problem", "nan")
) else 1)
PY
}

audit_geometry() {
  local gro=$1 output=$2 stage=$3
  "$PYTHON" - "$candidate_root/audit/build_audit.json" "$gro" "$output" "$stage" <<'PY'
import json
import math
import sys
audit_path, gro_path, output_path, stage = sys.argv[1:]
audit = json.load(open(audit_path))
atom_ids = audit["reactive_global_atoms"]
lines = open(gro_path).read().splitlines()
xyz = {}
for name, atom_id in atom_ids.items():
    line = lines[atom_id + 1]
    xyz[name] = tuple(float(line[start:end]) for start, end in ((20, 28), (28, 36), (36, 44)))
def sub(a, b):
    return tuple(x - y for x, y in zip(a, b))
def norm(v):
    return math.sqrt(sum(x * x for x in v))
attack = sub(xyz["thr_og1"], xyz["carbonyl_c"])
carbonyl = sub(xyz["carbonyl_o"], xyz["carbonyl_c"])
distance = norm(attack)
cosine = sum(x * y for x, y in zip(attack, carbonyl)) / (norm(attack) * norm(carbonyl))
angle = math.degrees(math.acos(max(-1.0, min(1.0, cosine))))
result = {
    "stage": stage,
    "distance_thr267_og1_to_carbonyl_c_nm": distance,
    "angle_thr267_og1_carbonyl_c_carbonyl_o_deg": angle,
    "nac_final_frame": distance <= 0.35 and 95.0 <= angle <= 115.0,
    "distance_gate_nm": 0.35,
    "angle_gate_deg": [95.0, 115.0],
}
json.dump(result, open(output_path, "w"), indent=2)
open(output_path, "a").write("\n")
PY
}

run_stage() {
  local stage=$1 coord=$2 checkpoint=${3:-}
  local stage_root="$run_root/$stage"
  mkdir -p "$stage_root"
  local args=(
    grompp -f "$FLOW/mdp/$stage.mdp"
    -c "$coord" -r "$build/rebuilt.gro" -p "$build/topol.top"
    -o "$stage_root/run.tpr" -po "$stage_root/run.processed.mdp" -maxwarn 0
  )
  if [[ -n $checkpoint ]]; then
    args+=( -t "$checkpoint" )
  fi
  record "$stage" START "cpus=${SLURM_CPUS_PER_TASK:-8};dcu=1"
  (
    cd "$build"
    "$GMX" "${args[@]}" >"$stage_root/grompp.stdout" 2>"$stage_root/grompp.stderr"
  )
  local mdrun_args=(
    mdrun -deffnm "$stage_root/run"
    -ntomp "${SLURM_CPUS_PER_TASK:-8}" -nb gpu -pin on
  )
  if [[ $stage == em_nac_contact ]]; then
    mdrun_args+=( -pme cpu -bonded cpu )
  else
    mdrun_args+=( -pme gpu -bonded gpu )
  fi
  mpirun -np 1 "$GMX" "${mdrun_args[@]}" >"$stage_root/mdrun.stdout" 2>"$stage_root/mdrun.stderr"
  test -s "$stage_root/run.gro"
  test -s "$stage_root/run.log"
  scan_log "$stage_root/run.log" "$stage_root/numerical_audit.json"
  audit_geometry "$stage_root/run.gro" "$stage_root/nac_final_frame.json" "$stage"
  record "$stage" PASS_TECHNICAL "restrained=yes"
}

failure_status=1
current_stage=initialization
on_exit() {
  local rc=$?
  if [[ $failure_status -ne 0 ]]; then
    "$PYTHON" - "$run_root/FAIL.json" "$candidate" "$current_stage" "$rc" <<'PY'
import json
import sys
path, candidate, stage, rc = sys.argv[1:]
json.dump({
    "status": "NOT_EVALUATED_TECHNICAL_PILOT_FAIL",
    "candidate_id": candidate,
    "failed_stage": stage,
    "returncode": int(rc),
    "scientific_nac_evidence": False,
}, open(path, "w"), indent=2)
open(path, "a").write("\n")
PY
    record nylc_step1_pilot FAIL "stage=$current_stage;rc=$rc"
  fi
  exit "$rc"
}
trap on_exit EXIT

current_stage=em_nac_contact
run_stage em_nac_contact "$build/rebuilt.gro"
em_log="$run_root/em_nac_contact/run.log"
maximum_force="$(awk '/Maximum force/{v=$4} END{print v}' "$em_log")"
potential_energy="$(awk '/Potential Energy/{v=$4} END{print v}' "$em_log")"
"$PYTHON" - "$run_root/em_nac_contact/em_metrics.json" "$maximum_force" "$potential_energy" <<'PY'
import json
import sys
path, maximum_force, potential_energy = sys.argv[1:]
json.dump({
    "status": "EM_CONVERGENCE_NOT_REQUIRED",
    "maximum_force_kj_mol_nm": float(maximum_force),
    "potential_energy_kj_mol": float(potential_energy),
    "acceptance_basis": "finished_without_FATAL_LINCS_SETTLE_NaN_and_produced_coordinates",
}, open(path, "w"), indent=2)
open(path, "a").write("\n")
PY

current_stage=nvt10_nac
run_stage nvt10_nac "$run_root/em_nac_contact/run.gro"

current_stage=nvt50_nac
run_stage nvt50_nac "$run_root/nvt10_nac/run.gro" "$run_root/nvt10_nac/run.cpt"

"$PYTHON" - "$run_root/nvt50_nac/nac_final_frame.json" <<'PY'
import json
import sys
raise SystemExit(0 if json.load(open(sys.argv[1]))["nac_final_frame"] else 32)
PY

"$PYTHON" - "$run_root/TECHNICAL_PILOT_PASS.json" "$candidate" "$job_id" "$task_index" <<'PY'
import json
import sys
path, candidate, job, task = sys.argv[1:]
json.dump({
    "status": "TECHNICAL_PILOT_PASS",
    "candidate_id": candidate,
    "slurm_array_job_id": job,
    "slurm_array_task_id": int(task),
    "restrained": True,
    "scientific_nac_evidence": False,
    "next_gate": "review restrained NAC and numerical audit before any warmer or unrestrained stage",
}, open(path, "w"), indent=2)
open(path, "a").write("\n")
PY
record nylc_step1_pilot TECHNICAL_PILOT_PASS "scientific_nac_evidence=false"
failure_status=0
trap - EXIT
