#!/bin/bash
#SBATCH -J l4tol2_rebal
#SBATCH -p xahdnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --gres=dcu:1
#SBATCH --mem=28G
#SBATCH -t 24:00:00
#SBATCH --array=0-4%5
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/rebalance_%A_%a.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/rebalance_%A_%a.err

set -euo pipefail
TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
REPO="$TASK_ROOT/repo"
FLOW="$REPO/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
PYTHON=/work/home/acshdt1dks/anaconda3/bin/python3.9
GMX=/public/software/apps/Gromacs-DCU2/2022.1/mpi/bin/gmx_mpi
mkdir -p "$TASK_ROOT/logs"
source /work/home/acshdt1dks/opt/gromacs-fastest/env.sh
export GMX_MAXBACKUP=-1

mapfile -t candidates < <(
  "$PYTHON" - "$FLOW/manifests/candidates.json" <<'PY'
import json
import sys
for row in json.load(open(sys.argv[1]))["candidates"]:
    print(row["id"])
PY
)
task_index="${SLURM_ARRAY_TASK_ID:?array task id is required}"
candidate="${candidates[$task_index]}"
candidate_root="$TASK_ROOT/candidates/$candidate"
build="$candidate_root/build"
preflight="$candidate_root/preflight"
run_root="$candidate_root/runs"
mkdir -p "$run_root"

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
    "$PYTHON" - "$stamp" "$event" "$state" "$detail" "$candidate" \
      "${SLURM_ARRAY_JOB_ID:-NA}" "$task_index" >&8 <<'PY'
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
}, sort_keys=True))
PY
  ) 8>>"$TASK_ROOT/run_history.jsonl"
}

current_stage=initialization
finished=0
fail_trap() {
  rc=$?
  if [[ $finished -eq 0 ]]; then
    "$PYTHON" - "$candidate_root/FAIL.json" "$candidate" "$current_stage" "$rc" <<'PY'
import json
import sys
path, candidate, stage, rc = sys.argv[1:]
json.dump({
    "status": "NOT_EVALUATED_STAGE_FAIL",
    "candidate_id": candidate,
    "failed_stage": stage,
    "returncode": int(rc),
    "reason": "inspect stage logs; no automatic redocking or retry",
}, open(path, "w"), indent=2)
open(path, "a").write("\n")
PY
    record rebalance "FAIL_STAGE" "stage=$current_stage;rc=$rc"
  fi
  exit "$rc"
}
trap fail_trap EXIT

if ! "$PYTHON" - "$preflight/grompp_preflight.json" <<'PY'
import json
import sys
raise SystemExit(0 if json.load(open(sys.argv[1]))["state"] == "GROMPP_PASS" else 1)
PY
then
  current_stage=preflight
  record rebalance NOT_EVALUATED_GROMPP_FAIL "preflight did not pass"
  finished=1
  trap - EXIT
  exit 20
fi

scan_stage() {
  local log=$1
  local report=$2
  "$PYTHON" - "$log" "$report" <<'PY'
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
counts = {name: len(re.findall(pattern, text)) for name, pattern in patterns.items()}
counts["finished_mdrun"] = bool(re.search(r"Finished mdrun", text))
json.dump(counts, open(output_path, "w"), indent=2)
open(output_path, "a").write("\n")
raise SystemExit(0 if counts["finished_mdrun"] and not any(
    counts[name] for name in ("fatal", "lincs_warning", "settle_problem", "nan")
) else 1)
PY
}

run_stage() {
  local name=$1
  local coord=$2
  local checkpoint=${3:-}
  local restrained=${4:-yes}
  local stage="$run_root/$name"
  current_stage=$name
  mkdir -p "$stage"
  if [[ -s "$stage/PASS.json" ]]; then
    record "$name" SKIP_IMMUTABLE_PASS "existing PASS.json"
    return 0
  fi
  local args=(
    grompp -f "$FLOW/mdp/$name.mdp"
    -c "$coord" -p "$build/topol.top"
    -o "$stage/run.tpr" -po "$stage/run.processed.mdp" -maxwarn 0
  )
  if [[ $restrained == yes ]]; then
    args+=( -r "$build/rebuilt.gro" )
  fi
  if [[ -n $checkpoint ]]; then
    args+=( -t "$checkpoint" )
  fi
  record "$name" START "cpus=${SLURM_CPUS_PER_TASK:-8};dcu=1"
  (
    cd "$build"
    "$GMX" "${args[@]}" >"$stage/grompp.stdout" 2>"$stage/grompp.stderr"
  )
  local mdrun_args=(
    mdrun -deffnm "$stage/run"
    -ntomp "${SLURM_CPUS_PER_TASK:-8}" -nb gpu -pin on
  )
  if [[ $name == em ]]; then
    mdrun_args+=( -pme cpu -bonded cpu )
  else
    mdrun_args+=( -pme gpu -bonded gpu )
  fi
  mpirun -np 1 "$GMX" "${mdrun_args[@]}" >"$stage/mdrun.stdout" 2>"$stage/mdrun.stderr"
  test -s "$stage/run.gro"
  test -s "$stage/run.log"
  scan_stage "$stage/run.log" "$stage/numerical_audit.json"
  if [[ $name == em ]]; then
    fmax="$(awk '/Maximum force/{value=$4} END{print value}' "$stage/run.log")"
    test -n "$fmax"
    awk -v value="$fmax" 'BEGIN { exit !(value <= 500.0) }'
  else
    fmax=""
  fi
  "$PYTHON" - "$stage/PASS.json" "$candidate" "$name" "$restrained" "$fmax" <<'PY'
import json
import sys
path, candidate, stage, restrained, fmax = sys.argv[1:]
json.dump({
    "status": "PASS_TECHNICAL_STAGE",
    "candidate_id": candidate,
    "stage": stage,
    "restrained": restrained == "yes",
    "scientific_nac_evidence": stage == "npt300free" and restrained == "no",
    "maximum_force_kj_mol_nm": float(fmax) if fmax else None,
}, open(path, "w"), indent=2)
open(path, "a").write("\n")
PY
  record "$name" PASS "restrained=$restrained"
}

run_stage em "$build/rebuilt.gro" "" yes
run_stage nvt50 "$run_root/em/run.gro" "" yes
run_stage nvt150 "$run_root/nvt50/run.gro" "$run_root/nvt50/run.cpt" yes
run_stage nvt300 "$run_root/nvt150/run.gro" "$run_root/nvt150/run.cpt" yes
run_stage npt300r "$run_root/nvt300/run.gro" "$run_root/nvt300/run.cpt" yes
run_stage npt300rel "$run_root/npt300r/run.gro" "$run_root/npt300r/run.cpt" yes
run_stage npt300free "$run_root/npt300rel/run.gro" "$run_root/npt300rel/run.cpt" no

"$PYTHON" - "$candidate_root/TECHNICAL_COMPLETE.json" "$candidate" <<'PY'
import json
import sys
path, candidate = sys.argv[1:]
json.dump({
    "status": "TECHNICAL_COMPLETE_PENDING_UNCONSTRAINED_SCIENTIFIC_AUDIT",
    "candidate_id": candidate,
    "required_scientific_window": "runs/npt300free",
}, open(path, "w"), indent=2)
open(path, "a").write("\n")
PY
record rebalance TECHNICAL_COMPLETE "free_1ns_finished;scientific_audit_pending"
finished=1
trap - EXIT
