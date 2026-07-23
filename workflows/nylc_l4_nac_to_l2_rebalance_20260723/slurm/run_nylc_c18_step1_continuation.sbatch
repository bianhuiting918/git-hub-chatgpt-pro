#!/bin/bash
#SBATCH -J nylc_c18_s1
#SBATCH -p xahdnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --gres=dcu:1
#SBATCH --mem-per-cpu=3500M
#SBATCH -t 08:00:00
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/nylc_c18_step1_%j.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/nylc_c18_step1_%j.err

set -euo pipefail
TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
FLOW="$TASK_ROOT/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
PYTHON=/work/home/acshdt1dks/anaconda3/bin/python3.9
GMX=/public/software/apps/Gromacs-DCU2/2022.1/mpi/bin/gmx_mpi
candidate=nylc_c18_11854ps
candidate_root="$TASK_ROOT/candidates/$candidate"
build="$candidate_root/build"
selected="$candidate_root/selected_step1_gs/c18_50k_late_lowest_potential_nac_9p55ps.gro"
expected_sha=c67cbb1d275863606be62628df6829e8ef15fbad39f3b5a51188c311f95ce235
job_id="${SLURM_JOB_ID:?job id is required}"
run_root="$candidate_root/runs_step1_gs/continuation_$job_id"
mkdir -p "$run_root" "$TASK_ROOT/logs"
source /work/home/acshdt1dks/opt/gromacs-fastest/env.sh
export GMX_MAXBACKUP=-1

test "$(sha256sum "$selected" | awk '{print $1}')" = "$expected_sha"

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
    "$PYTHON" - "$stamp" "$event" "$state" "$detail" "$candidate" "$job_id" >&8 <<'PY'
import json,sys
stamp,event,state,detail,candidate,job=sys.argv[1:]
print(json.dumps({
  "timestamp":stamp,"event":event,"state":state,"detail":detail,
  "candidate_id":candidate,"slurm_job_id":job,
  "workflow":"nylc_c18_step1_continuation"
},sort_keys=True))
PY
  ) 8>>"$TASK_ROOT/run_history.jsonl"
}

scan_log() {
  "$PYTHON" - "$1" "$2" <<'PY'
import json,re,sys
text=open(sys.argv[1],errors="replace").read()
result={
 "fatal":len(re.findall(r"(?i)fatal(?:\s+error)?",text)),
 "lincs_warning":len(re.findall(r"(?i)lincs\s+warning",text)),
 "settle_problem":len(re.findall(r"(?i)(?:cannot\s+be\s+settled|settle[^\n]*(?:warning|error|failed))",text)),
 "nan":len(re.findall(r"(?i)\bnan\b",text)),
 "finished_mdrun":bool(re.search(r"Finished mdrun",text)),
}
json.dump(result,open(sys.argv[2],"w"),indent=2);open(sys.argv[2],"a").write("\n")
raise SystemExit(0 if result["finished_mdrun"] and not any(result[k] for k in ("fatal","lincs_warning","settle_problem","nan")) else 1)
PY
}

run_stage() {
  local stage=$1 coord=$2 checkpoint=${3:-} unrestrained=${4:-no}
  local root="$run_root/$stage"
  mkdir -p "$root"
  local args=(grompp -f "$FLOW/mdp/$stage.mdp" -c "$coord" -p "$build/topol.top" -o "$root/run.tpr" -po "$root/run.processed.mdp" -maxwarn 0)
  if [[ $unrestrained == no ]]; then args+=(-r "$selected"); fi
  if [[ -n $checkpoint ]]; then args+=(-t "$checkpoint"); fi
  record "$stage" START "unrestrained=$unrestrained"
  (cd "$build" && "$GMX" "${args[@]}" >"$root/grompp.stdout" 2>"$root/grompp.stderr")
  mpirun -np 1 "$GMX" mdrun -deffnm "$root/run" -ntomp "${SLURM_CPUS_PER_TASK:-8}" -nb gpu -pme gpu -bonded gpu -pin on >"$root/mdrun.stdout" 2>"$root/mdrun.stderr"
  test -s "$root/run.gro"
  test -s "$root/run.cpt"
  scan_log "$root/run.log" "$root/numerical_audit.json"
  record "$stage" PASS_TECHNICAL "unrestrained=$unrestrained"
}

failure=1
stage=initialization
on_exit() {
  rc=$?
  if [[ $failure -ne 0 ]]; then
    "$PYTHON" - "$run_root/FAIL.json" "$candidate" "$stage" "$rc" <<'PY'
import json,sys
path,candidate,stage,rc=sys.argv[1:]
json.dump({"status":"NOT_EVALUATED_STAGE_FAIL","candidate_id":candidate,"failed_stage":stage,"returncode":int(rc)},open(path,"w"),indent=2)
open(path,"a").write("\n")
PY
    record nylc_c18_step1 FAIL "stage=$stage;rc=$rc"
  fi
  exit "$rc"
}
trap on_exit EXIT

stage=nvt100_nac
run_stage nvt100_nac "$selected"
stage=nvt150_nac
run_stage nvt150_nac "$run_root/nvt100_nac/run.gro" "$run_root/nvt100_nac/run.cpt"
stage=nvt300_nac_release
run_stage nvt300_nac_release "$run_root/nvt150_nac/run.gro" "$run_root/nvt150_nac/run.cpt"
stage=npt300_nac_free_pilot
run_stage npt300_nac_free_pilot "$run_root/nvt300_nac_release/run.gro" "$run_root/nvt300_nac_release/run.cpt" yes

"$PYTHON" - "$run_root/UNRESTRAINED_100PS_PILOT_COMPLETE.json" "$candidate" "$job_id" <<'PY'
import json,sys
path,candidate,job=sys.argv[1:]
json.dump({
 "status":"UNRESTRAINED_100PS_PILOT_COMPLETE",
 "candidate_id":candidate,
 "slurm_job_id":job,
 "unrestrained_pilot_ps":100,
 "required_final_unrestrained_ps":1000,
 "scientific_status":"NOT_EVALUATED_PENDING_TRAJECTORY_AUDIT_AND_1NS_EXTENSION"
},open(path,"w"),indent=2)
open(path,"a").write("\n")
PY
record nylc_c18_step1 UNRESTRAINED_100PS_PILOT_COMPLETE "audit_required_before_1ns"
failure=0
trap - EXIT
