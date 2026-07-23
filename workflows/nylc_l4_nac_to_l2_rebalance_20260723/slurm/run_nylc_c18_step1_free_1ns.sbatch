#!/bin/bash
#SBATCH -J nylc_c18_1ns
#SBATCH -p xahdnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --gres=dcu:1
#SBATCH --mem-per-cpu=3500M
#SBATCH -t 04:00:00
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/nylc_c18_free_1ns_%j.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/nylc_c18_free_1ns_%j.err

set -euo pipefail
echo "SUPERSEDED_WRONG_NUCLEOPHILE_IDENTITY: atom 8896 is Thr262 OG1, not catalytic Thr267 OG1" >&2
exit 42
TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
FLOW="$TASK_ROOT/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
PYTHON=/work/home/acshdt1dks/anaconda3/bin/python3.9
GMX=/public/software/apps/Gromacs-DCU2/2022.1/mpi/bin/gmx_mpi
candidate=nylc_c18_11854ps
parent_job=61687591
candidate_root="$TASK_ROOT/candidates/$candidate"
build="$candidate_root/build"
input_root="$candidate_root/runs_step1_gs/continuation_$parent_job/npt300_nac_free_pilot"
input_gro="$input_root/run.gro"
input_cpt="$input_root/run.cpt"
expected_gro_sha=bab252dff6f2ce377e023c2cf1d2441da32fb0f2e6795622b782109cc0e26305
expected_cpt_sha=0faa7e20ea5a4f571801c5bb4790777009556f18525df2cfdebf808f37363f57
job_id="${SLURM_JOB_ID:?job id is required}"
run_root="$candidate_root/runs_step1_gs/free_1ns_$job_id"
mkdir -p "$run_root" "$TASK_ROOT/logs"
source /work/home/acshdt1dks/opt/gromacs-fastest/env.sh
export GMX_MAXBACKUP=-1

test "$(sha256sum "$input_gro" | awk '{print $1}')" = "$expected_gro_sha"
test "$(sha256sum "$input_cpt" | awk '{print $1}')" = "$expected_cpt_sha"

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
  "workflow":"nylc_c18_step1_free_1ns"
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
    record nylc_c18_step1_free_1ns FAIL "stage=$stage;rc=$rc"
  fi
  exit "$rc"
}
trap on_exit EXIT

stage=grompp
record nylc_c18_step1_free_1ns START "parent_job=$parent_job;unrestrained=yes;duration_ps=1000"
(
  cd "$build"
  "$GMX" grompp -f "$FLOW/mdp/npt300free.mdp" -c "$input_gro" -t "$input_cpt" -p "$build/topol.top" -o "$run_root/run.tpr" -po "$run_root/run.processed.mdp" -maxwarn 0 >"$run_root/grompp.stdout" 2>"$run_root/grompp.stderr"
)
stage=mdrun
mpirun -np 1 "$GMX" mdrun -deffnm "$run_root/run" -ntomp "${SLURM_CPUS_PER_TASK:-8}" -nb gpu -pme gpu -bonded gpu -pin on >"$run_root/mdrun.stdout" 2>"$run_root/mdrun.stderr"
test -s "$run_root/run.gro"
test -s "$run_root/run.cpt"
scan_log "$run_root/run.log" "$run_root/numerical_audit.json"

"$PYTHON" - "$run_root/UNRESTRAINED_1NS_COMPLETE.json" "$candidate" "$job_id" "$parent_job" <<'PY'
import json,sys
path,candidate,job,parent=sys.argv[1:]
json.dump({
 "status":"UNRESTRAINED_1NS_COMPLETE",
 "candidate_id":candidate,
 "slurm_job_id":job,
 "parent_job_id":parent,
 "unrestrained_extension_ps":1000,
 "scientific_status":"NOT_EVALUATED_PENDING_TRAJECTORY_AUDIT"
},open(path,"w"),indent=2)
open(path,"a").write("\n")
PY
record nylc_c18_step1_free_1ns UNRESTRAINED_1NS_COMPLETE "scientific_audit_required"
failure=0
trap - EXIT
