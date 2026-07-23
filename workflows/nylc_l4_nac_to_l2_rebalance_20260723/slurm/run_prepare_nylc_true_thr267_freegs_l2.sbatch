#!/bin/bash
#SBATCH -J nylc_t267_l2
#SBATCH -p xahdnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --gres=dcu:1
#SBATCH --mem-per-cpu=3500M
#SBATCH -t 00:30:00
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/true_thr267_l2_build_%j.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/true_thr267_l2_build_%j.err

set -euo pipefail
TASK=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
FLOW="$TASK/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
PY=/work/home/acshdt1dks/anaconda3/bin/python3.9
GMX=/public/software/apps/Gromacs-DCU2/2022.1/mpi/bin/gmx_mpi
candidate=nylc_C18_trueT267_freeGS
root="$TASK/candidates/$candidate"
build="$root/build"
job="${SLURM_JOB_ID:?}"
source /work/home/acshdt1dks/opt/gromacs-fastest/env.sh
export GMX_MAXBACKUP=-1

record() {
  local state=$1 detail=$2 stamp
  stamp="$(date -Is)"
  (
    flock -x 9
    printf '%s\tnylc_true_thr267_freegs_l2_build\t%s\tjob=%s;%s\n' "$stamp" "$state" "$job" "$detail" >&9
  ) 9>>"$TASK/run_history.tsv"
  (
    flock -x 8
    "$PY" - "$stamp" "$state" "$detail" "$job" >&8 <<'PY'
import json,sys
stamp,state,detail,job=sys.argv[1:]
print(json.dumps({"timestamp":stamp,"event":"nylc_true_thr267_freegs_l2_build","state":state,"detail":detail,"slurm_job_id":job,"candidate_id":"nylc_C18_trueT267_freeGS"},sort_keys=True))
PY
  ) 8>>"$TASK/run_history.jsonl"
}
failure=1
stage=initialization
trap 'rc=$?; if [[ $failure -ne 0 ]]; then record FAIL "stage=$stage;rc=$rc"; fi; exit $rc' EXIT
record START "source=rare_unrestrained_true_Thr267_NAC"
stage=graph_build
"$PY" "$FLOW/scripts/prepare_nylc_true_thr267_freegs_l2.py" >"$root/build.stdout.json"
stage=grompp_preflight
(
  cd "$build"
  "$GMX" grompp -f "$FLOW/mdp/em.mdp" -c rebuilt.gro -r rebuilt.gro -p topol.top \
    -o em_preflight.tpr -po em_preflight.processed.mdp -maxwarn 0 \
    >grompp_preflight.stdout 2>grompp_preflight.stderr
)
"$PY" - "$root/audit/build_audit.json" "$root/audit/preflight.json" "$job" <<'PY'
import hashlib,json,sys
build_path,out_path,job=sys.argv[1:]
audit=json.load(open(build_path))
def sha(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for chunk in iter(lambda:f.read(1024*1024),b""): h.update(chunk)
    return h.hexdigest()
out={
 "schema_version":1,
 "candidate_id":"nylc_C18_trueT267_freeGS",
 "slurm_job_id":job,
 "build_state":audit["state"],
 "grompp_returncode":0,
 "maxwarn":0,
 "status":"PASS_BUILD_AND_GROMPP",
 "tpr_sha256":sha(build_path.rsplit("/audit/",1)[0]+"/build/em_preflight.tpr"),
 "scientific_status":"NOT_EVALUATED_L2_EM_AND_UNRESTRAINED_REBALANCE_REQUIRED",
}
open(out_path,"w").write(json.dumps(out,indent=2,sort_keys=True)+"\n")
PY
record PASS_BUILD_AND_GROMPP "audit=$root/audit/preflight.json"
failure=0
trap - EXIT
