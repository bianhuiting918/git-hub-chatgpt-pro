#!/bin/bash
#SBATCH -J nylc_m1_em
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem-per-cpu=2500M
#SBATCH -t 02:00:00
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/nylc_m1_em_%j.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/nylc_m1_em_%j.err

set -euo pipefail
TASK=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
FLOW="$TASK/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
GMX=/work/home/acshdt1dks/opt/gmx-cp2k/gromacs-2023.1-cp2k-2024.1/bin/gmx_mpi_d
candidate=nylc_C18_trueT267_freeGS
build="$TASK/candidates/$candidate/microstates/M1_nalphaH2_Asp306H/build_job_61801089/build"
build_pass="$TASK/candidates/$candidate/microstates/M1_nalphaH2_Asp306H/build_job_61801089/PASS.json"
job="${SLURM_JOB_ID:?}"
out="$TASK/candidates/$candidate/microstates/M1_nalphaH2_Asp306H/runs/em_job_$job"
start="$build/system_M1.gro"
mkdir -p "$TASK/logs"
test ! -e "$out"
mkdir -p "$out"
set +u
source /work/home/acshdt1dks/opt/gmx-cp2k/env.sh
set -u
export GMX_MAXBACKUP=-1

record() {
  local state=$1 detail=$2 stamp
  stamp="$(date -Is)"
  ( flock -x 9; printf '%s\tnylc_m1_em\t%s\tcandidate=%s;job=%s;%s\n' "$stamp" "$state" "$candidate" "$job" "$detail" >&9 ) 9>>"$TASK/run_history.tsv"
  ( flock -x 8; "$PY" - "$stamp" "$state" "$detail" "$job" >&8 <<'PY'
import json,sys
stamp,state,detail,job=sys.argv[1:]
print(json.dumps({"timestamp":stamp,"event":"nylc_m1_em","state":state,"detail":detail,"candidate_id":"nylc_C18_trueT267_freeGS","microstate":"M1_NalphaH2_Asp306H","slurm_job_id":job},sort_keys=True))
PY
  ) 8>>"$TASK/run_history.jsonl"
}
fail() { local rc=$?; trap - ERR; set +e; record FAIL_TECHNICAL "exit_code=$rc;out=$out"; exit "$rc"; }
trap fail ERR

"$PY" - "$build_pass" <<'PY'
import json,sys
d=json.load(open(sys.argv[1]))
assert d["status"]=="PASS_TECHNICAL_M1_BUILD_GROMPP"
m=d["build_audit"]["microstate"]
assert m=={"Asp306":"ASH_HD2_neutral","Asp308":"ASP_minus","Thr267_Nalpha":"NH2_neutral_proxy","Thr267_Ogamma":"OH"}
assert d["build_audit"]["old_chain_charge_e"]==d["build_audit"]["new_chain_charge_e"]==-4.0
PY
test -s "$start"
record START "source_build_job=61801089;microstate=NalphaH2_Asp306H;protein_and_L2_heavy_restraints=yes;scientific_status=NOT_EVALUATED"

(
  cd "$build"
  "$GMX" grompp -f "$FLOW/mdp/em_cg_flexible_m1.mdp" -c "$start" -r "$start" -p "$build/topol.top" -o "$out/run.tpr" -po "$out/run.processed.mdp" -maxwarn 0 >"$out/grompp.stdout" 2>"$out/grompp.stderr"
)
mpirun -np 1 "$GMX" mdrun -s "$out/run.tpr" -deffnm "$out/run" -ntomp "${SLURM_CPUS_PER_TASK:-8}" -nb cpu -pme cpu -bonded cpu -pin on >"$out/mdrun.stdout" 2>"$out/mdrun.stderr"
test -s "$out/run.gro"
fmax="$(awk '/Maximum force/{value=$4} END{print value}' "$out/run.log")"
test -n "$fmax"
hard=0
if grep -Eqi 'FATAL ERROR|LINCS WARNING|cannot be settled|SETTLE.*(warning|error|failed)|(^|[^[:alnum:]_])NaN([^[:alnum:]_]|$)' "$out/run.log" "$out/mdrun.stderr"; then hard=1; fi
status=FAIL
if [[ $hard -eq 0 ]] && awk -v value="$fmax" 'BEGIN{exit !(value<=500.0)}'; then status=PASS; fi
"$PY" - "$out/$status.json" "$candidate" "$fmax" "$hard" "$job" <<'PY'
import hashlib,json,sys
path,candidate,fmax,hard,job=sys.argv[1:]
gro=path.rsplit("/",1)[0]+"/run.gro"
payload={
 "schema_version":1,
 "candidate_id":candidate,
 "microstate":"M1_NalphaH2_Asp306H",
 "slurm_job_id":job,
 "source_build_job_id":"61801089",
 "stage":"m1_em_cg_flexible_double",
 "status":"PASS_TECHNICAL_EM" if path.endswith("/PASS.json") else "NOT_EVALUATED_EM_FAIL",
 "integrator":"cg",
 "flexible_water":True,
 "precision":"double",
 "gromacs_version":"2023.1",
 "protein_and_l2_heavy_atom_restraints":True,
 "fmax_kj_mol_nm":float(fmax),
 "threshold_kj_mol_nm":500.0,
 "hard_warning_count":int(hard),
 "scientific_nac_evidence":False,
 "run_gro_sha256":hashlib.sha256(open(gro,"rb").read()).hexdigest(),
}
open(path,"w").write(json.dumps(payload,indent=2,sort_keys=True)+"\n")
PY
if [[ "$status" != PASS ]]; then record FAIL_GATE "fmax=$fmax;hard=$hard"; exit 31; fi
record PASS_TECHNICAL "fmax=$fmax;hard=$hard;output=$out/PASS.json;scientific_status=NOT_EVALUATED"
trap - ERR
