#!/bin/bash
#SBATCH -J nylc_t267_emcgf
#SBATCH -p xahdnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --gres=dcu:1
#SBATCH --mem-per-cpu=2500M
#SBATCH -t 02:00:00
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/true_thr267_emcgffgs_%j.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/true_thr267_emcgffgs_%j.err

set -euo pipefail
TASK=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
FLOW="$TASK/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
PY=/work/home/acshdt1dks/anaconda3/bin/python3.9
GMX=/public/software/apps/Gromacs-DCU2/2022.1/mpi/bin/gmx_mpi
candidate=nylc_C18_trueT267_freeGS
root="$TASK/candidates/$candidate"
build="$root/build"
out="$root/runs/em_cg_flexible"
start="$build/rebuilt.gro"
job="${SLURM_JOB_ID:?}"
mkdir -p "$out"
source /work/home/acshdt1dks/opt/gromacs-fastest/env.sh
export GMX_MAXBACKUP=-1
for f in run.tpr run.gro PASS.json FAIL.json; do test ! -e "$out/$f"; done

record() {
  local state=$1 detail=$2 stamp
  stamp="$(date -Is)"
  printf '%s\ttrue_thr267_l2_em_cg_flexible\t%s\tcandidate=%s;job=%s;%s\n' "$stamp" "$state" "$candidate" "$job" "$detail" >>"$TASK/run_history.tsv"
  "$PY" - "$stamp" "$state" "$detail" "$job" >>"$TASK/run_history.jsonl" <<'PY'
import json,sys
stamp,state,detail,job=sys.argv[1:]
print(json.dumps({"timestamp":stamp,"event":"true_thr267_l2_em_cg_flexible","state":state,"detail":detail,"candidate_id":"nylc_C18_trueT267_freeGS","slurm_job_id":job},sort_keys=True))
PY
}
record START "reason=steep_EM_61707615_water_pair_collapse;start=immutable_rebuilt.gro;flexible_water=yes"
(
  cd "$build"
  "$GMX" grompp -f "$FLOW/mdp/em_cg_flexible.mdp" -c "$start" -r "$start" -p "$build/topol.top" \
    -o "$out/run.tpr" -po "$out/run.processed.mdp" -maxwarn 0 \
    >"$out/grompp.stdout" 2>"$out/grompp.stderr"
)
mpirun -np 1 "$GMX" mdrun -s "$out/run.tpr" -deffnm "$out/run" -ntomp "${SLURM_CPUS_PER_TASK:-8}" \
  -nb gpu -pme cpu -bonded cpu -pin on >"$out/mdrun.stdout" 2>"$out/mdrun.stderr"
fmax="$(awk '/Maximum force/{value=$4} END{print value}' "$out/run.log")"
test -n "$fmax"
hard=0
if grep -Eqi 'FATAL ERROR|LINCS WARNING|cannot be settled|SETTLE.*(warning|error|failed)|(^|[^[:alnum:]_])NaN([^[:alnum:]_]|$)' "$out/run.log" "$out/mdrun.stderr"; then hard=1; fi
status=FAIL
if [[ $hard -eq 0 ]] && awk -v value="$fmax" 'BEGIN{exit !(value<=500.0)}'; then status=PASS; fi
"$PY" - "$out/$status.json" "$candidate" "$fmax" "$hard" "$job" <<'PY'
import json,sys
path,candidate,fmax,hard,job=sys.argv[1:]
open(path,"w").write(json.dumps({
 "schema_version":1,
 "candidate_id":candidate,
 "slurm_job_id":job,
 "stage":"em_cg_flexible",
 "status":"PASS_TECHNICAL_EM" if path.endswith("/PASS.json") else "NOT_EVALUATED_EM_FAIL",
 "integrator":"cg",
 "flexible_water":True,
 "start":"immutable build/rebuilt.gro",
 "fmax_kj_mol_nm":float(fmax),
 "threshold_kj_mol_nm":500.0,
 "hard_warning_count":int(hard),
 "scientific_nac_evidence":False,
},indent=2,sort_keys=True)+"\n")
PY
if [[ "$status" != PASS ]]; then record FAIL_GATE "fmax=$fmax;hard=$hard"; exit 31; fi
record PASS_TECHNICAL "fmax=$fmax;hard=$hard;output=$out/PASS.json"
