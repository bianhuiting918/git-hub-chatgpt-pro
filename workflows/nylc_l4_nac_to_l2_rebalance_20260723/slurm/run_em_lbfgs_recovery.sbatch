#!/bin/bash
#SBATCH -J l4tol2_emlbfgs
#SBATCH -p xahdnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --gres=dcu:1
#SBATCH --mem-per-cpu=3500M
#SBATCH -t 02:00:00
#SBATCH --array=0-4%5
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/emlbfgs_%A_%a.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/emcg_%A_%a.err

set -euo pipefail
TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
FLOW="$TASK_ROOT/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
PYTHON=/work/home/acshdt1dks/anaconda3/bin/python3.9
GMX=/public/software/apps/Gromacs-DCU2/2022.1/mpi/bin/gmx_mpi
source /work/home/acshdt1dks/opt/gromacs-fastest/env.sh
export GMX_MAXBACKUP=-1
mapfile -t candidates < <("$PYTHON" - "$FLOW/manifests/candidates.json" <<'PY'
import json,sys
for row in json.load(open(sys.argv[1]))["candidates"]:
    print(row["id"])
PY
)
candidate="${candidates[${SLURM_ARRAY_TASK_ID}]}"
root="$TASK_ROOT/candidates/$candidate"
build="$root/build"
out="$root/runs/em_lbfgs"
mkdir -p "$out"
cd "$build"
"$GMX" grompp -f "$FLOW/mdp/em_lbfgs.mdp" -c "$build/rebuilt.gro" -r "$build/rebuilt.gro" \
  -p "$build/topol.top" -o "$out/run.tpr" -po "$out/run.processed.mdp" -maxwarn 0 \
  >"$out/grompp.stdout" 2>"$out/grompp.stderr"
mpirun -np 1 "$GMX" mdrun -deffnm "$out/run" -ntomp 8 -nb gpu -pme cpu -bonded cpu -pin on \
  >"$out/mdrun.stdout" 2>"$out/mdrun.stderr"
fmax="$(awk '/Maximum force/{value=$4} END{print value}' "$out/run.log")"
test -n "$fmax"
hard=0
if grep -Eqi 'FATAL ERROR|LINCS WARNING|cannot be settled|SETTLE.*(warning|error|failed)|(^|[^[:alnum:]_])NAN([^[:alnum:]_]|$)' "$out/run.log" "$out/mdrun.stderr"; then hard=1; fi
status=FAIL
if [[ $hard -eq 0 ]] && awk -v v="$fmax" 'BEGIN{exit !(v<=500.0)}'; then status=PASS; fi
"$PYTHON" - "$out/$status.json" "$candidate" "$fmax" "$hard" <<'PY'
import json,sys
path,candidate,fmax,hard=sys.argv[1:]
json.dump({
 "candidate_id":candidate,
 "stage":"em_lbfgs",
 "status":"PASS_TECHNICAL_EM" if path.endswith("/PASS.json") else "NOT_EVALUATED_EM_FAIL",
 "integrator":"l-bfgs",
 "start":"rebuilt.gro",
 "fmax_kj_mol_nm":float(fmax),
 "threshold_kj_mol_nm":500.0,
 "hard_warning_count":int(hard),
 "scientific_nac_evidence":False,
},open(path,"w"),indent=2)
open(path,"a").write("\n")
PY
if [[ $status != PASS ]]; then exit 31; fi
