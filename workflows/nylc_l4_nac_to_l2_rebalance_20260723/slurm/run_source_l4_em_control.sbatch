#!/bin/bash
#SBATCH -J l4_source_em
#SBATCH -p xahdnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --gres=dcu:1
#SBATCH --mem-per-cpu=3500M
#SBATCH -t 00:30:00
#SBATCH --array=3-4%2
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/source_em_%A_%a.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/source_em_%A_%a.err

set -euo pipefail
TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
FLOW="$TASK_ROOT/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
PYTHON=/work/home/acshdt1dks/anaconda3/bin/python3.9
GMX=/public/software/apps/Gromacs-DCU2/2022.1/mpi/bin/gmx_mpi
source /work/home/acshdt1dks/opt/gromacs-fastest/env.sh
export GMX_MAXBACKUP=-1
mapfile -t rows < <("$PYTHON" - "$FLOW/manifests/candidates.json" <<'PY'
import json,sys
for row in json.load(open(sys.argv[1]))["candidates"]:
    print(row["id"] + "\t" + row["branch_dir"])
PY
)
IFS=$'\t' read -r candidate branch <<< "${rows[${SLURM_ARRAY_TASK_ID}]}"
out="$TASK_ROOT/candidates/$candidate/controls/source_l4_em"
mkdir -p "$out"
cd "$branch"
"$GMX" grompp -f "$FLOW/mdp/em.mdp" \
  -c "$TASK_ROOT/inputs/$candidate/source.gro" \
  -r "$TASK_ROOT/inputs/$candidate/source.gro" \
  -p "$branch/topol.top" -o "$out/run.tpr" -po "$out/run.processed.mdp" -maxwarn 0 \
  >"$out/grompp.stdout" 2>"$out/grompp.stderr"
mpirun -np 1 "$GMX" mdrun -deffnm "$out/run" -ntomp 8 -nb gpu -pme cpu -bonded cpu -pin on \
  >"$out/mdrun.stdout" 2>"$out/mdrun.stderr"
fmax="$(awk '/Maximum force/{value=$4} END{print value}' "$out/run.log")"
"$PYTHON" - "$out/result.json" "$candidate" "$fmax" <<'PY'
import json,sys
path,candidate,fmax=sys.argv[1:]
value=float(fmax)
json.dump({
 "candidate_id":candidate,
 "control":"unchanged source L4 at selected XTC frame",
 "fmax_kj_mol_nm":value,
 "em_threshold":500.0,
 "status":"PASS" if value <= 500 else "FAIL",
},open(path,"w"),indent=2)
open(path,"a").write("\n")
PY
