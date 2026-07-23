#!/bin/bash
#SBATCH --job-name=iccg_step1_pair
#SBATCH --output=iccg_step1_pair.%j.out
#SBATCH --error=iccg_step1_pair.%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=02:00:00
#SBATCH --partition=xahcnormal
set -euo pipefail
module load amber/2018-hpcx-gcc-7.3.1
: "${AMBERHOME:?Amber18 AMBERHOME must be loaded before running}"
PROJECT=/work/home/acshdt1dks/iccg_qmmm_step1_pair_20260723
cd "$PROJECT"
python audit_iccg_step1_pair.py preflight_report.json --submit-check
for state in LG1 LG2; do
  sander.MPI -O -i "${state}.in" -o "${state}.out" -p pair.prmtop -c "${state}.inpcrd" -r "${state}.rst" -inf "${state}.mdinfo"
  python audit_iccg_step1_pair.py "generated_${state}_pair.pdb" --geometry-json "${state}_postrun_geometry.json"
  python collect_iccg_step1_pair.py --log "${state}.out" --out "${state}.PASS.json" --process-exit 0 --geometry-report "${state}_postrun_geometry.json"
done
