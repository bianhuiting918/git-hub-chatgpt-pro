#!/bin/bash
#SBATCH -J nylc_s1_network
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem-per-cpu=1500M
#SBATCH -t 00:30:00
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/step1_network_%j.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/step1_network_%j.err

set -euo pipefail
TASK=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
FLOW="$TASK/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
# Immutable inputs: npt300free/run.xtc and npt300free/run.tpr
FREE="$TASK/candidates/nylc_C18_trueT267_freeGS/runs/rebalance_after_em_double/npt300free"
POST="$TASK/candidates/nylc_C18_trueT267_freeGS/runs/rebalance_after_em_double/postprocess_job_61710861"
job="${SLURM_JOB_ID:?}"
OUT="$TASK/candidates/nylc_C18_trueT267_freeGS/qmmm/step1_network_job_${job}"
PY=/work/home/acshdt1dks/anaconda3/bin/python3.9
GMX=/public/software/apps/Gromacs-DCU2/2022.1/mpi/bin/gmx_mpi

record() {
  local state=$1 detail=$2 stamp
  stamp="$(date -Is)"
  printf '%s\tstep1_network_screen\t%s\tcandidate=nylc_C18_trueT267_freeGS;job=%s;%s\n' "$stamp" "$state" "$job" "$detail" >> "$TASK/run_history.tsv"
  "$PY" - "$stamp" "$state" "$detail" "$job" >> "$TASK/run_history.jsonl" <<'PY'
import json,sys
stamp,state,detail,job=sys.argv[1:]
print(json.dumps({"timestamp":stamp,"event":"step1_network_screen","state":state,"detail":detail,"candidate_id":"nylc_C18_trueT267_freeGS","slurm_job_id":job},sort_keys=True))
PY
}

record START "source=fully_unrestrained_NPT_1ns;postprocess=61710861"
test -s "$FREE/run.xtc"
test -s "$FREE/run.tpr"
test -s "$POST/nac_distance.xvg"
test -s "$POST/nac_angle.xvg"
test -s "$POST/potential_energy.xvg"
mkdir "$OUT"
set +u
source /work/home/acshdt1dks/opt/gromacs-fastest/env.sh
set -u
export GMX_MAXBACKUP=-1

"$GMX" distance -s "$FREE/run.tpr" -f "$FREE/run.xtc" -xvg none -oall "$OUT/network_distances.xvg" \
  -select 'atomnr 7156 plus atomnr 10288 plus atomnr 7157 plus atomnr 10288 plus atomnr 8240 plus atomnr 10288 plus atomnr 8241 plus atomnr 10288 plus atomnr 8242 plus atomnr 10288 plus atomnr 9572 plus atomnr 8949 plus atomnr 9573 plus atomnr 8949 plus atomnr 9572 plus atomnr 9591 plus atomnr 9572 plus atomnr 9592 plus atomnr 9573 plus atomnr 9591 plus atomnr 9573 plus atomnr 9592 plus atomnr 7768 plus atomnr 8961 plus atomnr 8239 plus atomnr 7768' \
  > "$OUT/distance.stdout" 2> "$OUT/distance.stderr"

cat > "$OUT/oxyanion_angles.ndx" <<'NDX'
[ OxyanionAngles ]
7156 7157 10288
8240 8241 10288
8240 8242 10288
NDX
printf 'OxyanionAngles\n' | "$GMX" angle -s "$FREE/run.tpr" -f "$FREE/run.xtc" \
  -n "$OUT/oxyanion_angles.ndx" -type angle -all -xvg none -ov "$OUT/network_angles.xvg" \
  > "$OUT/angle.stdout" 2> "$OUT/angle.stderr"

"$PY" "$FLOW/scripts/analyze_nylc_step1_network_nac.py" \
  --distance "$POST/nac_distance.xvg" --attack-angle "$POST/nac_angle.xvg" \
  --potential "$POST/potential_energy.xvg" --network-distances "$OUT/network_distances.xvg" \
  --network-angles "$OUT/network_angles.xvg" --output-dir "$OUT/result" \
  > "$OUT/analyze.stdout" 2> "$OUT/analyze.stderr"

AUDIT="$OUT/result/step1_network_audit.json"
status=$("$PY" -c 'import json,sys;print(json.load(open(sys.argv[1]))["status"])' "$AUDIT")
if [[ "$status" == NOT_EVALUATED_NO_FULL_CATALYTIC_PREORGANIZATION ]]; then
  record NOT_EVALUATED "status=$status;audit=$AUDIT;next=bounded_oxyanion_recapture_not_TS"
  exit 0
fi
if [[ "$status" != PASS_FULL_CATALYTIC_PREORGANIZATION_PRESENT ]]; then
  record FAIL_TECHNICAL "status=$status;audit=$AUDIT"
  exit 30
fi
selected_time=$("$PY" -c 'import json,sys;print(json.load(open(sys.argv[1]))["lowest_potential_full_preorganization_frame"]["time_ps"])' "$AUDIT")
printf 'System\n' | "$GMX" trjconv -s "$FREE/run.tpr" -f "$FREE/run.xtc" -dump "$selected_time" \
  -pbc mol -o "$OUT/result/source.tmp.gro" > "$OUT/trjconv.stdout" 2> "$OUT/trjconv.stderr"
"$PY" "$FLOW/scripts/analyze_nylc_step1_network_nac.py" \
  --validate-gro "$OUT/result/source.tmp.gro" --audit "$AUDIT" \
  > "$OUT/validate.stdout" 2> "$OUT/validate.stderr"
record PASS "status=$status;audit=$AUDIT;selected=$OUT/result/selected_step1_network_gs.gro"
