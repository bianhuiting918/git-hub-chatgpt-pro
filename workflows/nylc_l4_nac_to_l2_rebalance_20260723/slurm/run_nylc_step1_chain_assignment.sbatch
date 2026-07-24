#!/bin/bash
#SBATCH -J nylc_s1_chain
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem-per-cpu=1500M
#SBATCH -t 00:20:00
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/step1_chain_%j.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/step1_chain_%j.err

set -euo pipefail
TASK=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
FLOW="$TASK/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
FREE="$TASK/candidates/nylc_C18_trueT267_freeGS/runs/rebalance_after_em_double/npt300free"
BASE="$TASK/candidates/nylc_C18_trueT267_freeGS/qmmm/step1_network_job_61713551/result/step1_network_frames.tsv"
job="${SLURM_JOB_ID:?}"
OUT="$TASK/candidates/nylc_C18_trueT267_freeGS/qmmm/step1_chain_assignment_job_${job}"
PY=/work/home/acshdt1dks/anaconda3/bin/python3.9
GMX=/public/software/apps/Gromacs-DCU2/2022.1/mpi/bin/gmx_mpi

record() {
  local state=$1 detail=$2 stamp
  stamp="$(date -Is)"
  printf '%s\tstep1_chain_assignment\t%s\tcandidate=nylc_C18_trueT267_freeGS;job=%s;%s\n' "$stamp" "$state" "$job" "$detail" >> "$TASK/run_history.tsv"
  "$PY" - "$stamp" "$state" "$detail" "$job" >> "$TASK/run_history.jsonl" <<'PY'
import json,sys
stamp,state,detail,job=sys.argv[1:]
print(json.dumps({"timestamp":stamp,"event":"step1_chain_assignment","state":state,"detail":detail,"candidate_id":"nylc_C18_trueT267_freeGS","slurm_job_id":job},sort_keys=True))
PY
}

on_error() {
  local rc=$?
  trap - ERR
  set +e
  record FAIL_TECHNICAL "exit_code=$rc;out=$OUT"
  exit "$rc"
}
trap on_error ERR

record START "source=step1_network_job_61713551;scope=all_Tyr146_Asn219_copy_combinations"
test -s "$FREE/run.xtc"
test -s "$FREE/run.tpr"
test -s "$BASE"
mkdir "$OUT"
set +u
source /work/home/acshdt1dks/opt/gromacs-fastest/env.sh
set -u
export GMX_MAXBACKUP=-1

"$GMX" distance -s "$FREE/run.tpr" -f "$FREE/run.xtc" -xvg none -oall "$OUT/copy1_distances.xvg" \
  -select 'atomnr 2020 plus atomnr 10288 plus atomnr 2021 plus atomnr 10288 plus atomnr 3104 plus atomnr 10288 plus atomnr 3105 plus atomnr 10288 plus atomnr 3106 plus atomnr 10288' \
  > "$OUT/distance.stdout" 2> "$OUT/distance.stderr"

cat > "$OUT/copy1_angles.ndx" <<'NDX'
[ OxyanionAnglesCopy1 ]
2020 2021 10288
3104 3105 10288
3104 3106 10288
NDX
printf 'OxyanionAnglesCopy1\n' | "$GMX" angle -f "$FREE/run.xtc" \
  -n "$OUT/copy1_angles.ndx" -type angle -all -xvg none -ov "$OUT/copy1_angles.xvg" \
  > "$OUT/angle.stdout" 2> "$OUT/angle.stderr"

"$PY" "$FLOW/scripts/analyze_nylc_step1_chain_assignment.py" \
  --base-frames "$BASE" --copy1-distances "$OUT/copy1_distances.xvg" \
  --copy1-angles "$OUT/copy1_angles.xvg" --output-dir "$OUT/result" \
  > "$OUT/analyze.stdout" 2> "$OUT/analyze.stderr"

AUDIT="$OUT/result/step1_chain_assignment_audit.json"
status=$("$PY" -c 'import json,sys;print(json.load(open(sys.argv[1]))["status"])' "$AUDIT")
if [[ "$status" == PASS_CROSS_SUBUNIT_OXYANION_PREORGANIZATION_PRESENT ]]; then
  record PASS_AUDIT "status=$status;audit=$AUDIT;next=extract_and_validate_selected_GS"
elif [[ "$status" == NOT_EVALUATED_NO_CROSS_SUBUNIT_OXYANION_PREORGANIZATION ]]; then
  record NOT_EVALUATED "status=$status;audit=$AUDIT;next=bounded_oxyanion_recapture_not_TS"
else
  record FAIL_TECHNICAL "status=$status;audit=$AUDIT"
  exit 31
fi
