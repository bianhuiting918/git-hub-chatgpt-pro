#!/bin/bash
#SBATCH -J nylc_true267
#SBATCH -p xahdnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --gres=dcu:1
#SBATCH --mem-per-cpu=2000M
#SBATCH -t 02:00:00
#SBATCH --array=0-1
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/true_thr267_%A_%a.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/true_thr267_%A_%a.err

set -euo pipefail
TASK=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
FLOW="$TASK/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
SOURCE=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/apo_gate_l4_three_carbonyl_20260715/cyclic_gate_nac_20260719/branches/nylc_gyaq
PY=/work/home/acshdt1dks/anaconda3/bin/python3.9
GMX=/public/software/apps/Gromacs-DCU2/2022.1/mpi/bin/gmx_mpi
mkdir -p "$TASK/logs" "$TASK/audit/true_thr267_residence"
source /work/home/acshdt1dks/opt/gromacs-fastest/env.sh

specs=(
  "C18 10297 10298"
  "C23 10303 10304"
)
read -r branch carbonyl_c carbonyl_o <<<"${specs[${SLURM_ARRAY_TASK_ID:?}]}"
stamp="$(date -Is)"
printf '%s\ttrue_thr267_residence\tSTART\tjob=%s;branch=%s;legacy=SUPERSEDED_WRONG_NUCLEOPHILE\n' "$stamp" "$SLURM_ARRAY_JOB_ID" "$branch" >>"$TASK/run_history.tsv"

"$PY" "$FLOW/scripts/audit_true_nylc_thr267_residence.py" \
  --branch-root "$SOURCE/$branch" \
  --carbonyl-c "$carbonyl_c" \
  --carbonyl-o "$carbonyl_o" \
  --gmx "$GMX" \
  --output-dir "$TASK/audit/true_thr267_residence/$branch" \
  --workers 4

stamp="$(date -Is)"
status="$("$PY" -c "import json;print(json.load(open('$TASK/audit/true_thr267_residence/$branch/true_thr267_residence_audit.json'))['scientific_status'])")"
printf '%s\ttrue_thr267_residence\tCOMPLETE\tjob=%s;branch=%s;status=%s\n' "$stamp" "$SLURM_ARRAY_JOB_ID" "$branch" "$status" >>"$TASK/run_history.tsv"
"$PY" - "$TASK/run_history.jsonl" "$stamp" "$SLURM_ARRAY_JOB_ID" "$branch" "$status" <<'PY'
import json,sys
path,stamp,job,branch,status=sys.argv[1:]
with open(path,"a") as f:
    f.write(json.dumps({"timestamp":stamp,"event":"true_thr267_residence","state":"COMPLETE","slurm_array_job_id":job,"branch":branch,"scientific_status":status},sort_keys=True)+"\n")
PY
