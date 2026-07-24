#!/bin/bash
#SBATCH -J nylc_t267_dftb3
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem-per-cpu=8000M
#SBATCH -t 02:00:00
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/true_thr267_dftb3_%j.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/true_thr267_dftb3_%j.err

# Numerical entry gate only; this is not a TS search or PMF.
set -eo pipefail
TASK=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
FLOW="$TASK/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
POSTJOB=61709209
POST="$TASK/candidates/nylc_C18_trueT267_freeGS/runs/rebalance_after_em_double/postprocess_job_$POSTJOB"
OUT="$TASK/candidates/nylc_C18_trueT267_freeGS/qmmm/step1_dftb3_preflight_post$POSTJOB"
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
GMX=/public/software/apps/Gromacs-DCU2/2022.1/mpi/bin/gmx_mpi
AMBER_RUNTIME=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/nylc_gyaq_pa66_l2_nac_qmmm_20260723/amber_runtime
job="${SLURM_JOB_ID:?}"
set +u
source /work/home/acshdt1dks/opt/gromacs-fastest/env.sh
set -u
export GMX_MAXBACKUP=-1

record() {
  local state=$1 detail=$2 stamp
  stamp="$(date -Is)"
  printf '%s\tstep1_dftb3_preflight\t%s\tcandidate=nylc_C18_trueT267_freeGS;job=%s;%s\n' "$stamp" "$state" "$job" "$detail" >> "$TASK/run_history.tsv"
  "$PY" - "$stamp" "$state" "$detail" "$job" >> "$TASK/run_history.jsonl" <<'PY'
import json,sys
stamp,state,detail,job=sys.argv[1:]
print(json.dumps({"timestamp":stamp,"event":"step1_dftb3_preflight","state":state,"detail":detail,"candidate_id":"nylc_C18_trueT267_freeGS","slurm_job_id":job},sort_keys=True))
PY
}

record START "dependency_postprocess=$POSTJOB;numerical_preflight_only"
test -s "$POST/postprocess_summary.json"
"$PY" "$FLOW/scripts/prepare_nylc_true_thr267_dftb3_preflight.py" \
  --task-root "$TASK" --post-root "$POST" --output "$OUT" --gmx "$GMX" \
  > "$TASK/logs/true_thr267_dftb3_prepare_$job.stdout" \
  2> "$TASK/logs/true_thr267_dftb3_prepare_$job.stderr"
if [[ -s "$OUT/NOT_EVALUATED.json" ]]; then
  record NOT_EVALUATED "reason=NO_FREE_L2_NAC;output=$OUT/NOT_EVALUATED.json"
  exit 0
fi

module purge >/dev/null 2>&1 || true
module load amber/2018-hpcx-gcc-7.3.1
export AMBERHOME="$AMBER_RUNTIME"
test -s "$AMBERHOME/dat/slko/3ob-3-1/C-C.skf"
test -s "$OUT/qmmm_preflight_audit.json"
sander -O -i "$OUT/01_qmmm_one_step.in" -o "$OUT/01_qmmm_one_step.out" \
  -p "$OUT/system.prmtop" -c "$OUT/selected_free_l2_nac.rst7" -r "$OUT/01_qmmm_one_step.rst7"
if grep -Eqi 'Convergence could not be achieved|SANDER BOMB|segmentation|forrtl|vlimit[[:space:]]+exceeded|NaN|FATAL' "$OUT/01_qmmm_one_step.out"; then
  record FAIL_TECHNICAL "stage=one_step;output=$OUT/01_qmmm_one_step.out"
  exit 30
fi
sander -O -i "$OUT/02_qmmm_20_step.in" -o "$OUT/02_qmmm_20_step.out" \
  -p "$OUT/system.prmtop" -c "$OUT/01_qmmm_one_step.rst7" -r "$OUT/02_qmmm_20_step.rst7"
"$PY" "$FLOW/scripts/audit_nylc_dftb3_smoke.py" \
  --one-step "$OUT/01_qmmm_one_step.out" --segment "$OUT/02_qmmm_20_step.out" \
  --audit "$OUT/qmmm_preflight_audit.json" --output "$OUT/PASS.json"
record PASS_TECHNICAL "output=$OUT/PASS.json;interpretation=not_TS_search_or_PMF"
