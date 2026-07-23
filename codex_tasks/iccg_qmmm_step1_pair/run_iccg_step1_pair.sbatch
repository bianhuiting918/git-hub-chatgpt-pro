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
PROJECT=/work/home/acshdt1dks/iccg_qmmm_step1_pair_20260723
cd "$PROJECT"
job_id="${SLURM_JOB_ID:-NO_SLURM_JOB_ID}"
append_history() {
  printf '%s\tSBATCH\tjob=%s\tstate=%s\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$job_id" "$1" "$2" >> RUN_HISTORY.tsv
}
if ! command -v sander.MPI >/dev/null; then
  append_history RUNTIME_ENV NOT_SUBMITTED_RUNTIME_ENV_SANDER_MPI_MISSING
  printf '%s\n' NOT_SUBMITTED_RUNTIME_ENV_SANDER_MPI_MISSING > SUBMISSION_STATUS.txt
  exit 127
fi
python audit_iccg_step1_pair.py preflight_report.json --submit-check
for state in LG1 LG2; do
  append_history "$state" STATE_BEGIN
  set +e
  sander.MPI -O -i "${state}.in" -o "${state}.out" -p pair.prmtop -c "${state}.inpcrd" -r "${state}.rst" -inf "${state}.mdinfo"
  rc=$?
  set -e
  if [[ "$rc" -ne 0 ]]; then
    append_history "$state" "TECHNICAL_FAIL_SANDER_MPI_EXIT_${rc}"
    printf '%s\n' "NOT_EVALUATED_SANDER_MPI_EXIT_${rc}" > SUBMISSION_STATUS.txt
    exit "$rc"
  fi
  python audit_iccg_step1_pair.py "generated_${state}_pair.pdb" --geometry-json "${state}_postrun_geometry.json"
  python collect_iccg_step1_pair.py --log "${state}.out" --out "${state}.PASS.json" --process-exit 0 --geometry-report "${state}_postrun_geometry.json"
  append_history "$state" STATE_END
done
