#!/bin/bash
#SBATCH -J nylc_a1_acdiag
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem-per-cpu=1000M
#SBATCH -t 00:10:00
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-%x-%j.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-%x-%j.err

set -uo pipefail
TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
INPUT="$TASK_ROOT/a1_activated_nac_20260726/parameterization/job_61887449/model/NTA1_CAP.input.mol2"
OUT="$TASK_ROOT/a1_activated_nac_20260726/diagnostics/antechamber_job_${SLURM_JOB_ID:?}"
mkdir -p "$OUT"
module purge
module load amber/2018-hpcx-gcc-7.3.1
antechamber -h >"$OUT/help.stdout" 2>"$OUT/help.stderr"
help_rc=$?
antechamber -i "$INPUT" -fi mol2 -o "$OUT/read_charge.mol2" -fo mol2 \
  -c rc -nc 0 -m 1 -at gaff2 -rn NTA1 -s 2 \
  >"$OUT/read_charge.stdout" 2>"$OUT/read_charge.stderr"
parse_rc=$?
printf 'help_rc\t%s\nparse_rc\t%s\n' "$help_rc" "$parse_rc" >"$OUT/result.tsv"
now="$(date '+%Y-%m-%dT%H:%M:%S%z')"
state=PASS_DIAGNOSTIC
if [[ $parse_rc -ne 0 ]]; then state=FAIL_DIAGNOSTIC_INPUT_PARSE; fi
printf '%s\t%s\t%s\t%s\t%s\n' "$now" nylc_a1_antechamber_diagnostic "antechamber -c rc current_input" "$state" "job=${SLURM_JOB_ID};parse_rc=$parse_rc;output=$OUT" >>"$TASK_ROOT/run_history.tsv"
"$PY" - "$now" "$state" "$parse_rc" "$OUT" "${SLURM_JOB_ID}" >>"$TASK_ROOT/run_history.jsonl" <<'PY'
import json
import sys
now, state, rc, output, job = sys.argv[1:]
print(json.dumps({
    "time": now,
    "event": "nylc_a1_antechamber_diagnostic",
    "command": "antechamber -c rc current_input",
    "state": state,
    "parse_rc": int(rc),
    "output": output,
    "slurm_job_id": job,
}, sort_keys=True))
PY
exit "$parse_rc"
