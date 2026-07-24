#!/bin/bash
#SBATCH -J nylc_s1_oxyresp
#SBATCH -p xahdnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --gres=dcu:1
#SBATCH --mem-per-cpu=2000M
#SBATCH -t 01:00:00
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/step1_oxyresp_%j.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/step1_oxyresp_%j.err

set -euo pipefail
TASK=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
FLOW="$TASK/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
ROOT="$TASK/candidates/nylc_C18_trueT267_freeGS"
BUILD="$ROOT/build"
# Immutable source trajectory: runs/rebalance_after_em_double/npt300free/run.xtc
FREE="$ROOT/runs/rebalance_after_em_double/npt300free"
job="${SLURM_JOB_ID:?}"
OUT="$ROOT/qmmm/step1_oxyanion_recapture_response_job_${job}"
PY=/work/home/acshdt1dks/anaconda3/bin/python3.9
GMX=/public/software/apps/Gromacs-DCU2/2022.1/mpi/bin/gmx_mpi

record() {
  local state=$1 detail=$2 stamp
  stamp="$(date -Is)"
  printf '%s\tstep1_oxyanion_recapture_response\t%s\tcandidate=nylc_C18_trueT267_freeGS;job=%s;%s\n' "$stamp" "$state" "$job" "$detail" >> "$TASK/run_history.tsv"
  "$PY" - "$stamp" "$state" "$detail" "$job" >> "$TASK/run_history.jsonl" <<'PY'
import json,sys
stamp,state,detail,job=sys.argv[1:]
print(json.dumps({"timestamp":stamp,"event":"step1_oxyanion_recapture_response","state":state,"detail":detail,"candidate_id":"nylc_C18_trueT267_freeGS","slurm_job_id":job},sort_keys=True))
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

record START "source=fully_unrestrained_210ps;purpose=bounded_response_not_scientific_sampling"
mkdir "$OUT"
set +u
source /work/home/acshdt1dks/opt/gromacs-fastest/env.sh
set -u
export GMX_MAXBACKUP=-1

printf 'System\n' | "$GMX" trjconv -s "$FREE/run.tpr" -f "$FREE/run.xtc" -dump 210 -pbc mol \
  -o "$OUT/source.tmp.gro" > "$OUT/trjconv.stdout" 2> "$OUT/trjconv.stderr"
"$PY" "$FLOW/scripts/audit_nylc_step1_oxyanion_recapture.py" \
  --start-gro "$OUT/source.tmp.gro" --source-only --output "$OUT/source_validation.json" \
  > "$OUT/source_validation.stdout" 2> "$OUT/source_validation.stderr"
source_status=$("$PY" -c 'import json,sys;print(json.load(open(sys.argv[1]))["status"])' "$OUT/source_validation.json")
[[ "$source_status" == PASS_SOURCE_STRICT_NAC ]]
mv "$OUT/source.tmp.gro" "$OUT/source.gro"

printf 'q\\n' | "$GMX" make_ndx -f "$OUT/source.gro" -o "$OUT/recapture.ndx" > "$OUT/make_ndx.stdout" 2> "$OUT/make_ndx.stderr"\ncat >> "$OUT/recapture.ndx" <<'NDX'\n[ Thr267_OG1 ]
8961
[ L2_C12 ]
10287
[ L2_O2 ]
10288
[ Tyr146_N ]
7156
[ Tyr146_H ]
7157
[ Asn219_ND2 ]
8240
[ Asn219_HD21 ]
8241
NDX

(
  cd "$BUILD"
  "$GMX" grompp -f "$FLOW/mdp/nvt100_step1_oxyanion_recapture_response.mdp" \
    -c "$OUT/source.gro" -p "$BUILD/topol.top" -n "$OUT/recapture.ndx" \
    -o "$OUT/run.tpr" -po "$OUT/run.processed.mdp" -maxwarn 0 \
    > "$OUT/grompp.stdout" 2> "$OUT/grompp.stderr"
)
mpirun -np 1 "$GMX" mdrun -s "$OUT/run.tpr" -deffnm "$OUT/run" \
  -ntomp "${SLURM_CPUS_PER_TASK:-8}" -nb gpu -pme gpu -bonded gpu -pin on \
  > "$OUT/mdrun.stdout" 2> "$OUT/mdrun.stderr"
test -s "$OUT/run.gro"

"$PY" "$FLOW/scripts/audit_nylc_step1_oxyanion_recapture.py" \
  --start-gro "$OUT/source.gro" --end-gro "$OUT/run.gro" --log "$OUT/run.log" \
  --output "$OUT/recapture_response_audit.json" \
  > "$OUT/audit.stdout" 2> "$OUT/audit.stderr"
status=$("$PY" -c 'import json,sys;print(json.load(open(sys.argv[1]))["status"])' "$OUT/recapture_response_audit.json")
if [[ "$status" == PASS_TECHNICAL_RESPONSE ]]; then
  record PASS_TECHNICAL "status=$status;audit=$OUT/recapture_response_audit.json;scientific=NOT_EVALUATED_RESTRAINED"
else
  record NOT_EVALUATED "status=$status;audit=$OUT/recapture_response_audit.json;next=inspect_response_without_TS"
fi
