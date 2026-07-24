#!/bin/bash
#SBATCH -J nylc_ash306_probe
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem-per-cpu=1500M
#SBATCH -t 00:30:00
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/ash306_probe_%j.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/ash306_probe_%j.err

set -euo pipefail
TASK=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
FLOW="$TASK/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
SOURCE="$TASK/candidates/nylc_C18_trueT267_freeGS/qmmm/step1_dftb3_preflight_post61710861_job61712692/selected_free_l2_nac.gro"
job="${SLURM_JOB_ID:?}"
OUT="$TASK/candidates/nylc_C18_trueT267_freeGS/qmmm/ash306_chain_probe_job_${job}"
PY=/work/home/acshdt1dks/anaconda3/bin/python3.9
GMX=/public/software/apps/Gromacs-DCU2/2022.1/mpi/bin/gmx_mpi

record() {
  local state=$1 detail=$2 stamp
  stamp="$(date -Is)"
  printf '%s\tash306_chain_probe\t%s\tcandidate=nylc_C18_trueT267_freeGS;job=%s;%s\n' "$stamp" "$state" "$job" "$detail" >> "$TASK/run_history.tsv"
  "$PY" - "$stamp" "$state" "$detail" "$job" >> "$TASK/run_history.jsonl" <<'PY'
import json,sys
stamp,state,detail,job=sys.argv[1:]
print(json.dumps({"timestamp":stamp,"event":"ash306_chain_probe","state":state,"detail":detail,"candidate_id":"nylc_C18_trueT267_freeGS","slurm_job_id":job},sort_keys=True))
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

record START "source=434ps_GS_hypothesis;scope=isolated_chain_H_probe"
mkdir "$OUT"
"$PY" "$FLOW/scripts/prepare_nylc_ash306_chain_probe.py" --source-gro "$SOURCE" --export-pdb "$OUT/chainH_ASH306.pdb"
set +u
source /work/home/acshdt1dks/opt/gromacs-fastest/env.sh
set -u
(
  cd "$OUT"
  "$GMX" pdb2gmx -f chainH_ASH306.pdb -o chainH_ASH306.gro -p chainH.top -i posre_chainH.itp \
    -ff amber99sb-ildn -water tip3p -ignh > pdb2gmx.stdout 2> pdb2gmx.stderr
)
itp=$(find "$OUT" -maxdepth 1 -type f -name 'topol_Protein_chain_*.itp' -print -quit)
test -n "$itp"
set +e
"$PY" "$FLOW/scripts/prepare_nylc_ash306_chain_probe.py" --source-gro "$SOURCE" \
  --generated-gro "$OUT/chainH_ASH306.gro" --generated-itp "$itp" \
  --audit-output "$OUT/ash306_chain_probe_audit.json"
rc=$?
set -e
status=$("$PY" -c 'import json,sys;print(json.load(open(sys.argv[1]))["status"])' "$OUT/ash306_chain_probe_audit.json")
if [[ "$rc" -ne 0 || "$status" != PROBE_PASS_ASH306_CHAIN ]]; then
  record FAIL_TECHNICAL "status=$status;audit=$OUT/ash306_chain_probe_audit.json"
  exit 31
fi
cp "$OUT/ash306_chain_probe_audit.json" "$OUT/PASS.json"
record PASS_TECHNICAL "status=$status;audit=$OUT/PASS.json;interpretation=isolated_chain_probe_only"
