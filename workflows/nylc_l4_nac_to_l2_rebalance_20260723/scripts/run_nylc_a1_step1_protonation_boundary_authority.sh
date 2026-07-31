#!/usr/bin/env bash
set -euo pipefail

: "${A1_PB_CODE_ROOT:?set immutable code root}"
: "${A1_PB_GITHUB_COMMIT:?set immutable GitHub commit}"
INDEX="${SLURM_ARRAY_TASK_ID:?run as array task 0..3}"
if (( INDEX < 0 || INDEX > 3 )); then exit 2; fi

TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
CODE_ROOT="$A1_PB_CODE_ROOT"
COMMIT="$A1_PB_GITHUB_COMMIT"
ARRAY_JOB="${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-manual}}"
OUT_ROOT="$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_step1_protonation_boundary_authority"
OUT="$OUT_ROOT/attempt_${ARRAY_JOB}_${INDEX}"
SCRATCH_BASE="${SLURM_TMPDIR:-/tmp}/nylc_a1pb_${ARRAY_JOB}_${INDEX}"
DRIVER="$CODE_ROOT/scripts/prepare_audit_nylc_a1_step1_protonation_boundary_authority.py"
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
SUCCESS=0
FINAL_RC=0

append_history() {
    local state="$1" detail="$2" now
    now="$(date '+%Y-%m-%dT%H:%M:%S%z')"
    (
        flock -x 9
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$now" "nylc_a1_step1_protonation_boundary_authority" "run_nylc_a1_step1_protonation_boundary_authority.sh" "$state" "$COMMIT" "$detail" >> "$TASK_ROOT/run_history.tsv"
        printf '{"time":"%s","event":"nylc_a1_step1_protonation_boundary_authority","state":"%s","git_commit":"%s","detail":"%s","output":"%s"}\n' "$now" "$state" "$COMMIT" "$detail" "$OUT" >> "$TASK_ROOT/run_history.jsonl"
    ) 9>> "$TASK_ROOT/.run_history.lock"
}

on_exit() {
    local rc=$?
    if (( SUCCESS == 0 )); then
        append_history FAILED "job=$ARRAY_JOB;task=$INDEX;rc=$rc"
    fi
}
trap on_exit EXIT

test -s "$DRIVER"
test ! -e "$OUT"
test ! -e "$SCRATCH_BASE"
append_history STARTED "job=$ARRAY_JOB;task=$INDEX;blocks=4;unrestrained=true"

module purge >/dev/null 2>&1 || true
module load amber/2018-hpcx-gcc-7.3.1
export AMBERHOME="$TASK_ROOT/../nylc_gyaq_pa66_l2_nac_qmmm_20260723/amber_runtime"
export OMP_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
test -s "$AMBERHOME/dat/slko/3ob-3-1/C-C.skf"

"$PY" "$DRIVER" --mode initialize --task-index "$INDEX" --root "$OUT" --github-commit "$COMMIT"
READY=$("$PY" -c 'import json,sys; print(json.load(open(sys.argv[1]))["ready"])' "$OUT/SOURCE_MANIFEST.json")
if [[ "$READY" == True ]]; then
    for BLOCK in 1 2 3 4; do
        BOUT="$OUT/block_$(printf '%02d' "$BLOCK")"
        BSCRATCH="$SCRATCH_BASE/block_$(printf '%02d' "$BLOCK")"
        "$PY" "$DRIVER" --mode prepare-block --root "$OUT" --output "$BOUT" --scratch "$BSCRATCH" --block-index "$BLOCK"
        INPUT_RST=$("$PY" -c 'import json,sys; print(json.load(open(sys.argv[1]))["input_restart"])' "$BOUT/BLOCK_MANIFEST.json")
        PRMTOP=$("$PY" -c 'import json,sys; print(json.load(open(sys.argv[1]))["prmtop"])' "$OUT/SOURCE_MANIFEST.json")
        set +e
        (
            cd "$BSCRATCH"
            mpirun --bind-to none -np 8 sander.MPI -O -i block.in -o block.out -p "$PRMTOP" -c "$INPUT_RST" -r block.rst7 -inf block.mdinfo
        )
        ENGINE_RC=$?
        set -e
        "$PY" "$DRIVER" --mode audit-block --root "$OUT" --output "$BOUT" --scratch "$BSCRATCH"
        ACCEPTED=$("$PY" -c 'import json,sys; print(json.load(open(sys.argv[1]))["accepted"])' "$BOUT/RESULT.json")
        if (( ENGINE_RC != 0 )); then FINAL_RC=$ENGINE_RC; fi
        if [[ "$ACCEPTED" != True ]]; then break; fi
    done
fi

"$PY" "$DRIVER" --mode finalize --root "$OUT"
STATUS=$("$PY" -c 'import json,sys; print(json.load(open(sys.argv[1]))["status"])' "$OUT/RESULT.json")
append_history TERMINAL "job=$ARRAY_JOB;task=$INDEX;status=$STATUS;automatic_downstream=NONE"
if (( FINAL_RC != 0 )); then exit "$FINAL_RC"; fi
SUCCESS=1
printf '%s task=%s output=%s\n' "$STATUS" "$INDEX" "$OUT"
