#!/usr/bin/env bash
set -euo pipefail

: "${A1_ACTIVATED_ADDITION_CODE_ROOT:?set immutable code root}"
: "${A1_ACTIVATED_ADDITION_GITHUB_COMMIT:?set immutable GitHub commit}"
INDEX="${SLURM_ARRAY_TASK_ID:?run as array task 0 or 1}"
if (( INDEX < 0 || INDEX > 1 )); then exit 2; fi

TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
CODE_ROOT="$A1_ACTIVATED_ADDITION_CODE_ROOT"
COMMIT="$A1_ACTIVATED_ADDITION_GITHUB_COMMIT"
ARRAY_JOB="${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-manual}}"
OUT_ROOT="$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_step1_activated_mmframe_addition"
OUT="$OUT_ROOT/attempt_${ARRAY_JOB}_${INDEX}"
SCRATCH_BASE="${SLURM_TMPDIR:-/tmp}/nylc_a1_actadd_${ARRAY_JOB}_${INDEX}"
DRIVER="$CODE_ROOT/scripts/prepare_audit_nylc_a1_step1_activated_mmframe_addition.py"
PRMTOP="$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_unified_core_dftb3_preflight/attempt_62011285/prepared/system.prmtop"
EXPECTED_PRMTOP_SHA=a61d15bf0bf78675be93275d45f274e808ed6ae450fc1ca21a8e14aee8c12ca0
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
SUCCESS=0
FINAL_RC=0

append_history() {
    local state="$1" detail="$2" now
    now="$(date '+%Y-%m-%dT%H:%M:%S%z')"
    (
        flock -x 9
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$now" "nylc_a1_step1_activated_mmframe_addition" "run_nylc_a1_step1_activated_mmframe_addition.sh" "$state" "$COMMIT" "$detail" >> "$TASK_ROOT/run_history.tsv"
        printf '{"time":"%s","event":"nylc_a1_step1_activated_mmframe_addition","state":"%s","git_commit":"%s","detail":"%s","output":"%s"}\n' "$now" "$state" "$COMMIT" "$detail" "$OUT" >> "$TASK_ROOT/run_history.jsonl"
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
test -s "$PRMTOP"
test "$(sha256sum "$PRMTOP" | awk '{print $1}')" = "$EXPECTED_PRMTOP_SHA"
test ! -e "$OUT"
test ! -e "$SCRATCH_BASE"
append_history STARTED "job=$ARRAY_JOB;task=$INDEX;chemistry=NalphaH3plus_OG1minus;authority=unrestrained;windows=6"

module purge >/dev/null 2>&1 || true
module load amber/2018-hpcx-gcc-7.3.1
export AMBERHOME="$TASK_ROOT/../nylc_gyaq_pa66_l2_nac_qmmm_20260723/amber_runtime"
export OMP_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
test -s "$AMBERHOME/dat/slko/3ob-3-1/C-C.skf"

"$PY" "$DRIVER" --mode initialize --task-index "$INDEX" --root "$OUT" --github-commit "$COMMIT"
READY=$("$PY" -c 'import json,sys; print(json.load(open(sys.argv[1]))["ready_for_authority"])' "$OUT/SOURCE_MANIFEST.json")
if [[ "$READY" != True ]]; then
    "$PY" "$DRIVER" --mode finalize --root "$OUT"
    append_history TERMINAL "job=$ARRAY_JOB;task=$INDEX;gate=SOURCE_INTEGRITY;automatic_downstream=NONE"
    SUCCESS=1
    exit 0
fi

AUTH_SCRATCH="$SCRATCH_BASE/authority"
"$PY" "$DRIVER" --mode prepare-authority --root "$OUT" --scratch "$AUTH_SCRATCH"
SOURCE_RST=$("$PY" -c 'import json,sys; print(json.load(open(sys.argv[1]))["source"]["restart"])' "$OUT/SOURCE_MANIFEST.json")
set +e
(
    cd "$AUTH_SCRATCH"
    mpirun --bind-to none -np 8 sander.MPI -O -i authority.in -o authority.out -p "$PRMTOP" -c "$SOURCE_RST" -ref "$SOURCE_RST" -r authority.rst7 -inf authority.mdinfo
)
AUTH_RC=$?
set -e
"$PY" "$DRIVER" --mode audit-authority --root "$OUT" --scratch "$AUTH_SCRATCH"
ELIGIBLE=$("$PY" -c 'import json,sys; print(json.load(open(sys.argv[1]))["eligible_for_windows"])' "$OUT/AUTHORITY_RESULT.json")
if (( AUTH_RC != 0 )); then FINAL_RC=$AUTH_RC; fi

if [[ "$ELIGIBLE" == True ]]; then
    for WINDOW in 1 2 3 4 5 6; do
        WOUT="$OUT/window_$(printf '%02d' "$WINDOW")"
        WSCRATCH="$SCRATCH_BASE/window_$(printf '%02d' "$WINDOW")"
        "$PY" "$DRIVER" --mode prepare-window --root "$OUT" --output "$WOUT" --scratch "$WSCRATCH" --window-index "$WINDOW"
        INPUT_RST=$("$PY" -c 'import json,sys; print(json.load(open(sys.argv[1]))["input_restart"])' "$WOUT/WINDOW_MANIFEST.json")
        set +e
        (
            cd "$WSCRATCH"
            mpirun --bind-to none -np 8 sander.MPI -O -i stage.in -o stage.out -p "$PRMTOP" -c "$INPUT_RST" -ref "$INPUT_RST" -r stage.rst7 -inf stage.mdinfo
        )
        ENGINE_RC=$?
        set -e
        "$PY" "$DRIVER" --mode audit-window --root "$OUT" --output "$WOUT" --scratch "$WSCRATCH"
        ACCEPTED=$("$PY" -c 'import json,sys; print(json.load(open(sys.argv[1]))["accepted"])' "$WOUT/RESULT.json")
        if (( ENGINE_RC != 0 )); then FINAL_RC=$ENGINE_RC; fi
        if [[ "$ACCEPTED" != True ]]; then break; fi
    done
fi

"$PY" "$DRIVER" --mode finalize --root "$OUT"
STATUS=$("$PY" -c 'import json,sys; print(json.load(open(sys.argv[1]))["status"])' "$OUT/RESULT.json")
append_history TERMINAL "job=$ARRAY_JOB;task=$INDEX;authority_rc=$AUTH_RC;status=$STATUS;automatic_downstream=NONE"
if (( FINAL_RC != 0 )); then exit "$FINAL_RC"; fi
SUCCESS=1
printf '%s task=%s output=%s\n' "$STATUS" "$INDEX" "$OUT"
