#!/usr/bin/env bash
set -euo pipefail

: "${A1_BIDIRECTIONAL_2CV_CODE_ROOT:?set immutable code root}"
: "${A1_BIDIRECTIONAL_2CV_GITHUB_COMMIT:?set immutable GitHub commit}"

TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
CODE_ROOT="$A1_BIDIRECTIONAL_2CV_CODE_ROOT"
COMMIT="$A1_BIDIRECTIONAL_2CV_GITHUB_COMMIT"
INDEX="${SLURM_ARRAY_TASK_ID:?run as array task 0..7}"
ARRAY_JOB="${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-manual}}"
if (( INDEX < 0 || INDEX > 7 )); then exit 2; fi

OUT_ROOT="$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_step1_bidirectional_2cv"
OUT="$OUT_ROOT/attempt_${ARRAY_JOB}_${INDEX}"
SCRATCH_ROOT="${SLURM_TMPDIR:-/tmp}/nylc_a1_step1_bidirectional_2cv_${ARRAY_JOB}_${INDEX}"
DRIVER="$CODE_ROOT/scripts/prepare_audit_nylc_a1_step1_bidirectional_2cv.py"
PRMTOP="$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_unified_core_dftb3_preflight/attempt_62011285/prepared/system.prmtop"
EXPECTED_PRMTOP_SHA=a61d15bf0bf78675be93275d45f274e808ed6ae450fc1ca21a8e14aee8c12ca0
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
SUCCESS=0

append_history() {
    local state="$1" detail="$2" now
    now="$(date '+%Y-%m-%dT%H:%M:%S%z')"
    (
        flock -x 9
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$now" "nylc_a1_step1_bidirectional_2cv" "run_nylc_a1_step1_bidirectional_2cv.sh" "$state" "$COMMIT" "$detail" >> "$TASK_ROOT/run_history.tsv"
        printf '{"time":"%s","event":"nylc_a1_step1_bidirectional_2cv","state":"%s","git_commit":"%s","detail":"%s","output":"%s"}\n' "$now" "$state" "$COMMIT" "$detail" "$OUT" >> "$TASK_ROOT/run_history.jsonl"
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
test ! -e "$SCRATCH_ROOT"
append_history STARTED "job=$ARRAY_JOB;task=$INDEX;chains=8;model=bidirectional_2cv"

module purge >/dev/null 2>&1 || true
module load amber/2018-hpcx-gcc-7.3.1
export AMBERHOME="$TASK_ROOT/../nylc_gyaq_pa66_l2_nac_qmmm_20260723/amber_runtime"
export OMP_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
test -s "$AMBERHOME/dat/slko/3ob-3-1/C-C.skf"

"$PY" "$DRIVER" --mode initialize --task-index "$INDEX" --root "$OUT" --github-commit "$COMMIT"
ACCEPTED_RST7="$OUT/accepted.rst7"
test -s "$ACCEPTED_RST7"

for WINDOW in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19; do
    PERSIST=$(printf '%s/window_%02d' "$OUT" "$WINDOW")
    SCRATCH=$(printf '%s/window_%02d' "$SCRATCH_ROOT" "$WINDOW")
    "$PY" "$DRIVER" --mode prepare-window --root "$OUT" --output "$PERSIST" --scratch "$SCRATCH" --input-rst7 "$ACCEPTED_RST7" --window-index "$WINDOW"
    set +e
    (
        cd "$SCRATCH"
        mpirun --bind-to none -np 8 sander.MPI -O -i stage.in -o stage.out -p "$PRMTOP" -c "$ACCEPTED_RST7" -ref "$ACCEPTED_RST7" -r stage.rst7 -inf stage.mdinfo
    )
    ENGINE_RC=$?
    set -e
    "$PY" "$DRIVER" --mode audit-window --root "$OUT" --output "$PERSIST" --scratch "$SCRATCH"
    TERMINAL=$("$PY" -c 'import json,sys; print(int(json.load(open(sys.argv[1], encoding="utf-8")).get("terminal", False)))' "$OUT/MANIFEST.json")
    if (( ENGINE_RC != 0 )) || [[ "$TERMINAL" = 1 ]]; then
        break
    fi
done

"$PY" "$DRIVER" --mode finalize --root "$OUT"
test -s "$OUT/SHA256.tsv"
(
    flock -x 9
    "$PY" "$DRIVER" --mode merge-if-ready --output-root "$OUT_ROOT" --array-job "$ARRAY_JOB"
) 9>> "$OUT_ROOT/.merge.lock"
STATUS=$("$PY" -c 'import json,sys; print(json.load(open(sys.argv[1], encoding="utf-8"))["status"])' "$OUT/RESULT.json")
GATE=$("$PY" -c 'import json,sys; print(json.load(open(sys.argv[1], encoding="utf-8"))["scientific_gate"])' "$OUT/RESULT.json")
append_history TERMINAL "job=$ARRAY_JOB;task=$INDEX;status=$STATUS;gate=$GATE;automatic_downstream=NONE"
SUCCESS=1
printf '%s task=%s gate=%s output=%s\n' "$STATUS" "$INDEX" "$GATE" "$OUT"
