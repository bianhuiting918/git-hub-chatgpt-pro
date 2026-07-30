#!/usr/bin/env bash
set -euo pipefail

: "${A1_RAW_NAC_FORWARD_CODE_ROOT:?set immutable code root}"
: "${A1_RAW_NAC_FORWARD_GITHUB_COMMIT:?set immutable GitHub commit}"

TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
CODE_ROOT="$A1_RAW_NAC_FORWARD_CODE_ROOT"
COMMIT="$A1_RAW_NAC_FORWARD_GITHUB_COMMIT"
INDEX="${SLURM_ARRAY_TASK_ID:?run as array task 0..17}"
ARRAY_JOB="${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-manual}}"
if (( INDEX < 0 || INDEX > 17 )); then exit 2; fi

OUT_ROOT="$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_step1_raw_nac_forward_calibration"
OUT="$OUT_ROOT/attempt_${ARRAY_JOB}_${INDEX}"
SCRATCH="${SLURM_TMPDIR:-/tmp}/nylc_a1_raw_nac_forward_${ARRAY_JOB}_${INDEX}"
DRIVER="$CODE_ROOT/scripts/prepare_audit_nylc_a1_step1_raw_nac_forward_calibration.py"
PRMTOP="$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_unified_core_dftb3_preflight/attempt_62011285/prepared/system.prmtop"
EXPECTED_PRMTOP_SHA=a61d15bf0bf78675be93275d45f274e808ed6ae450fc1ca21a8e14aee8c12ca0
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
SUCCESS=0

append_history() {
    local state="$1" detail="$2" now
    now="$(date '+%Y-%m-%dT%H:%M:%S%z')"
    (
        flock -x 9
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$now" "nylc_a1_step1_raw_nac_forward_calibration" "run_nylc_a1_step1_raw_nac_forward_calibration.sh" "$state" "$COMMIT" "$detail" >> "$TASK_ROOT/run_history.tsv"
        printf '{"time":"%s","event":"nylc_a1_step1_raw_nac_forward_calibration","state":"%s","git_commit":"%s","detail":"%s","output":"%s"}\n' "$now" "$state" "$COMMIT" "$detail" "$OUT" >> "$TASK_ROOT/run_history.jsonl"
    ) 9>> "$TASK_ROOT/.run_history.lock"
}

on_exit() {
    local rc=$?
    if (( SUCCESS == 0 )); then append_history FAILED "job=$ARRAY_JOB;task=$INDEX;rc=$rc"; fi
}
trap on_exit EXIT

test -s "$DRIVER"
test -s "$PRMTOP"
test "$(sha256sum "$PRMTOP" | awk '{print $1}')" = "$EXPECTED_PRMTOP_SHA"
test ! -e "$OUT"
test ! -e "$SCRATCH"
append_history STARTED "job=$ARRAY_JOB;task=$INDEX;matrix=2seeds_x_baseline_plus_2mechanisms_x_4force;raw_nac=true"

module purge >/dev/null 2>&1 || true
module load amber/2018-hpcx-gcc-7.3.1
export AMBERHOME="$TASK_ROOT/../nylc_gyaq_pa66_l2_nac_qmmm_20260723/amber_runtime"
export OMP_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
test -s "$AMBERHOME/dat/slko/3ob-3-1/C-C.skf"

"$PY" "$DRIVER" --mode initialize --task-index "$INDEX" --root "$OUT" --github-commit "$COMMIT"
"$PY" "$DRIVER" --mode prepare --root "$OUT" --scratch "$SCRATCH"
SOURCE_RST7=$("$PY" -c 'import json,sys; print(json.load(open(sys.argv[1]))["source"]["restart"])' "$OUT/MANIFEST.json")
test -s "$SOURCE_RST7"

set +e
(
    cd "$SCRATCH"
    mpirun --bind-to none -np 8 sander.MPI -O -i stage.in -o stage.out -p "$PRMTOP" -c "$SOURCE_RST7" -ref "$SOURCE_RST7" -r stage.rst7 -inf stage.mdinfo
)
ENGINE_RC=$?
set -e

"$PY" "$DRIVER" --mode audit --root "$OUT" --scratch "$SCRATCH"
(
    flock -x 9
    "$PY" "$DRIVER" --mode merge-if-ready --output-root "$OUT_ROOT" --array-job "$ARRAY_JOB"
) 9>> "$OUT_ROOT/.merge.lock"
STATUS=$("$PY" -c 'import json,sys; print(json.load(open(sys.argv[1]))["status"])' "$OUT/RESULT.json")
GATE=$("$PY" -c 'import json,sys; print(json.load(open(sys.argv[1]))["scientific_gate"])' "$OUT/RESULT.json")
append_history TERMINAL "job=$ARRAY_JOB;task=$INDEX;engine_rc=$ENGINE_RC;status=$STATUS;gate=$GATE;automatic_downstream=NONE"
if (( ENGINE_RC != 0 )); then exit "$ENGINE_RC"; fi
SUCCESS=1
printf '%s task=%s gate=%s output=%s\n' "$STATUS" "$INDEX" "$GATE" "$OUT"
