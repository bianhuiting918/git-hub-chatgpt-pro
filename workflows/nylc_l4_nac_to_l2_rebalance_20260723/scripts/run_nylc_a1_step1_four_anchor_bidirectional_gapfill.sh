#!/usr/bin/env bash
set -euo pipefail
: "${A1_FOUR_ANCHOR_CODE_ROOT:?set immutable code root}"
: "${A1_FOUR_ANCHOR_GITHUB_COMMIT:?set immutable GitHub commit}"

TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
CODE_ROOT=$A1_FOUR_ANCHOR_CODE_ROOT
COMMIT=$A1_FOUR_ANCHOR_GITHUB_COMMIT
INDEX=$SLURM_ARRAY_TASK_ID
ARRAY_JOB=$SLURM_ARRAY_JOB_ID
OUT_ROOT=$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_step1_four_anchor_bidirectional_gapfill
OUT=$OUT_ROOT/attempt_$ARRAY_JOB'_'$INDEX
SCRATCH=${SLURM_TMPDIR:-/tmp}/nylc_a1_four_anchor_$ARRAY_JOB'_'$INDEX
DRIVER=$CODE_ROOT/scripts/prepare_audit_nylc_a1_step1_four_anchor_bidirectional_gapfill.py
PRMTOP=$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_unified_core_dftb3_preflight/attempt_62011285/prepared/system.prmtop
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
STOP=TECHNICAL_FAILURE
FINALIZED=0

finish_failure() {
    rc=$?
    if [[ $FINALIZED == 0 && -s $OUT/MANIFEST.json ]]; then
        "$PY" "$DRIVER" --mode finalize --root "$OUT" --stop-reason "$STOP" --technical-failure || true
    fi
    exit "$rc"
}
trap finish_failure ERR

test "$INDEX" -ge 0
test "$INDEX" -le 15
test -s "$DRIVER"
test -s "$PRMTOP"
test ! -e "$OUT"
mkdir -p "$OUT_ROOT" "$SCRATCH"
START=$SCRATCH/start.rst7

"$PY" "$DRIVER" --mode initialize --task-index "$INDEX" --root "$OUT"     --start-rst7 "$START" --github-commit "$COMMIT"

module purge >/dev/null 2>&1 || true
module load amber/2018-hpcx-gcc-7.3.1
export AMBERHOME=$TASK_ROOT/../nylc_gyaq_pa66_l2_nac_qmmm_20260723/amber_runtime
export OMP_NUM_THREADS=1
test -s "$AMBERHOME/dat/slko/3ob-3-1/C-C.skf"

PLAN=$("$PY" "$DRIVER" --mode describe-task --task-index "$INDEX")
COUNT=$("$PY" - "$PLAN" <<'PY'
import json,sys
print(len(json.loads(sys.argv[1])["stage_specs"]))
PY
)
CURRENT=$START
STOP=COMPLETED_ALL_STAGES

for ((STAGE_INDEX=0; STAGE_INDEX<COUNT; STAGE_INDEX++)); do
    STAGE=$(printf 'stage_%02d' "$STAGE_INDEX")
    PERSIST=$OUT/$STAGE
    WORK=$SCRATCH/$STAGE
    "$PY" "$DRIVER" --mode prepare-stage --root "$OUT" --output "$PERSIST"         --scratch "$WORK" --input-rst7 "$CURRENT" --stage-index "$STAGE_INDEX"
    (
        cd "$WORK"
        mpirun --bind-to none -np 8 sander.MPI -O             -i stage.in -o stage.out -p "$PRMTOP" -c "$CURRENT" -ref "$CURRENT"             -r stage.rst7 -inf stage.mdinfo
    )
    "$PY" "$DRIVER" --mode audit-stage --root "$OUT" --output "$PERSIST" --scratch "$WORK"
    CURRENT=$PERSIST/stage.rst7
    test -s "$CURRENT"
    GUARD=$("$PY" - "$PERSIST/RESULT.json" <<'PY'
import json,sys
print(int(json.load(open(sys.argv[1]))["guard_stop"]))
PY
)
    if [[ "$GUARD" == 1 ]]; then
        STOP=SCIENTIFIC_GUARD_STOP
        break
    fi
done

"$PY" "$DRIVER" --mode finalize --root "$OUT" --stop-reason "$STOP"
test -s "$OUT/SHA256.tsv"
FINALIZED=1
(
    flock -x 9
    "$PY" "$DRIVER" --mode merge-if-ready --output-root "$OUT_ROOT" --array-job "$ARRAY_JOB"
) 9>"$OUT_ROOT/.merge.lock"

HISTORY=$TASK_ROOT/run_history.tsv
(
    flock -x 9
    printf '%s\t%s\t%s\t%s\t%s\n' "$(date -Iseconds)" "$ARRAY_JOB" "$INDEX" "$COMMIT" "$STOP" >> "$HISTORY"
) 9>"$TASK_ROOT/.run_history.lock"
