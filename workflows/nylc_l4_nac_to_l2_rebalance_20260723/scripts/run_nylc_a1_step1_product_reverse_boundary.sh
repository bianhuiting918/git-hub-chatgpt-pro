#!/usr/bin/env bash
set -euo pipefail

TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
CODE_ROOT="${A1_PRODUCT_REVERSE_CODE_ROOT:?set immutable code root}"
GITHUB_COMMIT="${A1_PRODUCT_REVERSE_GITHUB_COMMIT:?set immutable GitHub commit}"
INDEX="${SLURM_ARRAY_TASK_ID:?run as array task 0 or 1}"
ARRAY_JOB="${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-manual}}"
ATTEMPT="${ARRAY_JOB}_${INDEX}"
OUTPUT_ROOT="$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_step1_product_reverse_boundary"
OUT="$OUTPUT_ROOT/attempt_$ATTEMPT"
SCRATCH_ROOT="${SLURM_TMPDIR:-/tmp}/nylc_a1_product_reverse_$ATTEMPT"
DRIVER="$CODE_ROOT/scripts/prepare_audit_nylc_a1_step1_product_reverse_boundary.py"
PRMTOP="$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_unified_core_dftb3_preflight/attempt_62011285/prepared/system.prmtop"
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python

case "$INDEX" in
    0)
        SOURCE_RST7="$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_step1_acyl_release_md_continuation/attempt_62216380_0/release_md_endpoint.rst7"
        EXPECTED_SOURCE_SHA=5d8f76d2c90e3e8c707b640c55a93938f53e18dc30d6f93d54adda467da25f41
        ;;
    1)
        SOURCE_RST7="$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_step1_acyl_release_md_continuation/attempt_62216380_1/release_md_endpoint.rst7"
        EXPECTED_SOURCE_SHA=4cc60ad4d7be099b3f76040f1ee8c49b91deecebb20ed98570bc192b3d511132
        ;;
    *)
        printf 'invalid array index: %s\n' "$INDEX" >&2
        exit 2
        ;;
esac

EVENT=nylc_a1_step1_product_reverse_boundary
COMMAND=run_nylc_a1_step1_product_reverse_boundary.sh
STATE=STARTED
CURRENT=preflight
DETAIL="array_index=$INDEX;source=$SOURCE_RST7;source_sha=$EXPECTED_SOURCE_SHA;output=$OUT"
OWNED=0
FINALIZED=0

append_history() {
    local now
    now="$(date '+%Y-%m-%dT%H:%M:%S%z')"
    (
        flock -x 9
        "$PY" - "$TASK_ROOT/run_history.tsv" "$TASK_ROOT/run_history.jsonl" \
            "$now" "$EVENT" "$COMMAND" "$STATE" "$GITHUB_COMMIT" "$DETAIL" "$OUT" <<'PY'
import json
import os
import sys

tsv_path, jsonl_path, now, event, command, state, commit, detail, output = sys.argv[1:]
tsv_row = "\t".join((now, event, command, state, commit, detail)) + "\n"
jsonl_row = json.dumps(
    {
        "time": now,
        "event": event,
        "command": command,
        "state": state,
        "git_commit": commit,
        "detail": detail,
        "output": output,
    },
    sort_keys=True,
) + "\n"
with open(tsv_path, "a+", encoding="utf-8") as tsv, open(
    jsonl_path, "a+", encoding="utf-8"
) as jsonl:
    tsv.seek(0, os.SEEK_END)
    jsonl.seek(0, os.SEEK_END)
    tsv_offset, jsonl_offset = tsv.tell(), jsonl.tell()
    try:
        tsv.write(tsv_row)
        tsv.flush()
        os.fsync(tsv.fileno())
        jsonl.write(jsonl_row)
        jsonl.flush()
        os.fsync(jsonl.fileno())
    except BaseException:
        tsv.seek(tsv_offset)
        tsv.truncate()
        jsonl.seek(jsonl_offset)
        jsonl.truncate()
        raise
PY
    ) 9>>"$TASK_ROOT/run_history.tsv"
}

merge_if_ready_locked() {
    mkdir -p "$OUTPUT_ROOT"
    (
        flock -x 9
        "$PY" "$DRIVER" --mode merge-if-ready \
            --output-root "$OUTPUT_ROOT" --array-job "$ARRAY_JOB"
    ) 9>"$OUTPUT_ROOT/.merge.lock"
}

result_is_complete() {
    test -s "$OUT/SOURCE_MANIFEST.json"
    test -s "$OUT/RESULT.json"
    test -s "$OUT/SHA256.tsv"
    local sentinels=0
    if [[ -s "$OUT/PASS.json" ]]; then
        sentinels=$((sentinels + 1))
    fi
    if [[ -s "$OUT/NOT_EVALUATED.json" ]]; then
        sentinels=$((sentinels + 1))
    fi
    (( sentinels == 1 ))
}

finish() {
    local code=$?
    trap - EXIT
    if (( OWNED && ! FINALIZED )); then
        CURRENT=failure_finalize
        if "$PY" "$DRIVER" --mode finalize --root "$OUT" \
            --stop-reason "TECHNICAL_FAILURE" --technical-failure; then
            if result_is_complete; then
                FINALIZED=1
                merge_if_ready_locked || true
            fi
        fi
    fi
    STATE=NOT_EVALUATED
    DETAIL="array_index=$INDEX;stage=$CURRENT;exit_code=$code;status=NOT_EVALUATED_TECHNICAL_A1_PRODUCT_REVERSE_BOUNDARY_SCOUT;output=$OUT"
    append_history || true
    exit "$code"
}
trap finish EXIT

append_history

test -s "$DRIVER"
test -s "$PRMTOP"
test -s "$SOURCE_RST7"
ACTUAL_SOURCE_SHA="$(sha256sum "$SOURCE_RST7" | awk '{print $1}')"
test "$ACTUAL_SOURCE_SHA" = "$EXPECTED_SOURCE_SHA"
test ! -e "$OUT"

mkdir -p "$OUTPUT_ROOT" "$SCRATCH_ROOT"
START_RST7="$SCRATCH_ROOT/start.rst7"

CURRENT=initialize
"$PY" "$DRIVER" --mode initialize --seed-index "$INDEX" \
    --root "$OUT" --start-rst7 "$START_RST7" \
    --github-commit "$GITHUB_COMMIT"
OWNED=1

module purge >/dev/null 2>&1 || true
module load amber/2018-hpcx-gcc-7.3.1
export AMBERHOME="$TASK_ROOT/../nylc_gyaq_pa66_l2_nac_qmmm_20260723/amber_runtime"
export OMP_NUM_THREADS=1
test -s "$AMBERHOME/dat/slko/3ob-3-1/C-C.skf"

CURRENT_RST7="$START_RST7"
STOP_REASON=COMPLETED_ALL_WINDOWS

for WINDOW_INDEX in 0 1 2 3; do
    WINDOW_TAG="$(printf '%02d' "$WINDOW_INDEX")"
    WINDOW_OUT="$OUT/window_$WINDOW_TAG"
    WINDOW_SCRATCH="$SCRATCH_ROOT/window_$WINDOW_TAG"
    NEXT_RST7="$WINDOW_SCRATCH/stage.rst7"

    CURRENT="prepare_window_$WINDOW_TAG"
    "$PY" "$DRIVER" --mode prepare-window \
        --window-index "$WINDOW_INDEX" --root "$OUT" \
        --output "$WINDOW_OUT" --scratch "$WINDOW_SCRATCH" \
        --input-rst7 "$CURRENT_RST7"

    CURRENT="run_window_$WINDOW_TAG"
    (
        cd "$WINDOW_SCRATCH"
        mpirun --bind-to none -np "${SLURM_NTASKS:-8}" sander.MPI -O \
            -i stage.in -o stage.out -p "$PRMTOP" \
            -c "$CURRENT_RST7" -ref "$CURRENT_RST7" \
            -r stage.rst7 -inf stage.mdinfo
    )
    test -s "$NEXT_RST7"

    CURRENT="audit_window_$WINDOW_TAG"
    "$PY" "$DRIVER" --mode audit-window \
        --root "$OUT" --output "$WINDOW_OUT" --scratch "$WINDOW_SCRATCH"

    read -r WINDOW_TECHNICAL GUARD_STOP RESTRAINED_HINT < <(
        "$PY" - "$WINDOW_OUT/RESULT.json" <<'PY'
import json
import sys

result = json.load(open(sys.argv[1], encoding="utf-8"))
print(
    int(result["technical_pass"]),
    int(result["guard_stop"]),
    int(result["restrained_hint"]),
)
PY
    )
    test "$WINDOW_TECHNICAL" = 1
    CURRENT_RST7="$NEXT_RST7"

    if [[ "$GUARD_STOP" == 1 ]]; then
        STOP_REASON=GUARD_STOP
        break
    fi

    if [[ "$RESTRAINED_HINT" == 1 ]]; then
        RELEASE_OUT="$OUT/release_window_$WINDOW_TAG"
        RELEASE_SCRATCH="$SCRATCH_ROOT/release_window_$WINDOW_TAG"

        CURRENT="prepare_release_$WINDOW_TAG"
        "$PY" "$DRIVER" --mode prepare-release \
            --window-index "$WINDOW_INDEX" --root "$OUT" \
            --output "$RELEASE_OUT" --scratch "$RELEASE_SCRATCH" \
            --input-rst7 "$CURRENT_RST7"

        CURRENT="run_release_$WINDOW_TAG"
        (
            cd "$RELEASE_SCRATCH"
            mpirun --bind-to none -np "${SLURM_NTASKS:-8}" sander.MPI -O \
                -i stage.in -o stage.out -p "$PRMTOP" \
                -c "$CURRENT_RST7" -r stage.rst7 \
                -x release.mdcrd -inf stage.mdinfo
        )
        test -s "$RELEASE_SCRATCH/stage.rst7"
        test -s "$RELEASE_SCRATCH/release.mdcrd"

        CURRENT="audit_release_$WINDOW_TAG"
        "$PY" "$DRIVER" --mode audit-release \
            --root "$OUT" --output "$RELEASE_OUT" \
            --scratch "$RELEASE_SCRATCH"

        read -r RELEASE_TECHNICAL CANDIDATE_TOTAL < <(
            "$PY" - "$RELEASE_OUT/RESULT.json" "$OUT/SOURCE_MANIFEST.json" <<'PY'
import json
import sys

release = json.load(open(sys.argv[1], encoding="utf-8"))
manifest = json.load(open(sys.argv[2], encoding="utf-8"))
print(int(release["technical_pass"]), len(manifest["candidates"]))
PY
        )
        test "$RELEASE_TECHNICAL" = 1
        if (( CANDIDATE_TOTAL >= 2 )); then
            STOP_REASON=CANDIDATE_LIMIT_REACHED
            break
        fi
    fi
done

CURRENT=finalize
"$PY" "$DRIVER" --mode finalize --root "$OUT" \
    --stop-reason "$STOP_REASON"
result_is_complete
FINALIZED=1

read -r TECHNICAL_COMPLETE SCIENTIFIC_GATE CANDIDATE_COUNT < <(
    "$PY" - "$OUT/RESULT.json" <<'PY'
import json
import sys

result = json.load(open(sys.argv[1], encoding="utf-8"))
print(
    int(result["technical_complete"]),
    result["scientific_gate"],
    result["candidate_count"],
)
PY
)
test "$TECHNICAL_COMPLETE" = 1

CURRENT=merge_if_ready
merge_if_ready_locked

STATE=PASS_TECHNICAL
CURRENT=terminal_history
DETAIL="array_index=$INDEX;status=PASS_TECHNICAL_A1_PRODUCT_REVERSE_BOUNDARY_SCOUT;scientific_gate=$SCIENTIFIC_GATE;candidate_count=$CANDIDATE_COUNT;stop_reason=$STOP_REASON;source_sha=$EXPECTED_SOURCE_SHA;output=$OUT/RESULT.json"
append_history
trap - EXIT
printf 'PASS_TECHNICAL_A1_PRODUCT_REVERSE_BOUNDARY_SCOUT seed_index=%s gate=%s candidates=%s output=%s\n' \
    "$INDEX" "$SCIENTIFIC_GATE" "$CANDIDATE_COUNT" "$OUT"
