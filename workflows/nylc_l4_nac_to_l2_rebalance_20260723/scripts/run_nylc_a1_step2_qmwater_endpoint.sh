#!/usr/bin/env bash
# Minimal NylC Step2 A2 QM-water preflight.  Two independent eight-rank legs
# run concurrently inside each sixteen-rank array task; this never constructs product.
set -euo pipefail

TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
CODE_ROOT="${A1_STEP2_QMWATER_CODE_ROOT:?set immutable workflow code root}"
GITHUB_COMMIT="${A1_STEP2_QMWATER_GITHUB_COMMIT:?set runtime code snapshot commit}"
INDEX="${SLURM_ARRAY_TASK_ID:?run as Slurm array task 0 or 1}"
ARRAY_JOB="${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-manual}}"
ATTEMPT="${ARRAY_JOB}_${INDEX}"
OUTPUT_ROOT="$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_step2_qmwater_a2_preflight"
OUT="$OUTPUT_ROOT/attempt_$ATTEMPT"
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
DRIVER="$CODE_ROOT/scripts/prepare_audit_nylc_a1_step2_qmwater_endpoint.py"
PRMTOP="$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_unified_core_dftb3_preflight/attempt_62011285/prepared/system.prmtop"
SCRATCH_ROOT="${SLURM_TMPDIR:-/tmp}/nylc_a1_step2_qmwater_$ATTEMPT"

case "$INDEX" in
    0)
        SEED=seed26723
        SOURCE_RST7="$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_step1_acyl_release_md_continuation/attempt_62216380_0/release_md_endpoint.rst7"
        EXPECTED_SOURCE_SHA=5d8f76d2c90e3e8c707b640c55a93938f53e18dc30d6f93d54adda467da25f41
        ;;
    1)
        SEED=seed26737
        SOURCE_RST7="$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_step1_acyl_release_md_continuation/attempt_62216380_1/release_md_endpoint.rst7"
        EXPECTED_SOURCE_SHA=4cc60ad4d7be099b3f76040f1ee8c49b91deecebb20ed98570bc192b3d511132
        ;;
    *)
        printf 'invalid array index: %s\n' "$INDEX" >&2
        exit 2
        ;;
esac

EVENT=nylc_a1_step2_qmwater_a2_preflight
COMMAND=run_nylc_a1_step2_qmwater_endpoint.sh
STATE=STARTED
CURRENT=initialization
DETAIL="array_index=$INDEX;seed=$SEED;source=$SOURCE_RST7;source_sha=$EXPECTED_SOURCE_SHA;execution=parallel_two_8rank_legs;output=$OUT"
OWNED=0
TERMINAL=0

append_history() {
    local now
    now="$(date '+%Y-%m-%dT%H:%M:%S%z')"
    (
        flock -x 9
        "$PY" - "$TASK_ROOT/run_history.tsv" "$TASK_ROOT/run_history.jsonl"             "$now" "$EVENT" "$COMMAND" "$STATE" "$GITHUB_COMMIT"             "$DETAIL" "$OUT" <<'PY'
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
    ) 9>>"$TASK_ROOT/.run_history.lock"
}

write_hashes() {
    [[ -d "$OUT" ]] || return 0
    (
        cd "$OUT"
        for name in             A2_MANIFEST.json READY.json WATER_SELECTION.json WATER_CANDIDATES.tsv             A2_LEG_0.json A2_LEG_1.json RESULT.json PASS.json FAIL.json             NOT_EVALUATED.json a2_leg0_endpoint.rst7 a2_leg1_endpoint.rst7
        do
            [[ -f "$name" ]] && sha256sum "$name"
        done | sort -k2
    ) >"$OUT/SHA256.tsv"
}

finish() {
    local code=$?
    trap - EXIT
    if (( code != 0 && ! TERMINAL )); then
        CURRENT=unexpected_failure_finalize
        "$PY" - "$OUT" "$SEED" "$CURRENT" "$code" <<'PY' || true
import json
import pathlib
import sys

out, seed, stage, code = sys.argv[1:]
root = pathlib.Path(out)
root.mkdir(parents=True, exist_ok=True)
payload = {
    "schema_version": 1,
    "status": "NOT_EVALUATED_STEP2_A2_RUNTIME_FAILURE",
    "seed": seed,
    "classification": "NOT_EVALUATED_STEP2_A2_RUNTIME_FAILURE",
    "failed_stage": stage,
    "exit_code": int(code),
    "NEXT": "STEP2_PRODUCT_ENDPOINT_BLOCKED_PENDING_A2_PASS",
    "product_endpoint_implemented": False,
}
for name in ("RESULT.json", "NOT_EVALUATED.json"):
    (root / name).write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
PY
        write_hashes || true
        STATE=NOT_EVALUATED
        DETAIL="array_index=$INDEX;seed=$SEED;stage=$CURRENT;exit_code=$code;status=NOT_EVALUATED_STEP2_A2_RUNTIME_FAILURE;output=$OUT"
        append_history || true
    fi
    exit "$code"
}
trap finish EXIT

append_history
test ! -e "$OUT"
test -s "$DRIVER"
test -s "$PRMTOP"
test -s "$SOURCE_RST7"
mkdir -p "$OUTPUT_ROOT"

module purge >/dev/null 2>&1 || true
module load amber/2018-hpcx-gcc-7.3.1
export AMBERHOME="$TASK_ROOT/../nylc_gyaq_pa66_l2_nac_qmmm_20260723/amber_runtime"
export OMP_NUM_THREADS=1
test -s "$AMBERHOME/dat/slko/3ob-3-1/C-C.skf"

CURRENT=prepare_water_and_contract
"$PY" "$DRIVER"     --mode prepare     --seed-index "$INDEX"     --output "$OUT"     --scratch "$SCRATCH_ROOT"     --code-root "$CODE_ROOT"     --github-commit "$GITHUB_COMMIT"
OWNED=1

if [[ ! -s "$OUT/READY.json" ]]; then
    test -s "$OUT/RESULT.json"
    write_hashes
    CLASSIFICATION="$("$PY" - "$OUT/RESULT.json" <<'PY'
import json, sys
print(json.load(open(sys.argv[1], encoding="utf-8"))["classification"])
PY
)"
    STATE=NOT_EVALUATED
    CURRENT=terminal_non_evaluated_preparation
    DETAIL="array_index=$INDEX;seed=$SEED;classification=$CLASSIFICATION;exit_code=0;next=STEP2_PRODUCT_ENDPOINT_BLOCKED_PENDING_A2_PASS;output=$OUT/RESULT.json"
    append_history
    TERMINAL=1
    trap - EXIT
    exit 0
fi

# The task owns sixteen ranks. Launch two exclusive eight-rank Slurm steps
# concurrently; each writes only to its own scratch directory.
declare -a LEG_PIDS=()
for LEG in 0 1; do
    CURRENT="launch_a2_leg_$LEG"
    STAGE="$SCRATCH_ROOT/a2_leg$LEG"
    test -s "$STAGE/stage.in"
    (
        cd "$STAGE"
        set +e
        srun --exclusive -N 1 -n 8 sander.MPI -O \
            -i stage.in -o stage.out -p "$PRMTOP" -c "$SOURCE_RST7" \
            -r stage.rst7 -x a2.mdcrd -inf stage.mdinfo
        ENGINE_RC=$?
        printf '%s\n' "$ENGINE_RC" >engine.rc
        exit 0
    ) &
    LEG_PIDS[$LEG]=$!
done

CURRENT=wait_for_parallel_a2_legs
for LEG in 0 1; do
    wait "${LEG_PIDS[$LEG]}"
done

for LEG in 0 1; do
    CURRENT="audit_a2_leg_$LEG"
    STAGE="$SCRATCH_ROOT/a2_leg$LEG"
    "$PY" "$DRIVER" \
        --mode audit-leg \
        --output "$OUT" \
        --scratch "$SCRATCH_ROOT" \
        --leg "$LEG"

    TECHNICAL="$("$PY" - "$OUT/A2_LEG_$LEG.json" <<'PY'
import json, sys
print(int(json.load(open(sys.argv[1], encoding="utf-8"))["technical_pass"]))
PY
)"
    if [[ "$TECHNICAL" = 1 ]]; then
        test -s "$STAGE/stage.rst7"
        cp "$STAGE/stage.rst7" "$OUT/a2_leg${LEG}_endpoint.rst7"
    fi
done

CURRENT=finalize_a2_denominator
"$PY" "$DRIVER" --mode finalize --output "$OUT"
test -s "$OUT/RESULT.json"
write_hashes

read -r STATUS CLASSIFICATION < <("$PY" - "$OUT/RESULT.json" <<'PY'
import json, sys
result = json.load(open(sys.argv[1], encoding="utf-8"))
print(result["status"], result["classification"])
PY
)
case "$STATUS" in
    PASS_EXPLORATORY_STEP2_A2_ACYL_BASIN_REVALIDATED)
        STATE=PASS_EXPLORATORY_A2
        ;;
    FAIL_STEP2_A2_NOT_REVALIDATED)
        STATE=FAIL_SCIENTIFIC_A2_ENDPOINT
        ;;
    *)
        STATE=NOT_EVALUATED
        ;;
esac
CURRENT=terminal_history
DETAIL="array_index=$INDEX;seed=$SEED;status=$STATUS;classification=$CLASSIFICATION;legs=2;execution=parallel_2x8ranks;trajectory=scratch_only;next=STEP2_PRODUCT_ENDPOINT_BLOCKED_PENDING_A2_PASS;output=$OUT/RESULT.json"
append_history
TERMINAL=1
trap - EXIT
printf '%s seed=%s output=%s\n' "$STATUS" "$SEED" "$OUT"
