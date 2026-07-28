#!/usr/bin/env bash
set -euo pipefail

TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
CODE_ROOT="${A1_ACYL_ENDPOINT_CODE_ROOT:?set immutable code root}"
GITHUB_COMMIT="${A1_ACYL_ENDPOINT_GITHUB_COMMIT:?set immutable GitHub commit}"
INDEX="${SLURM_ARRAY_TASK_ID:?run as array task 0 or 1}"
ARRAY_JOB="${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-manual}}"
ATTEMPT="${ARRAY_JOB}_${INDEX}"
OUTPUT_ROOT="$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_step1_acyl_endpoint_stability"
OUT="$OUTPUT_ROOT/attempt_$ATTEMPT"
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
DRIVER="$CODE_ROOT/scripts/prepare_audit_nylc_a1_step1_acyl_endpoint_stability.py"
PRMTOP="$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_unified_core_dftb3_preflight/attempt_62011285/prepared/system.prmtop"
SCRATCH_ROOT="${SLURM_TMPDIR:-/tmp}/nylc_a1_acyl_endpoint_$ATTEMPT"
START_RST7="$SCRATCH_ROOT/start.rst7"

EVENT=nylc_a1_acyl_endpoint_stability
COMMAND=run_nylc_a1_step1_acyl_endpoint_stability.sh
STATE=STARTED
CURRENT=initialization
DETAIL="array_index=$INDEX;scope=acyl_endpoint_stability_not_barrier_or_mechanism;output=$OUT"
OWNED=0
SEED_FINALIZED=0

append_history() {
    local now
    now="$(date '+%Y-%m-%dT%H:%M:%S%z')"
    (
        flock -x 9
        "$PY" - "$TASK_ROOT/run_history.tsv" "$TASK_ROOT/run_history.jsonl"             "$now" "$EVENT" "$COMMAND" "$STATE" "$GITHUB_COMMIT" "$DETAIL" "$OUT" <<'PY'
import json,os,sys
tsv_path,jsonl_path,now,event,command,state,commit,detail,output=sys.argv[1:]
tsv_row="\t".join((now,event,command,state,commit,detail))+"\n"
jsonl_row=json.dumps({"time":now,"event":event,"command":command,"state":state,
 "git_commit":commit,"detail":detail,"output":output},sort_keys=True)+"\n"
with open(tsv_path,"a+",encoding="utf-8") as tsv,open(jsonl_path,"a+",encoding="utf-8") as jsonl:
 tsv.seek(0,os.SEEK_END); jsonl.seek(0,os.SEEK_END); toff,joff=tsv.tell(),jsonl.tell()
 try:
  tsv.write(tsv_row); tsv.flush(); os.fsync(tsv.fileno())
  jsonl.write(jsonl_row); jsonl.flush(); os.fsync(jsonl.fileno())
 except BaseException:
  tsv.seek(toff); tsv.truncate(); jsonl.seek(joff); jsonl.truncate(); raise
PY
    ) 9>>"$TASK_ROOT/run_history.tsv"
}

write_hashes() {
    (
        cd "$OUT"
        for name in ENDPOINT_MANIFEST.json RESULT.json PASS.json NOT_EVALUATED.json                     constructed_product.rst7 released_endpoint.rst7; do
            [[ -f "$name" ]] && sha256sum "$name"
        done | sort -k2
    ) >"$OUT/SHA256.tsv"
}

seed_result_complete() {
    test -s "$OUT/ENDPOINT_MANIFEST.json"
    test -s "$OUT/RESULT.json"
    test -s "$OUT/SHA256.tsv"
    local sentinels=0
    if [[ -s "$OUT/PASS.json" ]]; then sentinels=$((sentinels + 1)); fi
    if [[ -s "$OUT/NOT_EVALUATED.json" ]]; then sentinels=$((sentinels + 1)); fi
    (( sentinels == 1 ))
}

merge_if_ready_locked() {
    (
        flock -x 9
        "$PY" "$DRIVER" --mode merge-if-ready --output-root "$OUTPUT_ROOT" --array-job "$ARRAY_JOB"
    ) 9<"$OUTPUT_ROOT"
}

finish() {
    local code=$?
    trap - EXIT
    if (( OWNED && ! SEED_FINALIZED )); then
        CURRENT=failure_finalize_seed
        if "$PY" "$DRIVER" --mode finalize-seed --output "$OUT" --technical-failure; then
            if write_hashes && seed_result_complete; then
                SEED_FINALIZED=1
                merge_if_ready_locked || true
            fi
        fi
    fi
    STATE=NOT_EVALUATED
    DETAIL="array_index=$INDEX;stage=$CURRENT;exit_code=$code;status=NOT_EVALUATED_A1_ACYL_ENDPOINT_STABILITY;output=$OUT"
    append_history || true
    exit "$code"
}
trap finish EXIT
append_history

mkdir -p "$SCRATCH_ROOT"
CURRENT=initialize
"$PY" "$DRIVER" --mode initialize --seed-index "$INDEX" --output "$OUT"     --start-rst7 "$START_RST7" --github-commit "$GITHUB_COMMIT"
OWNED=1

module purge >/dev/null 2>&1 || true
module load amber/2018-hpcx-gcc-7.3.1
export AMBERHOME="$TASK_ROOT/../nylc_gyaq_pa66_l2_nac_qmmm_20260723/amber_runtime"
export OMP_NUM_THREADS=1
test -s "$AMBERHOME/dat/slko/3ob-3-1/C-C.skf"

CURRENT_RST7="$START_RST7"
PREVIOUS_RESULT=""
for STAGE in intermediate product local_release full_release release_md; do
    STAGE_SCRATCH="$SCRATCH_ROOT/$STAGE"
    NEXT_RST7="$STAGE_SCRATCH/stage.rst7"
    CURRENT="prepare_$STAGE"
    PREVIOUS_ARGS=()
    if [[ "$STAGE" == release_md ]]; then
        PREVIOUS_ARGS=(--previous-result "$PREVIOUS_RESULT")
    fi
    "$PY" "$DRIVER" --mode prepare --stage "$STAGE" --output "$OUT"         --scratch "$STAGE_SCRATCH" --input-rst7 "$CURRENT_RST7" "${PREVIOUS_ARGS[@]}"

    CURRENT="run_$STAGE"
    if [[ "$STAGE" == release_md ]]; then
        (
            cd "$STAGE_SCRATCH"
            mpirun --bind-to none -np "${SLURM_NTASKS:-8}" sander.MPI -O                 -i stage.in -o stage.out -p "$PRMTOP" -c "$CURRENT_RST7"                 -r stage.rst7 -x release.mdcrd -inf stage.mdinfo
        )
    elif [[ "$STAGE" == intermediate || "$STAGE" == product || "$STAGE" == local_release ]]; then
        (
            cd "$STAGE_SCRATCH"
            mpirun --bind-to none -np "${SLURM_NTASKS:-8}" sander.MPI -O                 -i stage.in -o stage.out -p "$PRMTOP" -c "$CURRENT_RST7"                 -ref "$CURRENT_RST7" -r stage.rst7 -inf stage.mdinfo
        )
    else
        (
            cd "$STAGE_SCRATCH"
            mpirun --bind-to none -np "${SLURM_NTASKS:-8}" sander.MPI -O                 -i stage.in -o stage.out -p "$PRMTOP" -c "$CURRENT_RST7"                 -r stage.rst7 -inf stage.mdinfo
        )
    fi
    test -s "$NEXT_RST7"

    CURRENT="audit_$STAGE"
    "$PY" "$DRIVER" --mode audit-stage --stage "$STAGE" --output "$OUT"         --scratch "$STAGE_SCRATCH"
    PREVIOUS_RESULT="$STAGE_SCRATCH/STAGE_RESULT.json"
    read -r TECHNICAL INPUT_SHA OUTPUT_SHA < <("$PY" - "$PREVIOUS_RESULT" <<'PY'
import json,sys
result=json.load(open(sys.argv[1],encoding="utf-8"))
print(int(result["technical_pass"]),result["input_restart_sha256"],result["restart_sha256"])
PY
)
    [[ "$TECHNICAL" == 1 ]]
    [[ "$(sha256sum "$CURRENT_RST7" | awk '{print $1}')" == "$INPUT_SHA" ]]
    [[ "$(sha256sum "$NEXT_RST7" | awk '{print $1}')" == "$OUTPUT_SHA" ]]

    if [[ "$STAGE" == product ]]; then
        cp "$NEXT_RST7" "$OUT/constructed_product.rst7"
        [[ "$(sha256sum "$OUT/constructed_product.rst7" | awk '{print $1}')" == "$OUTPUT_SHA" ]]
    elif [[ "$STAGE" == full_release ]]; then
        cp "$NEXT_RST7" "$OUT/released_endpoint.rst7"
        [[ "$(sha256sum "$OUT/released_endpoint.rst7" | awk '{print $1}')" == "$OUTPUT_SHA" ]]
    fi
    CURRENT_RST7="$NEXT_RST7"
done

CURRENT=finalize_seed
"$PY" "$DRIVER" --mode finalize-seed --output "$OUT"
write_hashes
seed_result_complete
SEED_FINALIZED=1

CURRENT=merge_if_ready
merge_if_ready_locked

STATE=PASS_TECHNICAL
CURRENT=terminal_history
DETAIL="array_index=$INDEX;status=PASS_TECHNICAL_A1_ACYL_ENDPOINT_STABILITY;scientific_status=NOT_EVALUATED_TS_PMF_BARRIER_MECHANISM;output=$OUT/RESULT.json"
append_history
trap - EXIT
printf 'PASS_TECHNICAL_A1_ACYL_ENDPOINT_STABILITY seed_index=%s output=%s\n' "$INDEX" "$OUT"
