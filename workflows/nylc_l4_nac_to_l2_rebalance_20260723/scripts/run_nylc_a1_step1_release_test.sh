#!/usr/bin/env bash
# Two-stage release of the A1 constrained attack seed under the frozen 146-atom Hamiltonian.
set -euo pipefail
TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
CODE_ROOT="${A1_STEP1_RELEASE_CODE_ROOT:?set immutable code snapshot}"
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
QMMM_SOURCE="$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_unified_core_dftb3_preflight/attempt_62011285"
SEED_SOURCE="$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_unified_step1_qattack_extension/attempt_62021985"
START_RESTART="$SEED_SOURCE/scan/q06_1p45A/window.rst7"
ATTEMPT="${A1_STEP1_RELEASE_ATTEMPT:-${SLURM_JOB_ID:-manual_$(date -u '+%Y%m%dT%H%M%SZ')}}"
OUT="$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_step1_release_test/attempt_$ATTEMPT"
GITHUB_COMMIT="${A1_STEP1_RELEASE_GITHUB_COMMIT:-unknown}"
[[ ! -e "$OUT" ]] || { printf 'refusing to overwrite %s\n' "$OUT" >&2; exit 2; }
mkdir -p "$OUT"
EVENT=nylc_a1_step1_release_test
COMMAND=run_nylc_a1_step1_release_test.sh
STATE=STARTED
CURRENT=initialization
DETAIL="source=$START_RESTART;stages=local_release,full_release;scope=attack_contact_release_not_TS_TI_PMF_barrier;output=$OUT"

append_history() {
    local now
    now="$(date '+%Y-%m-%dT%H:%M:%S%z')"
    (
        flock -x 9
        "$PY" - "$TASK_ROOT/run_history.tsv" "$TASK_ROOT/run_history.jsonl" \
            "$now" "$EVENT" "$COMMAND" "$STATE" "$GITHUB_COMMIT" "$DETAIL" "$OUT" <<'PY'
import json,os,sys
tsv_path,jsonl_path,now,event,command,state,commit,detail,output=sys.argv[1:]
tsv_row="\t".join((now,event,command,state,commit,detail))+"\n"
jsonl_row=json.dumps({"time":now,"event":event,"command":command,"state":state,"git_commit":commit,"detail":detail,"output":output},sort_keys=True)+"\n"
with open(tsv_path,"a+",encoding="utf-8") as tsv, open(jsonl_path,"a+",encoding="utf-8") as jsonl:
    tsv.seek(0,os.SEEK_END); jsonl.seek(0,os.SEEK_END)
    tsv_offset=tsv.tell(); jsonl_offset=jsonl.tell()
    try:
        tsv.write(tsv_row); tsv.flush(); os.fsync(tsv.fileno())
        jsonl.write(jsonl_row); jsonl.flush(); os.fsync(jsonl.fileno())
    except BaseException:
        tsv.seek(tsv_offset); tsv.truncate()
        jsonl.seek(jsonl_offset); jsonl.truncate()
        raise
PY
    ) 9>"$TASK_ROOT/.run_history.lock"
}

finish() {
    local code=$?
    trap - EXIT
    STATE=NOT_EVALUATED
    DETAIL="stage=$CURRENT;exit_code=$code;scientific_status=NOT_EVALUATED;output=$OUT"
    "$PY" - "$OUT/NOT_EVALUATED.json" "$CURRENT" "$code" <<'PY'
import json,pathlib,sys
path,stage,code=sys.argv[1:]
pathlib.Path(path).write_text(json.dumps({
    "schema_version":1,
    "status":"NOT_EVALUATED_A1_STEP1_RELEASE_TEST_TECHNICAL_FAILURE",
    "scientific_status":"NOT_EVALUATED",
    "failed_stage":stage,
    "exit_code":int(code),
},indent=2,sort_keys=True)+"\n",encoding="utf-8")
PY
    append_history || true
    exit "$code"
}
trap finish EXIT
append_history

CURRENT=source_gate
[[ -s "$QMMM_SOURCE/PASS.json" && -s "$QMMM_SOURCE/prepared/system.prmtop" ]]
[[ -s "$SEED_SOURCE/PASS.json" && -s "$START_RESTART" ]]

CURRENT=preparation
"$PY" "$CODE_ROOT/scripts/prepare_nylc_a1_step1_release_test.py" --output "$OUT/release"

module purge >/dev/null 2>&1 || true
module load amber/2018-hpcx-gcc-7.3.1
export AMBERHOME="$TASK_ROOT/../nylc_gyaq_pa66_l2_nac_qmmm_20260723/amber_runtime"
test -s "$AMBERHOME/dat/slko/3ob-3-1/C-C.skf"

CURRENT=local_release
local_dir="$OUT/release/local_release"
(
    cd "$local_dir"
    mpirun --bind-to none -np "${SLURM_NTASKS:-1}" sander.MPI -O \
        -i stage.in -o stage.out \
        -p "$QMMM_SOURCE/prepared/system.prmtop" \
        -c "$START_RESTART" -ref "$START_RESTART" \
        -r stage.rst7 -inf stage.mdinfo
)
test -s "$local_dir/stage.rst7"

CURRENT=full_release
full_dir="$OUT/release/full_release"
(
    cd "$full_dir"
    mpirun --bind-to none -np "${SLURM_NTASKS:-1}" sander.MPI -O \
        -i stage.in -o stage.out \
        -p "$QMMM_SOURCE/prepared/system.prmtop" \
        -c "$local_dir/stage.rst7" \
        -r stage.rst7 -inf stage.mdinfo
)
test -s "$full_dir/stage.rst7"

CURRENT=independent_audit
"$PY" "$CODE_ROOT/scripts/audit_nylc_a1_step1_release_test.py" \
    --release-root "$OUT/release" \
    --source-prmtop "$QMMM_SOURCE/prepared/system.prmtop" \
    --start-restart "$START_RESTART" \
    --output "$OUT/A1_STEP1_RELEASE_TEST_RESULT.json"
"$PY" - "$OUT/A1_STEP1_RELEASE_TEST_RESULT.json" <<'PY'
import json,sys
if json.load(open(sys.argv[1],encoding="utf-8")).get("status")!="PASS_TECHNICAL_A1_STEP1_RELEASE_TEST":
    raise SystemExit("Step1 release test did not technically PASS")
PY

CURRENT=hashing
find "$OUT/release" -type f \( -name 'stage.in' -o -name 'stage.out' -o -name 'stage.rst7' \) -print0 | sort -z | xargs -0 sha256sum >"$OUT/SHA256.tsv"
sha256sum "$OUT/release/RELEASE_MANIFEST.json" "$OUT/A1_STEP1_RELEASE_TEST_RESULT.json" >>"$OUT/SHA256.tsv"
cp "$OUT/A1_STEP1_RELEASE_TEST_RESULT.json" "$OUT/PASS.json"
RELEASE_GATE="$("$PY" - "$OUT/PASS.json" <<'PY'
import json,sys
print(json.load(open(sys.argv[1],encoding="utf-8"))["release_gate"])
PY
)"
STATE=PASS_TECHNICAL
CURRENT=terminal_history
DETAIL="status=PASS_TECHNICAL_A1_STEP1_RELEASE_TEST;release_gate=$RELEASE_GATE;scientific_status=NOT_EVALUATED_TS_PMF_BARRIER;output=$OUT/PASS.json"
append_history
trap - EXIT
