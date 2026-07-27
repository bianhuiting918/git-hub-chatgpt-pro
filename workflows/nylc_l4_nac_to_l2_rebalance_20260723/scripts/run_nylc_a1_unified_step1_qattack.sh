#!/usr/bin/env bash
# Sequential 146-QM-atom Step1 q_attack scouting; constrained seed generation only.
set -euo pipefail
TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
CODE_ROOT="${A1_STEP1_QATTACK_CODE_ROOT:?set immutable code snapshot}"
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
SOURCE="$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_unified_core_dftb3_preflight/attempt_62011285"
ATTEMPT="${A1_STEP1_QATTACK_ATTEMPT:-${SLURM_JOB_ID:-manual_$(date -u '+%Y%m%dT%H%M%SZ')}}"
OUT="$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_unified_step1_qattack/attempt_$ATTEMPT"
GITHUB_COMMIT="${A1_STEP1_QATTACK_GITHUB_COMMIT:-unknown}"
[[ ! -e "$OUT" ]] || { printf 'refusing to overwrite %s\n' "$OUT" >&2; exit 2; }
mkdir -p "$OUT"
EVENT=nylc_a1_unified_step1_qattack
COMMAND=run_nylc_a1_unified_step1_qattack.sh
STATE=STARTED
CURRENT=initialization
DETAIL="source=$SOURCE;scope=constrained_attack_scout_not_TS_PMF_barrier;output=$OUT"
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
pathlib.Path(path).write_text(json.dumps({"schema_version":1,"status":"NOT_EVALUATED_A1_UNIFIED_STEP1_QATTACK_TECHNICAL_FAILURE","scientific_status":"NOT_EVALUATED","failed_stage":stage,"exit_code":int(code)},indent=2,sort_keys=True)+"\n",encoding="utf-8")
PY
    append_history || true
    exit "$code"
}
trap finish EXIT
append_history
CURRENT=source_gate
[[ -s "$SOURCE/PASS.json" && -s "$SOURCE/SHA256.tsv" ]]
"$PY" - "$SOURCE/PASS.json" <<'PY'
import json,sys
if json.load(open(sys.argv[1],encoding="utf-8")).get("status")!="PASS_A1_UNIFIED_CORE_DFTB3_NUMERICAL_PREFLIGHT":
    raise SystemExit("unified-core preflight is not PASS")
PY
CURRENT=preparation
"$PY" "$CODE_ROOT/scripts/prepare_nylc_a1_unified_step1_qattack.py" --output "$OUT/scan"
module purge >/dev/null 2>&1 || true
module load amber/2018-hpcx-gcc-7.3.1
export AMBERHOME="$TASK_ROOT/../nylc_gyaq_pa66_l2_nac_qmmm_20260723/amber_runtime"
test -s "$AMBERHOME/dat/slko/3ob-3-1/C-C.skf"
mapfile -t WINDOWS < <("$PY" - "$OUT/scan/SCAN_MANIFEST.json" <<'PY'
import json,sys
for item in json.load(open(sys.argv[1],encoding="utf-8"))["windows"]:
    print(item["name"])
PY
)
[[ "${#WINDOWS[@]}" -eq 9 ]]
current="$SOURCE/prepared/02_qmmm_20_step.rst7"
for name in "${WINDOWS[@]}"; do
    CURRENT="$name"
    directory="$OUT/scan/$name"
    (
        cd "$directory"
        mpirun --bind-to none -np "${SLURM_NTASKS:-1}" sander.MPI -O \
            -i window.in -o window.out \
            -p "$SOURCE/prepared/system.prmtop" \
            -c "$current" -ref "$current" \
            -r window.rst7 -inf window.mdinfo
    )
    test -s "$directory/window.rst7"
    current="$directory/window.rst7"
done
CURRENT=independent_audit
"$PY" "$CODE_ROOT/scripts/audit_nylc_a1_unified_step1_qattack.py" \
    --scan-root "$OUT/scan" \
    --source-prmtop "$SOURCE/prepared/system.prmtop" \
    --output "$OUT/A1_UNIFIED_STEP1_QATTACK_RESULT.json"
"$PY" - "$OUT/A1_UNIFIED_STEP1_QATTACK_RESULT.json" <<'PY'
import json,sys
if json.load(open(sys.argv[1],encoding="utf-8")).get("status")!="PASS_TECHNICAL_A1_UNIFIED_STEP1_QATTACK_SCAN":
    raise SystemExit("Step1 q_attack audit did not technically PASS")
PY
CURRENT=hashing
find "$OUT/scan" -type f \( -name 'window.in' -o -name 'window.RST' -o -name 'window.out' -o -name 'window.rst7' \) -print0 | sort -z | xargs -0 sha256sum >"$OUT/SHA256.tsv"
sha256sum "$OUT/scan/SCAN_MANIFEST.json" "$OUT/A1_UNIFIED_STEP1_QATTACK_RESULT.json" >>"$OUT/SHA256.tsv"
cp "$OUT/A1_UNIFIED_STEP1_QATTACK_RESULT.json" "$OUT/PASS.json"
STATE=PASS_TECHNICAL
CURRENT=terminal_history
DETAIL="status=PASS_TECHNICAL_A1_UNIFIED_STEP1_QATTACK_SCAN;scientific_status=NOT_EVALUATED_TS_PMF_BARRIER;output=$OUT/PASS.json"
append_history
trap - EXIT
