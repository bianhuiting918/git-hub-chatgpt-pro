#!/usr/bin/env bash
# Analyze all corrected nine-replica A1 trajectories for NAC-conditioned PT2 geometry.
set -euo pipefail
TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
CODE_ROOT="${A1_PT2_NAC_CODE_ROOT:?set immutable code snapshot}"
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
AUDIT_ROOT="$TASK_ROOT/a1_activated_nac_20260726/nac_audit"
ATTEMPT="${A1_PT2_NAC_ATTEMPT:-${SLURM_JOB_ID:-manual_$(date -u '+%Y%m%dT%H%M%SZ')}}"
OUT="$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_pt2_nac_ensemble_scan/attempt_$ATTEMPT"
GITHUB_COMMIT="${A1_PT2_NAC_GITHUB_COMMIT:-unknown}"
[[ ! -e "$OUT" ]] || { printf 'refusing to overwrite %s\n' "$OUT" >&2; exit 2; }
mkdir -p "$OUT"
EVENT=nylc_a1_pt2_nac_ensemble_scan
COMMAND=run_nylc_a1_pt2_nac_ensemble.sh
STATE=STARTED
CURRENT=initialization
DETAIL="source_audit_array=61976600;replica_denominator=9;scope=fixed_topology_geometry_not_proton_transfer_TS_PMF_barrier;output=$OUT"

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
    "status":"NOT_EVALUATED_A1_PT2_NAC_ENSEMBLE_TECHNICAL_FAILURE",
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

CURRENT=analysis
"$PY" "$CODE_ROOT/scripts/analyze_nylc_a1_pt2_nac_ensemble.py" \
    --audit-root "$AUDIT_ROOT" \
    --output "$OUT/analysis"
test -s "$OUT/analysis/RESULT.json"
test -s "$OUT/analysis/candidate_frames.tsv"

CURRENT=terminal_gate
"$PY" - "$OUT/analysis/RESULT.json" <<'PY'
import json,sys
result=json.load(open(sys.argv[1],encoding="utf-8"))
if result.get("status")!="PASS_TECHNICAL_A1_PT2_NAC_ENSEMBLE_SCAN":
    raise SystemExit("PT2 NAC ensemble scan did not technically PASS")
if result.get("actual_replica_denominator")!=9:
    raise SystemExit("replica denominator is not nine")
PY
sha256sum "$OUT/analysis/RESULT.json" "$OUT/analysis/candidate_frames.tsv" >"$OUT/SHA256.tsv"
cp "$OUT/analysis/RESULT.json" "$OUT/PASS.json"
GATE="$("$PY" - "$OUT/PASS.json" <<'PY'
import json,sys
print(json.load(open(sys.argv[1],encoding="utf-8"))["pt2_preorganization_gate"])
PY
)"
STATE=PASS_TECHNICAL
CURRENT=terminal_history
DETAIL="status=PASS_TECHNICAL_A1_PT2_NAC_ENSEMBLE_SCAN;pt2_gate=$GATE;scientific_status=NOT_EVALUATED_PROTON_TRANSFER_TS_PMF_BARRIER_MECHANISM;output=$OUT/PASS.json"
append_history
trap - EXIT
