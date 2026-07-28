#!/usr/bin/env bash
set -euo pipefail

TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
CODE_ROOT="${A1_PT2_CN_CODE_ROOT:?set immutable code root}"
GITHUB_COMMIT="${A1_PT2_CN_GITHUB_COMMIT:?set immutable GitHub commit}"
INDEX="${SLURM_ARRAY_TASK_ID:?run as array task 0-15}"
ARRAY_JOB="${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-manual}}"
ATTEMPT="${ARRAY_JOB}_${INDEX}"
OUT="$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_pt2_cn_3d_scout/attempt_$ATTEMPT"
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
PREPARE="$CODE_ROOT/scripts/prepare_audit_nylc_a1_step1_pt2_cn_scout.py"
PRMTOP="$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_unified_core_dftb3_preflight/attempt_62011285/prepared/system.prmtop"
SCRATCH_ROOT="${SLURM_TMPDIR:-/tmp}/nylc_a1_pt2_cn_${ATTEMPT}"
START_RST7="$SCRATCH_ROOT/start.rst7"
GUIDE_RST7="$SCRATCH_ROOT/guide.rst7"

EVENT=nylc_a1_pt2_cn_3d_scout
COMMAND=run_nylc_a1_step1_pt2_cn_scout.sh
STATE=STARTED
CURRENT=initialization
DETAIL="array_index=$INDEX;scope=bounded_scout_not_TS_PMF_barrier;output=$OUT"
OWNED=0

append_history() {
    local now
    now="$(date '+%Y-%m-%dT%H:%M:%S%z')"
    (
        flock -x 9
        "$PY" - "$TASK_ROOT/run_history.tsv" "$TASK_ROOT/run_history.jsonl" \
            "$now" "$EVENT" "$COMMAND" "$STATE" "$GITHUB_COMMIT" \
            "$DETAIL" "$OUT" <<'PY'
import json
import os
import sys
tsv_path,jsonl_path,now,event,command,state,commit,detail,output=sys.argv[1:]
tsv_row="\t".join((now,event,command,state,commit,detail))+"\n"
jsonl_row=json.dumps({
    "time":now,"event":event,"command":command,"state":state,
    "git_commit":commit,"detail":detail,"output":output,
},sort_keys=True)+"\n"
with open(tsv_path,"a+",encoding="utf-8") as tsv, open(
    jsonl_path,"a+",encoding="utf-8"
) as jsonl:
    tsv.seek(0,os.SEEK_END); jsonl.seek(0,os.SEEK_END)
    toff,joff=tsv.tell(),jsonl.tell()
    try:
        tsv.write(tsv_row); tsv.flush(); os.fsync(tsv.fileno())
        jsonl.write(jsonl_row); jsonl.flush(); os.fsync(jsonl.fileno())
    except BaseException:
        tsv.seek(toff); tsv.truncate()
        jsonl.seek(joff); jsonl.truncate()
        raise
PY
    ) 9>"$TASK_ROOT/.run_history.lock"
}

finish() {
    local code=$?
    trap - EXIT
    if ((OWNED)); then
        "$PY" - "$OUT/NOT_EVALUATED.json" "$CURRENT" "$code" "$INDEX" <<'PY'
import json,pathlib,sys
path,stage,code,index=sys.argv[1:]
target=pathlib.Path(path)
if not target.exists() and not target.with_name("PASS.json").exists():
    target.write_text(json.dumps({
        "schema_version":1,
        "status":"NOT_EVALUATED_A1_PT2_CN_SCOUT",
        "scientific_status":"NOT_EVALUATED",
        "array_index":int(index),
        "failed_stage":stage,
        "exit_code":int(code),
    },indent=2,sort_keys=True)+"\n",encoding="utf-8")
PY
    fi
    STATE=NOT_EVALUATED
    DETAIL="array_index=$INDEX;stage=$CURRENT;exit_code=$code;scientific_status=NOT_EVALUATED;output=$OUT"
    append_history || true
    exit "$code"
}
trap finish EXIT
append_history

CURRENT=prepare
mkdir -p "$SCRATCH_ROOT"
"$PY" "$PREPARE" \
    --mode prepare \
    --index "$INDEX" \
    --output "$OUT" \
    --start-rst7 "$START_RST7" \
    --github-commit "$GITHUB_COMMIT"
OWNED=1

module purge >/dev/null 2>&1 || true
module load amber/2018-hpcx-gcc-7.3.1
export AMBERHOME="$TASK_ROOT/../nylc_gyaq_pa66_l2_nac_qmmm_20260723/amber_runtime"
export OMP_NUM_THREADS=1
test -s "$AMBERHOME/dat/slko/3ob-3-1/C-C.skf"

CURRENT=guide
(
    cd "$OUT/guide"
    mpirun --bind-to none -np "${SLURM_NTASKS:-8}" sander.MPI -O \
        -i stage.in -o stage.out \
        -p "$PRMTOP" \
        -c "$START_RST7" -ref "$START_RST7" \
        -r "$GUIDE_RST7" -inf stage.mdinfo
)
test -s "$GUIDE_RST7"

CURRENT=target
(
    cd "$OUT/target"
    mpirun --bind-to none -np "${SLURM_NTASKS:-8}" sander.MPI -O \
        -i stage.in -o stage.out \
        -p "$PRMTOP" \
        -c "$GUIDE_RST7" -ref "$GUIDE_RST7" \
        -r stage.rst7 -inf stage.mdinfo
)
test -s "$OUT/target/stage.rst7"

CURRENT=audit
"$PY" "$PREPARE" \
    --mode audit \
    --index "$INDEX" \
    --output "$OUT"
test -s "$OUT/PASS.json"

CURRENT=hashing
sha256sum \
    "$OUT/SCOUT_MANIFEST.json" \
    "$OUT/guide/stage.in" \
    "$OUT/guide/restraints.RST" \
    "$OUT/guide/stage.out" \
    "$OUT/target/stage.in" \
    "$OUT/target/restraints.RST" \
    "$OUT/target/stage.out" \
    "$OUT/target/stage.rst7" \
    "$OUT/RESULT.json" \
    "$OUT/PASS.json" >"$OUT/SHA256.tsv"

CANDIDATE_GATE="$("$PY" - "$OUT/PASS.json" <<'PY'
import json,sys
print(json.load(open(sys.argv[1],encoding="utf-8"))["candidate_gate"])
PY
)"
STATE=PASS_TECHNICAL
CURRENT=terminal_history
DETAIL="array_index=$INDEX;status=PASS_TECHNICAL_A1_PT2_CN_SCOUT;candidate_gate=$CANDIDATE_GATE;scientific_status=NOT_EVALUATED_TS_PMF_BARRIER;output=$OUT/PASS.json"
append_history
trap - EXIT
printf 'PASS_TECHNICAL_A1_PT2_CN_SCOUT index=%s gate=%s output=%s\n' \
    "$INDEX" "$CANDIDATE_GATE" "$OUT"
