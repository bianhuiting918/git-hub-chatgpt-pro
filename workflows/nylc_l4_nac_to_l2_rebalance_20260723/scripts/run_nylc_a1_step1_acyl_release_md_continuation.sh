#!/usr/bin/env bash
set -euo pipefail

TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
CODE_ROOT="${A1_ACYL_RELEASE_CONT_CODE_ROOT:?set immutable code root}"
GITHUB_COMMIT="${A1_ACYL_RELEASE_CONT_GITHUB_COMMIT:?set immutable GitHub commit}"
INDEX="${SLURM_ARRAY_TASK_ID:?run as array task 0 or 1}"
ARRAY_JOB="${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-manual}}"
ATTEMPT="${ARRAY_JOB}_${INDEX}"
OUTPUT_ROOT="$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_step1_acyl_release_md_continuation"
OUT="$OUTPUT_ROOT/attempt_$ATTEMPT"
SOURCE_OUT="$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_step1_acyl_endpoint_stability/attempt_62186505_$INDEX"
SOURCE_RST7="$SOURCE_OUT/released_endpoint.rst7"
SOURCE_MANIFEST="$SOURCE_OUT/ENDPOINT_MANIFEST.json"
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
DRIVER="$CODE_ROOT/scripts/prepare_audit_nylc_a1_step1_acyl_endpoint_stability.py"
PRMTOP="$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_unified_core_dftb3_preflight/attempt_62011285/prepared/system.prmtop"
SCRATCH_ROOT="${SLURM_TMPDIR:-/tmp}/nylc_a1_acyl_release_md_$ATTEMPT"
STAGE_SCRATCH="$SCRATCH_ROOT/release_md"

case "$INDEX" in
    0) EXPECTED_INPUT_SHA=c667198dc1609890377dfad3921ea79e9b5dfa58400fb429b787313c9e2b2ddb ;;
    1) EXPECTED_INPUT_SHA=d1941473daf8542cb723dc73ca1682dc11fb50fc7d23c3b1a5fcffd747bd55c3 ;;
    *) printf 'invalid array index: %s\n' "$INDEX" >&2; exit 2 ;;
esac

EVENT=nylc_a1_acyl_release_md_continuation
COMMAND=run_nylc_a1_step1_acyl_release_md_continuation.sh
STATE=STARTED
CURRENT=initialization
DETAIL="array_index=$INDEX;source=$SOURCE_RST7;source_sha=$EXPECTED_INPUT_SHA;output=$OUT"
OWNED=0
SEED_FINALIZED=0

append_history() {
    local now
    now="$(date '+%Y-%m-%dT%H:%M:%S%z')"
    (
        flock -x 9
        "$PY" - "$TASK_ROOT/run_history.tsv" "$TASK_ROOT/run_history.jsonl" "$now" \
            "$EVENT" "$COMMAND" "$STATE" "$GITHUB_COMMIT" "$DETAIL" "$OUT" <<'PY'
import json, os, sys
tsv_path,jsonl_path,now,event,command,state,commit,detail,output=sys.argv[1:]
tsv_row="\t".join((now,event,command,state,commit,detail))+"\n"
jsonl_row=json.dumps({"time":now,"event":event,"command":command,"state":state,
 "git_commit":commit,"detail":detail,"output":output},sort_keys=True)+"\n"
with open(tsv_path,"a+",encoding="utf-8") as tsv, open(jsonl_path,"a+",encoding="utf-8") as jsonl:
    tsv.seek(0,os.SEEK_END); jsonl.seek(0,os.SEEK_END)
    toff,joff=tsv.tell(),jsonl.tell()
    try:
        tsv.write(tsv_row); tsv.flush(); os.fsync(tsv.fileno())
        jsonl.write(jsonl_row); jsonl.flush(); os.fsync(jsonl.fileno())
    except BaseException:
        tsv.seek(toff); tsv.truncate(); jsonl.seek(joff); jsonl.truncate(); raise
PY
    ) 9>>"$TASK_ROOT/.run_history.lock"
}

write_hashes() {
    (
        cd "$OUT"
        for name in ENDPOINT_MANIFEST.json FULL_RELEASE_SOURCE_RESULT.json RESULT.json PASS.json \
            NOT_EVALUATED.json release_md_endpoint.rst7; do
            [[ -f "$name" ]] && sha256sum "$name"
        done | sort -k2
    ) >"$OUT/SHA256.tsv"
}

merge_if_ready_locked() {
    (
        flock -x 9
        "$PY" "$DRIVER" --mode merge-if-ready --output-root "$OUTPUT_ROOT" --array-job "$ARRAY_JOB"
    ) 9>>"$OUTPUT_ROOT/.merge.lock"
}

finish() {
    local code=$?
    trap - EXIT
    if (( OWNED && ! SEED_FINALIZED )); then
        CURRENT=failure_finalize_seed
        if "$PY" "$DRIVER" --mode finalize-seed --output "$OUT" --technical-failure; then
            write_hashes || true
            SEED_FINALIZED=1
            merge_if_ready_locked || true
        fi
    fi
    STATE=NOT_EVALUATED
    DETAIL="array_index=$INDEX;stage=$CURRENT;exit_code=$code;status=NOT_EVALUATED_A1_ACYL_RELEASE_MD_CONTINUATION;output=$OUT"
    append_history || true
    exit "$code"
}
trap finish EXIT
append_history

test ! -e "$OUT"
test -s "$SOURCE_RST7"
test -s "$SOURCE_MANIFEST"
ACTUAL_INPUT_SHA="$(sha256sum "$SOURCE_RST7" | awk '{print $1}')"
test "$ACTUAL_INPUT_SHA" = "$EXPECTED_INPUT_SHA"
mkdir -p "$OUTPUT_ROOT" "$SCRATCH_ROOT"

CURRENT=initialize_continuation
"$PY" - "$SOURCE_MANIFEST" "$SOURCE_RST7" "$EXPECTED_INPUT_SHA" "$OUT" "$GITHUB_COMMIT" <<'PY'
import json, pathlib, sys
source_manifest, source_rst7, expected_sha, output, commit = sys.argv[1:]
source = json.load(open(source_manifest, encoding="utf-8"))
stages = source.get("stages", [])
expected_order = ["intermediate", "product", "local_release", "full_release"]
if len(stages) < 4 or [s.get("stage") for s in stages[:4]] != expected_order:
    raise SystemExit("source four-stage order mismatch")
if not all(s.get("technical_pass") is True for s in stages[:4]):
    raise SystemExit("source includes non-PASS pre-release stage")
full_release = stages[3]
if full_release.get("restart_sha256") != expected_sha:
    raise SystemExit("source full-release SHA mismatch")
if source.get("qm_contract", {}).get("qm_atom_count") != 146:
    raise SystemExit("source QM atom count changed")
if source.get("qm_contract", {}).get("qmcharge") != 0:
    raise SystemExit("source QM charge changed")
if source.get("qm_contract", {}).get("electron_count_including_link_h") != 510:
    raise SystemExit("source electron count changed")
if source.get("qm_contract", {}).get("link_atom_count") != 6:
    raise SystemExit("source link count changed")
if source.get("qm_contract", {}).get("step1_qm_water_count") != 0:
    raise SystemExit("Step1 unexpectedly contains QM water")
out = pathlib.Path(output)
out.mkdir(parents=True, exist_ok=False)
manifest = dict(source)
manifest["status"] = "READY_A1_ACYL_RELEASE_MD_CONTINUATION"
manifest["github_commit"] = commit
manifest["stages"] = stages[:4]
manifest["continuation_source"] = {
    "source_manifest": source_manifest,
    "source_restart": source_rst7,
    "source_restart_sha256": expected_sha,
    "source_array_job": "62186505",
}
(out / "ENDPOINT_MANIFEST.json").write_text(
    json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
)
(out / "FULL_RELEASE_SOURCE_RESULT.json").write_text(
    json.dumps(full_release, indent=2, sort_keys=True) + "\n", encoding="utf-8"
)
PY
OWNED=1

module purge >/dev/null 2>&1 || true
module load amber/2018-hpcx-gcc-7.3.1
export AMBERHOME="$TASK_ROOT/../nylc_gyaq_pa66_l2_nac_qmmm_20260723/amber_runtime"
export OMP_NUM_THREADS=1
test -s "$AMBERHOME/dat/slko/3ob-3-1/C-C.skf"

CURRENT=prepare_release_md
"$PY" "$DRIVER" --mode prepare --stage release_md --output "$OUT" \
    --scratch "$STAGE_SCRATCH" --input-rst7 "$SOURCE_RST7" \
    --previous-result "$OUT/FULL_RELEASE_SOURCE_RESULT.json"

CURRENT=run_release_md
(
    cd "$STAGE_SCRATCH"
    mpirun --bind-to none -np "${SLURM_NTASKS:-8}" sander.MPI -O \
        -i stage.in -o stage.out -p "$PRMTOP" -c "$SOURCE_RST7" \
        -r stage.rst7 -x release.mdcrd -inf stage.mdinfo
)
test -s "$STAGE_SCRATCH/stage.rst7"

CURRENT=audit_release_md
"$PY" "$DRIVER" --mode audit-stage --stage release_md --output "$OUT" --scratch "$STAGE_SCRATCH"
read -r TECHNICAL INPUT_SHA OUTPUT_SHA FRAME_COUNT < <("$PY" - "$STAGE_SCRATCH/STAGE_RESULT.json" <<'PY'
import json,sys
r=json.load(open(sys.argv[1],encoding="utf-8"))
print(int(r["technical_pass"]),r["input_restart_sha256"],r["restart_sha256"],r["frame_count"])
PY
)
test "$TECHNICAL" = 1
test "$INPUT_SHA" = "$EXPECTED_INPUT_SHA"
test "$FRAME_COUNT" = 50
test "$(sha256sum "$STAGE_SCRATCH/stage.rst7" | awk '{print $1}')" = "$OUTPUT_SHA"
cp "$STAGE_SCRATCH/stage.rst7" "$OUT/release_md_endpoint.rst7"

CURRENT=finalize_seed
"$PY" "$DRIVER" --mode finalize-seed --output "$OUT"
cp "$STAGE_SCRATCH/stage.rst7" "$OUT/release_md_endpoint.rst7"
write_hashes
SEED_FINALIZED=1

CURRENT=merge_if_ready
merge_if_ready_locked

STATE=PASS_TECHNICAL
CURRENT=terminal_history
DETAIL="array_index=$INDEX;status=PASS_TECHNICAL_A1_ACYL_RELEASE_MD_CONTINUATION;input_sha=$EXPECTED_INPUT_SHA;output=$OUT/RESULT.json"
append_history
trap - EXIT
printf 'PASS_TECHNICAL_A1_ACYL_RELEASE_MD_CONTINUATION seed_index=%s output=%s\n' "$INDEX" "$OUT"
