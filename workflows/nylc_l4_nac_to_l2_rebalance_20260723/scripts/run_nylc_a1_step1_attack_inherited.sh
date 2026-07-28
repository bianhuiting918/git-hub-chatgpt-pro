#!/usr/bin/env bash
set -euo pipefail

TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
CODE_ROOT="${A1_ATTACK_CODE_ROOT:?set immutable code root}"
GITHUB_COMMIT="${A1_ATTACK_GITHUB_COMMIT:?set immutable GitHub commit}"
INDEX="${SLURM_ARRAY_TASK_ID:?run as array task 0 or 1}"
ARRAY_JOB="${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-manual}}"
ATTEMPT="${ARRAY_JOB}_${INDEX}"
OUT="$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_attack_inherited_chain/attempt_$ATTEMPT"
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
PREPARE="$CODE_ROOT/scripts/prepare_audit_nylc_a1_step1_attack_inherited.py"
PRMTOP="$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_unified_core_dftb3_preflight/attempt_62011285/prepared/system.prmtop"
SCRATCH_ROOT="${SLURM_TMPDIR:-/tmp}/nylc_a1_attack_inherited_${ATTEMPT}"
CHAIN_STATE="$OUT/CHAIN_START.json"
START_RST7="$SCRATCH_ROOT/start.rst7"

EVENT=nylc_a1_attack_inherited_chain
COMMAND=run_nylc_a1_step1_attack_inherited.sh
STATE=STARTED
CURRENT=initialization
DETAIL="array_index=$INDEX;scope=inherited_attack_only_not_TS_PMF_barrier;output=$OUT"
OWNED=0
CHAIN_TECHNICAL=1
INHERITANCE_SHA256_VERIFIED=0

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
    ) 9>"$TASK_ROOT/.run_history.lock"
}

finish() {
    local code=$?
    trap - EXIT
    if (( OWNED )) && [[ ! -s "$OUT/CHAIN_MANIFEST.json" ]]; then
        "$PY" - "$OUT/CHAIN_MANIFEST.json" "$INDEX" "$CURRENT" "$code" "$GITHUB_COMMIT" <<'PY'
import json,pathlib,sys
path,index,stage,code,commit=sys.argv[1:]
pathlib.Path(path).write_text(json.dumps({"schema_version":1,
 "status":"NOT_EVALUATED_A1_ATTACK_INHERITED_CHAIN",
 "scientific_status":"NOT_EVALUATED_TS_PMF_BARRIER_MECHANISM","next_action":"DO_NOT_START_PMF",
 "seed_index":int(index),"failed_stage":stage,"exit_code":int(code),"github_commit":commit},
 indent=2,sort_keys=True)+"\n",encoding="utf-8")
PY
    fi
    STATE=NOT_EVALUATED
    DETAIL="array_index=$INDEX;stage=$CURRENT;exit_code=$code;scientific_status=NOT_EVALUATED_TS_PMF_BARRIER_MECHANISM;output=$OUT"
    append_history || true
    exit "$code"
}
trap finish EXIT
append_history

mkdir -p "$OUT" "$SCRATCH_ROOT"
OWNED=1
CURRENT=initialize
"$PY" "$PREPARE" --mode initialize --seed-index "$INDEX" --start-rst7 "$START_RST7" \
    --chain-state "$CHAIN_STATE" --github-commit "$GITHUB_COMMIT"

module purge >/dev/null 2>&1 || true
module load amber/2018-hpcx-gcc-7.3.1
export AMBERHOME="$TASK_ROOT/../nylc_gyaq_pa66_l2_nac_qmmm_20260723/amber_runtime"
export OMP_NUM_THREADS=1
test -s "$AMBERHOME/dat/slko/3ob-3-1/C-C.skf"

CURRENT_RST7="$START_RST7"
FIRST_CANDIDATE=""
SHORTEST_CANDIDATE=""
SHORTEST_ATTACK=""
TERMINAL_RST7=""
STOP_REASON="COMPLETED_NINE_WINDOWS"

for WINDOW_INDEX in 0 1 2 3 4 5 6 7 8; do
    WINDOW_OUT="$OUT/window_$(printf '%02d' "$WINDOW_INDEX")"
    WINDOW_SCRATCH="$SCRATCH_ROOT/window_$(printf '%02d' "$WINDOW_INDEX")"
    NEXT_RST7="$WINDOW_SCRATCH/stage.rst7"
    mkdir -p "$WINDOW_SCRATCH"
    CURRENT=prepare_window_$WINDOW_INDEX
    "$PY" "$PREPARE" --mode prepare --seed-index "$INDEX" --window-index "$WINDOW_INDEX" \
        --output "$WINDOW_OUT" --start-rst7 "$CURRENT_RST7" --chain-state "$CHAIN_STATE" \
        --github-commit "$GITHUB_COMMIT"
    # Amber resolves DISANG relative to its working directory, not to stage.in.
    cp "$WINDOW_OUT/restraints.RST" "$WINDOW_SCRATCH/restraints.RST"
    CURRENT=run_window_$WINDOW_INDEX
    (
        cd "$WINDOW_SCRATCH"
        mpirun --bind-to none -np "${SLURM_NTASKS:-8}" sander.MPI -O \
            -i "$WINDOW_OUT/stage.in" -o stage.out -p "$PRMTOP" \
            -c "$CURRENT_RST7" -ref "$CURRENT_RST7" -r stage.rst7 -inf stage.mdinfo
    )
    test -s "$NEXT_RST7"
    CURRENT=audit_window_$WINDOW_INDEX
    "$PY" "$PREPARE" --mode audit --seed-index "$INDEX" --output "$WINDOW_OUT" \
        --scratch "$WINDOW_SCRATCH" --chain-state "$CHAIN_STATE"
    test -s "$WINDOW_OUT/RESULT.json"
    TERMINAL_RST7="$NEXT_RST7"
    read -r TECHNICAL CANDIDATE GUARD_STOP < <("$PY" - "$WINDOW_OUT/RESULT.json" <<'PY'
import json,sys
r=json.load(open(sys.argv[1],encoding="utf-8"))
print(int(r["gates"]["technical"]),int(r["gates"]["TETRAHEDRAL_LIKE_RESTRAINED"]),
      int(r["stopped_by_reactant_guard"]))
PY
)
    if [[ "$TECHNICAL" != 1 ]]; then
        CHAIN_TECHNICAL=0
        STOP_REASON="TECHNICAL_NOT_EVALUATED"
        break
    fi
    if [[ "$CANDIDATE" == 1 ]]; then
        [[ -n "$FIRST_CANDIDATE" ]] || FIRST_CANDIDATE="$NEXT_RST7"
        ATTACK="$("$PY" - "$WINDOW_OUT/RESULT.json" <<'PY'
import json,sys
print(json.load(open(sys.argv[1],encoding="utf-8"))["observed"]["attack_A"])
PY
)"
        if [[ -z "$SHORTEST_ATTACK" ]] || "$PY" - "$ATTACK" "$SHORTEST_ATTACK" <<'PY'
import sys
raise SystemExit(0 if float(sys.argv[1]) < float(sys.argv[2]) else 1)
PY
        then
            SHORTEST_ATTACK="$ATTACK"
            SHORTEST_CANDIDATE="$NEXT_RST7"
        fi
    fi
    if [[ "$GUARD_STOP" == 1 ]]; then
        STOP_REASON="OBSERVED_COUPLED_RESPONSE_OUTSIDE_REACTANT_GUARD"
        break
    fi
    CURRENT_RST7="$NEXT_RST7"
done

CURRENT=finalize
# Verify provenance from the persisted manifests/results, not only shell variables.
INHERITANCE_SHA256_VERIFIED="$("$PY" - "$OUT" "$CHAIN_STATE" <<'PY'
import json,pathlib,sys
root,state_path=map(pathlib.Path,sys.argv[1:])
state=json.load(open(state_path,encoding="utf-8"))
previous=state["start_rst7_sha256"]
records=[]
for manifest_path in sorted(root.glob("window_*/WINDOW_MANIFEST.json")):
    result_path=manifest_path.with_name("RESULT.json")
    if not result_path.is_file():
        break
    manifest=json.load(open(manifest_path,encoding="utf-8"))
    result=json.load(open(result_path,encoding="utf-8"))
    observed_input=manifest.get("inherited_restart_sha256")
    observed_output=result.get("output_restart_sha256")
    if observed_input != previous or not observed_output:
        raise SystemExit("restart SHA256 inheritance mismatch")
    records.append({"window":manifest["window"]["window_index"],
      "input_restart":manifest["inherited_restart"],"input_sha256":observed_input,
      "output_restart":result.get("output_restart"),"output_sha256":observed_output})
    previous=observed_output
(root/"inheritance_verification.json").write_text(json.dumps({
 "inheritance_sha256_verified":True,"start_rst7":state["start_rst7"],
 "start_rst7_sha256":state["start_rst7_sha256"],"records":records},
 indent=2,sort_keys=True)+"\n",encoding="utf-8")
print(1)
PY
)"

if [[ -n "$FIRST_CANDIDATE" ]]; then
    cp "$FIRST_CANDIDATE" "$OUT/first_tetrahedral_like_candidate.rst7"
    if [[ "$SHORTEST_CANDIDATE" != "$FIRST_CANDIDATE" ]]; then
        cp "$SHORTEST_CANDIDATE" "$OUT/shortest_attack_tetrahedral_like_candidate.rst7"
    fi
elif [[ -n "$TERMINAL_RST7" ]]; then
    cp "$TERMINAL_RST7" "$OUT/terminal_endpoint.rst7"
fi
if [[ "$STOP_REASON" == "OBSERVED_COUPLED_RESPONSE_OUTSIDE_REACTANT_GUARD" ]]; then
    cp "$TERMINAL_RST7" "$OUT/guard_stop_endpoint.rst7"
fi

"$PY" - "$OUT" "$INDEX" "$GITHUB_COMMIT" "$STOP_REASON" "$CHAIN_TECHNICAL" \
    "$FIRST_CANDIDATE" "$SHORTEST_CANDIDATE" "$TERMINAL_RST7" "$INHERITANCE_SHA256_VERIFIED" <<'PY'
import hashlib,json,pathlib,sys
out,index,commit,reason,technical,first,shortest,terminal,verified=sys.argv[1:]
root=pathlib.Path(out)
def record(path):
 p=pathlib.Path(path)
 if not path or not p.is_file(): return None
 h=hashlib.sha256(p.read_bytes()).hexdigest()
 return {"path":str(p),"sha256":h}
windows=[json.load(open(p,encoding="utf-8")) for p in sorted(root.glob("window_*/RESULT.json"))]
payload={"schema_version":1,
 "status":"PASS_TECHNICAL_A1_ATTACK_INHERITED_CHAIN" if technical=="1" else "NOT_EVALUATED_A1_ATTACK_INHERITED_CHAIN",
 "scientific_status":"NOT_EVALUATED_TS_PMF_BARRIER_MECHANISM","next_action":"DO_NOT_START_PMF",
 "seed_index":int(index),"github_commit":commit,"stop_reason":reason,"window_count":len(windows),
 "inheritance_sha256_verified":verified=="1",
 "persistent_restarts":{"first_tetrahedral_like":record(root/"first_tetrahedral_like_candidate.rst7"),
  "shortest_attack_tetrahedral_like":record(root/"shortest_attack_tetrahedral_like_candidate.rst7"),
  "terminal_endpoint":record(root/"terminal_endpoint.rst7"),
  "guard_stop_endpoint":record(root/"guard_stop_endpoint.rst7")},
 "windows":[{"window":w["window"],"status":w["status"],"classification":w["classification"],
             "input_restart_sha256":json.load(open(root/f"window_{w['window']['window_index']:02d}/WINDOW_MANIFEST.json",encoding="utf-8"))["inherited_restart_sha256"],
             "output_restart_sha256":w["output_restart_sha256"]} for w in windows]}
(root/"CHAIN_MANIFEST.json").write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n",encoding="utf-8")
PY

CURRENT=hashing
find "$OUT" -type f \( -name 'WINDOW_MANIFEST.json' -o -name 'RESULT.json' -o -name 'PASS.json' -o -name 'NOT_EVALUATED.json' -o -name 'CHAIN_*.json' -o -name 'inheritance_verification.json' -o -name '*.rst7' \) -print0 \
    | sort -z | xargs -0r sha256sum >"$OUT/SHA256.tsv"

if [[ "$CHAIN_TECHNICAL" != 1 ]]; then
    STATE=NOT_EVALUATED
    CURRENT=terminal_history
    DETAIL="array_index=$INDEX;status=NOT_EVALUATED_A1_ATTACK_INHERITED_CHAIN;stop_reason=$STOP_REASON;output=$OUT/CHAIN_MANIFEST.json"
    append_history
    trap - EXIT
    exit 1
fi

STATE=PASS_TECHNICAL
CURRENT=terminal_history
DETAIL="array_index=$INDEX;status=PASS_TECHNICAL_A1_ATTACK_INHERITED_CHAIN;stop_reason=$STOP_REASON;scientific_status=NOT_EVALUATED_TS_PMF_BARRIER_MECHANISM;output=$OUT/CHAIN_MANIFEST.json"
append_history
trap - EXIT
printf 'PASS_TECHNICAL_A1_ATTACK_INHERITED_CHAIN seed_index=%s stop_reason=%s output=%s\n' "$INDEX" "$STOP_REASON" "$OUT"
