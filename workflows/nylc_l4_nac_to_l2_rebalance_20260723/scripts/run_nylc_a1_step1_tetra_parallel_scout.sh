#!/usr/bin/env bash
set -euo pipefail

TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
CODE_ROOT="${A1_TETRA_SCOUT_CODE_ROOT:?set immutable code root}"
GITHUB_COMMIT="${A1_TETRA_SCOUT_GITHUB_COMMIT:?set immutable GitHub commit}"
INDEX="${SLURM_ARRAY_TASK_ID:?run as array task 0..5}"
ARRAY_JOB="${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-manual}}"
ATTEMPT="${ARRAY_JOB}_${INDEX}"
OUT="$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_step1_tetra_parallel_scout/attempt_$ATTEMPT"
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
DRIVER="$CODE_ROOT/scripts/prepare_audit_nylc_a1_step1_tetra_parallel_scout.py"
PRMTOP="$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_unified_core_dftb3_preflight/attempt_62011285/prepared/system.prmtop"
SCRATCH_ROOT="${SLURM_TMPDIR:-/tmp}/nylc_a1_tetra_parallel_$ATTEMPT"
START_RST7="$SCRATCH_ROOT/start.rst7"

EVENT=nylc_a1_tetra_parallel_scout
COMMAND=run_nylc_a1_step1_tetra_parallel_scout.sh
STATE=STARTED
CURRENT=initialization
DETAIL="array_index=$INDEX;scope=bounded_tetrahedralization_not_path_barrier;output=$OUT"
OWNED=0
TERMINAL_RST7=""

append_history() {
    local now
    now="$(date '+%Y-%m-%dT%H:%M:%S%z')"
    (
        flock -x 9
        "$PY" - "$TASK_ROOT/run_history.tsv" "$TASK_ROOT/run_history.jsonl" "$now"             "$EVENT" "$COMMAND" "$STATE" "$GITHUB_COMMIT" "$DETAIL" "$OUT" <<'PY'
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
    ) 9>>"$TASK_ROOT/.run_history.lock"
}

finish() {
    local code=$?
    trap - EXIT
    if (( OWNED )); then
        args=(--mode finalize --task-index "$INDEX" --root "$OUT" --stop-reason "TECHNICAL_EXIT_$CURRENT")
        [[ -n "$TERMINAL_RST7" ]] && args+=(--terminal-rst7 "$TERMINAL_RST7")
        "$PY" "$DRIVER" "${args[@]}" || true
    fi
    STATE=NOT_EVALUATED
    DETAIL="array_index=$INDEX;stage=$CURRENT;exit_code=$code;scientific_status=NOT_EVALUATED_TS_PMF_BARRIER_MECHANISM;output=$OUT"
    append_history || true
    exit "$code"
}
trap finish EXIT
append_history

test ! -e "$OUT"
mkdir -p "$SCRATCH_ROOT"
CURRENT=initialize
"$PY" "$DRIVER" --mode initialize --task-index "$INDEX" --root "$OUT"     --start-rst7 "$START_RST7" --github-commit "$GITHUB_COMMIT"
OWNED=1

module purge >/dev/null 2>&1 || true
module load amber/2018-hpcx-gcc-7.3.1
export AMBERHOME="$TASK_ROOT/../nylc_gyaq_pa66_l2_nac_qmmm_20260723/amber_runtime"
export OMP_NUM_THREADS=1
test -s "$AMBERHOME/dat/slko/3ob-3-1/C-C.skf"

CURRENT_RST7="$START_RST7"
STOP_REASON=COMPLETED_ALL_WINDOWS
if (( INDEX < 2 )); then
    WINDOWS="0 1 2 3"
else
    WINDOWS="0 1"
fi

for WINDOW_INDEX in $WINDOWS; do
    WINDOW_OUT="$OUT/window_$(printf '%02d' "$WINDOW_INDEX")"
    WINDOW_SCRATCH="$SCRATCH_ROOT/window_$(printf '%02d' "$WINDOW_INDEX")"
    NEXT_RST7="$WINDOW_SCRATCH/stage.rst7"
    mkdir -p "$WINDOW_SCRATCH"
    CURRENT="prepare_window_$WINDOW_INDEX"
    "$PY" "$DRIVER" --mode prepare --task-index "$INDEX" --window-index "$WINDOW_INDEX"         --root "$OUT" --output "$WINDOW_OUT" --input-rst7 "$CURRENT_RST7"
    cp "$WINDOW_OUT/restraints.RST" "$WINDOW_SCRATCH/restraints.RST"
    CURRENT="run_window_$WINDOW_INDEX"
    (
        cd "$WINDOW_SCRATCH"
        mpirun --bind-to none -np "${SLURM_NTASKS:-8}" sander.MPI -O             -i "$WINDOW_OUT/stage.in" -o stage.out -p "$PRMTOP"             -c "$CURRENT_RST7" -ref "$CURRENT_RST7" -r stage.rst7 -inf stage.mdinfo
    )
    test -s "$NEXT_RST7"
    CURRENT="audit_window_$WINDOW_INDEX"
    "$PY" "$DRIVER" --mode audit --task-index "$INDEX" --root "$OUT"         --output "$WINDOW_OUT" --scratch "$WINDOW_SCRATCH"
    TERMINAL_RST7="$NEXT_RST7"
    read -r TECHNICAL CANDIDATE GUARD_STOP < <("$PY" - "$WINDOW_OUT/RESULT.json" <<'PY'
import json,sys
r=json.load(open(sys.argv[1],encoding="utf-8"))
print(int(r["technical_pass"]),int(r["gates"]["TETRAHEDRAL_LIKE_RESTRAINED"]),
      int(r["stopped_by_reactant_guard"]))
PY
)
    if [[ "$TECHNICAL" != 1 ]]; then
        STOP_REASON=TECHNICAL_NOT_EVALUATED
        break
    fi
    if [[ "$GUARD_STOP" == 1 ]]; then
        STOP_REASON=OUT_OF_SCOPE_COUPLED_PT_OR_CN_RESPONSE
        break
    fi
    if [[ "$CANDIDATE" == 1 ]]; then
        RELEASE_OUT="$OUT/full_release"
        RELEASE_SCRATCH="$SCRATCH_ROOT/full_release"
        mkdir -p "$RELEASE_SCRATCH"
        CURRENT=prepare_full_release
        "$PY" "$DRIVER" --mode prepare-release --task-index "$INDEX" --root "$OUT"             --output "$RELEASE_OUT" --input-rst7 "$NEXT_RST7"
        CURRENT=run_full_release
        (
            cd "$RELEASE_SCRATCH"
            mpirun --bind-to none -np "${SLURM_NTASKS:-8}" sander.MPI -O                 -i "$RELEASE_OUT/stage.in" -o stage.out -p "$PRMTOP"                 -c "$NEXT_RST7" -r stage.rst7 -inf stage.mdinfo
        )
        test -s "$RELEASE_SCRATCH/stage.rst7"
        CURRENT=audit_full_release
        "$PY" "$DRIVER" --mode audit-release --task-index "$INDEX" --root "$OUT"             --output "$RELEASE_OUT" --scratch "$RELEASE_SCRATCH"
        TERMINAL_RST7="$RELEASE_SCRATCH/stage.rst7"
        STOP_REASON=FIRST_RESTRAINED_CANDIDATE_FULLY_RELEASED
        break
    fi
    CURRENT_RST7="$NEXT_RST7"
done

CURRENT=finalize
final_args=(--mode finalize --task-index "$INDEX" --root "$OUT" --stop-reason "$STOP_REASON")
[[ -n "$TERMINAL_RST7" ]] && final_args+=(--terminal-rst7 "$TERMINAL_RST7")
"$PY" "$DRIVER" "${final_args[@]}"
test -s "$OUT/CHAIN_MANIFEST.json"
test -s "$OUT/SHA256.tsv"
TECHNICAL="$("$PY" - "$OUT/CHAIN_MANIFEST.json" <<'PY'
import json,sys
print(int(json.load(open(sys.argv[1],encoding="utf-8"))["status"]=="PASS_TECHNICAL_A1_TETRA_SCOUT"))
PY
)"
if [[ "$TECHNICAL" != 1 ]]; then
    STATE=NOT_EVALUATED
    CURRENT=terminal_history
    DETAIL="array_index=$INDEX;status=NOT_EVALUATED_TECHNICAL_A1_TETRA_SCOUT;stop_reason=$STOP_REASON;output=$OUT/CHAIN_MANIFEST.json"
    append_history
    trap - EXIT
    exit 1
fi

STATE=PASS_TECHNICAL
CURRENT=terminal_history
DETAIL="array_index=$INDEX;status=PASS_TECHNICAL_A1_TETRA_SCOUT;stop_reason=$STOP_REASON;scientific_status=NOT_EVALUATED_TS_PMF_BARRIER_MECHANISM;output=$OUT/CHAIN_MANIFEST.json"
append_history
trap - EXIT
printf 'PASS_TECHNICAL_A1_TETRA_SCOUT task=%s stop_reason=%s output=%s\n' "$INDEX" "$STOP_REASON" "$OUT"
