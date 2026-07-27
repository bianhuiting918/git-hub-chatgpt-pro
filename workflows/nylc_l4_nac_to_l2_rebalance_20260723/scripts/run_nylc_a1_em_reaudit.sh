#!/usr/bin/env bash
set -euo pipefail
TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
FLOW="$TASK_ROOT/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
SOURCE="${1:?usage: run_nylc_a1_em_reaudit.sh SOURCE_EM_CANDIDATE_DIR [CODE_ROOT]}"
CODE_ROOT="${2:-${A1_SCRIPT_ROOT:-$FLOW/scripts}}"
GITHUB_COMMIT="${A1_GITHUB_COMMIT:-unknown}"
[[ -d "$SOURCE/em_free" && -d "$SOURCE/input" ]] || { printf 'invalid EM candidate directory: %s\n' "$SOURCE" >&2; exit 2; }
for name in audit_nylc_a1_em.py audit_nylc_a1_full_system.py build_nylc_a1_full_system.py; do
    [[ -f "$CODE_ROOT/$name" ]] || { printf 'missing frozen audit source: %s\n' "$CODE_ROOT/$name" >&2; exit 2; }
done
CANDIDATE="$(basename "$SOURCE")"
ATTEMPT="${SLURM_ARRAY_JOB_ID:-manual}_${SLURM_ARRAY_TASK_ID:-0}_${SLURM_JOB_ID:-manual}"
OUT="$TASK_ROOT/a1_activated_nac_20260726/em_reaudit/attempt_$ATTEMPT/$CANDIDATE"
[[ ! -e "$OUT" ]] || { printf 'refusing to overwrite %s\n' "$OUT" >&2; exit 2; }
mkdir -p "$OUT"
EVENT=nylc_a1_em_reaudit
COMMAND="run_nylc_a1_em_reaudit.sh $SOURCE $CODE_ROOT"
STATE=FAIL_TECHNICAL
DETAIL="candidate=$CANDIDATE; initializing"
append_history(){
    local now
    now="$(date '+%Y-%m-%dT%H:%M:%S%z')"
    (
        flock -x 9
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$now" "$EVENT" "$COMMAND" "$STATE" "$GITHUB_COMMIT" "$DETAIL" >>"$TASK_ROOT/run_history.tsv"
        "$PY" - "$now" "$EVENT" "$COMMAND" "$STATE" "$GITHUB_COMMIT" "$DETAIL" "$OUT" >>"$TASK_ROOT/run_history.jsonl" <<'PY'
import json,sys
now,event,command,state,commit,detail,output=sys.argv[1:]
print(json.dumps({"time":now,"event":event,"command":command,"state":state,"git_commit":commit,"detail":detail,"output":output},sort_keys=True))
PY
    ) 9>"$TASK_ROOT/.run_history.lock"
}
finish(){
    code=$?
    trap - EXIT
    if [[ $code -ne 0 ]]; then DETAIL="candidate=$CANDIDATE; exit_code=$code; output=$OUT"; fi
    append_history
    exit "$code"
}
trap finish EXIT
"$PY" "$CODE_ROOT/audit_nylc_a1_em.py" \
    --gro "$SOURCE/em_free/run.gro" \
    --chain-itp "$SOURCE/input/topol_Protein_chain_H.itp" \
    --build-audit "$SOURCE/input/A1_FULL_SYSTEM_BUILD.json" \
    --log "$SOURCE/em_free/run.log" \
    --output "$OUT/A1_EM_REAUDIT.json"
sha256sum \
    "$SOURCE/em_free/run.gro" "$SOURCE/em_free/run.log" "$SOURCE/em_free/run.edr" \
    "$CODE_ROOT/audit_nylc_a1_em.py" "$CODE_ROOT/audit_nylc_a1_full_system.py" \
    "$CODE_ROOT/build_nylc_a1_full_system.py" "$OUT/A1_EM_REAUDIT.json" >"$OUT/sha256.tsv"
STATE=PASS_TECHNICAL
DETAIL="candidate=$CANDIDATE; status=PASS_A1_EM; source=$SOURCE; output=$OUT"
