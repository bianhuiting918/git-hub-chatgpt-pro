#!/bin/bash
#SBATCH -J nylc_a1_nac_merge
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem-per-cpu=2000M
#SBATCH -t 00:30:00
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-%x-%j.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-%x-%j.err
# Submit with --dependency=afterany:<AUDIT_ARRAY_JOB_ID>.
set -euo pipefail
TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
FLOW="$TASK_ROOT/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
AUDIT_ARRAY_JOB_ID="${AUDIT_ARRAY_JOB_ID:?submit with afterany dependency and AUDIT_ARRAY_JOB_ID}"
GITHUB_COMMIT="${A1_GITHUB_COMMIT:-unknown}"
JOB_ID="${SLURM_JOB_ID:?SLURM_JOB_ID required}"
SNAP="$TASK_ROOT/a1_activated_nac_20260726/code_snapshots/nac_merge_$JOB_ID"
MERGE_ROOT="$TASK_ROOT/a1_activated_nac_20260726/nac_audit_merge"
OUT="$MERGE_ROOT/job_$JOB_ID"
[[ ! -e "$OUT" ]] || { printf 'refusing to overwrite %s\n' "$OUT" >&2; exit 2; }
mkdir -p "$MERGE_ROOT"
mkdir -p "$SNAP/scripts" "$SNAP/manifests"
cp "$FLOW/scripts/merge_nylc_a1_nac_audits.py" "$SNAP/scripts/"
cp "$FLOW/manifests/nylc_a1_nac_audit_universe.json" "$SNAP/manifests/"
sha256sum "$SNAP"/scripts/* "$SNAP"/manifests/* >"$SNAP/SNAPSHOT_SHA256.tsv"
EVENT=nylc_a1_nac_merge
COMMAND="merge_nylc_a1_nac_audits.py --audit-array-job-id $AUDIT_ARRAY_JOB_ID"
STATE=FAIL_TECHNICAL
DETAIL="audit_array=$AUDIT_ARRAY_JOB_ID;initializing"
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
finish(){ code=$?; trap - EXIT; if [[ $code -ne 0 ]]; then DETAIL="audit_array=$AUDIT_ARRAY_JOB_ID;exit_code=$code;output=$OUT"; fi; append_history; exit "$code"; }
trap finish EXIT
"$PY" "$SNAP/scripts/merge_nylc_a1_nac_audits.py" \
    --universe "$SNAP/manifests/nylc_a1_nac_audit_universe.json" \
    --audit-root "$TASK_ROOT/a1_activated_nac_20260726/nac_audit" \
    --audit-array-job-id "$AUDIT_ARRAY_JOB_ID" --output-dir "$OUT" \
    >"$OUT.stdout" 2>"$OUT.stderr"
sha256sum "$OUT"/* >"$OUT/sha256.tsv"
read -r SCI ELIG MISS < <("$PY" - "$OUT/A1_NAC_AUDIT_SUMMARY.json" <<'PY'
import json,sys
p=json.load(open(sys.argv[1]))
print(p["scientific_status"],p["eligible_replica_count"],len(p["missing_slots"]))
PY
)
STATE=PASS_TECHNICAL
DETAIL="audit_array=$AUDIT_ARRAY_JOB_ID;scientific_status=$SCI;eligible_replica_count=$ELIG;missing_slots=$MISS;output=$OUT"
trap - EXIT
append_history
