#!/usr/bin/env bash
set -euo pipefail
TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
FLOW="$TASK_ROOT/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
CODE_ROOT="${A1_NAC_CODE_ROOT:-$FLOW}"
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
GMX=/public/software/apps/gromacs/2022.2/hpcx-gcc7.3.1/bin/gmx_mpi
SLOT="${1:-${SLURM_ARRAY_TASK_ID:?usage: run_nylc_a1_nac_audit.sh SLOT}}"
GITHUB_COMMIT="${A1_GITHUB_COMMIT:-unknown}"
[[ "$SLOT" =~ ^[0-8]$ ]] || exit 2
UNIVERSE="$CODE_ROOT/manifests/nylc_a1_nac_audit_universe.json"
mapfile -t VALUES < <("$PY" - "$UNIVERSE" "$SLOT" <<'PY'
import json,sys
p=json.load(open(sys.argv[1])); row=p["replicas"][int(sys.argv[2])]
assert row["slot"]==int(sys.argv[2])
for key in ("candidate_id","velocity_seed","free_run_root","source_cycle_ndx","completion_manifest"):
    print(row[key])
print(p["selection_manifest"])
PY
)
CANDIDATE="${VALUES[0]}"; SEED="${VALUES[1]}"; FREE="${VALUES[2]}"
SOURCE_NDX="${VALUES[3]}"; COMPLETE="${VALUES[4]}"; SELECTION="${VALUES[5]}"
ATTEMPT="${SLURM_ARRAY_JOB_ID:-manual}_${SLURM_ARRAY_TASK_ID:-$SLOT}_${SLURM_JOB_ID:-manual}"
OUT="$TASK_ROOT/a1_activated_nac_20260726/nac_audit/attempt_$ATTEMPT/$CANDIDATE/seed$SEED"
[[ ! -e "$OUT" ]] || { printf 'refusing to overwrite %s\n' "$OUT" >&2; exit 2; }
mkdir -p "$OUT"
EVENT=nylc_a1_nac_audit
COMMAND="run_nylc_a1_nac_audit.sh $SLOT"
STATE=FAIL_TECHNICAL
DETAIL="slot=$SLOT;candidate=$CANDIDATE;seed=$SEED;initializing"
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
CURRENT=validation
finish(){
    code=$?
    trap - EXIT
    if [[ $code -ne 0 ]]; then
        STATE=FAIL_TECHNICAL
        DETAIL="slot=$SLOT;candidate=$CANDIDATE;seed=$SEED;stage=$CURRENT;exit_code=$code;output=$OUT"
        if [[ ! -s "$OUT/replica_audit.json" ]]; then
            "$PY" - "$OUT/NOT_EVALUATED.json" "$SLOT" "$CANDIDATE" "$SEED" "$CURRENT" "$code" <<'PY'
import json,sys
path,slot,candidate,seed,stage,code=sys.argv[1:]
open(path,"w").write(json.dumps({"schema_version":1,"slot":int(slot),"candidate_id":candidate,"velocity_seed":int(seed),"technical_status":"FAIL","scientific_status":"NOT_EVALUATED_AUDIT_TECHNICAL_FAIL","failed_stage":stage,"exit_code":int(code)},indent=2,sort_keys=True)+"\n")
PY
        fi
    fi
    append_history
    exit "$code"
}
trap finish EXIT
"$PY" - "$COMPLETE" <<'PY'
import json,sys
p=json.load(open(sys.argv[1]))
if p.get("status")!="PASS_TECHNICAL_A1_EQUILIBRATION" or not p.get("fully_unrestrained"):
    raise SystemExit("parent equilibration is not a fully unrestrained technical PASS")
PY
for path in "$FREE/run.tpr" "$FREE/run.xtc" "$FREE/run.edr" "$FREE/run.log" "$SOURCE_NDX" "$SELECTION"; do
    [[ -s "$path" ]] || { printf 'missing audit input %s\n' "$path" >&2; exit 2; }
done
module purge
module load gromacs/2022.2-hpcx-gcc-7.3.1
export GMX_MAXBACKUP=-1 OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
STATE=STARTED
DETAIL="slot=$SLOT;candidate=$CANDIDATE;seed=$SEED;free_root=$FREE;scientific_status=NOT_EVALUATED"
append_history
CURRENT=primitive_generation
"$PY" "$CODE_ROOT/scripts/generate_nylc_m1_ensemble_primitives.py" \
    --tpr "$FREE/run.tpr" --xtc "$FREE/run.xtc" \
    --source-cycle-ndx "$SOURCE_NDX" --selection-manifest "$SELECTION" \
    --candidate-id "$CANDIDATE" --velocity-seed "$SEED" \
    --output-dir "$OUT/primitives" >"$OUT/generate.stdout" 2>"$OUT/generate.stderr"
mv "$OUT/primitives/replica_manifest.json" "$OUT/primitives/generator_manifest_m1_legacy.json"
"$PY" - "$UNIVERSE" "$SLOT" "$OUT/primitives/a1_replica_manifest.json" <<'PY'
import json,sys
u=json.load(open(sys.argv[1])); row=u["replicas"][int(sys.argv[2])]
p={"schema_version":1,"slot":row["slot"],"candidate_id":row["candidate_id"],"velocity_seed":row["velocity_seed"],"microstate":u["microstate"],"gate_definition":"NylC residues 261-266; Thr267 excluded","fully_unrestrained":True,"nac_contract":u["nac"],"analysis_window_ps":u["analysis_window_ps"],"source_free_run_root":row["free_run_root"],"proton_path_evidence_included":False}
open(sys.argv[3],"w").write(json.dumps(p,indent=2,sort_keys=True)+"\n")
PY
ln -s "$FREE/run.log" "$OUT/primitives/run.log"
CURRENT=energy_extraction
# gmx energy extracts Potential, Temperature, Pressure, and Volume from the same free-window EDR.
printf 'Potential\n0\n' | "$GMX" energy -f "$FREE/run.edr" -o "$OUT/primitives/potential_energy.xvg" >"$OUT/potential.stdout" 2>"$OUT/potential.stderr"
printf 'Temperature\nPressure\nVolume\n0\n' | "$GMX" energy -f "$FREE/run.edr" -o "$OUT/primitives/thermo.xvg" >"$OUT/thermo.stdout" 2>"$OUT/thermo.stderr"
CURRENT=replica_audit
"$PY" "$CODE_ROOT/scripts/audit_nylc_a1_nac_replica.py" \
    --run-root "$OUT/primitives" --manifest "$OUT/primitives/a1_replica_manifest.json" \
    --window-start-ps 0 --window-end-ps 1000 --output "$OUT/replica_audit.json" \
    >"$OUT/audit.stdout" 2>"$OUT/audit.stderr"
CURRENT=hashing
sha256sum "$FREE/run.tpr" "$FREE/run.xtc" "$FREE/run.edr" "$FREE/run.log" \
    "$OUT/primitives/nac_distance.xvg" "$OUT/primitives/nac_angle.xvg" \
    "$OUT/primitives/gate_opening.xvg" "$OUT/replica_audit.json" >"$OUT/sha256.tsv"
STATE=PASS_TECHNICAL
DETAIL="$("$PY" - "$OUT/replica_audit.json" <<'PY'
import json,sys
p=json.load(open(sys.argv[1]))
print(f"slot={p.get('slot','NA')};candidate={p['candidate_id']};seed={p['velocity_seed']};scientific_status={p['scientific_status']};nac_occupancy={p['nac']['nac_occupancy']}")
PY
)"
trap - EXIT
append_history
