#!/usr/bin/env bash
set -euo pipefail
TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
FLOW="$TASK_ROOT/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
CODE_ROOT="${A1_EQ_CODE_ROOT:-$FLOW}"
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
GMX=/public/software/apps/Gromacs-DCU2/2022.1/mpi/bin/gmx_mpi
TASK_ID="${1:-${SLURM_ARRAY_TASK_ID:?usage: run_nylc_a1_equilibration.sh TASK_ID}}"
GITHUB_COMMIT="${A1_GITHUB_COMMIT:-unknown}"
[[ "$TASK_ID" =~ ^[0-8]$ ]] || exit 2
candidate_index=$((TASK_ID / 3))
seed_index=$((TASK_ID % 3))
CANDIDATES=(nac_evt18_time1206ps nac_evt25_time1462ps nac_evt08_time1086ps)
SEEDS=(26711 26723 26737)
EM_ROOTS=(
"$TASK_ROOT/a1_activated_nac_20260726/em/attempt_61968026_0_61968027/nac_evt18_time1206ps"
"$TASK_ROOT/a1_activated_nac_20260726/em/attempt_61968026_1_61968026/nac_evt25_time1462ps"
"$TASK_ROOT/a1_activated_nac_20260726/em/attempt_61969851_2_61969851/nac_evt08_time1086ps"
)
AUDITS=(
"$TASK_ROOT/a1_activated_nac_20260726/em_reaudit/attempt_manual_0_61969849/nac_evt18_time1206ps/A1_EM_REAUDIT.json"
"$TASK_ROOT/a1_activated_nac_20260726/em_reaudit/attempt_manual_0_61969850/nac_evt25_time1462ps/A1_EM_REAUDIT.json"
"$TASK_ROOT/a1_activated_nac_20260726/em/attempt_61969851_2_61969851/nac_evt08_time1086ps/em_free/A1_EM_AUDIT.json"
)
CANDIDATE="${CANDIDATES[$candidate_index]}"
SEED="${SEEDS[$seed_index]}"
EM_ROOT="${EM_ROOTS[$candidate_index]}"
EM_AUDIT="${AUDITS[$candidate_index]}"
ATTEMPT="${SLURM_ARRAY_JOB_ID:-manual}_${SLURM_ARRAY_TASK_ID:-$TASK_ID}_${SLURM_JOB_ID:-manual}"
OUT="$TASK_ROOT/a1_activated_nac_20260726/equilibration/$CANDIDATE/seed$SEED/attempt_$ATTEMPT"
[[ ! -e "$OUT" ]] || { printf 'refusing to overwrite %s\n' "$OUT" >&2; exit 2; }
mkdir -p "$OUT"
EVENT=nylc_a1_equilibration
COMMAND="run_nylc_a1_equilibration.sh $TASK_ID"
STATE=FAIL_TECHNICAL
DETAIL="candidate=$CANDIDATE;seed=$SEED;initializing"
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
    if [[ $code -ne 0 ]]; then DETAIL="candidate=$CANDIDATE;seed=$SEED;stage=${CURRENT:-initialization};exit_code=$code;output=$OUT"; fi
    append_history
    exit "$code"
}
trap finish EXIT
if [[ ! -s "$EM_AUDIT" ]]; then
    STATE=NOT_EVALUATED_EM_PENDING
    DETAIL="candidate=$CANDIDATE;seed=$SEED;missing_em_audit=$EM_AUDIT"
    "$PY" - "$OUT/NOT_EVALUATED.json" "$DETAIL" <<'PY'
import json,sys
open(sys.argv[1],"w").write(json.dumps({"status":"NOT_EVALUATED_EM_PENDING","reason":sys.argv[2]},indent=2)+"\n")
PY
    exit 0
fi
"$PY" - "$EM_AUDIT" <<'PY'
import json,sys
if json.load(open(sys.argv[1])).get("status")!="PASS_A1_EM":
    raise SystemExit("input is not PASS_A1_EM")
PY
for path in "$EM_ROOT/em_free/run.gro" "$EM_ROOT/input/topol.top"; do [[ -s "$path" ]] || { printf 'missing input %s\n' "$path" >&2; exit 2; }; done
set +u
source /work/home/acshdt1dks/opt/gromacs-fastest/env.sh
set -u
export GMX_MAXBACKUP=-1 OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
"$PY" - "$CODE_ROOT/mdp/nvt50_m1.mdp" "$OUT/nvt50.rendered.mdp" "$SEED" <<'PY'
import re,sys
text=open(sys.argv[1]).read()
text=re.sub(r"(?m)^gen-seed\s*=.*$",f"gen-seed                 = {sys.argv[3]}",text)
open(sys.argv[2],"w").write(text)
PY
scan_log(){
    "$PY" - "$1" "$2" <<'PY'
import json,re,sys
text=open(sys.argv[1],errors="replace").read()
audit={
 "fatal":len(re.findall(r"(?i)fatal(?:\s+error)?",text)),
 "lincs_warning":len(re.findall(r"(?i)lincs\s+warning",text)),
 "settle_problem":len(re.findall(r"(?i)(?:cannot\s+be\s+settled|settle[^\n]*(?:warning|error|failed))",text)),
 "nan":len(re.findall(r"(?i)\bnan\b",text)),
 "finished_mdrun":bool(re.search(r"Finished mdrun",text)),
}
open(sys.argv[2],"w").write(json.dumps(audit,indent=2,sort_keys=True)+"\n")
if not audit["finished_mdrun"] or any(audit[k] for k in ("fatal","lincs_warning","settle_problem","nan")): raise SystemExit(1)
PY
}
run_stage(){
    local name="$1" mdp="$2" coord="$3" checkpoint="$4" restrained="$5"
    local stage="$OUT/$name"
    mkdir -p "$stage"
    local args=(grompp -f "$mdp" -c "$coord" -p "$EM_ROOT/input/topol.top" -o "$stage/run.tpr" -po "$stage/run.processed.mdp" -maxwarn 0)
    [[ "$restrained" == yes ]] && args+=(-r "$EM_ROOT/em_free/run.gro")
    [[ -n "$checkpoint" ]] && args+=(-t "$checkpoint")
    "$GMX" "${args[@]}" >"$stage/grompp.stdout" 2>"$stage/grompp.stderr"
    mpirun -np 1 "$GMX" mdrun -s "$stage/run.tpr" -deffnm "$stage/run" -ntomp "${SLURM_CPUS_PER_TASK:-8}" -nb gpu -pme gpu -bonded gpu -pin on >"$stage/mdrun.stdout" 2>"$stage/mdrun.stderr"
    test -s "$stage/run.gro"; test -s "$stage/run.cpt"
    scan_log "$stage/run.log" "$stage/numerical_audit.json"
}
STATE=STARTED
DETAIL="candidate=$CANDIDATE;seed=$SEED;em_audit=$EM_AUDIT;scientific_status=NOT_EVALUATED"
append_history
CURRENT=nvt50
run_stage nvt50 "$OUT/nvt50.rendered.mdp" "$EM_ROOT/em_free/run.gro" "" yes
CURRENT=nvt150
run_stage nvt150 "$CODE_ROOT/mdp/nvt150_m1.mdp" "$OUT/nvt50/run.gro" "$OUT/nvt50/run.cpt" yes
CURRENT=nvt300
run_stage nvt300 "$CODE_ROOT/mdp/nvt300_m1.mdp" "$OUT/nvt150/run.gro" "$OUT/nvt150/run.cpt" yes
CURRENT=npt300r
run_stage npt300r "$CODE_ROOT/mdp/npt300r_m1.mdp" "$OUT/nvt300/run.gro" "$OUT/nvt300/run.cpt" yes
CURRENT=npt300rel
run_stage npt300rel "$CODE_ROOT/mdp/npt300rel_m1.mdp" "$OUT/npt300r/run.gro" "$OUT/npt300r/run.cpt" yes
CURRENT=npt300free
FREE="$OUT/npt300free"
mkdir -p "$FREE"
"$GMX" grompp -f "$CODE_ROOT/mdp/npt300free_m1.mdp" -c "$OUT/npt300rel/run.gro" -t "$OUT/npt300rel/run.cpt" -p "$EM_ROOT/input/topol.top" -o "$FREE/run.tpr" -po "$FREE/run.processed.mdp" -maxwarn 0 >"$FREE/grompp.stdout" 2>"$FREE/grompp.stderr"
"$GMX" dump -s "$FREE/run.tpr" >"$FREE/tpr.dump" 2>"$FREE/tpr_dump.stderr"
"$PY" - "$FREE/tpr.dump" "$FREE/run.processed.mdp" "$FREE/free_tpr_contract.json" <<'PY'
import json,re,sys
dump=open(sys.argv[1],errors="replace").read(); mdp=open(sys.argv[2],errors="replace").read()
pos=sum(map(int,re.findall(r"#posres_xA\s*=\s*(\d+)",dump)))
dis=len(re.findall(r"(?i)\bDISRES\b",dump))
defines=[]
for line in mdp.splitlines():
    body=line.split(";",1)[0].strip()
    if "=" in body and body.split("=",1)[0].strip().lower()=="define" and body.split("=",1)[1].strip():
        defines.append(body.split("=",1)[1].strip())
payload={"position_restraints":pos,"distance_restraints":dis,"nonempty_defines":defines,"scientific_window":"fully_unrestrained_NPT_1ns"}
open(sys.argv[3],"w").write(json.dumps(payload,indent=2,sort_keys=True)+"\n")
if pos or dis or defines: raise SystemExit("free TPR contains restraints")
PY
mpirun -np 1 "$GMX" mdrun -s "$FREE/run.tpr" -deffnm "$FREE/run" -ntomp "${SLURM_CPUS_PER_TASK:-8}" -nb gpu -pme gpu -bonded gpu -pin on >"$FREE/mdrun.stdout" 2>"$FREE/mdrun.stderr"
test -s "$FREE/run.gro"; test -s "$FREE/run.cpt"
scan_log "$FREE/run.log" "$FREE/numerical_audit.json"
"$PY" - "$OUT/EQUILIBRATION_COMPLETE.json" "$CANDIDATE" "$SEED" "$EM_AUDIT" <<'PY'
import hashlib,json,pathlib,sys
out=pathlib.Path(sys.argv[1]); free=out.parent/"npt300free"
payload={"schema_version":1,"status":"PASS_TECHNICAL_A1_EQUILIBRATION","scientific_status":"NOT_EVALUATED_PENDING_NAC_AUDIT","candidate_id":sys.argv[2],"velocity_seed":int(sys.argv[3]),"parent_em_audit":sys.argv[4],"fully_unrestrained":True,"free_window_ps":[0.0,1000.0],"final_gro_sha256":hashlib.sha256((free/"run.gro").read_bytes()).hexdigest()}
out.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
PY
STATE=PASS_TECHNICAL
DETAIL="candidate=$CANDIDATE;seed=$SEED;fully_unrestrained_1ns_complete=yes;scientific_status=NOT_EVALUATED_PENDING_NAC_AUDIT"
trap - EXIT
append_history
