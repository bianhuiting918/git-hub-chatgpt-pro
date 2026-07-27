#!/usr/bin/env bash
# Amber18 A1 unified Step1/Step2 protein core DFTB3 numerical preflight; no MM or frame selection occurs here.
set -euo pipefail
TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
CODE_ROOT="${A1_UNIFIED_CORE_DFTB3_CODE_ROOT:?set immutable code snapshot}"
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
export GMXDATA=/public/software/apps/Gromacs-DCU2/2022.1/mpi/share/gromacs
AMBER_RUNTIME="$TASK_ROOT/../nylc_gyaq_pa66_l2_nac_qmmm_20260723/amber_runtime"
REP="$TASK_ROOT/a1_activated_nac_20260726/representative_frame/attempt_61990814"
ATTEMPT="${A1_UNIFIED_CORE_DFTB3_ATTEMPT:-${SLURM_JOB_ID:-manual_$(date -u '+%Y%m%dT%H%M%SZ')}}"
OUT="$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_unified_core_dftb3_preflight/attempt_$ATTEMPT"
GITHUB_COMMIT="${A1_UNIFIED_CORE_GITHUB_COMMIT:-unknown}"
[[ ! -e "$OUT" ]] || { printf 'refusing to overwrite %s\n' "$OUT" >&2; exit 2; }
mkdir -p "$OUT"
EVENT=nylc_a1_unified_core_dftb3_numerical_preflight
COMMAND=run_nylc_a1_unified_core_dftb3_preflight.sh
STATE=STARTED
CURRENT=initialization
DETAIL="source=$REP;scope=unified_protein_core_numerical_preflight_not_TS_PMF_barrier;output=$OUT"
append_history() {
    local now
    now="$(date '+%Y-%m-%dT%H:%M:%S%z')"
    (
        flock -x 9
        "$PY" - "$TASK_ROOT/run_history.tsv" "$TASK_ROOT/run_history.jsonl" \
            "$now" "$EVENT" "$COMMAND" "$STATE" "$GITHUB_COMMIT" "$DETAIL" "$OUT" <<'PY'
import json
import os
import sys
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
        tsv.seek(tsv_offset); tsv.truncate(); tsv.flush(); os.fsync(tsv.fileno())
        jsonl.seek(jsonl_offset); jsonl.truncate(); jsonl.flush(); os.fsync(jsonl.fileno())
        raise
PY
    ) 9>"$TASK_ROOT/.run_history.lock"
}
PROMOTED=0
demote_promoted_outputs() {
    if [[ -e "$OUT/PASS.json" ]]; then
        mv "$OUT/PASS.json" "$OUT/FAILED_NOT_PROMOTED_PASS.json"
    fi
    PROMOTED=0
}
finish() {
    local code=$?
    trap - EXIT
    if ((PROMOTED)); then
        demote_promoted_outputs
    fi
    STATE=NOT_EVALUATED
    DETAIL="stage=$CURRENT;exit_code=$code;scientific_status=NOT_EVALUATED;output=$OUT"
    "$PY" - "$OUT/NOT_EVALUATED.json" "$CURRENT" "$code" <<'PY'
import json,pathlib,sys
path,stage,code=sys.argv[1:]
pathlib.Path(path).write_text(json.dumps({"schema_version":1,"status":"NOT_EVALUATED_A1_UNIFIED_CORE_DFTB3_TECHNICAL_FAILURE","scientific_status":"NOT_EVALUATED","failed_stage":stage,"exit_code":int(code),"promoted":False},indent=2,sort_keys=True)+"\n",encoding="utf-8")
PY
    append_history || printf 'failed to append terminal history\n' >&2
    exit "$code"
}
trap finish EXIT
append_history
CURRENT=representative_gate
[[ -s "$REP/PASS.json" && -s "$REP/representative_354ps.gro" && -s "$REP/preflight.tpr" ]]
CURRENT=pbc_whole
set +u
source /work/home/acshdt1dks/opt/gromacs-fastest/env.sh
set -u
echo 0 | /public/software/apps/Gromacs-DCU2/2022.1/mpi/bin/gmx_mpi trjconv \
    -s "$REP/preflight.tpr" \
    -f "$REP/representative_354ps.gro" \
    -o "$OUT/representative_whole.gro" \
    -pbc mol -ur compact >"$OUT/trjconv.stdout" 2>"$OUT/trjconv.stderr"
test -s "$OUT/representative_whole.gro"
CURRENT=preparation
"$PY" "$CODE_ROOT/scripts/prepare_nylc_a1_unified_core_dftb3_preflight.py" \
    --coordinate "$OUT/representative_whole.gro" \
    --output "$OUT/prepared" >"$OUT/prepare.stdout" 2>"$OUT/prepare.stderr"
test -s "$OUT/prepared/qmmm_preflight_audit.json"
# Amber18 runtime and the 3ob-3-1 Slater-Koster set are mandatory.
module purge >/dev/null 2>&1 || true
module load amber/2018-hpcx-gcc-7.3.1
export AMBERHOME="$AMBER_RUNTIME"
test -s "$AMBERHOME/dat/slko/3ob-3-1/C-C.skf"
CURRENT=one_step
mpirun --bind-to none -np "${SLURM_NTASKS:-1}" sander.MPI -O -i "$OUT/prepared/01_qmmm_one_step.in" -o "$OUT/prepared/01_qmmm_one_step.out" -p "$OUT/prepared/system.prmtop" -c "$OUT/prepared/representative_354ps.rst7" -r "$OUT/prepared/01_qmmm_one_step.rst7" -inf "$OUT/prepared/01_qmmm_one_step.mdinfo"
CURRENT=twenty_step
mpirun --bind-to none -np "${SLURM_NTASKS:-1}" sander.MPI -O -i "$OUT/prepared/02_qmmm_20_step.in" -o "$OUT/prepared/02_qmmm_20_step.out" -p "$OUT/prepared/system.prmtop" -c "$OUT/prepared/01_qmmm_one_step.rst7" -r "$OUT/prepared/02_qmmm_20_step.rst7" -inf "$OUT/prepared/02_qmmm_20_step.mdinfo"
CURRENT=numerical_audit
set +e
"$PY" "$CODE_ROOT/scripts/audit_nylc_a1_unified_core_dftb3_smoke.py" --one-step "$OUT/prepared/01_qmmm_one_step.out" --segment "$OUT/prepared/02_qmmm_20_step.out" --audit "$OUT/prepared/qmmm_preflight_audit.json" --output "$OUT/A1_UNIFIED_CORE_DFTB3_PREFLIGHT_RESULT.json"
audit_rc=$?
set -e
[[ "$audit_rc" -eq 0 ]]
"$PY" - "$OUT/A1_UNIFIED_CORE_DFTB3_PREFLIGHT_RESULT.json" <<'PY'
import json,sys
if json.load(open(sys.argv[1],encoding="utf-8")).get("status")!="PASS_A1_UNIFIED_CORE_DFTB3_NUMERICAL_PREFLIGHT":
    raise SystemExit("A1 DFTB3 audit did not PASS")
PY
CURRENT=hashing
sha256sum "$REP/PASS.json" "$REP/representative_354ps.gro" "$OUT/representative_whole.gro" "$OUT/prepared/qmmm_preflight_audit.json" "$OUT/prepared/system.prmtop" "$OUT/prepared/representative_354ps.rst7" "$OUT/prepared/01_qmmm_one_step.in" "$OUT/prepared/01_qmmm_one_step.out" "$OUT/prepared/02_qmmm_20_step.in" "$OUT/prepared/02_qmmm_20_step.out" "$OUT/A1_UNIFIED_CORE_DFTB3_PREFLIGHT_RESULT.json" >"$OUT/SHA256.tsv"
PROMOTED=1
cp "$OUT/A1_UNIFIED_CORE_DFTB3_PREFLIGHT_RESULT.json" "$OUT/PASS.json"
STATE=PASS_TECHNICAL
CURRENT=terminal_history
DETAIL="status=PASS_A1_UNIFIED_CORE_DFTB3_NUMERICAL_PREFLIGHT;scientific_status=NOT_EVALUATED_TS_PMF_BARRIER;output=$OUT/PASS.json"
append_history
trap - EXIT
