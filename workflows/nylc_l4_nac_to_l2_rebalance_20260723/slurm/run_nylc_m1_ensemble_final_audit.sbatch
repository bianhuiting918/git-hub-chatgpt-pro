#!/bin/bash
#SBATCH -J nylc_m1_final
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH -t 04:00:00
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/final_audit_%j.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/final_audit_%j.err

set -euo pipefail

TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
FLOW="$TASK_ROOT/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
SELECTION="$TASK_ROOT/ensemble/selection_job_61813011/ensemble_selection.json"
PYTHON=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
GMX=/public/software/apps/Gromacs-DCU2/2022.1/mpi/bin/gmx_mpi
STAGEA_RANK_JOB_ID="${STAGEA_RANK_JOB_ID:?STAGEA_RANK_JOB_ID is required}"
STAGEB_ARRAY_JOB_ID="${STAGEB_ARRAY_JOB_ID:?STAGEB_ARRAY_JOB_ID is required}"
STAGEB_COUNT="${STAGEB_COUNT:?STAGEB_COUNT is required}"
JOB_ID="${SLURM_JOB_ID:?SLURM_JOB_ID is required}"
TOP6="$TASK_ROOT/ensemble/stageA_rank_job_${STAGEA_RANK_JOB_ID}/results/stageA_top6_manifest.json"
STAGEA_AUDITS="$TASK_ROOT/ensemble/stageA_rank_job_${STAGEA_RANK_JOB_ID}/results/stageA_replica_audits.json"
OUT="$TASK_ROOT/ensemble/final_audit_job_${JOB_ID}"
EVENT=nylc_m1_real_nac_final_audit
COMMAND="audit Stage B raw primitives and nominate QM/MM preflight medoids"

set +u
source /work/home/acshdt1dks/opt/gromacs-fastest/env.sh
set -u
export GMX_MAXBACKUP=-1
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"

append_history() {
    local state="$1"
    local detail="$2"
    local timestamp
    timestamp="$(date '+%Y-%m-%dT%H:%M:%S%z')"
    (
        flock -x 9
        printf '%s\t%s\t%s\t%s\t%s\n' "$timestamp" "$EVENT" "$COMMAND" "$state" "$detail" >>"$TASK_ROOT/run_history.tsv"
        "$PYTHON" - "$timestamp" "$state" "$detail" "$JOB_ID" "$STAGEA_RANK_JOB_ID" "$STAGEB_ARRAY_JOB_ID" "$OUT" >>"$TASK_ROOT/run_history.jsonl" <<'PY'
import json
import sys
timestamp, state, detail, job_id, rank_job, stageb_job, output = sys.argv[1:]
print(json.dumps({
    "timestamp": timestamp,
    "event": "nylc_m1_real_nac_final_audit",
    "state": state,
    "detail": detail,
    "slurm_job_id": job_id,
    "stageA_rank_job_id": rank_job,
    "stageB_array_job_id": stageb_job,
    "output_root": output,
}, sort_keys=True))
PY
    ) 9>>"$TASK_ROOT/.run_history.lock"
}

write_record() {
    local candidate_id="$1"
    local seed="$2"
    local run_root="$3"
    local status="$4"
    shift 4
    "$PYTHON" - "$RECORDS_JSONL" "$candidate_id" "$seed" "$run_root" "$status" "$@" <<'PY'
import json
import sys
output, candidate, seed, run_root, status, *paths = sys.argv[1:]
record = {
    "candidate_id": candidate,
    "velocity_seed": int(seed),
    "path": run_root,
    "completion_status": status,
}
if status == "PASS_TECHNICAL_STAGE_B":
    names = [
        "complete_json_path",
        "free_tpr_contract_path",
        "numerical_audit_path",
        "replica_audit_path",
        "primitives_dir",
        "tpr",
        "xtc",
    ]
    if len(paths) != len(names):
        raise SystemExit("PASS Stage B record lacks required path fields")
    record.update(dict(zip(names, paths)))
with open(output, "a") as handle:
    handle.write(json.dumps(record, sort_keys=True) + "\n")
PY
}

if [[ -e "$OUT" ]]; then
    append_history FAIL_TECHNICAL_AUDIT "reason=output_exists;refuse_overwrite=yes"
    exit 2
fi
mkdir -p "$OUT"
append_history START "stageA_rank=$STAGEA_RANK_JOB_ID;stageB_array=$STAGEB_ARRAY_JOB_ID;stageB_count=$STAGEB_COUNT"
failed=1
current=initialization
trap 'rc=$?; if [[ $failed -ne 0 ]]; then append_history FAIL_TECHNICAL_AUDIT "stage=$current;exit_code=$rc;outputs_preserved=yes"; fi; exit $rc' EXIT

read -r selected_count expected_count < <(
    "$PYTHON" - "$TOP6" <<'PY'
import json
import sys
selected = json.load(open(sys.argv[1]))["selected_candidate_ids"]
if not 1 <= len(selected) <= 6 or len(selected) != len(set(selected)):
    raise SystemExit("Stage A selected universe must contain 1-6 unique candidates")
print(len(selected), 3 * len(selected))
PY
)
if [[ "$STAGEB_COUNT" -ne "$expected_count" ]]; then
    echo "Stage B denominator mismatch: got $STAGEB_COUNT expected $expected_count" >&2
    exit 2
fi

RECORDS_JSONL="$OUT/stageB_records.jsonl"
: >"$RECORDS_JSONL"
corrupt=0
seeds=(26711 26723 26737)

for task_id in $(seq 0 $((STAGEB_COUNT - 1))); do
    current="stageB_slot_$task_id"
    read -r candidate_id seed candidate_index < <(
        "$PYTHON" - "$TOP6" "$SELECTION" "$task_id" <<'PY'
import json
import sys
top_path, selection_path, task_id = sys.argv[1:]
selected = json.load(open(top_path))["selected_candidate_ids"]
authority = json.load(open(selection_path))["selected"]
task_id = int(task_id)
candidate_slot, seed_slot = divmod(task_id, 3)
candidate = selected[candidate_slot]
indices = [
    index for index, item in enumerate(authority)
    if item["candidate_id"] == candidate
]
if len(indices) != 1:
    raise SystemExit("selected Stage B candidate is absent or duplicated in authority")
print(candidate, (26711, 26723, 26737)[seed_slot], indices[0])
PY
    )
    CANDIDATE_ROOT="$TASK_ROOT/ensemble/candidates/$candidate_id"
    BUILD="$CANDIDATE_ROOT/build_job_61813799_$candidate_index/build"
    RUN_ROOT="$CANDIDATE_ROOT/replicas/seed$seed/stageB_attempt_${STAGEB_ARRAY_JOB_ID}_${task_id}"

    if [[ -s "$RUN_ROOT/STAGE_B_COMPLETE.json" ]]; then
        required=(run.tpr run.xtc run.edr run.gro run.log free_tpr_contract.json numerical_audit.json)
        missing=()
        for name in "${required[@]}"; do
            [[ -s "$RUN_ROOT/$name" ]] || missing+=("$name")
        done
        if [[ "${#missing[@]}" -ne 0 ]]; then
            corrupt=1
            write_record "$candidate_id" "$seed" "$RUN_ROOT" NOT_EVALUATED_STAGEB_MISSING_ARTIFACT
            continue
        fi
        AUDIT_ROOT="$RUN_ROOT/stageB_audit_job_$JOB_ID"
        if [[ -e "$AUDIT_ROOT" ]]; then
            corrupt=1
            write_record "$candidate_id" "$seed" "$RUN_ROOT" NOT_EVALUATED_STAGEB_AUDIT_OUTPUT_EXISTS
            continue
        fi
        mkdir -p "$AUDIT_ROOT"
        PRIMITIVES="$AUDIT_ROOT/primitives"
        if ! "$PYTHON" "$FLOW/scripts/generate_nylc_m1_ensemble_primitives.py"             --tpr "$RUN_ROOT/run.tpr"             --xtc "$RUN_ROOT/run.xtc"             --source-cycle-ndx "$BUILD/source_cycle.ndx"             --selection-manifest "$SELECTION"             --candidate-id "$candidate_id"             --velocity-seed "$seed"             --output-dir "$PRIMITIVES"             >"$AUDIT_ROOT/generate.stdout" 2>"$AUDIT_ROOT/generate.stderr"; then
            corrupt=1
            write_record "$candidate_id" "$seed" "$RUN_ROOT" NOT_EVALUATED_STAGEB_PRIMITIVE_GENERATION_FAIL
            continue
        fi
        ln -s "$RUN_ROOT/run.log" "$PRIMITIVES/run.log"
        if ! printf 'Potential\n0\n' | "$GMX" energy -f "$RUN_ROOT/run.edr" -o "$PRIMITIVES/potential_energy.xvg"             >"$AUDIT_ROOT/potential.stdout" 2>"$AUDIT_ROOT/potential.stderr"; then
            corrupt=1
            write_record "$candidate_id" "$seed" "$RUN_ROOT" NOT_EVALUATED_STAGEB_ENERGY_EXTRACTION_FAIL
            continue
        fi
        if ! printf 'Temperature\nPressure\nVolume\n0\n' | "$GMX" energy -f "$RUN_ROOT/run.edr" -o "$PRIMITIVES/thermo.xvg"             >"$AUDIT_ROOT/thermo.stdout" 2>"$AUDIT_ROOT/thermo.stderr"; then
            corrupt=1
            write_record "$candidate_id" "$seed" "$RUN_ROOT" NOT_EVALUATED_STAGEB_THERMO_EXTRACTION_FAIL
            continue
        fi
        REPLICA_AUDIT="$AUDIT_ROOT/replica_audit.json"
        if ! "$PYTHON" "$FLOW/scripts/audit_nylc_m1_ensemble_replica.py"             --run-root "$PRIMITIVES"             --manifest "$PRIMITIVES/replica_manifest.json"             --window-start-ps 100             --window-end-ps 1000             --output "$REPLICA_AUDIT"             >"$AUDIT_ROOT/audit.stdout" 2>"$AUDIT_ROOT/audit.stderr"; then
            corrupt=1
            write_record "$candidate_id" "$seed" "$RUN_ROOT" NOT_EVALUATED_STAGEB_REPLICA_AUDIT_FAIL
            continue
        fi
        write_record "$candidate_id" "$seed" "$RUN_ROOT" PASS_TECHNICAL_STAGE_B             "$RUN_ROOT/STAGE_B_COMPLETE.json"             "$RUN_ROOT/free_tpr_contract.json"             "$RUN_ROOT/numerical_audit.json"             "$REPLICA_AUDIT"             "$PRIMITIVES"             "$RUN_ROOT/run.tpr"             "$RUN_ROOT/run.xtc"
    elif [[ -s "$RUN_ROOT/NOT_EVALUATED.json" ]]; then
        parent_status="$("$PYTHON" - "$RUN_ROOT/NOT_EVALUATED.json" <<'PY'
import json
import sys
print(json.load(open(sys.argv[1])).get("status", "NOT_EVALUATED_STAGEB_PARENT"))
PY
)"
        write_record "$candidate_id" "$seed" "$RUN_ROOT" "$parent_status"
    elif [[ -d "$RUN_ROOT" ]]; then
        write_record "$candidate_id" "$seed" "$RUN_ROOT" NOT_EVALUATED_STAGEB_INCOMPLETE
    else
        write_record "$candidate_id" "$seed" "$RUN_ROOT" NOT_EVALUATED_STAGEB_MISSING_SLOT
    fi
done

current=manifest_validation
"$PYTHON" - "$RECORDS_JSONL" "$OUT/stageB_records_manifest.json" "$STAGEB_COUNT" <<'PY'
import json
import sys
source, output, expected = sys.argv[1:]
records = [json.loads(line) for line in open(source) if line.strip()]
expected = int(expected)
pairs = [(item["candidate_id"], item["velocity_seed"]) for item in records]
if len(records) != expected or len(pairs) != len(set(pairs)):
    raise SystemExit("Stage B record denominator is incomplete or duplicated")
open(output, "w").write(json.dumps({
    "schema_version": 1,
    "expected_replica_count": expected,
    "records": records,
}, indent=2, sort_keys=True) + "\n")
PY

if [[ "$corrupt" -ne 0 ]]; then
    "$PYTHON" - "$OUT/AUDIT_CORRUPTION.json" "$STAGEB_ARRAY_JOB_ID" <<'PY'
import json
import sys
path, parent = sys.argv[1:]
open(path, "w").write(json.dumps({
    "schema_version": 1,
    "status": "FAIL_TECHNICAL_AUDIT_CORRUPTION",
    "parent_stageB_array_job_id": parent,
    "reason": "one or more Stage B COMPLETE records lacked auditable raw primitives",
}, indent=2, sort_keys=True) + "\n")
PY
    exit 3
fi

current=assemble_raw_input
"$PYTHON" "$FLOW/scripts/assemble_nylc_m1_final_audit_input.py"     --authority "$SELECTION"     --stage-a-audits "$STAGEA_AUDITS"     --stage-a-top6 "$TOP6"     --stage-b-records "$OUT/stageB_records_manifest.json"     --output-dir "$OUT/assembled"     --forbidden-job-id 61801874     --forbidden-job-id 61803121     --forbidden-job-id 61814393     >"$OUT/assemble.stdout" 2>"$OUT/assemble.stderr"

current=independent_final_gate
"$PYTHON" "$FLOW/scripts/independent_audit_nylc_m1_ensemble.py"     --input "$OUT/assembled/final_audit_input.json"     --output-dir "$OUT/results"     >"$OUT/final_audit.stdout" 2>"$OUT/final_audit.stderr"

scientific_status="$("$PYTHON" - "$OUT/results/FINAL_AUDIT.json" "$OUT/PASS.json" "$JOB_ID" "$STAGEB_ARRAY_JOB_ID" <<'PY'
import json
import sys
source, output, job_id, parent = sys.argv[1:]
audit = json.load(open(source))
payload = {
    "schema_version": 1,
    "status": "PASS_TECHNICAL_FINAL_AUDIT",
    "scientific_status": audit["scientific_status"],
    "slurm_job_id": job_id,
    "parent_stageB_array_job_id": parent,
    "qmmm_eligible_count": audit["qmmm_eligible_count"],
    "scientific_boundary": audit["scientific_boundary"],
}
open(output, "w").write(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(audit["scientific_status"])
PY
)"
append_history PASS_TECHNICAL "stageB_array=$STAGEB_ARRAY_JOB_ID;scientific_status=$scientific_status;results=$OUT/results/FINAL_AUDIT.json"
failed=0
trap - EXIT
exit 0
