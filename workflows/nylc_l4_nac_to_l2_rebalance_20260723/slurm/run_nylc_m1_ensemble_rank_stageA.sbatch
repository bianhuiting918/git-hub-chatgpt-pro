#!/bin/bash
#SBATCH -J nylc_m1_rankA
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem-per-cpu=2500M
#SBATCH -t 08:00:00
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-%x-%j.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-%x-%j.err

# Submit only with --dependency=afterany:<complete Stage A array job>.
set -euo pipefail

TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
FLOW="$TASK_ROOT/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
SELECTION="$TASK_ROOT/ensemble/selection_job_61813011/ensemble_selection.json"
PYTHON=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
GMX=/public/software/apps/Gromacs-DCU2/2022.1/mpi/bin/gmx_mpi
STAGEA_ARRAY_JOB_ID="${STAGEA_ARRAY_JOB_ID:?submit with --export=ALL,STAGEA_ARRAY_JOB_ID=<Stage A array>}"
JOB_ID="${SLURM_JOB_ID:?SLURM_JOB_ID is required}"
OUT="$TASK_ROOT/ensemble/stageA_rank_job_$JOB_ID"
expected_stageA_slots=36
seeds=(26711 26723 26737)
EVENT=nylc_m1_real_nac_stageA_rank
COMMAND="audit 36 Stage A slots and select top six conformations"

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
        "$PYTHON" - "$timestamp" "$state" "$detail" "$JOB_ID" "$STAGEA_ARRAY_JOB_ID" "$OUT" >>"$TASK_ROOT/run_history.jsonl" <<'PY'
import json
import sys
timestamp, state, detail, job_id, parent_job, output_root = sys.argv[1:]
print(json.dumps({
    "timestamp": timestamp,
    "event": "nylc_m1_real_nac_stageA_rank",
    "state": state,
    "slurm_job_id": job_id,
    "parent_stageA_array_job_id": parent_job,
    "output_root": output_root,
    "detail": detail,
}, sort_keys=True))
PY
    ) 9>>"$TASK_ROOT/.run_history.lock"
}

if [[ -e "$OUT" ]]; then
    append_history FAIL_TECHNICAL_AUDIT "reason=output_exists;refuse_overwrite=yes"
    exit 2
fi
mkdir -p "$OUT"
append_history START "parent_stageA_array=$STAGEA_ARRAY_JOB_ID;expected_stageA_slots=$expected_stageA_slots;dependency=afterany"

failed=1
current=initialization
trap 'rc=$?; if [[ $failed -ne 0 ]]; then append_history FAIL_TECHNICAL_AUDIT "stage=$current;exit_code=$rc;outputs_preserved=yes"; fi; exit $rc' EXIT

AUDITS_JSONL="$OUT/replica_audits.jsonl"
FINGERPRINT_TSV="$OUT/fingerprint_paths.tsv"
: >"$AUDITS_JSONL"
: >"$FINGERPRINT_TSV"
corrupt=0

write_not_evaluated() {
    local output="$1"
    local candidate_id="$2"
    local seed="$3"
    local status="$4"
    local reason="$5"
    "$PYTHON" - "$output" "$candidate_id" "$seed" "$status" "$reason" <<'PY'
import json
import sys
path, candidate_id, seed, status, reason = sys.argv[1:]
payload = {
    "schema_version": 1,
    "candidate_id": candidate_id,
    "velocity_seed": int(seed),
    "microstate": "M1_NalphaH2_Asp306H",
    "gate_definition": "NylC residues 261-266; Thr267 excluded",
    "fully_unrestrained": False,
    "technical_status": "FAIL",
    "scientific_status": status,
    "nac": {
        "nac_frame_count": 0,
        "nac_frame_count_after_20ps": 0,
        "analysis_frame_count_after_20ps": 0,
        "nac_occupancy_after_20ps": 0.0,
        "longest_event_after_20ps": {"duration_ps": 0.0, "frame_count": 0},
    },
    "reason": reason,
}
open(path, "w").write(json.dumps(payload, indent=2, sort_keys=True) + "\n")
PY
}

for task_id in $(seq 0 35); do
    current="slot_$task_id"
    candidate_index=$((task_id / 3))
    seed_index=$((task_id % 3))
    seed="${seeds[$seed_index]}"
    candidate_id="$("$PYTHON" - "$SELECTION" "$candidate_index" <<'PY'
import json
import sys
print(json.load(open(sys.argv[1]))["selected"][int(sys.argv[2])]["candidate_id"])
PY
)"
    CANDIDATE_ROOT="$TASK_ROOT/ensemble/candidates/$candidate_id"
    BUILD="$CANDIDATE_ROOT/build_job_61813799_$candidate_index/build"
    RUN_ROOT="$CANDIDATE_ROOT/replicas/seed$seed/stageA_job_${STAGEA_ARRAY_JOB_ID}_${task_id}"
    AUDIT_ROOT="$RUN_ROOT/stageA_audit_job_$JOB_ID"
    mkdir -p "$AUDIT_ROOT"
    AUDIT_JSON="$AUDIT_ROOT/replica_audit.json"

    if [[ -s "$RUN_ROOT/STAGE_A_COMPLETE.json" ]]; then
        FREE="$RUN_ROOT/npt300free_stageA"
        required=(run.tpr run.xtc run.edr run.gro run.log free_tpr_contract.json)
        missing=()
        for name in "${required[@]}"; do
            [[ -s "$FREE/$name" ]] || missing+=("$name")
        done
        if [[ "${#missing[@]}" -ne 0 ]]; then
            write_not_evaluated "$AUDIT_JSON" "$candidate_id" "$seed" "NOT_EVALUATED_MISSING_PRIMITIVE" "STAGE_A_COMPLETE exists but missing ${missing[*]}"
            corrupt=1
        else
            PRIMITIVES="$AUDIT_ROOT/primitives"
            if ! "$PYTHON" "$FLOW/scripts/generate_nylc_m1_ensemble_primitives.py" \
                --tpr "$FREE/run.tpr" \
                --xtc "$FREE/run.xtc" \
                --source-cycle-ndx "$BUILD/source_cycle.ndx" \
                --selection-manifest "$SELECTION" \
                --candidate-id "$candidate_id" \
                --velocity-seed "$seed" \
                --output-dir "$PRIMITIVES" \
                >"$AUDIT_ROOT/generate.stdout" 2>"$AUDIT_ROOT/generate.stderr"; then
                write_not_evaluated "$AUDIT_JSON" "$candidate_id" "$seed" "NOT_EVALUATED_AUDIT_GENERATION_FAIL" "primitive generator failed; see generate.stderr"
                corrupt=1
            else
                ln -s "$FREE/run.log" "$PRIMITIVES/run.log"
                printf 'Potential\n0\n' | "$GMX" energy -f "$FREE/run.edr" -o "$PRIMITIVES/potential_energy.xvg" \
                    >"$AUDIT_ROOT/potential.stdout" 2>"$AUDIT_ROOT/potential.stderr"
                printf 'Temperature\nPressure\nVolume\n0\n' | "$GMX" energy -f "$FREE/run.edr" -o "$PRIMITIVES/thermo.xvg" \
                    >"$AUDIT_ROOT/thermo.stdout" 2>"$AUDIT_ROOT/thermo.stderr"
                read -r window_start window_end < <(
                    "$PYTHON" - "$PRIMITIVES/nac_distance.xvg" <<'PY'
import sys
rows = [
    line.split() for line in open(sys.argv[1])
    if line.strip() and not line.lstrip().startswith(("#", "@"))
]
print(rows[0][0], rows[-1][0])
PY
                )
                if ! "$PYTHON" "$FLOW/scripts/audit_nylc_m1_ensemble_replica.py" \
                    --run-root "$PRIMITIVES" \
                    --manifest "$PRIMITIVES/replica_manifest.json" \
                    --window-start-ps "$window_start" \
                    --window-end-ps "$window_end" \
                    --output "$AUDIT_JSON" \
                    >"$AUDIT_ROOT/audit.stdout" 2>"$AUDIT_ROOT/audit.stderr"; then
                    [[ -s "$AUDIT_JSON" ]] || corrupt=1
                fi
                printf '%s\t%s\t%s\n' "$candidate_id" "$seed" "$PRIMITIVES/nac_fingerprints.npz" >>"$FINGERPRINT_TSV"
            fi
        fi
    elif [[ -s "$RUN_ROOT/NOT_EVALUATED.json" ]]; then
        parent_status="$("$PYTHON" - "$RUN_ROOT/NOT_EVALUATED.json" <<'PY'
import json
import sys
print(json.load(open(sys.argv[1])).get("status", "NOT_EVALUATED_STAGEA_PARENT"))
PY
)"
        write_not_evaluated "$AUDIT_JSON" "$candidate_id" "$seed" "$parent_status" "Stage A parent explicitly did not run"
    elif [[ -d "$RUN_ROOT" ]]; then
        write_not_evaluated "$AUDIT_JSON" "$candidate_id" "$seed" "NOT_EVALUATED_STAGEA_INCOMPLETE" "Stage A run root exists without completion record"
    else
        write_not_evaluated "$AUDIT_JSON" "$candidate_id" "$seed" "NOT_EVALUATED_MISSING_STAGEA_SLOT" "No Stage A run root or explicit failure record"
        corrupt=1
    fi
    if [[ ! -s "$AUDIT_JSON" ]]; then
        corrupt=1
        continue
    fi
    "$PYTHON" - "$AUDIT_JSON" "$task_id" >>"$AUDITS_JSONL" <<'PY'
import json
import sys
record = json.load(open(sys.argv[1]))
record["stageA_slot"] = int(sys.argv[2])
print(json.dumps(record, sort_keys=True))
PY
done

current=universe_validation
"$PYTHON" - "$AUDITS_JSONL" "$FINGERPRINT_TSV" "$OUT/stageA_replica_audits_input.json" "$OUT/fingerprint_map.json" "$expected_stageA_slots" <<'PY'
import json
import sys
from collections import defaultdict
jsonl, fingerprint_tsv, audit_output, map_output, expected = sys.argv[1:]
audits = [json.loads(line) for line in open(jsonl) if line.strip()]
expected = int(expected)
slots = [item["stageA_slot"] for item in audits]
if len(audits) != expected or sorted(slots) != list(range(expected)):
    raise SystemExit(f"Stage A universe mismatch: {len(audits)} records, slots={sorted(slots)}")
keys = [(item["candidate_id"], item["velocity_seed"]) for item in audits]
if len(keys) != len(set(keys)):
    raise SystemExit("duplicate candidate/seed audit")
fingerprints = defaultdict(dict)
for raw in open(fingerprint_tsv):
    if not raw.strip():
        continue
    candidate, seed, path = raw.rstrip("\n").split("\t")
    fingerprints[candidate][seed] = path
open(audit_output, "w").write(json.dumps(audits, indent=2, sort_keys=True) + "\n")
open(map_output, "w").write(json.dumps(fingerprints, indent=2, sort_keys=True) + "\n")
PY

if [[ "$corrupt" -ne 0 ]]; then
    "$PYTHON" - "$OUT/AUDIT_CORRUPTION.json" "$STAGEA_ARRAY_JOB_ID" <<'PY'
import json
import sys
path, parent = sys.argv[1:]
open(path, "w").write(json.dumps({
    "schema_version": 1,
    "status": "FAIL_TECHNICAL_AUDIT_CORRUPTION",
    "parent_stageA_array_job_id": parent,
    "reason": "one or more complete or expected Stage A slots lacked auditable primitives",
}, indent=2, sort_keys=True) + "\n")
PY
    exit 3
fi

current=rank_stageA
"$PYTHON" "$FLOW/scripts/rank_nylc_m1_ensemble.py" \
    --audits "$OUT/stageA_replica_audits_input.json" \
    --fingerprint-map "$OUT/fingerprint_map.json" \
    --output-dir "$OUT/results" \
    >"$OUT/rank.stdout" 2>"$OUT/rank.stderr"

"$PYTHON" - "$OUT/results/stageA_ranking_summary.json" "$OUT/PASS.json" "$JOB_ID" "$STAGEA_ARRAY_JOB_ID" <<'PY'
import json
import sys
summary_path, output_path, job_id, parent = sys.argv[1:]
summary = json.load(open(summary_path))
payload = {
    "schema_version": 1,
    "status": "PASS_TECHNICAL_STAGEA_RANK",
    "scientific_status": summary["selection_status"],
    "slurm_job_id": job_id,
    "parent_stageA_array_job_id": parent,
    "stageA_slot_count": summary["replica_audit_count"],
    "selected_candidate_ids": summary["selected_candidate_ids"],
    "scientific_boundary": "Stage A triage only; Stage B is required.",
}
open(output_path, "w").write(json.dumps(payload, indent=2, sort_keys=True) + "\n")
PY

append_history PASS_TECHNICAL "parent_stageA_array=$STAGEA_ARRAY_JOB_ID;expected_stageA_slots=$expected_stageA_slots;summary=$OUT/results/stageA_ranking_summary.json;scientific_status=STAGEA_TRIAGE_ONLY"
failed=0
trap - EXIT
