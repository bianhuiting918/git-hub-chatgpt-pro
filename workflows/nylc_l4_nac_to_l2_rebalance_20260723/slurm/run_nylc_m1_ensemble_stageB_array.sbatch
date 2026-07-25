#!/bin/bash
#SBATCH -J nylc_m1_stageB
#SBATCH -p xahdnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --gres=dcu:1
#SBATCH --mem-per-cpu=2500M
#SBATCH -t 12:00:00
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-%x-%A_%a.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/slurm-%x-%A_%a.err

set -euo pipefail

TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
FLOW="$TASK_ROOT/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
SELECTION="$TASK_ROOT/ensemble/selection_job_61813011/ensemble_selection.json"
PYTHON=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
GMX=/public/software/apps/Gromacs-DCU2/2022.1/mpi/bin/gmx_mpi
STAGEA_RANK_JOB_ID="${STAGEA_RANK_JOB_ID:?STAGEA_RANK_JOB_ID is required}"
STAGEA_ARRAY_JOB_ID="${STAGEA_ARRAY_JOB_ID:?STAGEA_ARRAY_JOB_ID is required}"
STAGEA_RUN_PREFIX="${STAGEA_RUN_PREFIX:-stageA_job}"
STAGEB_COUNT="${STAGEB_COUNT:?STAGEB_COUNT must equal 3 times selected conformation count}"
JOB_ID="${SLURM_ARRAY_JOB_ID:?submit this script with a dynamic --array range}"
TASK_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"
TASK_COUNT="${SLURM_ARRAY_TASK_COUNT:?SLURM_ARRAY_TASK_COUNT is required}"
case "$STAGEA_RUN_PREFIX" in
    stageA_job|stageA_resume_job) ;;
    *) echo "invalid STAGEA_RUN_PREFIX: $STAGEA_RUN_PREFIX" >&2; exit 2 ;;
esac
TOP6="$TASK_ROOT/ensemble/stageA_rank_job_${STAGEA_RANK_JOB_ID}/results/stageA_top6_manifest.json"
EVENT=nylc_m1_real_nac_stageB

set +u
source /work/home/acshdt1dks/opt/gromacs-fastest/env.sh
set -u
export GMX_MAXBACKUP=-1
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"

read -r selected_count candidate_id seed candidate_index < <(
    "$PYTHON" - "$TOP6" "$SELECTION" "$TASK_ID" <<'PY'
import json
import sys
top_path, selection_path, task_id = sys.argv[1:]
top = json.load(open(top_path))
selected = top["selected_candidate_ids"]
task_id = int(task_id)
if not 0 <= len(selected) <= 6:
    raise SystemExit("selected conformation count is outside 0-6")
if not selected:
    raise SystemExit("Stage B must not be submitted for a zero-candidate manifest")
candidate_slot, seed_slot = divmod(task_id, 3)
if candidate_slot >= len(selected):
    raise SystemExit("array task exceeds selected conformation universe")
seeds = (26711, 26723, 26737)
candidate_id = selected[candidate_slot]
authority = json.load(open(selection_path))["selected"]
indices = [index for index, item in enumerate(authority) if item["candidate_id"] == candidate_id]
if len(indices) != 1:
    raise SystemExit("selected candidate is absent or duplicated in source authority")
print(len(selected), candidate_id, seeds[seed_slot], indices[0])
PY
)
expected_count=$((3 * selected_count))
if [[ "$STAGEB_COUNT" -ne "$expected_count" || "$TASK_COUNT" -ne "$expected_count" ]]; then
    echo "Stage B array count mismatch: export=$STAGEB_COUNT slurm=$TASK_COUNT expected=$expected_count" >&2
    exit 2
fi

CANDIDATE_ROOT="$TASK_ROOT/ensemble/candidates/$candidate_id"
BUILD="$CANDIDATE_ROOT/build_job_61813799_$candidate_index/build"
STAGEA_ROOT="$CANDIDATE_ROOT/replicas/seed$seed/${STAGEA_RUN_PREFIX}_${STAGEA_ARRAY_JOB_ID}_$((candidate_index * 3 + TASK_ID % 3))"
STAGEA_FREE="$STAGEA_ROOT/npt300free_stageA"
OUT="$CANDIDATE_ROOT/replicas/seed$seed/stageB_attempt_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
COMMAND="extend Stage A checkpoint by 900ps candidate=$candidate_id seed=$seed"

append_history() {
    local state="$1"
    local detail="$2"
    local timestamp
    timestamp="$(date '+%Y-%m-%dT%H:%M:%S%z')"
    (
        flock -x 9
        printf '%s\t%s\t%s\t%s\t%s\n' "$timestamp" "$EVENT" "$COMMAND" "$state" "$detail" >>"$TASK_ROOT/run_history.tsv"
        "$PYTHON" - "$timestamp" "$state" "$detail" "$JOB_ID" "$TASK_ID" "$candidate_id" "$seed" "$STAGEA_ROOT" "$OUT" >>"$TASK_ROOT/run_history.jsonl" <<'PY'
import json
import sys
timestamp, state, detail, job_id, task_id, candidate_id, seed, parent, output = sys.argv[1:]
print(json.dumps({
    "timestamp": timestamp,
    "event": "nylc_m1_real_nac_stageB",
    "state": state,
    "candidate_id": candidate_id,
    "velocity_seed": int(seed),
    "slurm_job_id": job_id,
    "slurm_array_task_id": task_id,
    "parent_stageA_root": parent,
    "output_root": output,
    "detail": detail,
}, sort_keys=True))
PY
    ) 9>>"$TASK_ROOT/.run_history.lock"
}

scan_log() {
    "$PYTHON" - "$1" "$2" <<'PY'
import json
import re
import sys
text = open(sys.argv[1], errors="replace").read()
audit = {
    "fatal": len(re.findall(r"(?i)fatal(?:\s+error)?", text)),
    "lincs_warning": len(re.findall(r"(?i)lincs\s+warning", text)),
    "settle_problem": len(re.findall(r"(?i)(?:cannot\s+be\s+settled|settle[^\n]*(?:warning|error|failed))", text)),
    "nan": len(re.findall(r"(?i)\bnan\b", text)),
    "finished_mdrun": bool(re.search(r"Finished mdrun", text)),
}
open(sys.argv[2], "w").write(json.dumps(audit, indent=2, sort_keys=True) + "\n")
raise SystemExit(
    0 if audit["finished_mdrun"]
    and not any(audit[key] for key in ("fatal", "lincs_warning", "settle_problem", "nan"))
    else 1
)
PY
}

if [[ -e "$OUT" ]]; then
    append_history FAIL_TECHNICAL_STAGEB "reason=output_exists;refuse_overwrite=yes"
    exit 2
fi
mkdir -p "$OUT"
if [[ ! -s "$STAGEA_ROOT/STAGE_A_COMPLETE.json" ]]; then
    append_history NOT_EVALUATED_STAGEA_INCOMPLETE "parent_stageA=$STAGEA_ROOT"
    "$PYTHON" - "$OUT/NOT_EVALUATED.json" "$candidate_id" "$seed" "$STAGEA_ROOT" <<'PY'
import json
import sys
path, candidate_id, seed, parent = sys.argv[1:]
open(path, "w").write(json.dumps({
    "schema_version": 1,
    "status": "NOT_EVALUATED_STAGEA_INCOMPLETE",
    "candidate_id": candidate_id,
    "velocity_seed": int(seed),
    "parent_stageA_root": parent,
}, indent=2, sort_keys=True) + "\n")
PY
    exit 0
fi
for name in run.gro run.cpt run.tpr; do
    test -s "$STAGEA_FREE/$name"
done

parent_gro_sha="$(sha256sum "$STAGEA_FREE/run.gro" | awk '{print $1}')"
parent_cpt_sha="$(sha256sum "$STAGEA_FREE/run.cpt" | awk '{print $1}')"
append_history START "parent_stageA=$STAGEA_ROOT;parent_gro_sha256=$parent_gro_sha;parent_cpt_sha256=$parent_cpt_sha;extension_ps=900;velocities_regenerated=no"

failed=1
current=grompp
trap 'rc=$?; if [[ $failed -ne 0 ]]; then append_history FAIL_TECHNICAL_STAGEB "stage=$current;exit_code=$rc;outputs_preserved=yes"; fi; exit $rc' EXIT

(
    cd "$BUILD"
    "$GMX" grompp \
        -f "$FLOW/mdp/npt300free_m1_stageB_extend.mdp" \
        -c "$STAGEA_FREE/run.gro" \
        -t "$STAGEA_FREE/run.cpt" \
        -p "$BUILD/topol.top" \
        -o "$OUT/run.tpr" \
        -po "$OUT/run.processed.mdp" \
        -maxwarn 0 >"$OUT/grompp.stdout" 2>"$OUT/grompp.stderr"
)
"$GMX" dump -s "$OUT/run.tpr" >"$OUT/tpr.dump" 2>"$OUT/tpr_dump.stderr"
"$PYTHON" - "$OUT/tpr.dump" "$OUT/run.processed.mdp" "$OUT/free_tpr_contract.json" <<'PY'
import json
import re
import sys
dump_path, mdp_path, output_path = sys.argv[1:]
dump = open(dump_path, errors="replace").read()
mdp = open(mdp_path, errors="replace").read()
position = sum(int(value) for value in re.findall(r"#posres_xA\s*=\s*(\d+)", dump))
distance = len(re.findall(r"(?i)\bDISRES\b", dump))
define_values = []
for line in mdp.splitlines():
    body = line.split(";", 1)[0].strip()
    if "=" not in body:
        continue
    key, value = body.split("=", 1)
    if key.strip().lower() == "define" and value.strip():
        define_values.append(value.strip())
payload = {
    "position_restraints": position,
    "distance_restraints": distance,
    "processed_mdp_nonempty_define_values": define_values,
    "scientific_window": "fully_unrestrained_NPT_extension_900ps",
}
open(output_path, "w").write(json.dumps(payload, indent=2, sort_keys=True) + "\n")
if position or distance or define_values:
    raise SystemExit("Stage B TPR contains restraints")
PY

current=mdrun
mpirun -np 1 "$GMX" mdrun -s "$OUT/run.tpr" -deffnm "$OUT/run" \
    -ntomp "${SLURM_CPUS_PER_TASK:-8}" -nb gpu -pme gpu -bonded gpu -pin on \
    >"$OUT/mdrun.stdout" 2>"$OUT/mdrun.stderr"
test -s "$OUT/run.gro"
test -s "$OUT/run.cpt"
scan_log "$OUT/run.log" "$OUT/numerical_audit.json"
"$PYTHON" - "$OUT/STAGE_B_COMPLETE.json" "$candidate_id" "$seed" "$JOB_ID" "$TASK_ID" "$STAGEA_ROOT" "$parent_gro_sha" "$parent_cpt_sha" <<'PY'
import hashlib
import json
import pathlib
import sys
path, candidate_id, seed, job_id, task_id, parent, gro_sha, cpt_sha = sys.argv[1:]
root = pathlib.Path(path).parent
payload = {
    "schema_version": 1,
    "status": "PASS_TECHNICAL_STAGE_B",
    "scientific_status": "NOT_EVALUATED_PENDING_STAGE_B_AUDIT",
    "candidate_id": candidate_id,
    "velocity_seed": int(seed),
    "slurm_job_id": job_id,
    "slurm_array_task_id": task_id,
    "parent_stageA_root": parent,
    "parent_stageA_gro_sha256": gro_sha,
    "parent_stageA_cpt_sha256": cpt_sha,
    "extension_ps": 900.0,
    "fully_unrestrained": True,
    "velocities_regenerated": False,
    "final_gro_sha256": hashlib.sha256((root / "run.gro").read_bytes()).hexdigest(),
}
open(path, "w").write(json.dumps(payload, indent=2, sort_keys=True) + "\n")
PY

append_history PASS_TECHNICAL "extension_ps=900;fully_unrestrained=yes;velocities_regenerated=no;scientific_status=NOT_EVALUATED_PENDING_STAGE_B_AUDIT"
failed=0
trap - EXIT
