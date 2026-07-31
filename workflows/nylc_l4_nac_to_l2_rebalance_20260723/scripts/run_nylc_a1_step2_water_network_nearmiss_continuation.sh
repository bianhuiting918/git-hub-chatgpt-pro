#!/usr/bin/env bash
set -euo pipefail
: "${A1_STEP2_NEARMISS_CODE_ROOT:?set immutable code root}"
: "${A1_STEP2_NEARMISS_GITHUB_COMMIT:?set immutable GitHub commit}"

TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
CODE_ROOT="$A1_STEP2_NEARMISS_CODE_ROOT"
GITHUB_COMMIT="$A1_STEP2_NEARMISS_GITHUB_COMMIT"
INDEX="${SLURM_ARRAY_TASK_ID:?run as array task 0..7}"
ARRAY_JOB="${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-manual}}"
if (( INDEX < 0 || INDEX > 7 )); then exit 2; fi
SEED_INDEX=$((INDEX / 4))
REPLICA=$((INDEX % 4))

if (( SEED_INDEX == 0 )); then
    SEED=seed26723
    SOURCE_ATTEMPT=attempt_62425747_2
    SOURCE_RST7="$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_step2_water_network_sampling/$SOURCE_ATTEMPT/direct_near_miss_0.rst7"
    EXPECTED_SOURCE_SHA=309e419b7f14c73c9388fea7325d9e54e2434669019e24b209bf0052981058a4
    SOURCE_WATER_O=43796
    VELOCITY_SEED=$((26723401 + REPLICA))
else
    SEED=seed26737
    SOURCE_ATTEMPT=attempt_62425747_6
    SOURCE_RST7="$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_step2_water_network_sampling/$SOURCE_ATTEMPT/direct_near_miss_0.rst7"
    EXPECTED_SOURCE_SHA=42c46ebe61ad3016c86a91880ac1d6f94d7f7c7fdcda385bd89935277d378583
    SOURCE_WATER_O=13046
    VELOCITY_SEED=$((26737401 + REPLICA))
fi

PRMTOP="$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_unified_core_dftb3_preflight/attempt_62011285/prepared/system.prmtop"
EXPECTED_PRMTOP_SHA=a61d15bf0bf78675be93275d45f274e808ed6ae450fc1ca21a8e14aee8c12ca0
OUT_ROOT="$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_step2_water_network_nearmiss_continuation"
OUT="$OUT_ROOT/attempt_${ARRAY_JOB}_${INDEX}"
SCRATCH_ROOT="${SLURM_TMPDIR:-/tmp}/nylc_a1_step2_nearmiss_${ARRAY_JOB}_${INDEX}"
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
DRIVER="$CODE_ROOT/scripts/prepare_audit_nylc_a1_step2_water_network_sampling.py"
SUCCESS=0

append_history() {
    local state="$1" detail="$2" now
    now="$(date '+%Y-%m-%dT%H:%M:%S%z')"
    (
        flock -x 9
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$now" "nylc_a1_step2_nearmiss_continuation" "$(basename "$0")" "$state" "$GITHUB_COMMIT" "$detail" >> "$TASK_ROOT/run_history.tsv"
        printf '{"time":"%s","event":"nylc_a1_step2_nearmiss_continuation","state":"%s","git_commit":"%s","detail":"%s","output":"%s"}\n' "$now" "$state" "$GITHUB_COMMIT" "$detail" "$OUT" >> "$TASK_ROOT/run_history.jsonl"
    ) 9>> "$TASK_ROOT/.run_history.lock"
}

on_exit() {
    local rc=$?
    if (( SUCCESS == 0 )); then
        append_history FAILED "job=$ARRAY_JOB;task=$INDEX;seed=$SEED;replica=$REPLICA;rc=$rc"
    fi
}
trap on_exit EXIT

test -s "$DRIVER"
test -s "$SOURCE_RST7"
test -s "$PRMTOP"
test "$(sha256sum "$SOURCE_RST7" | awk '{print $1}')" = "$EXPECTED_SOURCE_SHA"
test "$(sha256sum "$PRMTOP" | awk '{print $1}')" = "$EXPECTED_PRMTOP_SHA"
test ! -e "$OUT"
test ! -e "$SCRATCH_ROOT"
append_history STARTED "job=$ARRAY_JOB;task=$INDEX;seed=$SEED;replica=$REPLICA;source=$SOURCE_ATTEMPT;source_sha=$EXPECTED_SOURCE_SHA;water_o=$SOURCE_WATER_O"

module purge >/dev/null 2>&1 || true
module load amber/2018-hpcx-gcc-7.3.1
export AMBERHOME="$TASK_ROOT/../nylc_gyaq_pa66_l2_nac_qmmm_20260723/amber_runtime"
export OMP_NUM_THREADS=1
test -s "$AMBERHOME/dat/slko/3ob-3-1/C-C.skf"

"$PY" "$DRIVER" --mode initialize --task-index "$INDEX" --output "$OUT" --scratch "$SCRATCH_ROOT" --code-root "$CODE_ROOT" --github-commit "$GITHUB_COMMIT"
sed -i "s/ig=[0-9][0-9]*/ig=$VELOCITY_SEED/" "$SCRATCH_ROOT/stage.in"
"$PY" - "$OUT/SOURCE_MANIFEST.json" "$SOURCE_ATTEMPT" "$SOURCE_RST7" "$EXPECTED_SOURCE_SHA" "$SOURCE_WATER_O" "$VELOCITY_SEED" <<'PY'
import json, os, pathlib, sys
path = pathlib.Path(sys.argv[1])
payload = json.loads(path.read_text(encoding="utf-8"))
payload["source"] = {
    "attempt": sys.argv[2],
    "restart": sys.argv[3],
    "restart_sha256": sys.argv[4],
    "route": "direct_near_miss",
    "water_oxygen_index1": int(sys.argv[5]),
}
payload["velocity_seed"] = int(sys.argv[6])
payload["lineage"] = {
    "parent_job": "62425747",
    "selection": "lowest direct near_miss_score per seed",
    "all_waters_remain_mm": True,
    "water_identity_fixed": False,
}
payload["protocol"]["near_miss_seeded_continuation"] = True
payload["protocol"]["reactive_restraints"] = False
payload["protocol"]["water_position_restraints"] = False
payload["protocol"]["water_identity_restraints"] = False
temporary = path.with_suffix(".json.tmp")
temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
os.replace(temporary, path)
PY

cd "$SCRATCH_ROOT"
set +e
mpirun --bind-to none -np 8 sander.MPI -O -i stage.in -o stage.out -p "$PRMTOP" -c "$SOURCE_RST7" -r stage.rst7 -x stage.nc -inf stage.mdinfo
ENGINE_RC=$?
set -e
printf '%s\n' "$ENGINE_RC" > engine.rc
"$PY" "$DRIVER" --mode audit --output "$OUT" --scratch "$SCRATCH_ROOT"
(
    flock -x 9
    "$PY" "$DRIVER" --mode merge-if-ready --output-root "$OUT_ROOT" --array-job "$ARRAY_JOB"
) 9>> "$OUT_ROOT/.merge.lock"
STATUS="$("$PY" - "$OUT/RESULT.json" <<'PY'
import json, sys
print(json.load(open(sys.argv[1], encoding="utf-8"))["status"])
PY
)"
append_history TERMINAL "job=$ARRAY_JOB;task=$INDEX;seed=$SEED;replica=$REPLICA;status=$STATUS;automatic_downstream=NONE"
SUCCESS=1
printf '%s seed=%s replica=%s source=%s output=%s\n' "$STATUS" "$SEED" "$REPLICA" "$SOURCE_ATTEMPT" "$OUT"
