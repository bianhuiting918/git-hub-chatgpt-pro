#!/usr/bin/env bash
set -euo pipefail
: "${A1_STEP2_WATER_PERSIST_CODE_ROOT:?set immutable code root}"
: "${A1_STEP2_WATER_PERSIST_GITHUB_COMMIT:?set immutable GitHub commit}"
TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
CODE_ROOT="$A1_STEP2_WATER_PERSIST_CODE_ROOT"
GITHUB_COMMIT="$A1_STEP2_WATER_PERSIST_GITHUB_COMMIT"
ARRAY_INDEX="${SLURM_ARRAY_TASK_ID:?run as combined array task 0..8}"
ARRAY_JOB="${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-manual}}"
if (( ARRAY_INDEX < 0 || ARRAY_INDEX > 8 )); then exit 2; fi

if (( ARRAY_INDEX == 0 )); then
    export A1_STEP2_WATER_REORG_CODE_ROOT="$CODE_ROOT"
    export A1_STEP2_WATER_REORG_GITHUB_COMMIT="$GITHUB_COMMIT"
    exec bash "$CODE_ROOT/scripts/run_nylc_a1_step2_water_reorganization_sampling.sh"
fi

PERSIST_INDEX=$((ARRAY_INDEX - 1))
SEED_INDEX=$((PERSIST_INDEX / 4))
REPLICA=$((PERSIST_INDEX % 4))
if (( SEED_INDEX == 0 )); then
    SEED=seed26723
    SOURCE_RST7="$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_step2_water_reorganization_sampling/attempt_62300430_3/selected_near_miss_0.rst7"
    EXPECTED_SOURCE_SHA=a65834ae6c0f82291fd0b8c38ece41983b501b2444deb98de1794b92d79622a7
else
    SEED=seed26737
    SOURCE_RST7="$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_step2_water_reorganization_sampling/attempt_62300430_7/selected_near_miss_0.rst7"
    EXPECTED_SOURCE_SHA=34ad5af01f5a2cd7e7ce0342590ab66f201902a31fd49236546441afa82f5ea6
fi
PRMTOP="$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_unified_core_dftb3_preflight/attempt_62011285/prepared/system.prmtop"
EXPECTED_PRMTOP_SHA=a61d15bf0bf78675be93275d45f274e808ed6ae450fc1ca21a8e14aee8c12ca0
OUT_ROOT="$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_step2_water_nearmiss_persistence"
OUT="$OUT_ROOT/attempt_${ARRAY_JOB}_${PERSIST_INDEX}"
SCRATCH_ROOT="${SLURM_TMPDIR:-/tmp}/nylc_a1_step2_water_persist_${ARRAY_JOB}_${PERSIST_INDEX}"
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
DRIVER="$CODE_ROOT/scripts/prepare_audit_nylc_a1_step2_water_nearmiss_persistence.py"
SUCCESS=0

append_history() {
    local state="$1" detail="$2" now
    now="$(date '+%Y-%m-%dT%H:%M:%S%z')"
    (
        flock -x 9
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$now" "nylc_a1_step2_water_nearmiss_persistence" "run_nylc_a1_step2_water_nearmiss_persistence.sh" "$state" "$GITHUB_COMMIT" "$detail" >> "$TASK_ROOT/run_history.tsv"
        printf '{"time":"%s","event":"nylc_a1_step2_water_nearmiss_persistence","state":"%s","git_commit":"%s","detail":"%s","output":"%s"}\n' "$now" "$state" "$GITHUB_COMMIT" "$detail" "$OUT" >> "$TASK_ROOT/run_history.jsonl"
    ) 9>> "$TASK_ROOT/.run_history.lock"
}

on_exit() {
    local rc=$?
    if (( SUCCESS == 0 )); then
        append_history FAILED "job=$ARRAY_JOB;array_task=$ARRAY_INDEX;persistence_task=$PERSIST_INDEX;seed=$SEED;replica=$REPLICA;rc=$rc"
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
append_history STARTED "job=$ARRAY_JOB;array_task=$ARRAY_INDEX;persistence_task=$PERSIST_INDEX;seed=$SEED;replica=$REPLICA;source_sha=$EXPECTED_SOURCE_SHA"

module purge >/dev/null 2>&1 || true
module load amber/2018-hpcx-gcc-7.3.1
export AMBERHOME="$TASK_ROOT/../nylc_gyaq_pa66_l2_nac_qmmm_20260723/amber_runtime"
export OMP_NUM_THREADS=1
test -s "$AMBERHOME/dat/slko/3ob-3-1/C-C.skf"

"$PY" "$DRIVER" --mode initialize --task-index "$PERSIST_INDEX" --output "$OUT" --scratch "$SCRATCH_ROOT" --code-root "$CODE_ROOT" --github-commit "$GITHUB_COMMIT"
cd "$SCRATCH_ROOT"
set +e
mpirun --bind-to none -np 8 sander.MPI -O -i stage.in -o stage.out -p "$PRMTOP" -c "$SOURCE_RST7" -r stage.rst7 -x stage.nc -inf stage.mdinfo
ENGINE_RC=$?
set -e
printf '%s\n' "$ENGINE_RC" > engine.rc
"$PY" "$DRIVER" --mode audit --output "$OUT" --scratch "$SCRATCH_ROOT"
"$PY" "$DRIVER" --mode merge-if-ready --output-root "$OUT_ROOT" --array-job "$ARRAY_JOB"
STATUS="$("$PY" - "$OUT/RESULT.json" <<'PY'
import json,sys
print(json.load(open(sys.argv[1],encoding="utf-8"))["status"])
PY
)"
append_history TERMINAL "job=$ARRAY_JOB;array_task=$ARRAY_INDEX;persistence_task=$PERSIST_INDEX;seed=$SEED;replica=$REPLICA;status=$STATUS;automatic_downstream=NONE"
SUCCESS=1
printf '%s seed=%s replica=%s output=%s\n' "$STATUS" "$SEED" "$REPLICA" "$OUT"
