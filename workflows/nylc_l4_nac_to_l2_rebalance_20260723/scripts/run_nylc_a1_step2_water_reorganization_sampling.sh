#!/usr/bin/env bash
set -euo pipefail
: "${A1_STEP2_WATER_REORG_CODE_ROOT:?set immutable code root}"
: "${A1_STEP2_WATER_REORG_GITHUB_COMMIT:?set immutable GitHub commit}"
TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
CODE_ROOT="$A1_STEP2_WATER_REORG_CODE_ROOT"
GITHUB_COMMIT="$A1_STEP2_WATER_REORG_GITHUB_COMMIT"
INDEX="${SLURM_ARRAY_TASK_ID:?run as array task 0..7}"
ARRAY_JOB="${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-manual}}"
if (( INDEX < 0 || INDEX > 7 )); then exit 2; fi
SEED_INDEX=$((INDEX / 4))
REPLICA=$((INDEX % 4))
if (( SEED_INDEX == 0 )); then
    SEED=seed26723
    SOURCE_RST7="$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_step1_acyl_release_md_continuation/attempt_62216380_0/release_md_endpoint.rst7"
    EXPECTED_SOURCE_SHA=5d8f76d2c90e3e8c707b640c55a93938f53e18dc30d6f93d54adda467da25f41
else
    SEED=seed26737
    SOURCE_RST7="$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_step1_acyl_release_md_continuation/attempt_62216380_1/release_md_endpoint.rst7"
    EXPECTED_SOURCE_SHA=4cc60ad4d7be099b3f76040f1ee8c49b91deecebb20ed98570bc192b3d511132
fi
PRMTOP="$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_unified_core_dftb3_preflight/attempt_62011285/prepared/system.prmtop"
EXPECTED_PRMTOP_SHA=a61d15bf0bf78675be93275d45f274e808ed6ae450fc1ca21a8e14aee8c12ca0
OUT_ROOT="$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_step2_water_reorganization_sampling"
OUT="$OUT_ROOT/attempt_${ARRAY_JOB}_${INDEX}"
SCRATCH_ROOT="${SLURM_TMPDIR:-/tmp}/nylc_a1_step2_water_reorg_${ARRAY_JOB}_${INDEX}"
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
DRIVER="$CODE_ROOT/scripts/prepare_audit_nylc_a1_step2_water_reorganization_sampling.py"
SUCCESS=0

append_history() {
    local state="$1" detail="$2" now
    now="$(date '+%Y-%m-%dT%H:%M:%S%z')"
    (
        flock -x 9
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$now" "nylc_a1_step2_water_reorganization_sampling" "run_nylc_a1_step2_water_reorganization_sampling.sh" "$state" "$GITHUB_COMMIT" "$detail" >> "$TASK_ROOT/run_history.tsv"
        printf '{"time":"%s","event":"nylc_a1_step2_water_reorganization_sampling","state":"%s","git_commit":"%s","detail":"%s","output":"%s"}\n' "$now" "$state" "$GITHUB_COMMIT" "$detail" "$OUT" >> "$TASK_ROOT/run_history.jsonl"
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
append_history STARTED "job=$ARRAY_JOB;task=$INDEX;seed=$SEED;replica=$REPLICA;source_sha=$EXPECTED_SOURCE_SHA"

module purge >/dev/null 2>&1 || true
module load amber/2018-hpcx-gcc-7.3.1
export AMBERHOME="$TASK_ROOT/../nylc_gyaq_pa66_l2_nac_qmmm_20260723/amber_runtime"
export OMP_NUM_THREADS=1
test -s "$AMBERHOME/dat/slko/3ob-3-1/C-C.skf"

"$PY" "$DRIVER" --mode initialize --task-index "$INDEX" --output "$OUT" --scratch "$SCRATCH_ROOT" --code-root "$CODE_ROOT" --github-commit "$GITHUB_COMMIT"
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
import json,sys
print(json.load(open(sys.argv[1],encoding="utf-8"))["status"])
PY
)"
append_history TERMINAL "job=$ARRAY_JOB;task=$INDEX;seed=$SEED;replica=$REPLICA;status=$STATUS;automatic_downstream=NONE"
SUCCESS=1
printf '%s seed=%s replica=%s output=%s\n' "$STATUS" "$SEED" "$REPLICA" "$OUT"
