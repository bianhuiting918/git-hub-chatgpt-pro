#!/usr/bin/env bash
set -euo pipefail
: "${A1_STEP2_WATER_RECRUIT_CODE_ROOT:?set immutable code root}"
: "${A1_STEP2_WATER_RECRUIT_GITHUB_COMMIT:?set immutable GitHub commit}"
TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
CODE_ROOT="$A1_STEP2_WATER_RECRUIT_CODE_ROOT"
GITHUB_COMMIT="$A1_STEP2_WATER_RECRUIT_GITHUB_COMMIT"
INDEX="${SLURM_ARRAY_TASK_ID:?run as array task 0..3}"
ARRAY_JOB="${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-manual}}"
if (( INDEX < 0 || INDEX > 3 )); then exit 2; fi
SEED_INDEX=$((INDEX / 2))
DONOR_BRANCH=$((INDEX % 2))
if (( SEED_INDEX == 0 )); then
    SEED=seed26723
    SOURCE_RST7="$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_step1_acyl_release_md_continuation/attempt_62216380_0/release_md_endpoint.rst7"
    EXPECTED_SHA=5d8f76d2c90e3e8c707b640c55a93938f53e18dc30d6f93d54adda467da25f41
else
    SEED=seed26737
    SOURCE_RST7="$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_step1_acyl_release_md_continuation/attempt_62216380_1/release_md_endpoint.rst7"
    EXPECTED_SHA=4cc60ad4d7be099b3f76040f1ee8c49b91deecebb20ed98570bc192b3d511132
fi
OUT_ROOT="$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_step2_water_recruitment"
OUT="$OUT_ROOT/attempt_${ARRAY_JOB}_${INDEX}"
SCRATCH_ROOT="${SLURM_TMPDIR:-/tmp}/nylc_a1_step2_water_recruit_${ARRAY_JOB}_${INDEX}"
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
DRIVER="$CODE_ROOT/scripts/prepare_audit_nylc_a1_step2_water_recruitment.py"
PRMTOP="$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_unified_core_dftb3_preflight/attempt_62011285/prepared/system.prmtop"

append_history() {
    local state="$1" detail="$2" now
    now="$(date '+%Y-%m-%dT%H:%M:%S%z')"
    (
        flock -x 9
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$now" "nylc_a1_step2_water_recruitment" "run_nylc_a1_step2_water_recruitment.sh" "$state" "$GITHUB_COMMIT" "$detail" >> "$TASK_ROOT/run_history.tsv"
        printf '{"time":"%s","event":"nylc_a1_step2_water_recruitment","state":"%s","git_commit":"%s","detail":"%s","output":"%s"}\n' "$now" "$state" "$GITHUB_COMMIT" "$detail" "$OUT" >> "$TASK_ROOT/run_history.jsonl"
    ) 9>> "$TASK_ROOT/.run_history.lock"
}

test -s "$DRIVER"
test -s "$PRMTOP"
test -s "$SOURCE_RST7"
test "$(sha256sum "$SOURCE_RST7" | awk '{print $1}')" = "$EXPECTED_SHA"
test ! -e "$OUT"
append_history STARTED "job=$ARRAY_JOB;task=$INDEX;seed=$SEED;donor_branch=$DONOR_BRANCH;source_sha=$EXPECTED_SHA"

module purge >/dev/null 2>&1 || true
module load amber/2018-hpcx-gcc-7.3.1
export AMBERHOME="$TASK_ROOT/../nylc_gyaq_pa66_l2_nac_qmmm_20260723/amber_runtime"
export OMP_NUM_THREADS=1
test -s "$AMBERHOME/dat/slko/3ob-3-1/C-C.skf"

"$PY" "$DRIVER" --mode initialize --task-index "$INDEX" --output "$OUT" --scratch "$SCRATCH_ROOT" --code-root "$CODE_ROOT" --github-commit "$GITHUB_COMMIT"
GUIDED="$SCRATCH_ROOT/guided"
(
    cd "$GUIDED"
    set +e
    mpirun --bind-to none -np 8 sander.MPI -O -i stage.in -o stage.out -p "$PRMTOP" -c "$SOURCE_RST7" -ref "$SOURCE_RST7" -r stage.rst7 -inf stage.mdinfo
    printf '%s\n' "$?" > engine.rc
    exit 0
)
"$PY" "$DRIVER" --mode audit-guided --output "$OUT" --scratch "$SCRATCH_ROOT"
GUIDED_PASS="$("$PY" - "$OUT/GUIDED_RESULT.json" <<'PY'
import json,sys
print(int(json.load(open(sys.argv[1],encoding="utf-8"))["guided_attack_pose_pass"]))
PY
)"
if [[ "$GUIDED_PASS" = 1 ]]; then
    declare -a PIDS=()
    for LEG in 0 1; do
        STAGE="$SCRATCH_ROOT/a2_leg$LEG"
        (
            cd "$STAGE"
            set +e
            mpirun --bind-to none -np 8 sander.MPI -O -i stage.in -o stage.out -p "$PRMTOP" -c "$GUIDED/stage.rst7" -r stage.rst7 -x a2.mdcrd -inf stage.mdinfo
            printf '%s\n' "$?" > engine.rc
            exit 0
        ) &
        PIDS[$LEG]=$!
    done
    for LEG in 0 1; do wait "${PIDS[$LEG]}"; done
    for LEG in 0 1; do
        "$PY" "$DRIVER" --mode audit-leg --output "$OUT" --scratch "$SCRATCH_ROOT" --leg "$LEG"
        TECH="$("$PY" - "$OUT/A2_LEG_$LEG.json" <<'PY'
import json,sys
print(int(json.load(open(sys.argv[1],encoding="utf-8"))["technical_pass"]))
PY
)"
        if [[ "$TECH" = 1 ]]; then cp "$SCRATCH_ROOT/a2_leg$LEG/stage.rst7" "$OUT/a2_leg${LEG}_endpoint.rst7"; fi
    done
fi
"$PY" "$DRIVER" --mode finalize --output "$OUT"
(
    flock -x 9
    "$PY" "$DRIVER" --mode merge-if-ready --output-root "$OUT_ROOT" --array-job "$ARRAY_JOB" --seed-index "$SEED_INDEX"
) 9>> "$OUT_ROOT/.merge_seed_$SEED_INDEX.lock"
STATUS="$("$PY" - "$OUT/RESULT.json" <<'PY'
import json,sys
print(json.load(open(sys.argv[1],encoding="utf-8"))["status"])
PY
)"
append_history TERMINAL "job=$ARRAY_JOB;task=$INDEX;seed=$SEED;donor_branch=$DONOR_BRANCH;status=$STATUS;automatic_downstream=NONE"
printf '%s seed=%s donor_branch=%s output=%s\n' "$STATUS" "$SEED" "$DONOR_BRANCH" "$OUT"
