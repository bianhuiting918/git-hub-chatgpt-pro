#!/usr/bin/env bash
set -euo pipefail
: "${A1_BALANCED_BRIDGE_CODE_ROOT:?set immutable code root}"
: "${A1_BALANCED_BRIDGE_GITHUB_COMMIT:?set immutable GitHub commit}"
TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
CODE_ROOT=$A1_BALANCED_BRIDGE_CODE_ROOT
GITHUB_COMMIT=$A1_BALANCED_BRIDGE_GITHUB_COMMIT
INDEX=$SLURM_ARRAY_TASK_ID
ARRAY_JOB=$SLURM_ARRAY_JOB_ID
OUT_ROOT=$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_step1_balanced_four_distance_bridge
OUT=$OUT_ROOT/attempt_$ARRAY_JOB'_'$INDEX
SCRATCH=$SLURM_TMPDIR/nylc_a1_balanced_bridge_$ARRAY_JOB'_'$INDEX
DRIVER=$CODE_ROOT/scripts/prepare_audit_nylc_a1_step1_balanced_four_distance_bridge.py
PRMTOP=$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_unified_core_dftb3_preflight/attempt_62011285/prepared/system.prmtop
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
case "$INDEX" in
    0) SEED_INDEX=0; BRANCH_INDEX=0; SEED=seed26723; SHA=5d8f76d2c90e3e8c707b640c55a93938f53e18dc30d6f93d54adda467da25f41 ;;
    1) SEED_INDEX=0; BRANCH_INDEX=1; SEED=seed26723; SHA=5d8f76d2c90e3e8c707b640c55a93938f53e18dc30d6f93d54adda467da25f41 ;;
    2) SEED_INDEX=0; BRANCH_INDEX=2; SEED=seed26723; SHA=5d8f76d2c90e3e8c707b640c55a93938f53e18dc30d6f93d54adda467da25f41 ;;
    3) SEED_INDEX=1; BRANCH_INDEX=0; SEED=seed26737; SHA=4cc60ad4d7be099b3f76040f1ee8c49b91deecebb20ed98570bc192b3d511132 ;;
    4) SEED_INDEX=1; BRANCH_INDEX=1; SEED=seed26737; SHA=4cc60ad4d7be099b3f76040f1ee8c49b91deecebb20ed98570bc192b3d511132 ;;
    5) SEED_INDEX=1; BRANCH_INDEX=2; SEED=seed26737; SHA=4cc60ad4d7be099b3f76040f1ee8c49b91deecebb20ed98570bc192b3d511132 ;;
    *) exit 2 ;;
esac
if (( SEED_INDEX == 0 )); then SOURCE_RST7=$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_step1_acyl_release_md_continuation/attempt_62216380_0/release_md_endpoint.rst7
else SOURCE_RST7=$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_step1_acyl_release_md_continuation/attempt_62216380_1/release_md_endpoint.rst7; fi
test -s "$DRIVER"; test -s "$PRMTOP"; test -s "$SOURCE_RST7"
test "$(sha256sum "$SOURCE_RST7"|awk '{print $1}')" = "$SHA"; test ! -e "$OUT"
mkdir -p "$OUT_ROOT" "$SCRATCH"
START=$SCRATCH/start.rst7
"$PY" "$DRIVER" --mode initialize --seed-index "$SEED_INDEX" --branch-index "$BRANCH_INDEX" --root "$OUT" --start-rst7 "$START" --github-commit "$GITHUB_COMMIT"
module purge >/dev/null 2>&1 || true; module load amber/2018-hpcx-gcc-7.3.1
export AMBERHOME=$TASK_ROOT/../nylc_gyaq_pa66_l2_nac_qmmm_20260723/amber_runtime OMP_NUM_THREADS=1
test -s "$AMBERHOME/dat/slko/3ob-3-1/C-C.skf"
WIN=$OUT/branch_$BRANCH_INDEX; WSC=$SCRATCH/branch_$BRANCH_INDEX
"$PY" "$DRIVER" --mode prepare-window --root "$OUT" --output "$WIN" --scratch "$WSC" --input-rst7 "$START"
( cd "$WSC"; mpirun --bind-to none -np "$SLURM_NTASKS" sander.MPI -O -i stage.in -o stage.out -p "$PRMTOP" -c "$START" -ref "$START" -r stage.rst7 -inf stage.mdinfo )
"$PY" "$DRIVER" --mode audit-window --root "$OUT" --output "$WIN" --scratch "$WSC"
read -r GUARD HINT < <("$PY" - "$WIN/RESULT.json" <<'PY'
import json,sys
x=json.load(open(sys.argv[1])); print(int(x["guard_stop"]),int(x["restrained_hint"]))
PY
)
STOP=COMPLETED_ONE_BRANCH
if [[ "$GUARD" == 1 ]]; then STOP=GUARD_STOP
elif [[ "$HINT" == 1 ]]; then
 REL=$OUT/release_branch_$BRANCH_INDEX; RSC=$SCRATCH/release_branch_$BRANCH_INDEX
 "$PY" "$DRIVER" --mode prepare-release --root "$OUT" --output "$REL" --scratch "$RSC" --input-rst7 "$WSC/stage.rst7"
 ( cd "$RSC"; mpirun --bind-to none -np "$SLURM_NTASKS" sander.MPI -O -i stage.in -o stage.out -p "$PRMTOP" -c "$WSC/stage.rst7" -r stage.rst7 -x release.mdcrd -inf stage.mdinfo )
 "$PY" "$DRIVER" --mode audit-release --root "$OUT" --output "$REL" --scratch "$RSC"
 if [[ "$("$PY" - "$OUT/SOURCE_MANIFEST.json" <<'PY'
import json,sys
print(len(json.load(open(sys.argv[1]))["candidates"]))
PY
)" == 1 ]]; then STOP=CANDIDATE_LIMIT_REACHED; fi
fi
"$PY" "$DRIVER" --mode finalize --root "$OUT" --stop-reason "$STOP"
( flock -x 9; "$PY" "$DRIVER" --mode merge-if-ready --output-root "$OUT_ROOT" --array-job "$ARRAY_JOB" --seed-index "$SEED_INDEX"; ) 9>"$OUT_ROOT/.merge_seed_$SEED_INDEX.lock"
