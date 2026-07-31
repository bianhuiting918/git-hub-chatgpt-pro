#!/usr/bin/env bash
set -euo pipefail

: "${A1_ANGLE_CODE_ROOT:?set immutable code root}"
: "${A1_ANGLE_GITHUB_COMMIT:?set immutable GitHub commit}"

TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
CODE_ROOT="$A1_ANGLE_CODE_ROOT"
COMMIT="$A1_ANGLE_GITHUB_COMMIT"
INDEX="${SLURM_ARRAY_TASK_ID:?run as array task 0..1}"
ARRAY_JOB="${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-manual}}"
if (( INDEX < 0 || INDEX > 1 )); then exit 2; fi

OUT_ROOT="$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_step1_angle_protected_attack"
OUT="$OUT_ROOT/attempt_${ARRAY_JOB}_${INDEX}"
SCRATCH="${SLURM_TMPDIR:-/tmp}/nylc_a1_angle_${ARRAY_JOB}_${INDEX}"
DRIVER="$CODE_ROOT/scripts/prepare_audit_nylc_a1_step1_angle_protected_attack.py"
PRMTOP="$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_unified_core_dftb3_preflight/attempt_62011285/prepared/system.prmtop"
EXPECTED_PRMTOP_SHA=a61d15bf0bf78675be93275d45f274e808ed6ae450fc1ca21a8e14aee8c12ca0
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python

test -s "$DRIVER"
test -s "$PRMTOP"
test "$(sha256sum "$PRMTOP" | awk '{print $1}')" = "$EXPECTED_PRMTOP_SHA"
test ! -e "$OUT"
test ! -e "$SCRATCH"

module purge >/dev/null 2>&1 || true
module load amber/2018-hpcx-gcc-7.3.1
export AMBERHOME="$TASK_ROOT/../nylc_gyaq_pa66_l2_nac_qmmm_20260723/amber_runtime"
export OMP_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
test -s "$AMBERHOME/dat/slko/3ob-3-1/C-C.skf"

"$PY" "$DRIVER" --mode initialize --task-index "$INDEX" --root "$OUT" --github-commit "$COMMIT"
"$PY" "$DRIVER" --mode prepare --root "$OUT" --scratch "$SCRATCH"
SOURCE_RST7=$("$PY" -c 'import json,sys; print(json.load(open(sys.argv[1]))["source"]["restart"])' "$OUT/MANIFEST.json")
test -s "$SOURCE_RST7"

set +e
(
    cd "$SCRATCH"
    mpirun --bind-to none -np 8 sander.MPI -O -i stage.in -o stage.out -p "$PRMTOP" -c "$SOURCE_RST7" -ref "$SOURCE_RST7" -r stage.rst7 -inf stage.mdinfo
)
ENGINE_RC=$?
set -e

printf '%s\n' "$ENGINE_RC" > "$SCRATCH/engine.rc"
for name in stage.in restraints.RST stage.out stage.mdinfo engine.rc; do
    if [[ -f "$SCRATCH/$name" ]]; then cp -p "$SCRATCH/$name" "$OUT/$name"; fi
done
"$PY" "$DRIVER" --mode audit --root "$OUT" --scratch "$SCRATCH"
STATUS=$("$PY" -c 'import json,sys; print(json.load(open(sys.argv[1]))["status"])' "$OUT/RESULT.json")
GATE=$("$PY" -c 'import json,sys; print(json.load(open(sys.argv[1]))["scientific_gate"])' "$OUT/RESULT.json")
if (( ENGINE_RC != 0 )); then exit "$ENGINE_RC"; fi
printf '%s task=%s gate=%s output=%s\n' "$STATUS" "$INDEX" "$GATE" "$OUT"
