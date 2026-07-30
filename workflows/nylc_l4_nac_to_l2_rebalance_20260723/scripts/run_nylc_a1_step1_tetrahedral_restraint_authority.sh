#!/usr/bin/env bash
set -euo pipefail

: "$A1_TETRA_AUTH_CODE_ROOT"
: "$A1_TETRA_AUTH_GITHUB_COMMIT"

TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
CODE_ROOT="$A1_TETRA_AUTH_CODE_ROOT"
COMMIT="$A1_TETRA_AUTH_GITHUB_COMMIT"
JOB_ID="$SLURM_JOB_ID"
OUT_ROOT="$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_step1_tetrahedral_restraint_authority"
OUT="$OUT_ROOT/attempt_"$JOB_ID"_task1"
set +u
SCRATCH_BASE="$SLURM_TMPDIR"
set -u
if [[ -z "$SCRATCH_BASE" ]]; then SCRATCH_BASE=/tmp; fi
SCRATCH="$SCRATCH_BASE/nylc_a1_tetra_auth_"$JOB_ID
DRIVER="$CODE_ROOT/scripts/prepare_audit_nylc_a1_step1_tetrahedral_bridge.py"
PRMTOP="$TASK_ROOT/a1_activated_nac_20260726/qmmm/a1_unified_core_dftb3_preflight/attempt_62011285/prepared/system.prmtop"
EXPECTED_PRMTOP_SHA=a61d15bf0bf78675be93275d45f274e808ed6ae450fc1ca21a8e14aee8c12ca0
EXPECTED_SOURCE_SHA=1ed5a5b6a72a43f9f56e41bf427bb95ffc96237c86811a34aec50fa3c89668f1
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
SUCCESS=0
ENGINE_RC=999
AUDIT_RC=999

append_history() {
    local state="$1" detail="$2" now
    now="$(date '+%Y-%m-%dT%H:%M:%S%z')"
    (
        flock -x 9
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$now" "nylc_a1_step1_tetrahedral_restraint_authority" "run_nylc_a1_step1_tetrahedral_restraint_authority.sh" "$state" "$COMMIT" "$detail" >> "$TASK_ROOT/run_history.tsv"
        printf '{"time":"%s","event":"nylc_a1_step1_tetrahedral_restraint_authority","state":"%s","git_commit":"%s","detail":"%s","output":"%s"}\n' "$now" "$state" "$COMMIT" "$detail" "$OUT" >> "$TASK_ROOT/run_history.jsonl"
    ) 9>> "$TASK_ROOT/.run_history.lock"
}

persist_evidence() {
    if [[ ! -d "$OUT" ]]; then return 0; fi
    mkdir -p "$OUT/evidence"
    for name in stage.in restraints.RST stage.out stage.mdinfo restraint.dat stage.rst7; do
        if [[ -f "$SCRATCH/$name" ]]; then cp -f "$SCRATCH/$name" "$OUT/evidence/$name"; fi
    done
    printf '%s\n' "$ENGINE_RC" > "$OUT/evidence/engine.rc"
    printf '%s\n' "$AUDIT_RC" > "$OUT/evidence/audit.rc"
}

on_exit() {
    local rc=$?
    persist_evidence || true
    if (( SUCCESS == 0 )); then
        append_history FAILED "job=$JOB_ID;task=1;rc=$rc;engine_rc=$ENGINE_RC;audit_rc=$AUDIT_RC" || true
    fi
}
trap on_exit EXIT

test -s "$DRIVER"
test -s "$PRMTOP"
test "$(sha256sum "$PRMTOP" | awk '{print $1}')" = "$EXPECTED_PRMTOP_SHA"
test ! -e "$OUT"
test ! -e "$SCRATCH"
mkdir -p "$OUT_ROOT"
append_history STARTED "job=$JOB_ID;task=1;schedule=PT_ASSISTED;window=0;diagnostic_only=1"

module purge >/dev/null 2>&1 || true
module load amber/2018-hpcx-gcc-7.3.1
export AMBERHOME="$TASK_ROOT/../nylc_gyaq_pa66_l2_nac_qmmm_20260723/amber_runtime"
export OMP_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
test -s "$AMBERHOME/dat/slko/3ob-3-1/C-C.skf"

"$PY" "$DRIVER" --mode initialize --task-index 1 --root "$OUT" --github-commit "$COMMIT"
ACCEPTED="$OUT/accepted.rst7"
test -s "$ACCEPTED"
test "$(sha256sum "$ACCEPTED" | awk '{print $1}')" = "$EXPECTED_SOURCE_SHA"
PERSIST="$OUT/window_00"
"$PY" "$DRIVER" --mode prepare-window --root "$OUT" --output "$PERSIST" --scratch "$SCRATCH" --input-rst7 "$ACCEPTED" --window-index 0
mkdir -p "$OUT/evidence"
cp "$SCRATCH/stage.in" "$OUT/evidence/stage.in"
cp "$SCRATCH/restraints.RST" "$OUT/evidence/restraints.RST"
cp "$PERSIST/WINDOW_MANIFEST.json" "$OUT/evidence/WINDOW_MANIFEST.json"

set +e
(
    cd "$SCRATCH"
    mpirun --bind-to none -np 8 sander.MPI -O -i stage.in -o stage.out -p "$PRMTOP" -c "$ACCEPTED" -ref "$ACCEPTED" -r stage.rst7 -inf stage.mdinfo
)
ENGINE_RC=$?
persist_evidence
"$PY" "$DRIVER" --mode audit-window --root "$OUT" --output "$PERSIST" --scratch "$SCRATCH"
AUDIT_RC=$?
set -e
persist_evidence

"$PY" - "$OUT" "$ENGINE_RC" "$AUDIT_RC" <<'PY'
import hashlib
import json
import pathlib
import sys

root = pathlib.Path(sys.argv[1])
engine_rc = int(sys.argv[2])
audit_rc = int(sys.argv[3])
ev = root / "evidence"
rst = (ev / "restraints.RST").read_text(encoding="utf-8")
stage_in = (ev / "stage.in").read_text(encoding="utf-8")
stage_out = (ev / "stage.out").read_text(encoding="utf-8", errors="replace") if (ev / "stage.out").is_file() else ""
tokens = ("NMR", "restraint", "RESTRAINT", "DFTB", "QM/MM", "QMMM", "nquant", "qmmm_int")
engine_authority_lines = [line for line in stage_out.splitlines() if any(token in line for token in tokens)][:500]
files = {}
for path in sorted(ev.iterdir()):
    if path.is_file():
        files[path.name] = {
            "bytes": path.stat().st_size,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        }
payload = {
    "schema_version": 1,
    "status": "PASS_EVIDENCE_CAPTURED" if engine_rc == 0 and audit_rc == 0 else "NOT_EVALUATED_TECHNICAL_RESTRAINT_AUTHORITY",
    "scientific_status": "NOT_EVALUATED_RESTRAINT_AUTHORITY_DIAGNOSTIC_ONLY",
    "github_commit": json.loads((root / "MANIFEST.json").read_text(encoding="utf-8"))["github_commit"],
    "source_task_index": 1,
    "schedule": "PT_ASSISTED",
    "window_index": 0,
    "engine_rc": engine_rc,
    "audit_rc": audit_rc,
    "input_contract": {
        "physical_rst_lines": len([line for line in rst.splitlines() if line.strip().startswith("&rst")]),
        "expected_rst_lines": 5,
        "nmropt_present": "nmropt=1" in stage_in.replace(" ", "").lower(),
        "disang_present": "DISANG=restraints.RST" in stage_in,
        "dumpave_present": "DUMPAVE=restraint.dat" in stage_in,
    },
    "engine_authority_lines": engine_authority_lines,
    "files": files,
    "automatic_downstream_action": "NONE",
}
(root / "DIAGNOSTIC_RESULT.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
PY

(
    cd "$OUT"
    find . -type f ! -name SHA256.tsv -print0 | sort -z | xargs -0 sha256sum > SHA256.tsv
)
test -s "$OUT/DIAGNOSTIC_RESULT.json"
test -s "$OUT/SHA256.tsv"
append_history TERMINAL "job=$JOB_ID;task=1;engine_rc=$ENGINE_RC;audit_rc=$AUDIT_RC;diagnostic_only=1;automatic_downstream=NONE"
if (( ENGINE_RC != 0 || AUDIT_RC != 0 )); then exit 1; fi
SUCCESS=1
printf 'PASS_EVIDENCE_CAPTURED job=%s output=%s\n' "$JOB_ID" "$OUT"
