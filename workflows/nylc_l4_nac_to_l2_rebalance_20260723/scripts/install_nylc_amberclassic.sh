#!/usr/bin/env bash
set -euo pipefail

TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
TOOLS_ROOT="$TASK_ROOT/tools"
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
AMBERCLASSIC_COMMIT=0b35bfeb96026ffa4e5876391a0828f39b3cfc8d
JOB_ID="${SLURM_JOB_ID:?SLURM_JOB_ID is required}"
CPUS="${SLURM_CPUS_PER_TASK:-8}"
PREFIX="$TOOLS_ROOT/amberclassic_job_${JOB_ID}_${AMBERCLASSIC_COMMIT:0:12}"
STAGED_ARCHIVE="$TOOLS_ROOT/source_archives/AmberClassic-${AMBERCLASSIC_COMMIT}.tar.gz"
ACTIVE="$TOOLS_ROOT/ACTIVE_AMBERCLASSIC.json"
EVENT=nylc_install_amberclassic
COMMAND="install_nylc_amberclassic.sh commit=$AMBERCLASSIC_COMMIT prefix=$PREFIX staged_archive=$STAGED_ARCHIVE"
STATE=FAIL_TECHNICAL
DETAIL=initializing

mkdir -p "$TOOLS_ROOT"
if [[ -e "$PREFIX" ]]; then
  printf 'Refusing to overwrite installation attempt: %s\n' "$PREFIX" >&2
  exit 2
fi
if [[ ! -f "$STAGED_ARCHIVE" ]]; then
  printf 'Required staged source archive is missing: %s\n' "$STAGED_ARCHIVE" >&2
  exit 3
fi
mkdir -p "$PREFIX"

append_history() {
  local now
  now="$(date '+%Y-%m-%dT%H:%M:%S%z')"
  (
    flock -x 9
    printf '%s\t%s\t%s\t%s\t%s\n' "$now" "$EVENT" "$COMMAND" "$STATE" "$DETAIL" >>"$TASK_ROOT/run_history.tsv"
    "$PY" - "$now" "$STATE" "$DETAIL" "$JOB_ID" "$PREFIX" "$AMBERCLASSIC_COMMIT" "$STAGED_ARCHIVE" >>"$TASK_ROOT/run_history.jsonl" <<'PY'
import json
import sys
now, state, detail, job, prefix, commit, staged_archive = sys.argv[1:]
print(json.dumps({
    "time": now,
    "event": "nylc_install_amberclassic",
    "command": "install_nylc_amberclassic.sh",
    "state": state,
    "detail": detail,
    "slurm_job_id": job,
    "prefix": prefix,
    "source_commit": commit,
    "staged_archive": staged_archive,
}, sort_keys=True))
PY
  ) 9>"$TASK_ROOT/.run_history.lock"
}
finish() {
  code=$?
  trap - EXIT
  if [[ $code -ne 0 ]]; then
    DETAIL="exit_code=$code; prefix=$PREFIX; staged_archive=$STAGED_ARCHIVE"
  fi
  append_history
  exit "$code"
}
trap finish EXIT

module purge
module load compiler/gcc/11.4.0
gcc --version >"$PREFIX/compiler_versions.txt"
gfortran --version >>"$PREFIX/compiler_versions.txt"

tar -tzf "$STAGED_ARCHIVE" >"$PREFIX/archive_members.txt"
SOURCE_SHA256="$(sha256sum "$STAGED_ARCHIVE" | awk '{print $1}')"
printf '%s  %s\n' "$SOURCE_SHA256" "$(basename "$STAGED_ARCHIVE")" >"$PREFIX/source.sha256"
tar -xzf "$STAGED_ARCHIVE" --strip-components=1 -C "$PREFIX"

cd "$PREFIX"
./configure --noopenmp --noboost >configure.stdout 2>configure.stderr
"$PY" - "$PREFIX/src" "$PREFIX" "$PREFIX/dependency_relocation.json" <<'PY'
import hashlib
import json
import pathlib
import sys
source_root = pathlib.Path(sys.argv[1])
prefix = sys.argv[2]
output = pathlib.Path(sys.argv[3])
stale = "/home/case/AmberClassic"
changes = []
for path in sorted(source_root.rglob("depend")):
    text = path.read_text(encoding="utf-8", errors="strict")
    count = text.count(stale)
    if not count:
        continue
    before = hashlib.sha256(text.encode("utf-8")).hexdigest()
    updated = text.replace(stale, prefix)
    path.write_text(updated, encoding="utf-8")
    changes.append({
        "path": str(path.relative_to(source_root.parent)),
        "occurrences": count,
        "before_sha256": before,
        "after_sha256": hashlib.sha256(updated.encode("utf-8")).hexdigest(),
    })
if not changes:
    raise SystemExit("expected at least one stale AmberClassic dependency path")
record = {
    "schema_version": 1,
    "status": "PASS_DEPENDENCY_PATH_RELOCATION",
    "stale_prefix": stale,
    "replacement_prefix": prefix,
    "files_changed": changes,
    "remaining_stale_occurrences": sum(
        path.read_text(encoding="utf-8", errors="strict").count(stale)
        for path in source_root.rglob("depend")
    ),
}
if record["remaining_stale_occurrences"] != 0:
    raise SystemExit("stale AmberClassic dependency paths remain")
output.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8")
PY
make -j 1 install >make_install.stdout 2>make_install.stderr

test -f "$PREFIX/dat/antechamber/CONNECT.TPL"
test -x "$PREFIX/bin/antechamber"
test -x "$PREFIX/bin/parmchk2"
test -x "$PREFIX/bin/tleap"
test -x "$PREFIX/bin/sqm"
test -x "$PREFIX/bin/msander"
set +u
source "$PREFIX/AmberClassic.sh"
set -u
antechamber -h >antechamber_help.stdout 2>antechamber_help.stderr

"$PY" - "$PREFIX" "$ACTIVE" "$AMBERCLASSIC_COMMIT" "$SOURCE_SHA256" "$JOB_ID" "$STAGED_ARCHIVE" <<'PY'
import hashlib
import json
import pathlib
import sys
prefix = pathlib.Path(sys.argv[1])
active = pathlib.Path(sys.argv[2])
commit, source_sha, job, staged_archive = sys.argv[3:]
required = [
    prefix / "dependency_relocation.json",
    prefix / "dat" / "antechamber" / "CONNECT.TPL",
    prefix / "bin" / "antechamber",
    prefix / "bin" / "parmchk2",
    prefix / "bin" / "tleap",
    prefix / "bin" / "sqm",
    prefix / "bin" / "msander",
]
record = {
    "schema_version": 1,
    "status": "PASS_AMBERCLASSIC_INSTALL",
    "source_repository": "Amber-MD/AmberClassic",
    "source_commit": commit,
    "source_archive": staged_archive,
    "source_archive_sha256": source_sha,
    "slurm_job_id": int(job),
    "prefix": str(prefix),
    "required_files": {
        str(path.relative_to(prefix)): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in required
    },
}
(prefix / "PASS.json").write_text(
    json.dumps(record, indent=2, sort_keys=True) + "\n",
    encoding="utf-8",
)
active.write_text(
    json.dumps(record, indent=2, sort_keys=True) + "\n",
    encoding="utf-8",
)
PY

STATE=PASS_TECHNICAL
DETAIL="status=PASS_AMBERCLASSIC_INSTALL; source_sha256=$SOURCE_SHA256; prefix=$PREFIX; PASS.json=$PREFIX/PASS.json"
