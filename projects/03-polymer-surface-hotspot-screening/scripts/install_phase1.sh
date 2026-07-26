#!/usr/bin/env bash
set -euo pipefail

ROOT=${ROOT:-/work/home/acshdt1dks/polymer_surface_hotspot_screen_20260725}
MAX_JOBS=${MAX_JOBS:-8}
PYTHON_BOOTSTRAP=${PYTHON_BOOTSTRAP:-/work/home/acshdt1dks/python3.11/bin/python3.11}
PROJECT_DIR=${PROJECT_DIR:-/work/home/acshdt1dks/polymer_surface_hotspot_screen_20260725/repo/projects/03-polymer-surface-hotspot-screening}
LOCK="$PROJECT_DIR/config/analysis-requirements.lock"
ADFR_BIN="$ROOT/software/ADFRsuite-1.0/bin"
ENV_DIR="$ROOT/envs/surface-screen-py311"
PIP_CACHE="$ROOT/cache/pip"
LOG_DIR="$ROOT/logs/install_phase1"
GATE_DIR="$ROOT/results/gates"
PASS_GATE="$GATE_DIR/INSTALL_PASS.json"
FAIL_GATE="$GATE_DIR/INSTALL_FAIL.json"

ROOT_REAL=$(readlink -f "$ROOT")
case "$ROOT_REAL/" in
  /work/home/acshdt1dks/*) ;;
  *)
    printf 'refusing unauthorized project root: %s\n' "$ROOT_REAL" >&2
    exit 2
    ;;
esac
if [ "$MAX_JOBS" -lt 1 ] || [ "$MAX_JOBS" -gt 64 ]; then
  printf 'MAX_JOBS must be between 1 and 64: %s\n' "$MAX_JOBS" >&2
  exit 2
fi
for guarded in "$ROOT/envs" "$ROOT/cache" "$ROOT/results" "$ROOT/logs"; do
  if [ -L "$guarded" ]; then
    guarded_real=$(readlink -f "$guarded")
    case "$guarded_real/" in
      "$ROOT_REAL"/*) ;;
      *)
        printf 'refusing symlink escaping project root: %s -> %s\n' "$guarded" "$guarded_real" >&2
        exit 2
        ;;
    esac
  fi
done

mkdir -p "$ROOT/envs" "$PIP_CACHE" "$LOG_DIR" "$GATE_DIR"
rm -f "$FAIL_GATE"
on_failure() {
  status=$?
  if [ "$status" -ne 0 ]; then
    printf '{"status":"INSTALL_FAIL","exit_code":%s,"timestamp":"%s"}\n' \
      "$status" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$FAIL_GATE"
  fi
  exit "$status"
}
trap on_failure EXIT

if [ ! -x "$PYTHON_BOOTSTRAP" ]; then
  printf 'bootstrap Python is not executable: %s\n' "$PYTHON_BOOTSTRAP" >&2
  exit 2
fi
if [ ! -d "$PROJECT_DIR" ]; then
  printf 'repository project directory is missing: %s\n' "$PROJECT_DIR" >&2
  exit 2
fi
PROJECT_REAL=$(readlink -f "$PROJECT_DIR")
case "$PROJECT_REAL/" in
  "$ROOT_REAL/repo/"*) ;;
  *)
    printf 'refusing project directory outside pinned repository: %s\n' "$PROJECT_REAL" >&2
    exit 2
    ;;
esac
if [ ! -f "$LOCK" ]; then
  printf 'dependency lock is missing: %s\n' "$LOCK" >&2
  exit 2
fi
for executable in autosite autogrid4 prepare_receptor; do
  if [ ! -x "$ADFR_BIN/$executable" ]; then
    printf 'required ADFRsuite executable is missing: %s\n' "$ADFR_BIN/$executable" >&2
    exit 2
  fi
done
if [ -e "$PASS_GATE" ]; then
  printf 'existing installation PASS retained without overwrite: %s\n' "$PASS_GATE"
  trap - EXIT
  exit 0
fi

if [ ! -x "$ENV_DIR/bin/python" ]; then
  "$PYTHON_BOOTSTRAP" -m venv "$ENV_DIR"
fi
"$ENV_DIR/bin/python" -m pip install --upgrade "pip==24.3.1" \
  --cache-dir "$PIP_CACHE" --no-input \
  >"$LOG_DIR/pip_bootstrap.stdout.log" 2>"$LOG_DIR/pip_bootstrap.stderr.log"
"$ENV_DIR/bin/python" -m pip install --no-input --cache-dir "$PIP_CACHE" \
  -r "$LOCK" \
  >"$LOG_DIR/pip_install.stdout.log" 2>"$LOG_DIR/pip_install.stderr.log"
"$ENV_DIR/bin/python" -m pip freeze > "$LOG_DIR/pip_freeze.txt"

"$ADFR_BIN/autosite" --version >"$LOG_DIR/autosite.version.log" 2>&1
"$ADFR_BIN/autogrid4" -h >"$LOG_DIR/autogrid4.version.log" 2>&1
"$ADFR_BIN/prepare_receptor" -h >"$LOG_DIR/prepare_receptor.version.log" 2>&1
"$ENV_DIR/bin/freesasa" --version >"$LOG_DIR/freesasa.version.log" 2>&1
"$ENV_DIR/bin/python" -c "import numpy, scipy, pandas, Bio, yaml" \
  >"$LOG_DIR/core_imports.log" 2>&1
"$ENV_DIR/bin/python" -c "import rdkit" >"$LOG_DIR/rdkit_import.log" 2>&1
"$ENV_DIR/bin/python" -c "import freesasa" >"$LOG_DIR/freesasa_import.log" 2>&1
"$ENV_DIR/bin/python" -m pytest --version >"$LOG_DIR/pytest.version.log" 2>&1

lock_sha=$(sha256sum "$LOCK" | awk '{print $1}')
adfr_archive_sha=$(sha256sum "$ROOT/software/downloads/ADFRsuite_x86_64Linux_1.0.tar.gz" | awk '{print $1}')
freeze_sha=$(sha256sum "$LOG_DIR/pip_freeze.txt" | awk '{print $1}')
"$ENV_DIR/bin/python" - "$PASS_GATE.tmp" "$lock_sha" "$adfr_archive_sha" "$freeze_sha" <<'PY'
import json
import pathlib
import sys
from datetime import datetime, timezone

path = pathlib.Path(sys.argv[1])
payload = {
    "status": "INSTALL_PASS",
    "timestamp_utc": datetime.now(timezone.utc).isoformat(),
    "analysis_lock_sha256": sys.argv[2],
    "adfrsuite_archive_sha256": sys.argv[3],
    "pip_freeze_sha256": sys.argv[4],
    "analysis_environment": "/work/home/acshdt1dks/polymer_surface_hotspot_screen_20260725/envs/surface-screen-py311",
    "adfrsuite_bin": "/work/home/acshdt1dks/polymer_surface_hotspot_screen_20260725/software/ADFRsuite-1.0/bin",
    "scientific_scope": "installation_only",
    "msms_status": "NOT_EVALUATED_OPTIONAL_FOR_CURRENT_AUTOSITE_MAP_STAGE",
}
path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
PY
mv "$PASS_GATE.tmp" "$GATE_DIR/INSTALL_PASS.json"
sha256sum "$GATE_DIR/INSTALL_PASS.json" "$LOG_DIR/pip_freeze.txt" "$LOCK" > "$GATE_DIR/INSTALL_SHA256SUMS"
rm -f "$FAIL_GATE"
trap - EXIT
printf 'INSTALL_PASS.json written: %s\n' "$GATE_DIR/INSTALL_PASS.json"
