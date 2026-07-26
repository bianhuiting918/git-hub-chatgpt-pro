#!/usr/bin/env bash
set -euo pipefail

TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
TOOLS_ROOT="$TASK_ROOT/tools"
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
AMBERCLASSIC_COMMIT=0b35bfeb96026ffa4e5876391a0828f39b3cfc8d
SOURCE_ROOT="$TOOLS_ROOT/source_archives"
FINAL_ARCHIVE="$SOURCE_ROOT/AmberClassic-${AMBERCLASSIC_COMMIT}.tar.gz"
URL="https://codeload.github.com/Amber-MD/AmberClassic/tar.gz/${AMBERCLASSIC_COMMIT}"
STAMP="$(date '+%Y%m%dT%H%M%S%z')"
AUDIT_DIR="$SOURCE_ROOT/audits/stage_${STAMP}_$$"
TMP=""
EVENT=nylc_stage_amberclassic_source
COMMAND="stage_nylc_amberclassic_source.sh commit=$AMBERCLASSIC_COMMIT"
STATE=FAIL_TECHNICAL
DETAIL=initializing

mkdir -p "$AUDIT_DIR"

append_history() {
  local now
  now="$(date '+%Y-%m-%dT%H:%M:%S%z')"
  (
    flock -x 9
    printf '%s\t%s\t%s\t%s\t%s\n' "$now" "$EVENT" "$COMMAND" "$STATE" "$DETAIL" >>"$TASK_ROOT/run_history.tsv"
    "$PY" - "$now" "$STATE" "$DETAIL" "$FINAL_ARCHIVE" "$AMBERCLASSIC_COMMIT" "$AUDIT_DIR" >>"$TASK_ROOT/run_history.jsonl" <<'PY'
import json
import sys
now, state, detail, archive, commit, audit_dir = sys.argv[1:]
print(json.dumps({
    "time": now,
    "event": "nylc_stage_amberclassic_source",
    "command": "stage_nylc_amberclassic_source.sh",
    "state": state,
    "detail": detail,
    "source_archive": archive,
    "source_commit": commit,
    "audit_dir": audit_dir,
}, sort_keys=True))
PY
  ) 9>"$TASK_ROOT/.run_history.lock"
}

finish() {
  code=$?
  trap - EXIT
  if [[ -n "$TMP" && -f "$TMP" ]]; then
    rm -f -- "$TMP"
  fi
  if [[ $code -ne 0 ]]; then
    DETAIL="exit_code=$code; audit_dir=$AUDIT_DIR"
  fi
  append_history
  exit "$code"
}
trap finish EXIT

if [[ ! -f "$FINAL_ARCHIVE" ]]; then
  TMP="$(mktemp "$AUDIT_DIR/AmberClassic-${AMBERCLASSIC_COMMIT}.tar.gz.tmp.XXXXXX")"
  curl -4 --retry 8 --retry-delay 5 -fL \
    --connect-timeout 30 --max-time 1800 "$URL" -o "$TMP"
  tar -tzf "$TMP" >"$AUDIT_DIR/archive_members.txt"
  mv -n "$TMP" "$FINAL_ARCHIVE"
  if [[ -f "$TMP" ]]; then
    rm -f -- "$TMP"
  fi
  TMP=""
else
  tar -tzf "$FINAL_ARCHIVE" >"$AUDIT_DIR/archive_members.txt"
fi

SOURCE_SHA256="$(sha256sum "$FINAL_ARCHIVE" | awk '{print $1}')"
SOURCE_BYTES="$(stat -c '%s' "$FINAL_ARCHIVE")"
printf '%s  %s\n' "$SOURCE_SHA256" "$(basename "$FINAL_ARCHIVE")" >"$AUDIT_DIR/source.sha256"

"$PY" - "$AUDIT_DIR/PASS_SOURCE.json" "$FINAL_ARCHIVE" "$AMBERCLASSIC_COMMIT" "$SOURCE_SHA256" "$SOURCE_BYTES" <<'PY'
import json
import pathlib
import sys
out = pathlib.Path(sys.argv[1])
archive = pathlib.Path(sys.argv[2])
commit, sha256, size = sys.argv[3:]
record = {
    "schema_version": 1,
    "status": "PASS_AMBERCLASSIC_SOURCE_STAGED",
    "source_repository": "Amber-MD/AmberClassic",
    "source_commit": commit,
    "source_archive": str(archive),
    "source_archive_sha256": sha256,
    "source_archive_bytes": int(size),
}
out.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8")
PY

STATE=PASS_TECHNICAL
DETAIL="status=PASS_AMBERCLASSIC_SOURCE_STAGED; source_sha256=$SOURCE_SHA256; bytes=$SOURCE_BYTES; archive=$FINAL_ARCHIVE"
