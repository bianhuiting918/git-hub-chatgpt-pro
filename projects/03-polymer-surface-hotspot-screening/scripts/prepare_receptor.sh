#!/usr/bin/env bash
set -euo pipefail

ROOT=${ROOT:-/work/home/acshdt1dks/polymer_surface_hotspot_screen_20260725}
ADFR_BIN=${ADFR_BIN:-$ROOT/software/ADFRsuite-1.0/bin}
ANALYSIS_PYTHON=${ANALYSIS_PYTHON:-$ROOT/envs/surface-screen-py311/bin/python}
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
CHECKER="$SCRIPT_DIR/check_receptor.py"

INPUT=
OUTPUT_DIR=
RECORD_ID=
CHAINS=
MODEL=1
CATALYTIC_RESIDUES=
RETAIN_HETERO=

usage() {
  cat <<'EOF'
Usage: prepare_receptor.sh --input INPUT.pdb --output-dir DIR --record-id ID \
  --chains A[,B] --catalytic-residues A:10,A:20 \
  [--model 1] [--retain-hetero A:ZN:401]
EOF
}

while [ "$#" -gt 0 ]; do
  case "$1" in
    --input) INPUT=$2; shift 2 ;;
    --output-dir) OUTPUT_DIR=$2; shift 2 ;;
    --record-id) RECORD_ID=$2; shift 2 ;;
    --chains) CHAINS=$2; shift 2 ;;
    --model) MODEL=$2; shift 2 ;;
    --catalytic-residues) CATALYTIC_RESIDUES=$2; shift 2 ;;
    --retain-hetero) RETAIN_HETERO=$2; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) printf 'unknown argument: %s\n' "$1" >&2; usage >&2; exit 2 ;;
  esac
done

for value_name in INPUT OUTPUT_DIR RECORD_ID CHAINS CATALYTIC_RESIDUES; do
  eval "value=\${$value_name}"
  if [ -z "$value" ]; then
    printf 'required argument is empty: %s\n' "$value_name" >&2
    exit 2
  fi
done
if [ ! -f "$INPUT" ]; then
  printf 'input PDB is missing: %s\n' "$INPUT" >&2
  exit 2
fi
for executable in "$ANALYSIS_PYTHON" "$ADFR_BIN/prepare_receptor"; do
  if [ ! -x "$executable" ]; then
    printf 'required executable is missing: %s\n' "$executable" >&2
    exit 2
  fi
done
if [ ! -f "$CHECKER" ]; then
  printf 'receptor checker is missing: %s\n' "$CHECKER" >&2
  exit 2
fi

ROOT_REAL=$(readlink -f "$ROOT")
mkdir -p "$(dirname "$OUTPUT_DIR")"
OUTPUT_REAL=$(readlink -m "$OUTPUT_DIR")
case "$OUTPUT_REAL/" in
  "$ROOT_REAL/"*) ;;
  *)
    printf 'refusing output outside project root: %s\n' "$OUTPUT_REAL" >&2
    exit 2
    ;;
esac

SELECTED="$OUTPUT_DIR/selected.pdb"
SELECTION_METADATA="$OUTPUT_DIR/selection_metadata.json"
PREPARED="$OUTPUT_DIR/receptor.pdbqt"
METADATA="$OUTPUT_DIR/metadata.json"
CHECKSUMS="$OUTPUT_DIR/SHA256SUMS"
FAILURE="$OUTPUT_DIR/RECEPTOR_FAIL.json"
STDOUT_LOG="$OUTPUT_DIR/prepare_receptor.stdout.log"
STDERR_LOG="$OUTPUT_DIR/prepare_receptor.stderr.log"

if [ -f "$METADATA" ] && "$ANALYSIS_PYTHON" - "$METADATA" <<'PY'
import json
import sys
with open(sys.argv[1], encoding="utf-8") as handle:
    payload = json.load(handle)
raise SystemExit(0 if payload.get("status") == "RECEPTOR_PASS" else 1)
PY
then
  printf 'existing receptor PASS retained without overwrite: %s\n' "$METADATA"
  exit 0
fi
if [ -e "$OUTPUT_DIR" ]; then
  if find "$OUTPUT_DIR" -mindepth 1 -maxdepth 1 -print -quit | grep -q .; then
    printf 'refusing partial or unverified receptor output: %s\n' "$OUTPUT_DIR" >&2
    exit 3
  fi
else
  mkdir -p "$OUTPUT_DIR"
fi

on_failure() {
  status=$?
  if [ "$status" -ne 0 ] && [ ! -e "$FAILURE" ]; then
    printf '{"status":"NOT_EVALUATED_RECEPTOR_PREP","exit_code":%s,"record_id":"%s"}\n' \
      "$status" "$RECORD_ID" > "$FAILURE"
  fi
  exit "$status"
}
trap on_failure EXIT

"$ANALYSIS_PYTHON" "$CHECKER" select \
  --input "$INPUT" \
  --output "$SELECTED" \
  --metadata "$SELECTION_METADATA" \
  --chains "$CHAINS" \
  --model "$MODEL" \
  --catalytic-residues "$CATALYTIC_RESIDUES" \
  --retain-hetero "$RETAIN_HETERO"

"$ADFR_BIN/prepare_receptor" \
  -r "$SELECTED" \
  -o "$PREPARED" \
  -A checkhydrogens \
  -U nphs_lps_waters \
  >"$STDOUT_LOG" 2>"$STDERR_LOG"

"$ANALYSIS_PYTHON" "$CHECKER" finalize \
  --selection-metadata "$SELECTION_METADATA" \
  --prepared-pdbqt "$PREPARED" \
  --output-metadata "$METADATA"

sha256sum \
  "$SELECTED" \
  "$SELECTION_METADATA" \
  "$PREPARED" \
  "$METADATA" \
  "$STDOUT_LOG" \
  "$STDERR_LOG" \
  > "$CHECKSUMS"

trap - EXIT
printf 'RECEPTOR_PASS metadata written: %s\n' "$METADATA"
