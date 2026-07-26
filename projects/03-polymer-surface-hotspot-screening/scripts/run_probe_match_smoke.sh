#!/usr/bin/env bash
set -euo pipefail

ROOT=/work/home/acshdt1dks/polymer_surface_hotspot_screen_20260725
COMMIT=""
MATERIAL=""
RECORD_ID=""
PROBE_IDS=""
SHELL_DIR=""
PATCH_DIR=""
MAPS_ROOT=""
RECEPTOR_PDBQT=""
OUTPUT_DIR=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --commit) COMMIT="$2"; shift 2 ;;
    --material-family) MATERIAL="$2"; shift 2 ;;
    --record-id) RECORD_ID="$2"; shift 2 ;;
    --probe-ids) PROBE_IDS="$2"; shift 2 ;;
    --shell-dir) SHELL_DIR="$2"; shift 2 ;;
    --patch-dir) PATCH_DIR="$2"; shift 2 ;;
    --maps-root) MAPS_ROOT="$2"; shift 2 ;;
    --receptor-pdbqt) RECEPTOR_PDBQT="$2"; shift 2 ;;
    --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done

[[ "$COMMIT" =~ ^[0-9a-f]{40}$ ]] || { echo "invalid commit" >&2; exit 2; }
[[ "$MATERIAL" == PET || "$MATERIAL" == NYLON ]] || { echo "invalid material family" >&2; exit 2; }
for value in "$RECORD_ID" "$PROBE_IDS" "$SHELL_DIR" "$PATCH_DIR" "$MAPS_ROOT" "$RECEPTOR_PDBQT" "$OUTPUT_DIR"; do
  [[ -n "$value" ]] || { echo "missing required argument" >&2; exit 2; }
done
case "$OUTPUT_DIR" in
  "$ROOT"/results/probe_match_smoke_*) ;;
  *) echo "output directory outside allowed smoke root" >&2; exit 2 ;;
esac
[[ ! -e "$OUTPUT_DIR" ]] || { echo "refusing existing output: $OUTPUT_DIR" >&2; exit 3; }

PYTHON="$ROOT/envs/surface-screen-py311/bin/python3.11"
MATCHER="$ROOT/deployments/$COMMIT/projects/03-polymer-surface-hotspot-screening/scripts/match_probe_fields.py"
[[ -x "$PYTHON" ]] || { echo "missing Python environment" >&2; exit 2; }
[[ -f "$MATCHER" ]] || { echo "missing commit-pinned matcher" >&2; exit 2; }

exec "$PYTHON" "$MATCHER"   --shell-dir "$SHELL_DIR"   --patch-dir "$PATCH_DIR"   --maps-root "$MAPS_ROOT"   --probe-sdf "$ROOT/results/probes_v1/probes.sdf"   --probe-manifest "$ROOT/results/probes_v1/validated_probes.tsv"   --receptor-pdbqt "$RECEPTOR_PDBQT"   --output-dir "$OUTPUT_DIR"   --record-id "$RECORD_ID"   --material-family "$MATERIAL"   --probe-ids "$PROBE_IDS"   --region catalytic_neighborhood
