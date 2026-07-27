#!/usr/bin/env bash
#SBATCH --job-name=surface_field
#SBATCH --partition=xahcnormal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=02:00:00
#SBATCH --output=/work/home/acshdt1dks/polymer_surface_hotspot_screen_20260725/logs/surface_fields/%A_%a.out
#SBATCH --error=/work/home/acshdt1dks/polymer_surface_hotspot_screen_20260725/logs/surface_fields/%A_%a.err

set -euo pipefail
ROOT=/work/home/acshdt1dks/polymer_surface_hotspot_screen_20260725
DEPLOY_COMMIT=${DEPLOY_COMMIT:?DEPLOY_COMMIT is required}
MANIFEST=${MANIFEST:?MANIFEST is required}
FAMILY=${FAMILY:?FAMILY is required}
OUTPUT_ROOT=${OUTPUT_ROOT:?OUTPUT_ROOT is required}
ROW_INDEX=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}
TOOL_ROOT="$ROOT/deployments/$DEPLOY_COMMIT/projects/03-polymer-surface-hotspot-screening"
PYTHON="$ROOT/envs/surface-screen-py311/bin/python3.11"
ADFR_ROOT="$ROOT/software/ADFRsuite-1.0"

case "$FAMILY" in PET|NYLON) ;; *) echo "invalid FAMILY=$FAMILY" >&2; exit 2 ;; esac
case "$OUTPUT_ROOT" in "$ROOT"/results/*) ;; *) echo "output outside results: $OUTPUT_ROOT" >&2; exit 2 ;; esac
test -f "$MANIFEST"
test -f "$TOOL_ROOT/scripts/run_surface_field_row.py"
mkdir -p "$ROOT/logs/surface_fields"

"$PYTHON" "$TOOL_ROOT/scripts/run_surface_field_row.py" \
  --manifest "$MANIFEST" \
  --row-index "$ROW_INDEX" \
  --family "$FAMILY" \
  --output-root "$OUTPUT_ROOT" \
  --project-root "$ROOT" \
  --tool-root "$TOOL_ROOT" \
  --adfr-root "$ADFR_ROOT"
