#!/usr/bin/env bash
#SBATCH --job-name=surface4
#SBATCH --partition=xahcnormal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2G
#SBATCH --time=24:00:00
#SBATCH --output=/work/home/acshdt1dks/polymer_surface_hotspot_screen_20260725/logs/surface_fields/%A_%a.out
#SBATCH --error=/work/home/acshdt1dks/polymer_surface_hotspot_screen_20260725/logs/surface_fields/%A_%a.err

set -euo pipefail
ROOT=/work/home/acshdt1dks/polymer_surface_hotspot_screen_20260725
DEPLOY_COMMIT=${DEPLOY_COMMIT:?DEPLOY_COMMIT is required}
MANIFEST=${MANIFEST:?MANIFEST is required}
FAMILY=${FAMILY:?FAMILY is required}
OUTPUT_ROOT=${OUTPUT_ROOT:?OUTPUT_ROOT is required}
EXISTING_PASS_ROOT=${EXISTING_PASS_ROOT:?EXISTING_PASS_ROOT is required}
BASE_START=${BASE_START:?BASE_START is required}
TOTAL_ROWS=${TOTAL_ROWS:?TOTAL_ROWS is required}
CHUNK_SIZE=${CHUNK_SIZE:?CHUNK_SIZE is required}
TASK_ID=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}
TOOL_ROOT="$ROOT/deployments/$DEPLOY_COMMIT/projects/03-polymer-surface-hotspot-screening"
PYTHON="$ROOT/envs/surface-screen-py311/bin/python3.11"
PARALLEL_RUNNER="$TOOL_ROOT/scripts/run_surface_field_parallel_chunk.py"
ROW_WORKER="$TOOL_ROOT/scripts/run_surface_field_row.py"
ADFR_ROOT="$ROOT/software/ADFRsuite-1.0"

case "$FAMILY" in PET|NYLON) ;; *) echo "invalid FAMILY=$FAMILY" >&2; exit 2 ;; esac
case "$OUTPUT_ROOT" in "$ROOT"/results/*) ;; *) echo "output outside results: $OUTPUT_ROOT" >&2; exit 2 ;; esac
case "$EXISTING_PASS_ROOT" in "$ROOT"/results/*) ;; *) echo "existing pass root outside results: $EXISTING_PASS_ROOT" >&2; exit 2 ;; esac
test -f "$MANIFEST"
test -f "$PARALLEL_RUNNER"
test -f "$ROW_WORKER"
mkdir -p "$ROOT/logs/surface_fields"

start=$((BASE_START + TASK_ID * CHUNK_SIZE))
end=$((start + CHUNK_SIZE))
if [ "$end" -gt "$TOTAL_ROWS" ]; then end=$TOTAL_ROWS; fi
if [ "$start" -ge "$TOTAL_ROWS" ]; then
  echo "EMPTY_CHUNK start=$start total=$TOTAL_ROWS"
  exit 0
fi

exec "$PYTHON" "$PARALLEL_RUNNER"   --manifest "$MANIFEST"   --family "$FAMILY"   --output-root "$OUTPUT_ROOT"   --existing-pass-root "$EXISTING_PASS_ROOT"   --project-root "$ROOT"   --tool-root "$TOOL_ROOT"   --adfr-root "$ADFR_ROOT"   --row-worker "$ROW_WORKER"   --start-row "$start"   --end-row "$end"   --workers "$SLURM_CPUS_PER_TASK"
