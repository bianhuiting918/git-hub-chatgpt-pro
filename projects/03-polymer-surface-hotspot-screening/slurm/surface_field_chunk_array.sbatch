#!/usr/bin/env bash
#SBATCH --job-name=surface_chunk
#SBATCH --partition=xahcnormal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
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
BASE_START=${BASE_START:?BASE_START is required}
TOTAL_ROWS=${TOTAL_ROWS:?TOTAL_ROWS is required}
CHUNK_SIZE=${CHUNK_SIZE:?CHUNK_SIZE is required}
TASK_ID=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}
TOOL_ROOT="$ROOT/deployments/$DEPLOY_COMMIT/projects/03-polymer-surface-hotspot-screening"
PYTHON="$ROOT/envs/surface-screen-py311/bin/python3.11"
WORKER="$TOOL_ROOT/scripts/run_surface_field_row.py"
ADFR_ROOT="$ROOT/software/ADFRsuite-1.0"

case "$FAMILY" in PET|NYLON) ;; *) echo "invalid FAMILY=$FAMILY" >&2; exit 2 ;; esac
case "$OUTPUT_ROOT" in "$ROOT"/results/*) ;; *) echo "output outside results: $OUTPUT_ROOT" >&2; exit 2 ;; esac
test -f "$MANIFEST"; test -f "$WORKER"; mkdir -p "$ROOT/logs/surface_fields"

start=$((BASE_START + TASK_ID * CHUNK_SIZE))
end=$((start + CHUNK_SIZE))
if [ "$end" -gt "$TOTAL_ROWS" ]; then end=$TOTAL_ROWS; fi
if [ "$start" -ge "$TOTAL_ROWS" ]; then
  echo "EMPTY_CHUNK start=$start total=$TOTAL_ROWS"
  exit 0
fi

failed=0
processed=0
for ((row=start; row<end; row++)); do
  set +e
  "$PYTHON" "$WORKER" \
    --manifest "$MANIFEST" --row-index "$row" --family "$FAMILY" \
    --output-root "$OUTPUT_ROOT" --project-root "$ROOT" \
    --tool-root "$TOOL_ROOT" --adfr-root "$ADFR_ROOT"
  rc=$?
  set -e
  processed=$((processed + 1))
  if [ "$rc" -ne 0 ]; then
    failed=$((failed + 1))
    echo "ROW_NOT_EVALUATED row=$row rc=$rc" >&2
  fi
done
echo "CHUNK_COMPLETE family=$FAMILY start=$start end=$end processed=$processed failed=$failed"
if [ "$failed" -ne 0 ]; then exit 10; fi
