#!/usr/bin/env bash
set -euo pipefail

ROOT="${1:-/data/bht/project01_baker_serhyd_routeB_20260701}"
REPO="${REPO:-/data/bht/git-hub-chatgpt-pro-codex}"
MAX_GPU_UTIL="${MAX_GPU_UTIL:-40}"
INTERVAL_SECONDS="${INTERVAL_SECONDS:-300}"
MAX_WAIT_MINUTES="${MAX_WAIT_MINUTES:-360}"

LOGDIR="$ROOT/logs"
MANIFEST_DIR="$ROOT/manifests"
MONITOR_STATUS="$MANIFEST_DIR/baker_theozyme_smoke_monitor_status.json"
MONITOR_LOG="$LOGDIR/baker_theozyme_smoke_monitor.log"
LOCKFILE="$LOGDIR/baker_theozyme_smoke_monitor.lock"

mkdir -p "$LOGDIR" "$MANIFEST_DIR"

exec 9>"$LOCKFILE"
if ! flock -n 9; then
  echo "Another monitor is already running: $LOCKFILE" >&2
  exit 0
fi

deadline=$(( $(date +%s) + MAX_WAIT_MINUTES * 60 ))
attempt=0

write_status() {
  local status="$1"
  local gpu_util="$2"
  local gpu_mem_used="$3"
  local gpu_mem_total="$4"
  local detail="$5"
  cat > "$MONITOR_STATUS" <<JSON
{
  "status": "$status",
  "route_label": "baker_theozyme_new_backbone",
  "root": "$ROOT",
  "attempt": $attempt,
  "max_gpu_util_percent": "$MAX_GPU_UTIL",
  "interval_seconds": "$INTERVAL_SECONDS",
  "max_wait_minutes": "$MAX_WAIT_MINUTES",
  "gpu_util_percent": "$gpu_util",
  "gpu_memory_used_mib": "$gpu_mem_used",
  "gpu_memory_total_mib": "$gpu_mem_total",
  "detail": "$detail",
  "monitor_log": "$MONITOR_LOG",
  "updated_at": "$(date -Is)"
}
JSON
}

echo "[$(date -Is)] monitor_start root=$ROOT max_gpu_util=$MAX_GPU_UTIL interval=$INTERVAL_SECONDS max_wait_minutes=$MAX_WAIT_MINUTES" >> "$MONITOR_LOG"

while [ "$(date +%s)" -le "$deadline" ]; do
  attempt=$((attempt + 1))
  gpu_util="NA"
  gpu_mem_used="NA"
  gpu_mem_total="NA"
  if command -v nvidia-smi >/dev/null 2>&1; then
    gpu_line="$(nvidia-smi --query-gpu=utilization.gpu,memory.used,memory.total --format=csv,noheader,nounits | head -1 | tr -d ' ')"
    IFS=',' read -r gpu_util gpu_mem_used gpu_mem_total <<< "$gpu_line"
  fi

  echo "[$(date -Is)] attempt=$attempt gpu_util=$gpu_util mem=${gpu_mem_used}/${gpu_mem_total}" >> "$MONITOR_LOG"

  if [ "$gpu_util" != "NA" ] && [ "$gpu_util" -le "$MAX_GPU_UTIL" ]; then
    write_status "GPU_AVAILABLE_LAUNCHING" "$gpu_util" "$gpu_mem_used" "$gpu_mem_total" "GPU below threshold; invoking launcher"
    cd "$REPO"
    git fetch --depth 1 origin phase1-reset-entrance-gate >> "$MONITOR_LOG" 2>&1 || true
    git reset --hard origin/phase1-reset-entrance-gate >> "$MONITOR_LOG" 2>&1 || true
    bash project01/phase1_20260701/scripts/launch_baker_theozyme_smoke.sh "$ROOT" >> "$MONITOR_LOG" 2>&1
    write_status "LAUNCHER_INVOKED" "$gpu_util" "$gpu_mem_used" "$gpu_mem_total" "See baker_theozyme_smoke_status.json for launch result"
    exit 0
  fi

  write_status "WAITING_GPU_BUSY" "$gpu_util" "$gpu_mem_used" "$gpu_mem_total" "GPU above threshold; not launching"
  sleep "$INTERVAL_SECONDS"
done

write_status "TIMEOUT_GPU_STILL_BUSY" "$gpu_util" "$gpu_mem_used" "$gpu_mem_total" "Max wait elapsed without launch"
echo "[$(date -Is)] monitor_timeout" >> "$MONITOR_LOG"
