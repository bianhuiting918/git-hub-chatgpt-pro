#!/usr/bin/env bash
set -euo pipefail

ROOT="${1:-/data/bht/project01_baker_serhyd_routeB_20260701}"
MAX_GPU_UTIL="${MAX_GPU_UTIL:-40}"
FORCE="${FORCE:-0}"
ALLOW_SHARED_GPU="${ALLOW_SHARED_GPU:-0}"

RFAA="/data/bht/design_tools/src/rf_diffusion_all_atom"
PY="/data/bht/design_tools/envs/rfaa_venv/bin/python"
CKPT="/data/bht/design_tools/src/rf_diffusion_all_atom/RFDiffusionAA_paper_weights.pt"
if [ ! -s "$CKPT" ]; then
  CKPT="/data/bht/design_tools/models/rfdiffusion_aa/RFDiffusionAA_paper_weights.pt"
fi

INPUT="$ROOT/inputs/baker_serhyd/design_pipeline/01_diffusion/inputs/theozyme.pdb"
OUTDIR="$ROOT/outputs/rfaa_baker_theozyme_smoke_n1"
LOGDIR="$ROOT/logs"
MANIFEST="$ROOT/manifests/baker_theozyme_smoke_manifest.tsv"
STATUS="$ROOT/manifests/baker_theozyme_smoke_status.json"
PIDFILE="$LOGDIR/rfaa_baker_theozyme_smoke.pid"
LOGFILE="$LOGDIR/rfaa_baker_theozyme_smoke.log"

mkdir -p "$OUTDIR" "$LOGDIR" "$ROOT/manifests"

if [ ! -s "$INPUT" ]; then
  echo "Missing input: $INPUT" >&2
  exit 2
fi
if [ ! -x "$PY" ]; then
  echo "Missing RFAA python: $PY" >&2
  exit 2
fi
if [ ! -s "$CKPT" ]; then
  echo "Missing RFAA checkpoint: $CKPT" >&2
  exit 2
fi

gpu_util="NA"
gpu_mem_used="NA"
gpu_mem_total="NA"
compute_process_count="NA"
compute_processes="NA"
if command -v nvidia-smi >/dev/null 2>&1; then
  gpu_line="$(nvidia-smi --query-gpu=utilization.gpu,memory.used,memory.total --format=csv,noheader,nounits | head -1 | tr -d ' ')"
  IFS=',' read -r gpu_util gpu_mem_used gpu_mem_total <<< "$gpu_line"
  compute_processes="$(nvidia-smi --query-compute-apps=pid,process_name,used_memory --format=csv,noheader,nounits 2>/dev/null | sed '/^[[:space:]]*$/d' || true)"
  if [ -z "$compute_processes" ]; then
    compute_process_count="0"
    compute_processes=""
  else
    compute_process_count="$(printf '%s\n' "$compute_processes" | wc -l | tr -d ' ')"
  fi
fi

if [ "$FORCE" != "1" ] && [ "$ALLOW_SHARED_GPU" != "1" ] && [ "$compute_process_count" != "NA" ] && [ "$compute_process_count" -gt 0 ]; then
  compute_processes_json="$(printf '%s' "$compute_processes" | python3 -c 'import json,sys; print(json.dumps(sys.stdin.read()))')"
  cat > "$STATUS" <<JSON
{
  "status": "BLOCKED_GPU_COMPUTE_PROCESS_PRESENT",
  "route_label": "baker_theozyme_new_backbone",
  "root": "$ROOT",
  "input": "$INPUT",
  "gpu_util_percent": "$gpu_util",
  "gpu_memory_used_mib": "$gpu_mem_used",
  "gpu_memory_total_mib": "$gpu_mem_total",
  "compute_process_count": "$compute_process_count",
  "compute_processes": $compute_processes_json,
  "allow_shared_gpu": "$ALLOW_SHARED_GPU",
  "next_action": "Wait until no other compute process is present, or intentionally set ALLOW_SHARED_GPU=1 / FORCE=1."
}
JSON
  echo "$STATUS"
  exit 0
fi

if [ "$FORCE" != "1" ] && [ "$gpu_util" != "NA" ] && [ "$gpu_util" -gt "$MAX_GPU_UTIL" ]; then
  cat > "$STATUS" <<JSON
{
  "status": "BLOCKED_GPU_BUSY_OR_UNAVAILABLE",
  "route_label": "baker_theozyme_new_backbone",
  "root": "$ROOT",
  "input": "$INPUT",
  "gpu_util_percent": "$gpu_util",
  "gpu_memory_used_mib": "$gpu_mem_used",
  "gpu_memory_total_mib": "$gpu_mem_total",
  "compute_process_count": "$compute_process_count",
  "max_gpu_util_percent": "$MAX_GPU_UTIL",
  "next_action": "Rerun this script when GPU util is below threshold, or set FORCE=1 intentionally."
}
JSON
  echo "$STATUS"
  exit 0
fi

printf 'run_id\troute_label\tinput_pdb\tligand\tcontig\tnum_designs\tdiffuser_T\tstatus\tpid\tlog\n' > "$MANIFEST"
printf 'baker_theozyme_smoke_n1\tbaker_theozyme_new_backbone\t%s\tbn1\t%s\t1\t20\tLAUNCHED\tPENDING\t%s\n' \
  "$INPUT" "12,A56-60,36,A83-85,15,A113-115,73,B145-147,10" "$LOGFILE" >> "$MANIFEST"

cd "$RFAA"
export CUDA_VISIBLE_DEVICES=0
export HYDRA_FULL_ERROR=1
nohup "$PY" run_inference.py \
  inference.deterministic=True \
  inference.input_pdb="$INPUT" \
  inference.ckpt_path="$CKPT" \
  inference.output_prefix="$OUTDIR/sample" \
  inference.ligand=bn1 \
  inference.num_designs=1 \
  inference.design_startnum=0 \
  diffuser.T=20 \
  contigmap.contigs=[\'12,A56-60,36,A83-85,15,A113-115,73,B145-147,10\'] \
  > "$LOGFILE" 2>&1 &

pid="$!"
echo "$pid" > "$PIDFILE"
python3 - "$MANIFEST" "$pid" <<'PY'
from pathlib import Path
import sys
p = Path(sys.argv[1])
pid = sys.argv[2]
text = p.read_text()
p.write_text(text.replace("\tPENDING\t", f"\t{pid}\t"))
PY

cat > "$STATUS" <<JSON
{
  "status": "LAUNCHED",
  "route_label": "baker_theozyme_new_backbone",
  "root": "$ROOT",
  "input": "$INPUT",
  "pid": "$pid",
  "pidfile": "$PIDFILE",
  "log": "$LOGFILE",
  "manifest": "$MANIFEST",
  "gpu_util_percent_at_launch": "$gpu_util",
  "gpu_memory_used_mib_at_launch": "$gpu_mem_used",
  "gpu_memory_total_mib": "$gpu_mem_total",
  "compute_process_count_at_launch": "$compute_process_count",
  "allow_shared_gpu": "$ALLOW_SHARED_GPU"
}
JSON

echo "$STATUS"
