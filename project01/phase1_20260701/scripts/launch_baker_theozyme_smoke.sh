#!/usr/bin/env bash
set -euo pipefail

ROOT="${1:-/data/bht/project01_baker_serhyd_routeB_20260701}"
MAX_GPU_UTIL="${MAX_GPU_UTIL:-40}"
FORCE="${FORCE:-0}"

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
if command -v nvidia-smi >/dev/null 2>&1; then
  gpu_line="$(nvidia-smi --query-gpu=utilization.gpu,memory.used,memory.total --format=csv,noheader,nounits | head -1 | tr -d ' ')"
  IFS=',' read -r gpu_util gpu_mem_used gpu_mem_total <<< "$gpu_line"
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
  "gpu_memory_total_mib": "$gpu_mem_total"
}
JSON

echo "$STATUS"
