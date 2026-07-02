#!/usr/bin/env bash
set -euo pipefail

ROOT="${1:-/data/bht/project01_baker_serhyd_routeB_20260701}"
CONTIG_SET="${CONTIG_SET:-pocket4_compact}"
RUN_ID="${RUN_ID:-ca_rfd_baker_pocket4_layered_${CONTIG_SET}_publicckpt_20260702}"
NUM_DESIGNS="${NUM_DESIGNS:-2}"
DESIGN_STARTNUM="${DESIGN_STARTNUM:-7000}"
MAX_GPU_UTIL="${MAX_GPU_UTIL:-40}"
ALLOW_SHARED_GPU="${ALLOW_SHARED_GPU:-0}"
FORCE="${FORCE:-0}"

CA_RFD="/data/bht/design_tools/src/CA_RFDiffusion"
PY="/data/bht/design_tools/envs/ca_rfd_release/bin/python"
PY_JSON="${PY_JSON:-/data/bht/design_tools/envs/rfaa_venv/bin/python}"
CKPT="/data/bht/design_tools/models/ca_rf_diffusion/ca_rfd_diffusion.pt"
INPUT="$ROOT/inputs/baker_serhyd/design_pipeline/01_diffusion/inputs/theozyme.pdb"
PREP="$ROOT/scripts/prepare_pocket4_layered_route.py"
POCKET_MANIFEST="$ROOT/manifests/pocket4_first_layered_route_manifest.json"
CONTIG_TSV="$ROOT/manifests/pocket4_first_layered_contigs.tsv"
OUTDIR="$ROOT/outputs/$RUN_ID"
LOGDIR="$ROOT/logs"
MANIFEST="$ROOT/manifests/${RUN_ID}_manifest.tsv"
STATUS="$ROOT/manifests/${RUN_ID}_status.json"
PIDFILE="$LOGDIR/${RUN_ID}.pid"
LOGFILE="$LOGDIR/${RUN_ID}.log"

mkdir -p "$OUTDIR" "$LOGDIR" "$ROOT/manifests"

if [ ! -s "$INPUT" ]; then echo "Missing input: $INPUT" >&2; exit 2; fi
if [ ! -x "$PY" ]; then echo "Missing python: $PY" >&2; exit 2; fi
if [ ! -s "$CKPT" ]; then echo "Missing checkpoint: $CKPT" >&2; exit 2; fi
if [ ! -f "$CA_RFD/rf_diffusion/run_inference.py" ]; then echo "Missing run_inference.py under $CA_RFD" >&2; exit 2; fi
if [ ! -f "$PREP" ]; then echo "Missing pocket4 prep script: $PREP" >&2; exit 2; fi

"$PY_JSON" "$PREP" \
  --input-pdb "$INPUT" \
  --ligand bn1 \
  --pocket-cutoff 4.0 \
  --out-json "$POCKET_MANIFEST" \
  --out-tsv "$ROOT/manifests/pocket4_first_layered_fixed_residues.tsv" \
  --contig-tsv "$CONTIG_TSV" >/dev/null

CONTIG_INFO="$("$PY_JSON" - "$CONTIG_TSV" "$CONTIG_SET" <<'PY'
from __future__ import annotations
import csv, sys
path, target = sys.argv[1], sys.argv[2]
for row in csv.DictReader(open(path), delimiter="\t"):
    if row["contig_set"] == target:
        print(row["contig"] + "\t" + row["fixed_segment_count"])
        raise SystemExit(0)
raise SystemExit(f"contig_set not found: {target}")
PY
)"
IFS=$'\t' read -r CONTIG FIXED_SEGMENT_COUNT <<< "$CONTIG_INFO"
IJ_VISIBLE_COUNT=$((FIXED_SEGMENT_COUNT + 1))
IJ_VISIBLE="$("$PY_JSON" - "$IJ_VISIBLE_COUNT" <<'PY'
import string, sys
n = int(sys.argv[1])
letters = string.ascii_lowercase + string.ascii_uppercase
if n > len(letters):
    raise SystemExit(f"Too many chunks for ij_visible alphabet: {n}")
print(letters[:n])
PY
)"

gpu_util="NA"; gpu_mem_used="NA"; gpu_mem_total="NA"; compute_process_count="NA"; compute_processes=""
if command -v nvidia-smi >/dev/null 2>&1; then
  gpu_line="$(nvidia-smi --query-gpu=utilization.gpu,memory.used,memory.total --format=csv,noheader,nounits | head -1 | tr -d ' ')"
  IFS=',' read -r gpu_util gpu_mem_used gpu_mem_total <<< "$gpu_line"
  compute_processes="$(nvidia-smi --query-compute-apps=pid,process_name,used_memory --format=csv,noheader,nounits 2>/dev/null | sed '/^[[:space:]]*$/d' || true)"
  if [ -z "$compute_processes" ]; then compute_process_count="0"; else compute_process_count="$(printf '%s\n' "$compute_processes" | wc -l | tr -d ' ')"; fi
fi

if [ "$FORCE" != "1" ] && [ "$ALLOW_SHARED_GPU" != "1" ] && [ "$compute_process_count" != "NA" ] && [ "$compute_process_count" -gt 0 ]; then
  compute_processes_json="$(printf '%s' "$compute_processes" | "$PY_JSON" -c 'import json,sys; print(json.dumps(sys.stdin.read()))')"
  cat > "$STATUS" <<JSON
{
  "status": "BLOCKED_GPU_COMPUTE_PROCESS_PRESENT",
  "run_id": "$RUN_ID",
  "route_label": "pocket4_first_layered_public_ca_rfdiffusion",
  "contig_set": "$CONTIG_SET",
  "contig": "$CONTIG",
  "fixed_segment_count": "$FIXED_SEGMENT_COUNT",
  "ij_visible": "$IJ_VISIBLE",
  "num_designs": "$NUM_DESIGNS",
  "design_startnum": "$DESIGN_STARTNUM",
  "gpu_util_percent": "$gpu_util",
  "gpu_memory_used_mib": "$gpu_mem_used",
  "gpu_memory_total_mib": "$gpu_mem_total",
  "compute_process_count": "$compute_process_count",
  "compute_processes": $compute_processes_json,
  "next_action": "Rerun when GPU compute process clears, or intentionally set ALLOW_SHARED_GPU=1."
}
JSON
  echo "$STATUS"
  exit 0
fi

if [ "$FORCE" != "1" ] && [ "$gpu_util" != "NA" ] && [ "$gpu_util" -gt "$MAX_GPU_UTIL" ]; then
  cat > "$STATUS" <<JSON
{
  "status": "BLOCKED_GPU_BUSY_OR_UNAVAILABLE",
  "run_id": "$RUN_ID",
  "route_label": "pocket4_first_layered_public_ca_rfdiffusion",
  "contig_set": "$CONTIG_SET",
  "contig": "$CONTIG",
  "fixed_segment_count": "$FIXED_SEGMENT_COUNT",
  "ij_visible": "$IJ_VISIBLE",
  "gpu_util_percent": "$gpu_util",
  "gpu_memory_used_mib": "$gpu_mem_used",
  "gpu_memory_total_mib": "$gpu_mem_total",
  "max_gpu_util_percent": "$MAX_GPU_UTIL",
  "next_action": "Rerun when GPU utilization is below threshold, or set FORCE=1 intentionally."
}
JSON
  echo "$STATUS"
  exit 0
fi

cat > "$MANIFEST" <<TSV
run_id	route_label	input_pdb	pocket_manifest	ckpt	ligand	contig_set	contig	fixed_segment_count	ij_visible	num_designs	design_startnum	diffuser_T	status	pid	log	output_prefix
$RUN_ID	pocket4_first_layered_public_ca_rfdiffusion	$INPUT	$POCKET_MANIFEST	$CKPT	bn1	$CONTIG_SET	$CONTIG	$FIXED_SEGMENT_COUNT	$IJ_VISIBLE	$NUM_DESIGNS	$DESIGN_STARTNUM	50	LAUNCHED	PENDING	$LOGFILE	$OUTDIR/sample
TSV

cd "$CA_RFD"
export CUDA_VISIBLE_DEVICES=0
export HYDRA_FULL_ERROR=1
export PYTHONPATH="$CA_RFD${PYTHONPATH:+:$PYTHONPATH}"

nohup "$PY" rf_diffusion/run_inference.py \
  --config-name=RFdiffusion_CA_inference \
  inference.deterministic=True \
  inference.input_pdb="$INPUT" \
  inference.ckpt_path="$CKPT" \
  inference.output_prefix="$OUTDIR/sample" \
  inference.ligand=bn1 \
  inference.ij_visible="$IJ_VISIBLE" \
  inference.num_designs="$NUM_DESIGNS" \
  inference.design_startnum="$DESIGN_STARTNUM" \
  inference.cautious=true \
  inference.write_trajectory=false \
  diffuser.T=50 \
  denoiser.noise_scale_frame=0.05 \
  denoiser.noise_scale_ca=0.0 \
  motif_only_2d=true \
  preprocess.eye_frames=true \
  "contigmap.contigs=['$CONTIG']" \
  > "$LOGFILE" 2>&1 &

pid="$!"
echo "$pid" > "$PIDFILE"
"$PY_JSON" - "$MANIFEST" "$pid" <<'PY'
from pathlib import Path
import sys
p = Path(sys.argv[1])
pid = sys.argv[2]
p.write_text(p.read_text().replace("\tPENDING\t", f"\t{pid}\t"))
PY

cat > "$STATUS" <<JSON
{
  "status": "LAUNCHED",
  "run_id": "$RUN_ID",
  "route_label": "pocket4_first_layered_public_ca_rfdiffusion",
  "contig_set": "$CONTIG_SET",
  "contig": "$CONTIG",
  "fixed_segment_count": "$FIXED_SEGMENT_COUNT",
  "ij_visible": "$IJ_VISIBLE",
  "input": "$INPUT",
  "pocket_manifest": "$POCKET_MANIFEST",
  "ckpt": "$CKPT",
  "num_designs": "$NUM_DESIGNS",
  "design_startnum": "$DESIGN_STARTNUM",
  "pid": "$pid",
  "pidfile": "$PIDFILE",
  "log": "$LOGFILE",
  "manifest": "$MANIFEST",
  "output_prefix": "$OUTDIR/sample",
  "checkpoint_caveat": "public ca_rfd_diffusion.pt, not Baker BFF_7.pt"
}
JSON

echo "$STATUS"