#!/usr/bin/env bash
set -euo pipefail

ROOT="${1:-/data/bht/project01_baker_serhyd_routeB_20260701}"
PARENT_SET="${PARENT_SET:-expanded}"
RUN_ID="${RUN_ID:-ca_rfd_baker_pocket4_l2_extend_${PARENT_SET}_n1_20260703}"
NUM_DESIGNS="${NUM_DESIGNS:-1}"
DESIGN_STARTNUM="${DESIGN_STARTNUM:-7500}"
N_FLANK="${N_FLANK:-24-36}"
C_FLANK="${C_FLANK:-24-36}"
MAX_GPU_UTIL="${MAX_GPU_UTIL:-40}"
ALLOW_SHARED_GPU="${ALLOW_SHARED_GPU:-0}"
FORCE="${FORCE:-0}"

CA_RFD="/data/bht/design_tools/src/CA_RFDiffusion"
PY="/data/bht/design_tools/envs/ca_rfd_release/bin/python"
PY_JSON="${PY_JSON:-/data/bht/design_tools/envs/rfaa_venv/bin/python}"
CKPT="/data/bht/design_tools/models/ca_rf_diffusion/ca_rfd_diffusion.pt"

case "$PARENT_SET" in
  compact)
    PARENT_PDB="${PARENT_PDB:-$ROOT/outputs/ca_rfd_baker_pocket4_layered_compact_n1_v2_20260702/sample_7200.pdb}"
    ;;
  medium)
    PARENT_PDB="${PARENT_PDB:-$ROOT/outputs/ca_rfd_baker_pocket4_layered_medium_n1_20260702/sample_7300.pdb}"
    ;;
  expanded)
    PARENT_PDB="${PARENT_PDB:-$ROOT/outputs/ca_rfd_baker_pocket4_layered_expanded_n1_20260702/sample_7400.pdb}"
    ;;
  *)
    if [ -z "${PARENT_PDB:-}" ]; then
      echo "Unknown PARENT_SET=$PARENT_SET and PARENT_PDB is not set" >&2
      exit 2
    fi
    ;;
esac

OUTDIR="$ROOT/outputs/$RUN_ID"
LOGDIR="$ROOT/logs"
MANIFEST="$ROOT/manifests/${RUN_ID}_manifest.tsv"
STATUS="$ROOT/manifests/${RUN_ID}_status.json"
PIDFILE="$LOGDIR/${RUN_ID}.pid"
LOGFILE="$LOGDIR/${RUN_ID}.log"

mkdir -p "$OUTDIR" "$LOGDIR" "$ROOT/manifests"

if [ ! -s "$PARENT_PDB" ]; then echo "Missing parent PDB: $PARENT_PDB" >&2; exit 2; fi
if [ ! -x "$PY" ]; then echo "Missing python: $PY" >&2; exit 2; fi
if [ ! -s "$CKPT" ]; then echo "Missing checkpoint: $CKPT" >&2; exit 2; fi
if [ ! -f "$CA_RFD/rf_diffusion/run_inference.py" ]; then echo "Missing run_inference.py under $CA_RFD" >&2; exit 2; fi

PARENT_INFO="$("$PY_JSON" - "$PARENT_PDB" <<'PY'
from pathlib import Path
import sys
p = Path(sys.argv[1])
ca = []
hetatm = 0
for line in p.read_text().splitlines():
    if line.startswith("ATOM") and line[12:16].strip() == "CA":
        ca.append((line[21].strip() or "_", int(line[22:26])))
    elif line.startswith("HETATM"):
        hetatm += 1
if not ca:
    raise SystemExit("No CA atoms found in parent PDB")
chains = sorted({c for c, _ in ca})
if chains != ["A"]:
    raise SystemExit(f"Expected single protein chain A, got {chains}")
print(f"{len(ca)}\t{ca[0][1]}\t{ca[-1][1]}\t{hetatm}")
PY
)"
IFS=$'\t' read -r PARENT_LEN FIRST_RES LAST_RES HETATM_COUNT <<< "$PARENT_INFO"
CONTIG="${N_FLANK},A${FIRST_RES}-${LAST_RES},${C_FLANK}"
IJ_VISIBLE="ab"

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
  "route_label": "pocket4_l2_extend_from_passed_l1",
  "parent_set": "$PARENT_SET",
  "parent_pdb": "$PARENT_PDB",
  "parent_len": "$PARENT_LEN",
  "contig": "$CONTIG",
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
  "route_label": "pocket4_l2_extend_from_passed_l1",
  "parent_set": "$PARENT_SET",
  "parent_pdb": "$PARENT_PDB",
  "parent_len": "$PARENT_LEN",
  "contig": "$CONTIG",
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
run_id	route_label	parent_set	parent_pdb	parent_len	parent_hetatm	ckpt	ligand	contig	ij_visible	num_designs	design_startnum	diffuser_T	status	pid	log	output_prefix
$RUN_ID	pocket4_l2_extend_from_passed_l1	$PARENT_SET	$PARENT_PDB	$PARENT_LEN	$HETATM_COUNT	$CKPT	bn1	$CONTIG	$IJ_VISIBLE	$NUM_DESIGNS	$DESIGN_STARTNUM	50	LAUNCHED	PENDING	$LOGFILE	$OUTDIR/sample
TSV

cd "$CA_RFD"
export CUDA_VISIBLE_DEVICES=0
export HYDRA_FULL_ERROR=1
export PYTHONPATH="$CA_RFD${PYTHONPATH:+:$PYTHONPATH}"

nohup "$PY" rf_diffusion/run_inference.py \
  --config-name=RFdiffusion_CA_inference \
  inference.deterministic=True \
  inference.input_pdb="$PARENT_PDB" \
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
  "route_label": "pocket4_l2_extend_from_passed_l1",
  "parent_set": "$PARENT_SET",
  "parent_pdb": "$PARENT_PDB",
  "parent_len": "$PARENT_LEN",
  "parent_hetatm": "$HETATM_COUNT",
  "contig": "$CONTIG",
  "ij_visible": "$IJ_VISIBLE",
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
