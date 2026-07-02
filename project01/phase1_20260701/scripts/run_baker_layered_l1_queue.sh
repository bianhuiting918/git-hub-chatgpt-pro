#!/usr/bin/env bash
set -euo pipefail

ROOT="${1:-/data/bht/project01_baker_serhyd_routeB_20260701}"
CONTIG_SETS="${CONTIG_SETS:-compact medium near_original}"
NUM_DESIGNS="${NUM_DESIGNS:-20}"
START_BASE="${START_BASE:-3000}"
PY_JSON="${PY_JSON:-/data/bht/design_tools/envs/rfaa_venv/bin/python}"
LAUNCHER="$ROOT/scripts/launch_baker_layered_l1_contig_sweep.sh"
GATE="$ROOT/scripts/gate_ca_rfdiffusion_theozyme_motif.py"
REF="$ROOT/inputs/baker_serhyd/design_pipeline/01_diffusion/inputs/theozyme.pdb"
QUEUE_STATUS="$ROOT/manifests/ca_rfd_baker_layered_l1_queue_status.json"

mkdir -p "$ROOT/manifests"

contig_for_set() {
  case "$1" in
    compact) echo "6,A56-60,18,A83-85,8,A113-115,32,B145-147,6" ;;
    medium) echo "8,A56-60,24,A83-85,10,A113-115,48,B145-147,8" ;;
    near_original) echo "10,A56-60,30,A83-85,12,A113-115,60,B145-147,8" ;;
    original) echo "12,A56-60,36,A83-85,15,A113-115,73,B145-147,10" ;;
    *) echo "Unknown contig set: $1" >&2; return 2 ;;
  esac
}

write_queue_status() {
  "$PY_JSON" - "$QUEUE_STATUS" "$@" <<'PY'
from __future__ import annotations
import json, sys, time
path = sys.argv[1]
items = sys.argv[2:]
data = {"updated": time.strftime("%Y-%m-%dT%H:%M:%S"), "route": "baker_layered_l1_contig_queue"}
for item in items:
    key, value = item.split("=", 1)
    data[key] = value
with open(path, "w") as handle:
    json.dump(data, handle, indent=2, sort_keys=True)
    handle.write("\n")
print(path)
PY
}

json_value() {
  local file="$1" key="$2"
  if [ ! -s "$file" ]; then return 1; fi
  "$PY_JSON" - "$file" "$key" <<'PY'
import json, sys
try:
    data = json.load(open(sys.argv[1]))
except Exception:
    sys.exit(1)
value = data.get(sys.argv[2], "")
print(value)
PY
}

run_gate_for_set() {
  local set_name="$1" run_id="$2" contig="$3" outdir="$4"
  local gate_tsv="$ROOT/manifests/${run_id}_motif_gate.tsv"
  local gate_json="$ROOT/manifests/${run_id}_motif_gate_summary.json"
  "$PY_JSON" "$GATE" \
    --reference-pdb "$REF" \
    --output-dir "$outdir" \
    --out-tsv "$gate_tsv" \
    --summary-json "$gate_json" \
    --ligand bn1 \
    --contig "$contig" \
    --rmsd-cutoff 1.0 \
    --pair-cutoff 1.0 >/dev/null
  "$PY_JSON" - "$gate_json" <<'PY'
import json, sys
p = sys.argv[1]
d = json.load(open(p))
print(f"evaluated={d.get('evaluated_pdb_count', 0)} counts={d.get('counts', {})} summary={p}")
PY
}

index=0
for set_name in $CONTIG_SETS; do
  contig="$(contig_for_set "$set_name")"
  run_id="ca_rfd_baker_layered_l1_${set_name}_publicckpt_20260702"
  status_file="$ROOT/manifests/${run_id}_status.json"
  outdir="$ROOT/outputs/$run_id"
  startnum=$((START_BASE + index * 1000))
  index=$((index + 1))

  status=""
  pid=""
  if [ -s "$status_file" ]; then
    status="$(json_value "$status_file" status || true)"
    pid="$(json_value "$status_file" pid || true)"
  fi

  pid_alive=0
  if [ "$status" = "LAUNCHED" ] && [ -n "$pid" ] && ps -p "$pid" >/dev/null 2>&1; then
    pid_alive=1
  fi

  pdb_count=0
  if [ -d "$outdir" ]; then
    pdb_count="$(find "$outdir" -maxdepth 1 -type f -name '*.pdb' | wc -l | tr -d ' ')"
  fi
  if [ "$pdb_count" -gt 0 ]; then
    gate_line="$(run_gate_for_set "$set_name" "$run_id" "$contig" "$outdir")"
    if [ "$pid_alive" -eq 1 ]; then
      write_queue_status status=GATED_RUNNING current_set="$set_name" current_run_id="$run_id" pid="$pid" pdb_count="$pdb_count" gate="$gate_line" status_file="$status_file" next_action="continue_monitoring_current_l1_and_regate_present_pdbs"
      exit 0
    fi
    write_queue_status status=GATED current_set="$set_name" current_run_id="$run_id" pdb_count="$pdb_count" gate="$gate_line" status_file="$status_file" next_action="continue_to_next_contig_or_sequence_generation_if_enough_PASS"
    continue
  fi

  if [ "$pid_alive" -eq 1 ]; then
    write_queue_status status=MONITORING current_set="$set_name" current_run_id="$run_id" pid="$pid" pdb_count="$pdb_count" status_file="$status_file" next_action="wait_for_current_l1_to_write_pdb_or_finish"
    exit 0
  fi

  CONTIG_SET="$set_name" NUM_DESIGNS="$NUM_DESIGNS" DESIGN_STARTNUM="$startnum" "$LAUNCHER" "$ROOT" >/tmp/${run_id}_launch_path.txt
  launched_status_file="$(tail -n 1 /tmp/${run_id}_launch_path.txt)"
  launched_status="$(json_value "$launched_status_file" status || true)"
  launched_pid="$(json_value "$launched_status_file" pid || true)"
  if [[ "$launched_status" == BLOCKED* ]]; then
    write_queue_status status="$launched_status" current_set="$set_name" current_run_id="$run_id" status_file="$launched_status_file" next_action="rerun_queue_when_gpu_is_available"
    exit 0
  fi
  if [ "$launched_status" = "LAUNCHED" ]; then
    write_queue_status status=LAUNCHED current_set="$set_name" current_run_id="$run_id" pid="$launched_pid" status_file="$launched_status_file" next_action="monitor_launched_l1"
    exit 0
  fi

  write_queue_status status=UNKNOWN current_set="$set_name" current_run_id="$run_id" status_file="$launched_status_file" observed_status="$launched_status" next_action="inspect_launcher_status"
  exit 1
done

write_queue_status status=DONE contig_sets="$CONTIG_SETS" next_action="inspect_gate_summaries_and_start_multiscaffold_sequence_generation"