#!/usr/bin/env bash
set -euo pipefail

ROOT=/data/bht/project01_phase1_reset_gpu
PLACER=/data/bht/design_tools/src/PLACER
PYTHON=/data/bht/design_tools/envs/placer_env/bin/python
IFILE="$ROOT/manifests/placer_n50_holo_inputs_capped10.ifile"
ODIR="$ROOT/placer_runs/capped10_n50_bonded"
LOGDIR="$ROOT/placer_runs/capped10_n50_bonded/per_input_logs"

mkdir -p "$ODIR" "$LOGDIR"
cd "$PLACER"

idx=0
while IFS= read -r pdb; do
  idx=$((idx+1))
  [[ -z "$pdb" ]] && continue

  base=$(basename "$pdb" .pdb)
  csv="$ODIR/${base}_bonded_bu2_n50.csv"

  if [[ -s "$csv" ]]; then
    echo "[$(date '+%F %T')] SKIP existing $idx $base"
    continue
  fi

  echo "[$(date '+%F %T')] START $idx $base"
  CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-0} "$PYTHON" run_PLACER.py \
    --ifile "$pdb" \
    --odir "$ODIR" \
    --suffix bonded_bu2_n50 \
    --nsamples 50 \
    --mutate 128A:75I \
    --predict_ligand X-bu2-1 \
    --bonds A-128-75I-OG:X-1-bu2-C1:1.533 \
    --residue_json /data/bht/design_tools/src/PLACER/examples/ligands/75I.json \
    --rerank prmsd \
    --cautious > "$LOGDIR/${base}.log" 2>&1
  echo "[$(date '+%F %T')] DONE $idx $base"
done < "$IFILE"
