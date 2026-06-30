#!/usr/bin/env bash
set -euo pipefail

cd /data/bht/design_tools/src/PLACER

CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-0} /data/bht/design_tools/envs/placer_env/bin/python run_PLACER.py \
  --ifile /data/bht/project01_phase1_reset_gpu/manifests/placer_n50_holo_inputs_capped10.ifile \
  --odir /data/bht/project01_phase1_reset_gpu/placer_runs/capped10_n50_bonded \
  --suffix bonded_bu2_n50 \
  --nsamples 50 \
  --mutate 128A:75I \
  --predict_ligand X-bu2-1 \
  --bonds A-128-75I-OG:X-1-bu2-C1:1.533 \
  --residue_json /data/bht/design_tools/src/PLACER/examples/ligands/75I.json \
  --rerank prmsd \
  --cautious
