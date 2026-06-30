#!/usr/bin/env bash
set -euo pipefail
ROOT=/data/bht/project01_phase1_reset_gpu/denovo_scaffold_routeB/serine_hydrolase
LMPNN=/data/bht/design_tools/src/LigandMPNN
PY=/data/bht/design_tools/envs/ligandmpnn_venv/bin/python
PDB=$ROOT/outputs/rfaa_serhyd_motif_H95_E102_S128_bu2_smoke/sample_0.pdb
OUT=$ROOT/outputs/ligandmpnn_serhyd_rfaa_smoke_n2
cd "$LMPNN"
export CUDA_VISIBLE_DEVICES=0
"$PY" run.py \
  --model_type "ligand_mpnn" \
  --checkpoint_ligand_mpnn "./model_params/ligandmpnn_v_32_005_25.pt" \
  --pdb_path "$PDB" \
  --out_folder "$OUT" \
  --seed 111 \
  --temperature 0.1 \
  --batch_size 2 \
  --number_of_batches 1 \
  --fixed_residues "A37 A50 A68" \
  --ligand_mpnn_use_side_chain_context 1 \
  --save_stats 1
