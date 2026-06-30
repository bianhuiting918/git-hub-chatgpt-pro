#!/usr/bin/env bash
set -euo pipefail
ROOT=/data/bht/project01_phase1_reset_gpu/denovo_scaffold_routeB/serine_hydrolase
RFAA=/data/bht/design_tools/src/rf_diffusion_all_atom
PY=/data/bht/design_tools/envs/rfaa_venv/bin/python
INPUT=/data/bht/project01_phase1_reset_gpu/manifests/denovo_SER_hydrolase_full_input.pdb
CKPT=/data/bht/design_tools/src/rf_diffusion_all_atom/RFDiffusionAA_paper_weights.pt
if [[ ! -s "$CKPT" ]]; then
  CKPT=/data/bht/design_tools/models/rfdiffusion_aa/RFDiffusionAA_paper_weights.pt
fi
OUT=$ROOT/outputs/rfaa_serhyd_motif_H95_E102_S128_bu2_smoke/sample
cd "$RFAA"
export CUDA_VISIBLE_DEVICES=0
export HYDRA_FULL_ERROR=1
"$PY" run_inference.py \
  inference.deterministic=True \
  inference.input_pdb="$INPUT" \
  inference.ckpt_path="$CKPT" \
  inference.output_prefix="$OUT" \
  inference.ligand=bu2 \
  inference.num_designs=1 \
  inference.design_startnum=0 \
  diffuser.T=20 \
  contigmap.contigs=[\'20-40,A95-95,5-15,A102-102,5-20,A128-128,20-40\'] \
  contigmap.length=90-140
