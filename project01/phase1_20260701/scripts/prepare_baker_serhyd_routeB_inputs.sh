#!/usr/bin/env bash
set -euo pipefail

ROOT="${1:-/data/bht/project01_baker_serhyd_routeB_20260701}"
SRC_URL="https://github.com/laukoag/serine-hydrolase-design.git"
SRC_DIR="$ROOT/external/serine-hydrolase-design"
INPUT_DIR="$ROOT/inputs/baker_serhyd"
MANIFEST="$ROOT/manifests/baker_serhyd_routeB_input_manifest.tsv"
STATUS="$ROOT/manifests/baker_serhyd_routeB_prepare_status.json"

mkdir -p "$ROOT/external" "$INPUT_DIR" "$ROOT/manifests" "$ROOT/logs" "$ROOT/scripts"

if [ ! -d "$SRC_DIR/.git" ]; then
  git clone --depth 1 "$SRC_URL" "$SRC_DIR"
else
  git -C "$SRC_DIR" fetch --depth 1 origin
  git -C "$SRC_DIR" reset --hard origin/main
fi

copy_input() {
  local role="$1"
  local rel="$2"
  local residue="$3"
  local dest="$INPUT_DIR/$rel"
  mkdir -p "$(dirname "$dest")"
  cp "$SRC_DIR/$rel" "$dest"
  printf '%s\t%s\t%s\t%s\t%s\n' "$role" "$residue" "$SRC_DIR/$rel" "$dest" "$(wc -c < "$dest")" >> "$MANIFEST"
}

printf 'role\tresidue_or_ligand\tsource_path\tstaged_path\tbytes\n' > "$MANIFEST"
copy_input "motif_generation_ts_like_params" "motif_gen/01_sampling_his_stub/inputs/mu1.params" "mu1"
copy_input "motif_generation_constraint" "motif_gen/01_sampling_his_stub/inputs/1LNS_mu1.cst" "mu1_SER_HIS"
copy_input "natural_ser_stub_reference" "motif_gen/01_sampling_his_stub/inputs/1LNS.pdb" "SER_stub"
copy_input "simple_theozyme_pdb" "motif_gen/02_diffusion/inputs/simple_theozyme.pdb" "mu1_SER_HIS"
copy_input "simple_theozyme_diffusion_tasks" "motif_gen/02_diffusion/inputs/tasks.json" "mu1"
copy_input "design_stage_substrate_params" "design_pipeline/03_design/inputs/bn1.params" "bn1"
copy_input "design_stage_theozyme_constraint" "design_pipeline/03_design/inputs/theozyme.cst" "bn1_SER_HIS_ASP_OXH"
copy_input "design_stage_theozyme_pdb" "design_pipeline/01_diffusion/inputs/theozyme.pdb" "bn1_SER_HIS_ASP_OXH"
copy_input "design_stage_diffusion_tasks" "design_pipeline/01_diffusion/inputs/tasks.json" "bn1"
copy_input "substrate_bound_gs_like_reference" "design_pipeline/06_PLACER/inputs/super_af2_bu2.pdb" "bu2"

gpu_status="NO_NVIDIA_SMI"
if command -v nvidia-smi >/dev/null 2>&1; then
  gpu_status="$(nvidia-smi --query-gpu=utilization.gpu,memory.used,memory.total --format=csv,noheader,nounits | head -1 | tr -d ' ')"
fi

cat > "$STATUS" <<JSON
{
  "status": "PREPARED_INPUTS_ONLY",
  "route_label": "baker_theozyme_new_backbone",
  "root": "$ROOT",
  "source_repo": "$SRC_URL",
  "source_commit": "$(git -C "$SRC_DIR" rev-parse HEAD)",
  "manifest": "$MANIFEST",
  "gpu_status_util_mem_used_mem_total": "$gpu_status",
  "next_action": "Launch a small RFAA/RFdiffusionAA smoke only when GPU capacity is available; do not use old fixed-reference-pocket route as completion evidence."
}
JSON

echo "$STATUS"
