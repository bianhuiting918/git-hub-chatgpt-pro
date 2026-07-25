#!/usr/bin/env bash
#SBATCH --job-name=polymer_manifest
#SBATCH --partition=xahcnormal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=/work/home/acshdt1dks/polymer_surface_hotspot_screen_20260725/logs/%x_%j.out
#SBATCH --error=/work/home/acshdt1dks/polymer_surface_hotspot_screen_20260725/logs/%x_%j.err

set -euo pipefail
ROOT=/work/home/acshdt1dks/polymer_surface_hotspot_screen_20260725
REPO="$ROOT/repo"
PYTHON=/work/home/acshdt1dks/python3.11/bin/python3.11
OUT="$ROOT/inputs/frozen_structure_manifests_20260726"

printf '%s START job_id=%s host=%s\n' "$(date -Is)" "${SLURM_JOB_ID:-NA}" "$(hostname)"
if [ -e "$OUT" ]; then
  printf '%s FAIL_EXISTING_OUTPUT path=%s\n' "$(date -Is)" "$OUT"
  exit 2
fi

"$PYTHON" "$REPO/projects/03-polymer-surface-hotspot-screening/scripts/build_structure_manifest.py" \
  --pet-manifest /work/home/acshdt1dks/petase_orbmol_lg1_lg4_layer8343_20260721/outputs/saprot_structure_conditional_likelihood_20260723/run_manifest.tsv \
  --nylon-manifest /work/home/acshdt1dks/petase_orbmol_lg1_lg4_layer8343_20260721/manifests/saprot_nylonase_4556_20260725_confscale_v2/run_manifest_all_3497.tsv \
  --nylon-recovered-manifest /work/home/acshdt1dks/petase_orbmol_lg1_lg4_layer8343_20260721/manifests/saprot_nylonase_recovered_20260725_v7/run_manifest_ready_735.tsv \
  --nylon-extra-manifest /work/home/acshdt1dks/petase_orbmol_lg1_lg4_layer8343_20260721/manifests/saprot_nylonase_recovered_20260725_v7/not_evaluated_52.tsv \
  --output-dir "$OUT"

sha256sum "$OUT"/*.tsv "$OUT"/*.json > "$OUT/SHA256SUMS"
printf '%s PASS output=%s\n' "$(date -Is)" "$OUT"
