#!/usr/bin/env bash
#SBATCH --job-name=polymer_manifest_v2
#SBATCH --partition=xahcnormal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=02:00:00
#SBATCH --output=/work/home/acshdt1dks/polymer_surface_hotspot_screen_20260725/logs/%x_%j.out
#SBATCH --error=/work/home/acshdt1dks/polymer_surface_hotspot_screen_20260725/logs/%x_%j.err

set -euo pipefail
ROOT=/work/home/acshdt1dks/polymer_surface_hotspot_screen_20260725
DEPLOY_COMMIT=${DEPLOY_COMMIT:?DEPLOY_COMMIT is required}
OUTPUT_DIR=${OUTPUT_DIR:?OUTPUT_DIR is required}
PROJECT="$ROOT/deployments/$DEPLOY_COMMIT/projects/03-polymer-surface-hotspot-screening"
PYTHON="$ROOT/envs/surface-screen-py311/bin/python3.11"
SCRIPT="$PROJECT/scripts/build_structure_manifest.py"

case "$OUTPUT_DIR" in
  "$ROOT"/inputs/*) ;;
  *)
    printf '%s FAIL_OUTPUT_OUTSIDE_INPUTS path=%s\n' "$(date -Is)" "$OUTPUT_DIR"
    exit 2
    ;;
esac

printf '%s START job_id=%s host=%s commit=%s output=%s\n' \
  "$(date -Is)" "${SLURM_JOB_ID:-NA}" "$(hostname)" "$DEPLOY_COMMIT" "$OUTPUT_DIR"

if [ ! -f "$SCRIPT" ]; then
  printf '%s FAIL_MISSING_SCRIPT path=%s\n' "$(date -Is)" "$SCRIPT"
  exit 2
fi
if [ -e "$OUTPUT_DIR" ]; then
  printf '%s FAIL_EXISTING_OUTPUT path=%s\n' "$(date -Is)" "$OUTPUT_DIR"
  exit 2
fi

"$PYTHON" "$SCRIPT" \
  --pet-manifest /work/home/acshdt1dks/petase_orbmol_lg1_lg4_layer8343_20260721/outputs/saprot_structure_conditional_likelihood_20260723/run_manifest.tsv \
  --nylon-authority-manifest /work/home/acshdt1dks/petase_orbmol_lg1_lg4_layer8343_20260721/manifests/saprot_nylonase_4556_20260725_confscale_v2/preflight_status.tsv \
  --nylon-manifest /work/home/acshdt1dks/petase_orbmol_lg1_lg4_layer8343_20260721/manifests/saprot_nylonase_4556_20260725_confscale_v2/run_manifest_all_3497.tsv \
  --nylon-recovered-manifest /work/home/acshdt1dks/petase_orbmol_lg1_lg4_layer8343_20260721/manifests/saprot_nylonase_recovered_20260725_v7/run_manifest_ready_735.tsv \
  --nylon-extra-manifest /work/home/acshdt1dks/petase_orbmol_lg1_lg4_layer8343_20260721/manifests/saprot_nylonase_recovered_20260725_v7/not_evaluated_52.tsv \
  --output-dir "$OUTPUT_DIR"

sha256sum "$OUTPUT_DIR"/*.tsv "$OUTPUT_DIR"/*.json > "$OUTPUT_DIR/SHA256SUMS"
printf '%s PASS_STRUCTURE_MANIFEST_V2 output=%s\n' "$(date -Is)" "$OUTPUT_DIR"
