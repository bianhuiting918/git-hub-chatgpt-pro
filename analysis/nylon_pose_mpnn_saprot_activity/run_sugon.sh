#!/usr/bin/env bash
set -euo pipefail
ROOT=/work/home/acshdt1dks/petase_orbmol_lg1_lg4_layer8343_20260721/outputs/nylon_pose_mpnn_saprot_activity_20260727_v1
PY=/work/home/acshdt1dks/anaconda3/envs/pytorch/bin/python
export PYTHONPATH="$ROOT/deps:/work/home/acshdt1dks/anaconda3/envs/matplotlib/lib/python3.9/site-packages:$ROOT/code:${PYTHONPATH:-}"
mkdir -p "$ROOT/logs" "$ROOT/figures" "$ROOT/clusters"
if [[ -f "$ROOT/PASS.json" ]]; then
  printf '%s\t%s\t%s\t%s\n' "$(date -Is)" "${SLURM_JOB_ID:-manual}" "SKIP_EXISTING_PASS" "$ROOT/PASS.json" >> "$ROOT/logs/run_history.tsv"
  exit 0
fi
printf '%s\t%s\t%s\t%s\n' "$(date -Is)" "${SLURM_JOB_ID:-manual}" "START" "run_analysis n_boot=10000" >> "$ROOT/logs/run_history.tsv"
set +e
"$PY" -m analysis.nylon_pose_mpnn_saprot_activity.run_analysis \
 --activity "$ROOT/inputs/chemcatal2025_nyl95_activity.tsv" \
 --pose "$ROOT/inputs/pose_enzyme_retention_summary.tsv" \
 --pose-strata "$ROOT/inputs/pose_enzyme_retention_by_stratum.tsv" \
 --pose-not-evaluated "$ROOT/inputs/pose_not_evaluated.tsv" \
 --mpnn /work/home/acshdt1dks/petase_orbmol_lg1_lg4_layer8343_20260721/outputs/proteinmpnn_compatibility_pet8329_nylon4167_20260727_v2_final/proteinmpnn_candidate_scores.tsv \
 --saprot /work/home/acshdt1dks/petase_orbmol_lg1_lg4_layer8343_20260721/manifests/proteinmpnn_compatibility_pet8329_nylon4167_20260727_v2_scopefix/run_manifest_all_12496_scopefixed.tsv \
 --output "$ROOT" --n-boot 10000 --seed 20260727
rc=$?
if [[ $rc -eq 0 ]]; then
  "$PY" -m analysis.nylon_pose_mpnn_saprot_activity.audit_run "$ROOT"
  rc=$?
fi
set -e
printf '%s\t%s\t%s\t%s\n' "$(date -Is)" "${SLURM_JOB_ID:-manual}" "END_RC_${rc}" "$ROOT" >> "$ROOT/logs/run_history.tsv"
exit "$rc"
