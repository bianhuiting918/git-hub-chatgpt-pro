# Nylonase pose-retention × ProteinMPNN × SaProt × activity RUNBOOK

Production root: `/work/home/acshdt1dks/petase_orbmol_lg1_lg4_layer8343_20260721/outputs/nylon_pose_mpnn_saprot_activity_20260727_v1`.

## Authority and frozen snapshot

The experimental authority is Nyl01–Nyl95 (95 exact sequence MD5s). The pose input is the 2026-07-22 Dell-derived frozen snapshot: 53 enzymes with a complete 180-pose AD-dimer Smina grid and 42 reason-specific technical exclusions (27 queued, 15 without a frozen D-D-T receptor). A live Dell refresh was attempted before this run; TCP/22 opened but the remote closed before SSH key exchange. Therefore the 53/42 split is explicitly snapshot provenance, not a current Dell completion claim.

## Scores and interpretation

ProteinMPNN uses exact canonical sequence rows with `family=Nylonase`, `proteinmpnn_status=PASS`, `scoring_sequence_scope=CANONICAL_EXACT`, and vanilla v_48_020 mean log-likelihood (soluble is sensitivity). SaProt uses exact canonical `FULL_PROTEIN` PASS mean log-likelihood. These are compatibility scores, not Tm, ΔG, activity, or barriers. Experimental spans are raster color-span proxies, not absolute activity.

## Run

```bash
cd /work/home/acshdt1dks/petase_orbmol_lg1_lg4_layer8343_20260721/outputs/nylon_pose_mpnn_saprot_activity_20260727_v1/code
export PYTHONPATH="$PWD/../deps:/work/home/acshdt1dks/anaconda3/envs/matplotlib/lib/python3.9/site-packages:$PWD"
pytest -q tests/nylon_pose_mpnn_saprot_activity
sbatch analysis/nylon_pose_mpnn_saprot_activity/slurm.sbatch
```

## Outputs and recovery

`authority_joined.tsv` is the enzyme-level exact-MD5 authority table. `not_evaluated.tsv` preserves technical reasons. Per-endpoint denominators, correlations, bootstrap CIs, standardized models, predictor-collinearity audit, sensitivity correlations and figures are written at the root. `input_sha256.tsv`, `logs/run_history.tsv`, `audit.json`, and `PASS.json` provide provenance. Re-run is idempotent: an existing root `PASS.json` causes a recorded skip. A failed run may be resumed only after inspecting Slurm logs; do not delete or overwrite the directory.
