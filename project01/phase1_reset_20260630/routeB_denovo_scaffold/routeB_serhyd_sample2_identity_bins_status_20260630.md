# Route B Serine Hydrolase Active-Pocket Identity Bins - 2026-06-30

## Status

target_complete: true
evaluated_total: 710
pass_total: 93
accepted_total: 50

## Gate

- Stage: sequence generation only; PLACER/QMMM are not acceptance gates here.
- Pocket definition for this sequence-stage run: active_site plus residues with CA-to-bu2 <= 8 A.
- Rationale: RFAA sidechain/heavy atoms in sample_2 contain diagnostic-only placeholder-like coordinates, so heavy-atom 4 A would incorrectly fix nearly the whole protein.
- Hard gate: ESMFold OK, generated fixed positions unchanged, Kabsch-aligned motif CA RMSD over A40/A49/A56 <= 1.0 A.
- Global CA RMSD and sidechain/ligand geometry are not hard gates at this stage.

## Accepted By Bin

| identity_bin | evaluated | pass | selected | best_motif_ca_rmsd_A | worst_selected_motif_ca_rmsd_A |
| --- | ---: | ---: | ---: | ---: | ---: |
| 90 | 70 | 14 | 10 | 0.1713 | 0.8496 |
| 80 | 200 | 22 | 10 | 0.4067 | 0.8388 |
| 70 | 140 | 16 | 10 | 0.2493 | 0.8006 |
| 60 | 140 | 21 | 10 | 0.2676 | 0.6517 |
| 50 | 160 | 20 | 10 | 0.4339 | 0.7155 |

## Lightweight Evidence Files

- serhyd_sample2_combined_summary.json
- serhyd_sample2_accepted_top10.tsv
- serhyd_sample2_accepted_top10.fasta
- serhyd_sample2_generate_round01.py
- serhyd_sample2_generate_round02_from_passes.py
- serhyd_sample2_gate.py

Large files intentionally not synced: ESMFold PDBs, raw logs, model caches, and trajectories.
