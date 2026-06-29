# Project 01 Phase1 Full-Ligand Mapping Status - 2026-06-30

## Current evaluated universe

PLACER pilot is complete for 5 rank-1 representatives, one per similarity bin: 90, 80, 70, 60, 50. Each bin has 3 PLACER conformers, for 15 crop models total.

- Pilot root: `/Dell/Dell14/bianht/project01_similarity_phase1_20260630/placer_pilot_rank1_20260630_0148`
- PLACER split manifest: `/Dell/Dell14/bianht/project01_similarity_phase1_20260630/placer_pilot_rank1_20260630_0148/placer_pilot_split_manifest.tsv`
- Full-ligand mapping root: `/Dell/Dell14/bianht/project01_similarity_phase1_20260630/placer_pilot_rank1_20260630_0148/full_ligand_mapped_no_min`
- Full-ligand mapping manifest: `/Dell/Dell14/bianht/project01_similarity_phase1_20260630/placer_pilot_rank1_20260630_0148/full_ligand_mapped_no_min/full_ligand_mapping_manifest.tsv`

## Completed

A new mapping script was created and run:

`/Dell/Dell14/bianht/project01_similarity_phase1_20260630/placer_pilot_rank1_20260630_0148/scripts/map_placer_crop_to_full_ligand.py`

It produced 15 mapped structures:

- Output file per model: `full_length_full_bu2_no_min.pdb`
- Full bu2 ligand atom count: 32 atoms
- PLACER ligand heavy atoms transferred: 18 atoms
- Ligand hydrogens restored by rigid fitting from the original ligand template: 14 atoms
- 75I/QEX extra atoms retained separately as `qex_extra_75i_from_placer.xyz`

## Strict geometry smoke

A strict all-atom pocket geometry smoke was run on one mapped model:

`/Dell/Dell14/bianht/project01_similarity_phase1_20260630/placer_pilot_rank1_20260630_0148/full_ligand_mapped_no_min/bin90_rank1_rec1_model_001/strict_geometry_smoke.json`

Result: `fail`, as expected at this stage.

Main reason: the current PLACER/LigandMPNN rank1 structures are still backbone-stage protein structures. Side-chain heavy atoms such as His95 ND1/NE2 and Ser126/Ser128 CB/OG are missing, so strict all-atom key-distance checks cannot yet pass.

## Next required step

Before QMMM-ready inputs can be claimed, run side-chain heavy-atom completion/packing for the 15 mapped structures, then rerun strict geometry filtering. OpenMM and Bio.PDB are available in the `orbmol_gpu` and `protrek_gpu` environments, but PDBFixer, SCWRL, PULCHRA, and Modeller were not found on the CPU server at this checkpoint.

## Audit labels

- PLACER output: complete, 15 conformers.
- Full bu2 ligand mapping: complete, 15 mapped structures.
- Strict all-atom pocket geometry: not passed yet because side-chain heavy atoms are missing.
- QMMM-ready status: blocked pending side-chain heavy completion/packing.
