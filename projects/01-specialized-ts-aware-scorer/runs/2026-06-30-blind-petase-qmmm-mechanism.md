# 2026-06-30 Blind PETase QM/MM Mechanism Task

## Status

A new blind first-principles PETase QM/MM mechanism task was created and Stage 1 structure-only setup has started.

## Scope correction

The task is not to reproduce the paper by consuming the paper's concrete trajectory/TS/RC results. The task is to reproduce the scientific discovery process from PETase structure and substrate chemistry, using only the paper's broad methodology as inspiration.

## Uploaded documents

- `projects/01-specialized-ts-aware-scorer/docs/petase_blind_qmmm_mechanism_plan.md`
- `project01/phase2_blind_petase_qmmm_20260630/README.md`
- `project01/phase2_blind_petase_qmmm_20260630/blind_work/01_system_setup/structure_selection.tsv`
- `project01/phase2_blind_petase_qmmm_20260630/blind_work/01_system_setup/stage1_system_setup_protocol.md`
- `project01/phase2_blind_petase_qmmm_20260630/blind_work/01_system_setup/gs_pose_manifest.tsv`
- `project01/phase2_blind_petase_qmmm_20260630/blind_work/01_system_setup/ligand_model_manifest.tsv`
- `project01/phase2_blind_petase_qmmm_20260630/blind_work/01_system_setup/ligand_model_definitions.md`
- `project01/phase2_blind_petase_qmmm_20260630/blind_work/01_system_setup/structure_download_manifest.csv`
- `project01/phase2_blind_petase_qmmm_20260630/blind_work/01_system_setup/ser_his_asp_triad_candidates.tsv`
- `project01/phase2_blind_petase_qmmm_20260630/blind_work/01_system_setup/pdb_preparation_log.md`
- `project01/phase2_blind_petase_qmmm_20260630/blind_work/01_system_setup/pdb_preparation_audit.tsv`
- `project01/phase2_blind_petase_qmmm_20260630/blind_work/01_system_setup/template_chain_decisions.tsv`
- `project01/phase2_blind_petase_qmmm_20260630/blind_work/01_system_setup/retained_water_candidates.tsv`
- `project01/phase2_blind_petase_qmmm_20260630/scripts/download_stage1_rcsb_structures.ps1`
- `project01/phase2_blind_petase_qmmm_20260630/scripts/identify_ser_his_asp_triads.py`
- `project01/phase2_blind_petase_qmmm_20260630/scripts/audit_stage1_pdb_preparation.py`

## Boundary

Before final validation, do not use:

- paper TS coordinates;
- paper reaction-coordinate formulas;
- paper selected CVs;
- paper umbrella windows;
- paper aimless-shooting trajectories;
- paper barrier values, rate constants, or rate-limiting-step assignment;
- paper conclusions about His motion, Asp role, Trp185 role, oxyanion-hole role, or tetrahedral TS/intermediate status.

Use only:

- PETase structures;
- substrate chemistry;
- generic serine-hydrolase chemistry;
- QM/MM, TS-search, PMF, committor, and trajectory-validation methods.

## Compute note

A previous Zenodo paper-trajectory download on the CPU server was stopped after the boundary correction, because paper shooting trajectories are concrete results and should not enter the blind workflow. The partial archive fragments were not deleted.

## Stage 1 progress

RCSB metadata was queried on 2026-06-30 for PETase structure candidates. The initial blind production template set is:

- primary: `6EQE`;
- secondary sensitivity templates: `5XJH`, `5YFE`, `6ILW`;
- backup WT-like templates: `6EQD`, `6EQF`, `6EQG`, `6EQH`, `6QGC`.

Mutant, non-PETase, and failed pre-query IDs were explicitly marked as excluded in `structure_selection.tsv`.

The selected WT-like/backup RCSB coordinate files were downloaded locally. A coordinate-only Ser-His-Asp geometric scan found the same active-site triad geometry across the WT-like templates, with `6EQE`, `5XJH`, and `6ILW` using chain A `SER160-HIS237-ASP206`, and `5YFE` using a numbering-shifted chain A `SER134-HIS211-ASP180`.

A PDB preparation audit was completed. It records selected chain decisions, missing residues/atoms from PDB headers, alternate conformers, non-water heterogens, geometric disulfide candidates, water counts, and crystallographic waters within 4 A of catalytic hetero atoms.

Blind substrate model definitions were added for:

- `PET_dimer_capped` as the primary PET-like acylation substrate;
- `BHET_like` as a small neutral docking/pose control;
- `MHET_like` as a product-side/deacylation reference fragment with explicit pH-7 carboxylate state;
- `MHET_like_acyl_enzyme_precursor` as a protein-covalent Ser160 acyl-enzyme precursor for deacylation setup.

Current structure-prep decisions:

- use chain A for the initial blind setup across selected templates;
- keep geometric disulfide candidates unless preparation software contradicts connectivity;
- retain listed catalytic-site water candidates for first local relaxation, then test sensitivity;
- resolve `6EQE` alternate conformers before substrate placement;
- preserve scissile ester atom labels through ligand 3D generation and topology conversion.

This is structure/substrate-derived preparation evidence and is not a paper-derived mechanism result.

## Next action

Continue Stage 1 by implementing repaired/protonation-ready coordinate files, assigning protonation states, generating 3D ligand conformers/parameters, and filling `gs_pose_manifest.tsv` with accepted and rejected Michaelis-complex candidates.
