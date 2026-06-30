# Blind PETase Stage 1 PDB Preparation Log

Date: 2026-06-30

## Boundary

This log records structure-only preparation steps for the blind PETase QM/MM mechanism workflow. No paper TS coordinates, RCs, selected CVs, umbrella windows, shooting trajectories, barrier values, or mechanistic conclusions were used.

## Completed

1. Queried RCSB metadata for candidate PETase structures.
2. Selected `6EQE` as the primary WT-like template, with `5XJH`, `5YFE`, and `6ILW` as secondary WT-like sensitivity templates.
3. Downloaded nine WT-like or backup PDB coordinate files from RCSB:
   - `6EQE`, `5XJH`, `5YFE`, `6ILW`, `6EQD`, `6EQF`, `6EQG`, `6EQH`, `6QGC`.
4. Ran a coordinate-only Ser-His-Asp geometric scan with these thresholds:
   - Ser OG to His ND1/NE2 <= 4.0 A;
   - His ND1/NE2 to Asp OD1/OD2 <= 4.0 A.
5. Ran PDB preparation audit for chains, missing residues, missing atoms, alternate conformers, non-water heterogens, crystallographic water, and geometric disulfide candidates.
6. Wrote explicit template-chain decisions and retained-water candidate tables.

## Triad Scan Result

The scan identifies a consistent catalytic triad geometry across the WT-like templates:

- `6EQE`: chain A `SER160-HIS237-ASP206`.
- `5XJH`: chain A `SER160-HIS237-ASP206`.
- `6ILW`: chain A `SER160-HIS237-ASP206`.
- `5YFE`: chain A `SER134-HIS211-ASP180`, likely a residue-numbering offset relative to the other templates.
- Multi-chain crystal forms report the same local geometry per chain.

This is a structure-derived active-site assignment, not a paper-derived mechanism result.

## PDB Preparation Audit Result

- Initial chain decision: use chain A for the primary blind setup because chain A contains the consistent Ser-His-Asp triad in every selected template.
- `6EQE` has a complete triad, 33 missing residues reported in the PDB header, no missing atoms reported, 26 alternate-conformer residues, two geometric Cys-Cys pairs, and two crystallographic water candidates within 4 A of catalytic hetero atoms.
- `5XJH` has a complete triad, 37 missing residues reported, no missing atoms or alternate-conformer residues, two geometric Cys-Cys pairs, and three retained-water candidates.
- `5YFE` uses shifted residue numbering, has four missing C-terminal His residues reported, glycerol and sulfate heterogens, two geometric Cys-Cys pairs, and no retained-water candidates within the 4 A triad cutoff.
- Backup multi-chain structures contain repeated chain copies and should remain sensitivity/backup templates rather than primary production templates.

## Current Decisions

1. Primary production template remains `6EQE`, chain A.
2. Secondary sensitivity templates remain `5XJH`, `5YFE`, and `6ILW`, chain A.
3. Geometric disulfide candidates should be kept unless later preparation software contradicts the connectivity.
4. Water molecules listed in `retained_water_candidates.tsv` should be retained for the first local relaxation pass, then tested by sensitivity runs.
5. Alternate conformers in `6EQE` must be resolved before ligand placement; default starting policy is highest occupancy, then active-site geometry inspection.

## Not Yet Done

- Actual repaired PDB/mmCIF generation.
- Alternate conformer selection in coordinate files.
- Crystal-water retention/removal implementation in coordinate files.
- Protonation-state assignment.
- Substrate model construction.
- Michaelis-complex generation and relaxation.

## Evidence Files

- `structure_selection.tsv`
- `structure_download_manifest.csv`
- `ser_his_asp_triad_candidates.tsv`
- `pdb_preparation_audit.tsv`
- `template_chain_decisions.tsv`
- `retained_water_candidates.tsv`
- `stage1_system_setup_protocol.md`
