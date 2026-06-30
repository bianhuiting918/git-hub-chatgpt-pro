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

## Triad Scan Result

The scan identifies a consistent catalytic triad geometry across the WT-like templates:

- `6EQE`: chain A `SER160-HIS237-ASP206`.
- `5XJH`: chain A `SER160-HIS237-ASP206`.
- `6ILW`: chain A `SER160-HIS237-ASP206`.
- `5YFE`: chain A `SER134-HIS211-ASP180`, likely a residue-numbering offset relative to the other templates.
- Multi-chain crystal forms report the same local geometry per chain.

This is a structure-derived active-site assignment, not a paper-derived mechanism result.

## Not Yet Done

- Missing residue/atom audit.
- Alternate conformer selection.
- Crystal-water retention/removal decisions.
- Protonation-state assignment.
- Disulfide verification.
- Substrate model construction.
- Michaelis-complex generation and relaxation.

## Evidence Files

- `structure_selection.tsv`
- `structure_download_manifest.csv`
- `ser_his_asp_triad_candidates.tsv`
- `stage1_system_setup_protocol.md`
