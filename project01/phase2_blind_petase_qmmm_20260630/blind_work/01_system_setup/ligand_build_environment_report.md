# Ligand Build Environment Report

Date: 2026-06-30

## Boundary

This report concerns only Stage 1 ligand construction environment readiness. It uses no paper-derived TS, RC, CV, trajectory, barrier, or mechanism result.

## Python Runtime

- Python: `3.12.13`
- RDKit import available: `no`
- Open Babel Python import available: `no`
- Biopython import available: `no`

## Status

`blocked_missing_chemistry_toolkit`

Ligand 3D conformer generation and atom-label validation require RDKit or Open Babel. The current bundled Python environment has neither, so this stage is intentionally stopped before producing unreliable ligand coordinates.

## Next Accepted Routes

1. Install or expose RDKit/Open Babel in the compute environment, then generate SDF/MOL2 from `ligand_model_definitions.md`.
2. Use an existing cluster chemistry stack, but record the exact executable path, version, command, and generated atom-label table.
3. If ligand coordinates are built manually in a GUI, export SDF/MOL2 plus an atom-label TSV and mark the build as manual in `ligand_model_manifest.tsv`.

## Grill Gate

Do not start docking or QM/MM with ligand coordinates unless the scissile ester labels in `ligand_model_manifest.tsv` are mapped to concrete atom names in the generated topology.
