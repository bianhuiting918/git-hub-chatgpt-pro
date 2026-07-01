# Blind Stage 1/2 Gate Next Actions

This runner records the current blind workflow state from structure-derived inputs.

Boundary: do not use article trajectories, article transition-state coordinates, article reaction coordinates, article selected CVs, article barriers, or article mechanism conclusions.

## Current Blocking Items

- None recorded by this runner.

## Required Order

1. Confirm compute-environment tools and versions.
2. Build ligand conformers and atom-label tables from blind substrate definitions.
3. Run protonation review on cleaned PETase coordinates.
4. Generate blind pose-generation queues from prepared structures and ligand atom labels.
5. Generate and score Stage 1 ground-state poses.
6. Queue Stage 2 classical MD only from accepted ground-state poses.
7. Do not activate Stage 4 QM/MM scans until productive Stage 2 conformers exist.
