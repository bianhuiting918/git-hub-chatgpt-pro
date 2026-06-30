# 2026-06-30 Blind PETase QM/MM Mechanism Task

## Status

A new blind first-principles PETase QM/MM mechanism task was created.

## Scope correction

The task is not to reproduce the paper by consuming the paper's concrete trajectory/TS/RC results. The task is to reproduce the scientific discovery process from PETase structure and substrate chemistry, using only the paper's broad methodology as inspiration.

## Uploaded documents

- `projects/01-specialized-ts-aware-scorer/docs/petase_blind_qmmm_mechanism_plan.md`
- `project01/phase2_blind_petase_qmmm_20260630/README.md`

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

## Next action

Start Stage 1: structure-only system setup.

Immediate tasks:

1. Select WT PETase structure templates from public structural databases.
2. Build PET dimer or BHET/MHET-like ester substrate models.
3. Generate candidate Michaelis complexes using only generic catalytic geometry.
4. Equilibrate and filter poses before any QM/MM reaction scans.
