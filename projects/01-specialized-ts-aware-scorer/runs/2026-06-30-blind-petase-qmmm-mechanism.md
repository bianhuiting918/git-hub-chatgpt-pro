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

## Next action

Continue Stage 1 by downloading the selected RCSB structures, preparing repaired/protonated models, building PET-like ligand fragments, and filling `gs_pose_manifest.tsv` with accepted and rejected Michaelis-complex candidates.
