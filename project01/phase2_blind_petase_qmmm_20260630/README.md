# Phase 2 - Blind PETase QM/MM Mechanism Reproduction

Created: 2026-06-30

## Purpose

This phase starts a blind first-principles PETase mechanism study. The goal is to derive the PETase acylation and deacylation mechanisms from structure and chemistry, using QM/MM calculations and dynamical validation, then compare against the paper only at the end.

## Boundary

Use the PETase paper only as methodological inspiration. Do not use the paper's concrete results as inputs.

Not allowed before final validation:

- paper TS coordinates;
- paper reaction-coordinate formulas;
- paper selected CVs;
- paper umbrella windows;
- paper aimless-shooting trajectories;
- paper barrier values, rate constants, or rate-limiting assignment;
- paper conclusions about His motion, Asp role, Trp185 role, oxyanion-hole role, or tetrahedral TS/intermediate status.

Allowed:

- PETase structures;
- generic serine hydrolase chemistry;
- QM/MM methodology;
- TS search, committor analysis, PMF, and unbiased/validated reaction-coordinate ideas.

## Main Plan

Primary plan document:

```text
projects/01-specialized-ts-aware-scorer/docs/petase_blind_qmmm_mechanism_plan.md
```

## Execution Stages

### Stage 1 - Structure-only setup

Deliverables:

```text
blind_work/01_system_setup/structure_selection.tsv
blind_work/01_system_setup/gs_pose_manifest.tsv
```

Tasks:

1. Select WT PETase structures from public structural databases.
2. Repair missing residues/atoms and assign disulfides.
3. Assign protonation states around pH 7.
4. Build PET dimer or BHET/MHET-like ester fragments.
5. Generate multiple substrate-bound Michaelis candidates.
6. Filter poses by generic catalytic geometry.

### Stage 2 - Classical ensemble

Deliverables:

```text
blind_work/02_classical_md/productive_conformer_manifest.tsv
blind_work/02_classical_md/rejected_pose_manifest.tsv
```

Tasks:

1. Solvate and equilibrate selected GS poses.
2. Run short MD replicates.
3. Cluster active-site conformations by catalytic geometry.
4. Select representative productive conformers for QM/MM exploration.

### Stage 3 - Mechanism hypothesis tree

Deliverables:

```text
blind_work/03_mechanism_tree/mechanism_hypotheses.yaml
blind_work/03_mechanism_tree/candidate_cv_sets.tsv
```

Tasks:

1. Enumerate acylation paths.
2. Enumerate deacylation paths.
3. Enumerate His/Asp proton-relay alternatives.
4. Define only generic CV candidates, not paper-derived CVs.

### Stage 4 - Low-cost QM/MM exploration

Deliverables:

```text
blind_work/04_qmmm_exploration/path_screening_table.tsv
blind_work/04_qmmm_exploration/ts_like_guess_manifest.tsv
```

Tasks:

1. Build initial QM regions from catalytic chemistry.
2. Run DFTB3/MM, xTB/MM, or low-cost DFT/MM relaxed scans.
3. Screen for chemically productive paths.
4. Extract TS-like guesses from scan maxima or string images.

### Stage 5 - TS refinement

Deliverables:

```text
blind_work/05_ts_refinement/ts_refinement_manifest.tsv
blind_work/05_ts_refinement/acylation_ts_candidates/
blind_work/05_ts_refinement/deacylation_ts_candidates/
```

Tasks:

1. Refine TS guesses with QM/MM TS optimization or constrained saddle refinement.
2. Expand QM region for top candidates.
3. Reject chemically wrong or unstable candidates.

### Stage 6 - TS ensemble construction

Deliverables:

```text
blind_work/06_ts_ensemble/acylation_ts_ensemble.tsv
blind_work/06_ts_ensemble/deacylation_ts_ensemble.tsv
```

Tasks:

1. Repeat TS search across multiple conformer clusters.
2. Cluster TS candidates.
3. Preserve ensemble diversity.
4. Select representatives for committor testing.

### Stage 7 - Dynamical validation

Deliverables:

```text
blind_work/07_committor/committor_summary.tsv
blind_work/07_committor/accepted_ts_manifest.tsv
blind_work/07_committor/rejected_ts_manifest.tsv
```

Tasks:

1. Run short forward/backward trajectories from TS candidates.
2. Estimate committors for representative TS clusters.
3. Verify product/reactant endpoint chemistry.

### Stage 8 - Free-energy barriers

Deliverables:

```text
blind_work/08_free_energy/acylation_pmf.tsv
blind_work/08_free_energy/deacylation_pmf.tsv
blind_work/08_free_energy/barrier_summary.md
```

Tasks:

1. Build reaction coordinates from our own successful trajectories.
2. Run umbrella/string/metadynamics sampling.
3. Analyze PMFs with MBAR/WHAM.
4. Estimate acylation and deacylation barriers.

### Stage 9 - Final paper validation

Deliverables:

```text
blind_work/09_paper_validation/blind_vs_paper_comparison.md
blind_work/09_paper_validation/discrepancy_audit.md
```

Tasks:

1. Open the paper results only after TS ensembles and preliminary PMFs exist.
2. Compare mechanism, TS geometry, barriers, and rate-limiting step.
3. Audit discrepancies by structure, substrate model, protonation, QM region, sampling, and RC choices.

## First compute action

Recommended first action:

```text
Build WT PETase + PET dimer/BHET-like ester GS candidates, then run classical equilibration and active-site geometry filtering.
```

Do not start mutant ranking until WT acylation/deacylation TS ensembles are chemically validated.
