# Codex tasks — Project 01: Enzyme-specific true-TS barrier scorer

Updated: 2026-06-28

This file gives concrete implementation tasks for Codex.

## Working rule

Project 01 is the **true-TS teacher** project.

It may use true/refined TS structures as inputs and should learn computed barriers from externally generated QM/MM or DFT labels.

Do not implement heavy QM/MM, DFT, PLACER, or TS-search internals in this repository. This repository should contain orchestration code, data schemas, feature extraction utilities, model code, tests, and documentation. Large raw datasets, PLACER ensemble dumps, QM/MM trajectories, DFT outputs, and model checkpoints should stay outside version control.

## Milestone 0 — repository scaffold

Create the following structure under `projects/01-specialized-ts-aware-scorer/`:

```text
configs/
data_schema/
docs/
scripts/
src/ts_aware_scorer/
tests/
notebooks/README.md
```

Do not commit large notebooks; only add a small `notebooks/README.md` explaining where analysis notebooks should live.

## Milestone 1 — data schemas

Create JSON schema files:

```text
data_schema/enzyme_variant.schema.json
data_schema/reaction_template.schema.json
data_schema/conformer_state_pair.schema.json
data_schema/qmmm_barrier_label.schema.json
data_schema/feature_table.schema.json
data_schema/project01_teacher_export.schema.json
```

### Required concepts

The schemas must represent:

- `enzyme_id`
- `variant_id`
- `catalytic_class`
- `reaction_template_id`
- `reaction_step`
- `sequence`
- `structure_pdb`
- `conformer_id`
- `conformer_source`
- `GS_complex_pdb`
- `true_TS_complex_pdb`
- `delta_G_dagger`, nullable
- `delta_E_TS_GS`, nullable
- `delta_delta_G_dagger_vs_reference`, nullable
- `reference_variant_id`, nullable
- `label_tier`
- `protocol_id`
- `geometry_features_path`
- `electrostatic_features_path`
- `ensemble_summary_path`

Allow missing labels with `null` because not every conformer/state pair will have QM/MM or DFT labels.

### Reaction template schema must include

- catalytic class, e.g. `serine_hydrolase`
- reaction step, e.g. `acylation_TS1`
- forming bonds
- breaking or weakening bonds
- proton transfers
- required catalytic residues
- required geometry features
- expected charge-shift or `delta_q` fields, if known

### Conformer state-pair schema must include

- raw PLACER or external conformer provenance
- geometry filter status
- cluster ID
- GS raw/relaxed/QM-MM structure paths
- true TS raw/relaxed/QM-MM structure paths
- link to external QM/MM or DFT label metadata

## Milestone 2 — feature extraction module

Create a Python package:

```text
src/ts_aware_scorer/
  __init__.py
  geometry.py
  electrostatics.py
  ensemble.py
  qmmm_labels.py
  featurize.py
```

Implement placeholder-safe functions with clear docstrings.

### `geometry.py`

Functions:

```python
compute_distance(coord_a, coord_b) -> float
compute_angle(coord_a, coord_b, coord_c) -> float
compute_dihedral(coord_a, coord_b, coord_c, coord_d) -> float
compute_reaction_geometry_features(gs_structure, ts_structure, reaction_template) -> dict
passes_geometry_filter(feature_row, thresholds) -> bool
```

Expected features:

- forming bond length
- breaking or weakening bond length
- nucleophilic attack angle
- proton transfer distance
- oxyanion-hole hydrogen-bond distances
- catalytic triad geometry
- leaving-group alignment

### `electrostatics.py`

Functions:

```python
compute_point_charge_potential(probe_points, atom_coords, atom_charges, dielectric=4.0, cutoff=None) -> list[float]
compute_reaction_projected_field(delta_q, potentials) -> float
compute_delta_q(charges_ts, charges_gs) -> dict
```

This is a lightweight classical baseline. Later QM/MM, APBS, or external electrostatic outputs can replace it.

### `ensemble.py`

Functions:

```python
summarize_placer_ensemble(conformer_feature_rows) -> dict
softmin(scores, temperature=0.593) -> float
productive_fraction(conformer_rows, geometry_thresholds) -> float
cluster_representative_summary(conformer_rows, cluster_col="cluster_id") -> dict
```

Temperature default `0.593` kcal/mol approximates RT at room temperature.

### `qmmm_labels.py`

Functions:

```python
load_qmmm_energy_label(path: str) -> dict
compute_delta_barrier(energy_gs: float, energy_ts: float) -> float
compute_delta_delta_barrier(delta_barrier: float, reference_delta_barrier: float) -> float
validate_protocol_id(label_row: dict) -> None
```

This module should parse external labels and validate metadata. It must not run QM/MM or DFT itself.

### `featurize.py`

Function:

```python
build_feature_row(sample_manifest: dict) -> dict
```

This should combine geometry, electrostatics, ensemble summaries, and QM/MM label metadata into one flat feature row.

## Milestone 3 — baseline models

Create:

```text
src/ts_aware_scorer/models.py
scripts/train_baseline.py
scripts/predict_rank.py
scripts/export_teacher_labels.py
configs/baseline.yaml
```

Baseline models:

1. `PhysicalFeatureRegressor`
2. `RidgeBarrierRegressor`
3. `RandomForestBarrierRegressor`
4. `TrueTSEmbeddingRegressor`, optional
5. `PairwiseRankingModel`, optional

Model targets:

- `delta_G_dagger`
- `delta_E_TS_GS`
- `delta_delta_G_dagger_vs_reference`
- within-template pairwise ranking

Do not require deep-learning dependencies at first. Start with `scikit-learn` if available.

## Milestone 4 — teacher export for Project 02

Create an export format and script:

```text
scripts/export_teacher_labels.py
```

The exported CSV/JSONL should include columns such as:

```text
sample_id
enzyme_id
variant_id
catalytic_class
reaction_template_id
reaction_step
protein_embedding_path
GS_embedding_path
true_TS_embedding_path
project01_delta_G_pred
project01_delta_delta_G_pred
project01_field_score
project01_geometry_score
project01_ensemble_score
project01_uncertainty
qmmm_delta_G_dagger
qmmm_delta_delta_G_dagger_vs_reference
protocol_id
```

This is the main interface from Project 01 to Project 02.

## Milestone 5 — example mini dataset

Create a tiny synthetic example in:

```text
examples/synthetic_project01/
```

Include only small JSON/CSV files, not real PDB, PLACER, QM/MM, or DFT data.

Purpose:

- let tests run
- demonstrate conformer-level state-pair format
- demonstrate true-TS feature table format
- demonstrate ranking output
- demonstrate Project 01 teacher export format

## Milestone 6 — tests

Create tests for:

```text
tests/test_geometry.py
tests/test_electrostatics.py
tests/test_ensemble.py
tests/test_qmmm_labels.py
tests/test_schema_examples.py
tests/test_teacher_export.py
```

Minimum expectations:

- distance is correct for simple coordinates
- angle is correct for a 90-degree toy case
- dihedral is deterministic for a simple toy case
- reaction-projected field equals dot product of `delta_q` and potential vector
- `delta_delta_barrier = delta_barrier - reference_delta_barrier`
- softmin behaves sensibly
- synthetic schema examples validate
- teacher export contains all fields needed by Project 02

## Milestone 7 — documentation for external compute

Create:

```text
docs/qmmm_barrier_protocol.md
docs/external_compute_interface.md
```

`qmmm_barrier_protocol.md` should define:

- reaction template ID
- state type: GS / TS / intermediate / product
- true TS construction or refinement source
- PLACER screening role
- QM region and MM environment
- protonation and charge state
- method, basis, functional, embedding, dispersion, solvent model
- constrained optimization or reaction-coordinate scan protocol
- label tiers, e.g. proxy / xTB / DFT cluster / QM-MM scan / QM-MM free energy
- required metadata for every label

`external_compute_interface.md` should document how external jobs return results:

```text
charges_GS.csv
charges_TS.csv
delta_q.csv
energies.json
optimized_GS.pdb
optimized_TS.pdb
protocol.yaml
```

Make clear that heavy QM/MM, DFT, and PLACER jobs are executed outside this repo.

## Milestone 8 — first runnable command

Add a minimal command in documentation:

```bash
python projects/01-specialized-ts-aware-scorer/scripts/train_baseline.py \
  --features examples/synthetic_project01/features.csv \
  --target delta_delta_G_dagger_vs_reference
```

Expected output:

```text
trained model artifact path
cross-validation metrics
ranked variants CSV
teacher export CSV/JSONL for Project 02
```

## Notes for Codex

- Keep functions deterministic and testable.
- Prefer small, typed utility functions.
- Do not add large binary files.
- Do not implement real quantum chemistry or PLACER internals.
- When unsure, create a small synthetic example rather than relying on real enzyme data.
