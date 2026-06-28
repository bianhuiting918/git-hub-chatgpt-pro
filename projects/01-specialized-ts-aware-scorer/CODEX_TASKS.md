# Codex tasks — Project 01: Specialized TS-aware enzyme scorer

This file gives concrete implementation tasks for Codex.

## Working rule

Do not implement heavy QM/MM, DFT, or PLACER internals in this repository. This repository should contain orchestration code, data schemas, feature extraction utilities, model code, and documentation. Large raw datasets and model checkpoints should stay outside version control.

## Milestone 0 — repository scaffold

Create the following structure under `projects/01-specialized-ts-aware-scorer/`:

```text
configs/
data_schema/
scripts/
src/ts_aware_scorer/
tests/
notebooks/README.md
```

Do not commit large notebooks; only add a small `notebooks/README.md` explaining where analysis notebooks should live.

## Milestone 1 — data schema

Create JSON schema files:

```text
data_schema/enzyme_variant.schema.json
data_schema/reaction_state_pair.schema.json
data_schema/feature_table.schema.json
```

Required fields:

- `enzyme_id`
- `variant_id`
- `reaction_id`
- `sequence`
- `structure_pdb`
- `GS_complex_pdb`
- `TS_complex_pdb`
- `delta_G_dagger`
- `delta_q_path`
- `geometry_features_path`
- `ensemble_summary_path`

Allow missing labels with `null` because not every variant will have DFT/QM/MM labels.

## Milestone 2 — feature extraction module

Create a Python package:

```text
src/ts_aware_scorer/
  __init__.py
  geometry.py
  electrostatics.py
  ensemble.py
  featurize.py
```

Implement placeholder-safe functions with clear docstrings:

### `geometry.py`

Functions:

```python
compute_distance(coord_a, coord_b) -> float
compute_angle(coord_a, coord_b, coord_c) -> float
compute_reaction_geometry_features(gs_structure, ts_structure, reaction_atoms) -> dict
```

Expected features:

- forming bond length
- breaking bond length
- nucleophilic attack angle
- proton transfer distance
- oxyanion-hole hydrogen-bond distances

### `electrostatics.py`

Functions:

```python
compute_point_charge_potential(probe_points, atom_coords, atom_charges, dielectric=4.0, cutoff=None) -> list[float]
compute_reaction_projected_field(delta_q, potentials) -> float
```

This is a lightweight classical baseline. Later QM/MM or APBS outputs can replace it.

### `ensemble.py`

Functions:

```python
summarize_placer_ensemble(conformer_feature_rows) -> dict
softmin(scores, temperature=0.593) -> float
productive_fraction(conformer_rows, geometry_thresholds) -> float
```

Temperature default 0.593 kcal/mol approximates RT at room temperature.

### `featurize.py`

Function:

```python
build_feature_row(sample_manifest: dict) -> dict
```

This should combine geometry, electrostatics, and ensemble summaries into one flat feature row.

## Milestone 3 — baseline models

Create:

```text
src/ts_aware_scorer/models.py
scripts/train_baseline.py
scripts/predict_rank.py
configs/baseline.yaml
```

Baseline models:

1. linear regression
2. ridge regression
3. small MLP, optional
4. pairwise ranking model, optional

Model targets:

- `delta_G_dagger`
- or `delta_delta_G_dagger_vs_WT`

Do not require deep-learning dependencies at first. Start with `scikit-learn` if available.

## Milestone 4 — example mini dataset

Create a tiny synthetic example in:

```text
examples/synthetic_project01/
```

Include only small JSON/CSV files, not real PDB or DFT data.

Purpose:

- let tests run
- demonstrate feature table format
- demonstrate ranking output

## Milestone 5 — tests

Create tests for:

```text
tests/test_geometry.py
tests/test_electrostatics.py
tests/test_ensemble.py
tests/test_schema_examples.py
```

Minimum expectations:

- distance is correct for simple coordinates
- angle is correct for a 90-degree toy case
- reaction-projected field equals dot product of `delta_q` and potential vector
- softmin behaves sensibly

## Milestone 6 — documentation for external compute

Create:

```text
docs/external_compute_interface.md
```

Document how external DFT/QM/MM jobs should return results:

```text
charges_GS.csv
charges_TS.csv
delta_q.csv
energies.json
optimized_GS.pdb
optimized_TS.pdb
```

Make clear that heavy DFT/QM/MM jobs are executed outside this repo.

## Milestone 7 — first runnable command

Add a minimal command in documentation:

```bash
python projects/01-specialized-ts-aware-scorer/scripts/train_baseline.py \
  --features examples/synthetic_project01/features.csv \
  --target delta_G_dagger
```

Expected output:

```text
trained model artifact path
cross-validation metrics
ranked variants CSV
```

## Notes for Codex

- Keep functions deterministic and testable.
- Prefer small, typed utility functions.
- Do not add large binary files.
- When unsure, create a small synthetic example rather than relying on real enzyme data.
