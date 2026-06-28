# Codex tasks — Project 02: Predicted-TS catalytic-class barrier predictor

Updated: 2026-06-28

This file gives concrete implementation tasks for Codex.

## Working rule

Project 02 is the **predicted-TS student** project.

It must not require manually supplied or true/refined TS complexes at inference time. It may use true-TS-derived labels, true-TS embeddings, and teacher outputs from Project 01 during training or distillation.

Do not commit large datasets, model checkpoints, downloaded protein models, raw PLACER ensembles, predicted TS bulk outputs, or raw QM/MM outputs.

## Milestone 0 — repository scaffold

Create the following structure under `projects/02-general-enzyme-prediction/`:

```text
configs/
data_schema/
docs/
scripts/
src/general_enzyme_predictor/
tests/
notebooks/README.md
```

## Milestone 1 — data schemas

Create JSON schema files:

```text
data_schema/general_enzyme_sample.schema.json
data_schema/reaction_prior.schema.json
data_schema/predicted_ts_prior.schema.json
data_schema/placer_screened_predicted_ts.schema.json
data_schema/project01_distilled_labels.schema.json
data_schema/evaluation_baseline_result.schema.json
```

Required fields:

- `sample_id`
- `enzyme_id`
- `variant_id`, nullable
- `catalytic_class`
- `reaction_template_id`
- `reaction_step`
- `substrate_id`
- `sequence`
- `GS_complex_pdb`, nullable
- `reaction_prior`, object
- `predicted_ts_prior`, nullable object
- `placer_screened_predicted_ts`, nullable object
- `distilled_project01_labels`, nullable object
- `activity_label`, nullable and optional only

The schema must allow optional true-TS supervision fields during training, but true TS fields must never be required for inference samples.

### `reaction_prior.schema.json` must include

- reaction template ID
- catalytic class
- reaction step
- atom mapping
- forming bonds
- breaking or weakening bonds
- proton transfers
- required catalytic residues
- expected geometry constraints

### `predicted_ts_prior.schema.json` must include

- `ts_prediction_model_id`
- predicted TS structure path, nullable
- predicted TS embedding path, nullable
- predicted `delta_q` path, nullable
- confidence or uncertainty score, nullable
- generation metadata

### `project01_distilled_labels.schema.json` must include

- Project 01 teacher barrier prediction
- QM/MM or DFT barrier label, if available
- `delta_delta_G_dagger_vs_reference`, if available
- Project 01 field / geometry / ensemble components
- true-TS embedding path, optional training-only field
- protocol ID

## Milestone 2 — reaction prior utilities

Create:

```text
src/general_enzyme_predictor/reaction_prior.py
```

Implement utilities:

```python
parse_reaction_smiles(reaction_smiles: str) -> dict
extract_bond_changes(atom_mapping: dict) -> dict
build_reaction_prior_features(prior: dict) -> dict
validate_reaction_template(prior: dict) -> None
```

Keep these functions lightweight. If RDKit is not installed, functions should fail gracefully with a clear error message.

## Milestone 3 — predicted TS interfaces

Create:

```text
src/general_enzyme_predictor/predicted_ts.py
```

Implement interface functions:

```python
load_predicted_ts_prior(path: str) -> dict
load_predicted_ts_embedding(path: str)
validate_predicted_ts_prior(prior: dict) -> None
summarize_predicted_ts_metadata(prior: dict) -> dict
```

This module should consume outputs from an external TS prediction model. It must not implement the TS prediction model itself.

## Milestone 4 — PLACER screening interfaces

Create:

```text
src/general_enzyme_predictor/placer_screening.py
```

Implement lightweight functions:

```python
load_placer_screening_summary(path: str) -> dict
summarize_predicted_ts_compatible_ensemble(rows) -> dict
productive_predicted_ts_fraction(rows, thresholds) -> float
```

This module should parse PLACER screening summaries and conformer features. It must not run PLACER itself.

## Milestone 5 — embedding input interfaces

Create:

```text
src/general_enzyme_predictor/embeddings.py
```

Functions:

```python
load_protein_embedding(path: str)
load_gs_embedding(path: str)
load_predicted_ts_embedding(path: str | None)
load_ensemble_embedding(path: str | None)
combine_embeddings(protein_embedding, gs_embedding, predicted_ts_embedding, ensemble_embedding, reaction_features)
```

Do not implement ESM inference or TS-model inference in this milestone. The first version should consume precomputed embeddings.

## Milestone 6 — model heads

Create:

```text
src/general_enzyme_predictor/models.py
```

Implement minimal model classes:

1. `GSOnlyRegressor`
2. `ReactionPriorRegressor`
3. `EmbeddingOnlyRegressor`
4. `PredictedTSStudentRegressor`
5. `MultiTaskPredictedTSRegressor`

Outputs:

```text
barrier_proxy
barrier_delta_proxy
catalytic_potential_score
field_proxy, optional
geometry_proxy, optional
ensemble_proxy, optional
uncertainty, optional
```

Start with simple PyTorch or scikit-learn models. Avoid unnecessary complexity.

## Milestone 7 — distillation from Project 01

Create:

```text
src/general_enzyme_predictor/distillation.py
scripts/train_from_project01_distillation.py
```

The script should train Project 02 using a CSV/JSONL exported from Project 01 with columns like:

```text
sample_id
enzyme_id
variant_id
catalytic_class
reaction_template_id
reaction_step
protein_embedding_path
GS_embedding_path
reaction_prior_path
predicted_ts_embedding_path
predicted_ts_prior_path
placer_screening_summary_path
project01_delta_G_pred
project01_delta_delta_G_pred
project01_field_score
project01_geometry_score
project01_ensemble_score
project01_uncertainty
qmmm_delta_G_dagger
qmmm_delta_delta_G_dagger_vs_reference
true_TS_embedding_path   # training-only, optional
protocol_id
```

The goal is to learn Project 01-style mechanism labels without requiring true TS at inference.

Distillation functions:

```python
load_project01_teacher_export(path: str)
build_student_training_table(export_df, predicted_ts_manifest_df)
compute_teacher_student_targets(row) -> dict
```

## Milestone 8 — losses

Create:

```text
src/general_enzyme_predictor/losses.py
```

Implement lightweight loss helpers:

```python
barrier_regression_loss(y_pred, y_true)
pairwise_ranking_loss(score_a, score_b, label, margin=1.0)
teacher_distillation_loss(student_outputs, teacher_outputs, weights=None)
ts_embedding_alignment_loss(z_pred, z_true)
```

The true TS embedding alignment loss is training-only. Inference must not depend on true TS embedding.

## Milestone 9 — evaluation splits

Create:

```text
src/general_enzyme_predictor/splits.py
```

Implement split helpers:

```python
split_within_enzyme(df, enzyme_id_col="enzyme_id")
split_leave_one_enzyme_out(df, enzyme_id_col="enzyme_id")
split_leave_one_family_out(df, family_col="family_id")
split_leave_one_catalytic_class_out(df, class_col="catalytic_class")
split_leave_one_reaction_template_out(df, template_col="reaction_template_id")
split_by_variant_not_conformer(df, variant_col="variant_id")
```

Avoid conformer leakage: conformers from the same variant/design must not be split across train and test in the main evaluation.

## Milestone 10 — required baselines and controls

Create:

```text
scripts/compare_baselines.py
```

Compare:

1. GS-only model
2. reaction-prior-only model
3. embedding-only model
4. predicted-TS student model
5. Project 01 true-TS teacher upper bound, imported as an external result table
6. wrong-TS control
7. shuffled-reaction-template control

Metrics:

- Spearman correlation with computed `ΔG‡` or `ΔΔG‡`
- pairwise ranking accuracy
- top-K enrichment
- leave-one-enzyme-out performance
- leave-one-reaction-template-out performance
- gap between predicted-TS student and true-TS teacher

The main diagnostic is whether predicted TS + PLACER screening closes the gap between GS-only and true-TS teacher.

## Milestone 11 — synthetic example

Create:

```text
examples/synthetic_project02/
```

Include tiny CSV/JSON files that demonstrate:

- an enzyme sequence identifier
- precomputed fake protein embedding path
- fake GS embedding path
- fake predicted TS embedding path
- reaction prior JSON
- PLACER screening summary JSON
- distilled Project 01 labels
- optional true TS embedding path for training-only alignment

Do not include large real embeddings, real PLACER outputs, or real TS prediction outputs.

## Milestone 12 — tests

Create tests for:

```text
tests/test_reaction_prior.py
tests/test_predicted_ts.py
tests/test_placer_screening.py
tests/test_embeddings.py
tests/test_distillation.py
tests/test_losses.py
tests/test_splits.py
tests/test_models_smoke.py
tests/test_baseline_controls.py
```

Tests should run without external databases, PLACER, TS prediction models, QM/MM, or heavy dependencies.

## Milestone 13 — documentation

Create:

```text
docs/predicted_ts_student_protocol.md
docs/inference_no_true_ts_boundary.md
```

`predicted_ts_student_protocol.md` should describe:

- expected TS prediction model outputs
- how predicted TS priors are passed to PLACER screening
- how Project 01 teacher labels are consumed
- how true-TS embeddings may be used during training only
- required baselines and controls

`inference_no_true_ts_boundary.md` should explicitly define what is forbidden at inference:

- no true TS coordinates
- no true-TS-derived embedding computed from target true TS
- no QM/MM barrier label for the target sample
- no leakage from Project 01 teacher outputs for test samples

## Milestone 14 — first runnable command

Add documentation for a minimal command:

```bash
python projects/02-general-enzyme-prediction/scripts/train_from_project01_distillation.py \
  --features examples/synthetic_project02/distillation_features.csv \
  --config projects/02-general-enzyme-prediction/configs/baseline.yaml
```

Expected output:

```text
trained student model artifact path
validation metrics
predicted ranking CSV
baseline comparison table
student-vs-teacher gap report
```

## Notes for Codex

- Keep Project 02 separate from Project 01.
- Project 02 may import exported CSV/JSONL files from Project 01, but should not import Project 01 internals unless a shared utility package is created later.
- Do not implement the actual TS prediction model, PLACER, QM/MM, or DFT engines.
- The first implementation should be a clean scaffold with tests, not a heavy end-to-end model.
