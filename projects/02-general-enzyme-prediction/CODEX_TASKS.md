# Codex tasks — Project 02: General enzyme prediction model

This file gives concrete implementation tasks for Codex.

## Working rule

Project 02 should not depend on manually supplied TS complexes at inference time. It may use TS-derived labels during training or distillation from Project 01.

Do not commit large datasets, model checkpoints, downloaded protein models, or raw QM/MM outputs.

## Milestone 0 — repository scaffold

Create the following structure under `projects/02-general-enzyme-prediction/`:

```text
configs/
data_schema/
scripts/
src/general_enzyme_predictor/
tests/
notebooks/README.md
```

## Milestone 1 — data schema

Create JSON schema files:

```text
data_schema/general_enzyme_sample.schema.json
data_schema/reaction_prior.schema.json
data_schema/distilled_labels.schema.json
```

Required fields:

- `enzyme_id`
- `reaction_id`
- `substrate_id`
- `sequence`
- `GS_complex_pdb`, nullable
- `reaction_prior`, object
- `activity_label`, nullable
- `distilled_project01_labels`, nullable

The schema must allow optional TS supervision fields but must not require them.

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
```

Keep these functions lightweight. If RDKit is not installed, functions should fail gracefully with a clear error message.

## Milestone 3 — embedding input interfaces

Create:

```text
src/general_enzyme_predictor/embeddings.py
```

Functions:

```python
load_protein_embedding(path: str)
load_ensemble_embedding(path: str)
load_reaction_embedding(path: str | None)
combine_embeddings(protein_embedding, ensemble_embedding, reaction_features)
```

Do not implement ESM inference in this milestone. The first version should consume precomputed embeddings.

## Milestone 4 — model heads

Create:

```text
src/general_enzyme_predictor/models.py
```

Implement minimal model classes:

1. `EmbeddingOnlyRegressor`
2. `ReactionConditionedRegressor`
3. `MultiTaskReactionRegressor`

Outputs:

```text
barrier_proxy
activity_proxy
field_proxy, optional
geometry_proxy, optional
ensemble_proxy, optional
```

Start with simple PyTorch or scikit-learn models. Avoid unnecessary complexity.

## Milestone 5 — distillation from Project 01

Create:

```text
src/general_enzyme_predictor/distillation.py
scripts/train_from_project01_distillation.py
```

The script should train Project 02 using a CSV exported from Project 01 with columns like:

```text
enzyme_id
variant_id
reaction_id
protein_embedding_path
GS_embedding_path
reaction_prior_path
project01_delta_G_pred
project01_field_score
project01_geometry_score
project01_ensemble_score
```

The goal is to learn Project 01-style mechanism labels without requiring TS at inference.

## Milestone 6 — evaluation splits

Create:

```text
src/general_enzyme_predictor/splits.py
```

Implement split helpers:

```python
split_within_enzyme(df, enzyme_id_col="enzyme_id")
split_leave_one_enzyme_out(df, enzyme_id_col="enzyme_id")
split_leave_one_reaction_family_out(df, family_col="reaction_family")
```

## Milestone 7 — synthetic example

Create:

```text
examples/synthetic_project02/
```

Include tiny CSV/JSON files that demonstrate:

- an enzyme sequence identifier
- precomputed fake protein embedding path
- fake GS ensemble embedding path
- reaction prior JSON
- activity label
- distilled Project 01 labels

Do not include large real embeddings.

## Milestone 8 — tests

Create tests for:

```text
tests/test_reaction_prior.py
tests/test_embeddings.py
tests/test_splits.py
tests/test_models_smoke.py
```

Tests should run without external databases or heavy dependencies.

## Milestone 9 — first runnable command

Add documentation for a minimal command:

```bash
python projects/02-general-enzyme-prediction/scripts/train_from_project01_distillation.py \
  --features examples/synthetic_project02/distillation_features.csv \
  --config projects/02-general-enzyme-prediction/configs/baseline.yaml
```

Expected output:

```text
trained model artifact path
validation metrics
predicted ranking CSV
```

## Milestone 10 — comparison baselines

Implement a comparison table script:

```text
scripts/compare_baselines.py
```

Compare:

1. embedding-only model
2. reaction-prior-only model
3. embedding + reaction prior
4. embedding + reaction prior + distilled Project 01 labels

Metrics:

- Spearman correlation
- pairwise ranking accuracy
- top-K enrichment

## Notes for Codex

- Keep Project 02 separate from Project 01.
- Project 02 may import exported CSVs from Project 01, but should not import Project 01 internals unless a shared utility package is created later.
- The first implementation should be a clean scaffold with tests, not a heavy end-to-end model.
