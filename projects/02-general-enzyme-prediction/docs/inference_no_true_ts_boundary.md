# No-true-TS inference boundary — Project 02

Created: 2026-06-28

Project 02 is designed for practical inference when true TS coordinates are unavailable.

This document defines what is allowed and forbidden at inference time.

## 1. Allowed inference inputs

Project 02 may use:

```text
enzyme sequence or structure
GS complex
substrate identity
reaction template
atom mapping
product structure
bond-breaking / bond-forming pattern
proton-transfer pattern
predicted TS geometry
predicted TS embedding
predicted Δq pattern
PLACER-screened predicted-TS-compatible ensemble
protein embeddings
GS embeddings
reaction-prior features
```

## 2. Forbidden inference inputs

Project 02 must not use:

```text
true/refined TS coordinates for the target sample
true-TS embedding computed from the target sample's true TS
QM/MM or DFT barrier label for the target sample
Project 01 teacher output computed with the target sample's true TS
manually curated target TS geometry when the benchmark claims no-true-TS inference
```

These fields may be available during training or for teacher-student analysis, but not for inference on held-out test samples.

## 3. Training-only fields

The following are allowed during training only:

```text
true_TS_complex_pdb
true_TS_embedding_path
qmmm_delta_G_dagger
qmmm_delta_delta_G_dagger_vs_reference
project01_teacher_delta_G_pred
project01_teacher_component_scores
```

They should be excluded from the inference feature builder and guarded by tests.

## 4. Evaluation rule

For a held-out test sample, the predicted-TS student should see only:

```text
protein + GS + reaction prior + predicted TS prior + PLACER screening summary
```

The true-TS teacher can be evaluated separately as an upper bound, but its true-TS inputs must not leak into the student prediction.

## 5. Required leakage controls

Project 02 tests should verify:

- true TS fields are absent from inference manifests.
- true TS embedding is not loaded by the inference feature builder.
- teacher labels are not included in student input features for test samples.
- conformers from the same variant/design are not split across train and test in main evaluations.
- wrong-TS and shuffled-reaction-template controls reduce performance if the model uses reaction-state information.

## 6. Recommended report language

Use:

```text
computed barrier proxy prediction
catalytic potential ranking
predicted-TS student inference
```

Avoid claiming:

```text
direct kcat prediction
universal enzyme activity prediction
true-TS-free quantum accuracy
```

Experimental activity can be reported later as optional external validation, not as the first-stage supervised objective.
