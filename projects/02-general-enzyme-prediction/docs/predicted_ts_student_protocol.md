# Predicted-TS student protocol — Project 02

Created: 2026-06-28

Project 02 trains a predicted-TS student model that approximates Project 01 true-TS teacher behavior without using true TS coordinates at inference time.

## 1. Core pipeline

```text
protein + GS + reaction prior
→ TS prediction model
→ predicted TS geometry / embedding / Δq prior
→ PLACER screening
→ predicted-TS-compatible ensemble features
→ predicted-TS student model
→ ΔG‡ proxy / ΔΔG‡ proxy / catalytic potential ranking
```

## 2. Inputs

Required inference-time inputs:

```yaml
enzyme_id:
sequence_or_structure:
GS_complex:
reaction_prior:
predicted_ts_prior:
placer_screening_summary:
```

Optional training-time inputs:

```yaml
project01_teacher_outputs:
true_TS_embedding:
qmmm_delta_G_dagger:
qmmm_delta_delta_G_dagger_vs_reference:
```

True TS fields are training-only and must not be required at inference.

## 3. TS prediction model output

Project 02 does not implement the TS prediction model. It consumes its output through a manifest:

```yaml
ts_prediction_model_id:
input_GS_structure:
input_reaction_template_id:
predicted_ts_structure:
predicted_ts_embedding:
predicted_delta_q:
confidence_score:
uncertainty_score:
generation_metadata:
```

The predicted TS output can be a structure, embedding, predicted charge-shift pattern, or any combination of these.

## 4. PLACER screening

The predicted TS prior should be passed to PLACER or a PLACER-compatible external screening workflow:

```text
protein + GS + predicted TS prior
→ PLACER-screened predicted-TS-compatible ensemble
```

The screening summary should include:

```yaml
placer_run_id:
input_predicted_ts_prior:
num_conformers_generated:
num_conformers_passing_geometry_filter:
productive_predicted_ts_fraction:
cluster_summaries:
ensemble_embedding_path:
geometry_feature_summary:
```

Raw PLACER ensembles should stay outside version control.

## 5. Teacher-student distillation

Project 01 exports teacher labels:

```text
project01_delta_G_pred
project01_delta_delta_G_pred
project01_field_score
project01_geometry_score
project01_ensemble_score
project01_uncertainty
qmmm_delta_G_dagger
qmmm_delta_delta_G_dagger_vs_reference
true_TS_embedding_path, optional training-only field
```

Project 02 learns from these labels using:

```text
barrier regression
relative barrier regression
pairwise ranking
teacher mechanism-component distillation
predicted-TS-to-true-TS embedding alignment, training only
```

## 6. Recommended baselines

Every Project 02 experiment should compare:

1. GS-only baseline.
2. reaction-prior-only baseline.
3. embedding-only baseline.
4. predicted-TS student.
5. Project 01 true-TS teacher upper bound.
6. wrong-TS control.
7. shuffled-reaction-template control.

The central diagnostic is:

```text
Does predicted TS + PLACER screening close the gap between GS-only and true-TS teacher?
```

## 7. Evaluation splits

Main splits:

```text
within-enzyme mutant split
leave-one-enzyme-out
leave-one-family-out
leave-one-reaction-template-out
within-catalytic-class transfer
```

Avoid conformer leakage. Conformers from the same variant/design should not be split across train and test in the main evaluation.

## 8. Output report

Each training run should produce:

```text
student predictions CSV
baseline comparison table
student-vs-teacher gap report
wrong-TS control report
shuffled-reaction-template control report
```
