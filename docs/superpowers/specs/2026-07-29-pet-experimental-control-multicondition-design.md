# PET experimental-control multicondition patch analysis design

## Objective

Test whether the previously observed 14 Å pattern—nominally higher relative external stickiness but lower catalytic-region stickiness in higher-activity PET enzymes—recurs across other Nature 2022 assay conditions.

This is a post hoc association analysis. It does not recompute structures, probe fields, shells, or patches, and it does not interpret endpoint product amounts as kinetics, barriers, or causal mechanisms.

## Evaluated universe

- Patch features: the 30 exact-canonical Nature 2022 controls with valid PET `PATCH_PASS`.
- Excluded: protein 202, retained as `NONCANONICAL_MAPPING_REQUIRED`.
- Activity authority: `activity_energy_long_authority.tsv`, restricted to `dataset=NATURAL2022`, `activity_metric=sum aromatic products`, and `activity_unit=mg/L`.
- PET and nylon denominators remain separate.

## Analysis tiers

### Primary replication: Table D3

Analyze every Table D3 temperature-buffer condition separately for the exact-canonical patch panel.

Primary 14 Å metrics:

1. `r14_relative_external_sticky_fraction`
2. `r14_composite_catalytic_candidate_percentile`
3. `r14_mean_catalytic_candidate_percentile_ACOA`

For each condition report the exact joined denominator, zero-product count, Spearman rho, permutation P value, bootstrap 95% interval, and Benjamini-Hochberg q value across the full Table D3 condition-by-metric family.

### Cross-condition robustness

Within each Table D3 condition, convert activity to tied ranks/percentiles. Average those condition-normalized ranks per protein, recording the number of available conditions. Correlate the cross-condition mean rank with the same three 14 Å metrics.

Also summarize direction consistency across Table D3 conditions:

- external metric: number and fraction of positive rhos;
- catalytic metrics: number and fraction of negative rhos;
- median and range of condition-wise rho.

These summaries are descriptive because assay endpoints are correlated repeated measurements.

### Radius sensitivity

Repeat the cross-condition aggregate analysis for 6 and 10 Å. Radius results are sensitivity analyses, not independent confirmations.

### Table D6 sensitivity

Analyze amorphous film, amorphous powder, and crystalline powder separately using within-condition rank normalization followed by per-protein aggregation. Report exact denominators. Because each raw condition contains only 11–15 records, label all D6 findings low-power exploratory.

### Exclusion

Table D4 is excluded from correlation analysis because each condition contains only 3–4 records.

## Statistical safeguards

- Never pool raw mg/L values across experimental conditions.
- Minimum per-condition joined denominator: 10.
- Use deterministic random seeds.
- Use permutation Spearman tests and bootstrap confidence intervals.
- Apply Benjamini-Hochberg correction to declared test families.
- Keep nominal P values and adjusted q values distinct.
- Missing activity rows are `NOT_EVALUATED_FOR_CONDITION`, not zero and not inactivity.

## Outputs and provenance

Create one versioned analysis directory containing:

- joined condition-level records;
- per-condition associations;
- cross-condition rank aggregates;
- aggregate associations and radius sensitivity;
- D6 substrate-form sensitivity;
- `ANALYSIS_PASS.json` with input/output hashes and denominator audits.

Add one re-runnable analysis script, update the remote RUNBOOK, and append one run-history record. Do not submit Slurm work because this is a small tabular analysis using existing artifacts.
