# Nylonase pose-retention, ProteinMPNN, SaProt, and activity correlation design

Date: 2026-07-27  
Branch: `codex/nylon-pose-mpnn-saprot-activity`

## Objective

Quantify whether PA66 AD-dimer pose retention, ProteinMPNN sequence-backbone compatibility, and SaProt sequence-structure compatibility associate with the Chem Catalysis 2025 Nyl01-Nyl95 product readouts. PA66-L1, PA66-L2, PA6, and PA66 L1/L2 product preference are separate endpoints. The analysis is descriptive and hypothesis-generating; it does not infer catalytic barriers or absolute enzyme activity.

## Authoritative inputs

All joins use exact `sequence_md5`.

1. Pose-retention tables from the frozen PA66 AD-dimer Smina analysis:
   - `enzyme_retention_summary.tsv`
   - `enzyme_retention_by_stratum.tsv`
   - Missing or incomplete grids remain `NOT_EVALUATED`.
2. ProteinMPNN scores:
   - vanilla `v_48_020` backbone-only score is primary;
   - soluble `v_48_020` is a sensitivity analysis.
3. SaProt:
   - full-protein mean log-likelihood over the exact canonical sequence is primary.
4. Chem Catalysis 2025 Nyl95 figure curation:
   - PA66 L1 and L2 color-span pixels;
   - PA6 dimer/trimer/tetramer color-span pixels and their maximum;
   - categorical screen labels only as a secondary analysis.

Raster color spans are image-derived proxies, not calibrated concentrations, rates, or absolute activity. PA66 and PA6 pixel spans are not combined numerically across figures.

## Evaluated universe

The activity authority contains Nyl01-Nyl95 plus the NylC control. Every output declares:

- activity-authority rows;
- pose-evaluated rows;
- ProteinMPNN-evaluated rows;
- SaProt-evaluated rows;
- exact four-way intersection;
- exclusions and reason-specific `NOT_EVALUATED` counts.

The enzyme is the statistical unit. Duplicate structures or repeated docking poses never become independent biological observations. The known pose table currently contains 53 evaluated enzymes; the exact final denominator is determined only after the four-way MD5 join.

## Variables

### Predictors

- Primary pose predictor: strict-NAC pose retention rate from the frozen enzyme-level summary.
- Pose sensitivities:
  - loose-gate retention rate;
  - weight=5 retention;
  - weight=10 retention.
- Primary ProteinMPNN predictor: vanilla mean log-likelihood.
- ProteinMPNN sensitivity: soluble-model mean log-likelihood.
- Primary SaProt predictor: full-protein mean log-likelihood.

No model score is called stability, activity, or activation energy.

### Experimental endpoints

1. `pa66_l1_color_span_px`
2. `pa66_l2_color_span_px`
3. `pa6_max_color_span_px`
4. PA66 product preference:
   `l1_fraction = L1 / (L1 + L2)`, evaluated only when `L1 + L2 > 0`.

The L1 fraction is conditional on a detected PA66 product signal and is not assigned to screen-negative rows.

Secondary endpoints use the curated categorical calls for L1, L2, any PA66 product, and PA6 where the authority supports them. Borderline PA6 raster calls remain explicit and are not silently forced into a binary class.

## Analyses

For each continuous endpoint:

- exact complete-case sample size;
- univariate Pearson and Spearman correlations for each predictor;
- raw and Holm-adjusted p-values within the declared endpoint family;
- bootstrap 95% confidence intervals;
- pairwise predictor correlation matrix;
- standardized multivariable linear model using pose, ProteinMPNN, and SaProt together;
- variance-inflation diagnostics, with multivariable coefficients withheld if the design is numerically unstable.

Because the Nyl95 panel contains homologous sequences, a sequence-cluster sensitivity analysis is required. Cluster-aware bootstrap or one-representative-per-cluster resampling will be reported separately from the enzyme-level primary analysis.

Categorical endpoints are sensitivity analyses using effect sizes and uncertainty, not replacements for the continuous proxy analysis.

## Figures

Produce four separate primary figures so unlike experimental scales are never conflated:

1. PA66-L1 proxy;
2. PA66-L2 proxy;
3. PA6 maximum proxy;
4. PA66 L1 fraction.

Each endpoint receives:

- a 3D scatter with axes pose retention, ProteinMPNN, and SaProt, colored by endpoint;
- three 2D predictor-versus-endpoint panels with fitted trend and uncertainty;
- a compact correlation/effect-size table;
- an interactive HTML version where practical.

Point labels identify enzymes on hover only; static plots label only predeclared extremes or influential observations to avoid cherry-picking.

## Audit and failure handling

- Failed joins, missing scores, incomplete docking grids, parse failures, and unsupported PA6 calls are `NOT_EVALUATED_<REASON>`.
- `NOT_EVALUATED` is never encoded as zero.
- A technically generated plot is not a scientific association.
- Correlation does not establish causality, catalytic competence, polymer binding, or a free-energy barrier.
- The full analysis must be reproducible from immutable input paths and SHA-256 hashes.
- Reruns use a new versioned output directory and never overwrite frozen inputs.

## Deliverables

GitHub stores only scripts, tests, this design, the implementation plan, RUNBOOK, compact TSV/JSON summaries, and small static figures. Large structure sets, trajectories, topologies, caches, model weights, credentials, and secrets are excluded.

The remote analysis directory stores:

- executable analysis and plotting scripts;
- `RUNBOOK.md`;
- `run_history.tsv` or JSONL;
- input manifest with paths and SHA-256 hashes;
- joined denominator table;
- `NOT_EVALUATED` table;
- correlation and model summaries;
- static and interactive figures.

## Verification gates

1. Exact MD5 joins and uniqueness tests pass.
2. Endpoint formulas and missing-value rules pass unit tests.
3. Reported denominators reconcile with the input manifest.
4. No `NOT_EVALUATED` value is converted to zero.
5. PA66 and PA6 scales remain separate.
6. Recomputed correlations match stored summaries within numerical tolerance.
7. Figures contain the same complete-case rows declared in their companion audit JSON.
