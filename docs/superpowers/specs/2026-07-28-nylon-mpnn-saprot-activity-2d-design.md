# Nylonase ProteinMPNN-SaProt 2D activity overlay design

Date: 2026-07-28  
Branch: `codex/nylon-pose-mpnn-saprot-activity`

## Objective

Create three two-dimensional nylonase compatibility maps using ProteinMPNN and SaProt, with separate experimental overlays for PA66-L1, PA66-L2, and PA6. The figures are descriptive compatibility-versus-endpoint views; they do not treat either score as activity, free energy, melting temperature, or a catalytic barrier.

## Authoritative universe and join

All joins use exact `sequence_md5`.

- Nylonase authority: 4,167 rows in the scope-fixed paired manifest.
- SaProt exact-canonical scope: `scoring_sequence_scope == CANONICAL_EXACT`.
- ProteinMPNN primary model: vanilla `v_48_020`.
- Background inclusion requires `proteinmpnn_status == PASS`, finite `vanilla_mean_log_likelihood`, finite `saprot_mean_log_likelihood`, and exact-canonical scoring.
- Live preflight on 2026-07-28 found 3,985 unique exact-MD5 background rows.
- Experimental authority is Chem Catalysis 2025 Nyl01-Nyl95, excluding NylC. Live preflight found 95 unique experimental sequences, 80 linked to the common background, and 15 `NOT_EVALUATED_COMPATIBILITY_MISSING_EXACT_SCORE`.
- Missing scores are never imputed or converted to zero.

## Axes and endpoints

Every figure uses the same axes and background:

- x: ProteinMPNN vanilla v_48_020 mean log-likelihood.
- y: SaProt full-protein mean log-likelihood.
- background: all 3,985 common exact-canonical nylonase sequences, uniform light gray.

Three endpoint-specific overlays are produced:

1. `pa66_l1_color_span_px`
2. `pa66_l2_color_span_px`
3. `pa6_max_color_span_px`

Each overlay uses only the 80 exact-MD5-linked experimental rows. Color encodes endpoint-specific percentile rank, where larger proxy values map to higher percentile. Hover text reports Nyl ID, raw raster color-span proxy, descending absolute rank, percentile rank, both compatibility scores, and sequence MD5.

The endpoint values are raster color-span proxies from the curated figure, not absolute enzyme activity or kinetic constants. Rankings are computed independently within each endpoint and must not be compared as a shared numerical activity scale across PA66 and PA6.

## Outputs

The versioned Sugon output directory contains:

- three standalone interactive HTML files;
- three standalone PNG and PDF files;
- one three-panel PNG and PDF summary;
- `nylon_mpnn_saprot_background.tsv`;
- `nyl95_experimental_overlay.tsv`;
- `nyl95_not_evaluated.tsv`;
- `endpoint_summary.tsv`;
- `input_sha256.tsv`;
- `PASS.json`, `RUNBOOK.md`, and append-only `logs/run_history.tsv`;
- reusable scripts and tests.

## Audit gates

A run passes only if:

- authority and input SHA256 values are recorded;
- background is exactly one row per scoring MD5;
- all background scores are finite and exact-canonical;
- experimental authority is 95 unique MD5 values;
- linked plus NOT_EVALUATED equals 95;
- all three linked endpoint values are finite;
- every linked row appears once in every endpoint figure;
- axes are identical across figures;
- gray background is not colored by length or another property;
- HTML, PNG, PDF, tables, and audit JSON all exist and are non-empty.

## Interpretation boundary

Higher (less negative) ProteinMPNN or SaProt log-likelihood means greater compatibility under that model. Neither model directly predicts PA66/PA6 activity. The plots expose associations and candidate regions only; scientific claims require effect sizes, uncertainty, sequence-cluster sensitivity, and experimental validation.
