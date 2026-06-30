# Project 01 Route B De Novo Scaffold Status - 2026-06-30

## Scope

This file tracks the true active-site constrained de novo/new-backbone route. It intentionally separates Route B from the already completed existing/reference-scaffold sequence panels.

## Current State

| Enzyme system | Route B status | Evaluated de novo backbones | Current decision |
| --- | --- | ---: | --- |
| serine hydrolase | NOT_STARTED | 0 | Ready to start with active-site motif extraction and motif-scaffolding smoke test. |
| KSI | DEFERRED_TS_GEOMETRY_NOT_DEFINED | 0 | Do not start until the steroid/enolate TS-like ligand, protonation, and motif constraints are explicit. |

## Already Completed But Not Route B

| Enzyme system | Completed current-stage panel | Why it is not Route B completion |
| --- | --- | --- |
| serine hydrolase | 10 accepted postseq entrance-pass sequences per 90/80/70/60/50 bin on the project-local reference scaffold | These are sequence variants/predictions around an existing reference, not newly generated backbones around a fixed motif. |
| KSI | 10 accepted postseq entrance-pass sequences per 90/80/70/60/50 bin on natural scaffold 6UBQ | This is the natural-scaffold route, not motif-scaffold generation. |

## Serine-Hydrolase Route B First Target

Reference/motif source:

```text
/data/bht/project01_phase1_reset_gpu/manifests/denovo_SER_hydrolase_full_input.pdb
```

Planned remote working root:

```text
/data/bht/project01_phase1_reset_gpu/denovo_scaffold_routeB/serine_hydrolase
```

First smoke target:

```text
Generate a small batch of active-site constrained new backbones around the fixed catalytic/reactive motif, then sequence, predict, and postseq-gate them. Do not expand to production bins until the smoke batch proves that motif extraction, scaffolding, sequencing, prediction, and gate accounting all work end to end.
```

## Acceptance Accounting

Evaluated universe for Route B is currently empty:

```text
de_novo_backbones_generated = 0
de_novo_backbones_sequenced = 0
postseq_structures_evaluated = 0
postseq_entrance_pass = 0
not_evaluated = 0
```

A candidate can only be counted as Route B accepted after all of these are true:

```text
1. It comes from a new-backbone/motif-scaffolding run, not from mutation of the existing reference scaffold.
2. The fixed catalytic/reactive motif is present and mapped.
3. A sequence is designed for that generated backbone.
4. The full protein structure is predicted or validated.
5. The postseq entrance gate passes under explicit thresholds.
```

Current gate thresholds to reuse unless deliberately changed:

```text
global_backbone_rmsd_max_A = 2.50
fixed_backbone_rmsd_max_A = 1.00
catalytic_heavy_rmsd_max_A = 0.75
protein_key_distance_abs_delta_max_A = 0.75
```

For true new backbones, these thresholds may need motif-local rather than whole-reference interpretation. Any change must be recorded as a project-policy threshold, not silently mixed with the existing/reference-scaffold gate.

## Immediate Remote Checks

When GPU access is available, run:

```bash
ls -lh /data/bht/project01_phase1_reset_gpu/manifests/denovo_SER_hydrolase_full_input.pdb
find /data/bht/design_tools -maxdepth 4 \( -iname '*rfdiff*' -o -iname '*RFdiff*' -o -iname '*rfaa*' \)
find /data/bht/project01_phase1_reset_gpu -maxdepth 5 -type f | grep -Ei 'rfdiff|motif|scaffold|denovo_SER|rfaa'
conda env list 2>/dev/null || true
```

## Non-Goals For This Stage

```text
PLACER_required_for_sequence_acceptance = false
QMMM_current_stage = false
large_outputs_to_github = false
```

Do not upload PDBs, trajectories, model weights, raw PLACER folders, or large logs. Sync only lightweight status markdown, manifests, and scripts.

## 2026-06-30 Serine-Hydrolase Smoke Update

A minimal GPU Route B smoke has now run for serine hydrolase:

```text
RFAA raw new backbone outputs: 1
LigandMPNN designed sequences: 2
Postseq ESMFold/gate evaluated: 0
```

See:

```text
project01/phase1_reset_20260630/routeB_denovo_scaffold/serhyd_routeB_smoke_status_20260630.md
project01/phase1_reset_20260630/routeB_denovo_scaffold/serhyd_routeB_smoke_public_manifest.tsv
```

This changes serine-hydrolase Route B from `NOT_STARTED` to `SMOKE_RFAA_LIGANDMPNN_COMPLETE__POSTSEQ_GATE_NOT_EVALUATED`. It is still not a completed Route B panel.
