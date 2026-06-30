# Project 01 Phase1 Reset Status - 2026-06-30

## Reset target

Run the sequence-similarity series through structure-aware gates before any downstream PLACER, TS-stability, or QMMM work. Fixed-template zero RMSD is not sufficient evidence that mutant folded structures preserve the active pocket.

The active control rule is now a sequence-generation loop:

```text
generate sequences/backbones -> predict full protein structure -> postseq entrance gate
-> add PASS rows to the sequence panel
-> if any bin lacks 10 postseq entrance-pass sequences, return to upstream generation
```

Current-stage acceptance rule:

```text
Each similarity bin must have 10 distinct accepted sequences.
Each accepted sequence must have:
  postseq_entrance_gate == PASS
```

PLACER and QMMM are deferred. The previous PLACER pilot is diagnostic only and is not used to reject sequences in this stage.

## Current audited denominator

- 50 selected sequences: 90/80/70/60/50 identity bins x 10 sequences per bin.
- All 50 have GPU ESMFold predicted structures under `/data/bht/project01_phase1_reset_gpu/postseq_structure_models/`.
- All 50 have postseq entrance-gate results under `/data/bht/project01_phase1_reset_gpu/postseq_structure_gate/`.

Postseq entrance-gate counts:

| Bin | Evaluated | PASS | FAIL | NOT_EVALUATED |
| --- | ---: | ---: | ---: | ---: |
| 90 | 10 | 9 | 1 | 0 |
| 80 | 10 | 8 | 2 | 0 |
| 70 | 10 | 4 | 6 | 0 |
| 60 | 10 | 0 | 10 | 0 |
| 50 | 10 | 0 | 10 | 0 |

Total:

- EVALUATED_WITH_OUTPUT: 50/50.
- Entrance PASS: 21/50.
- Entrance FAIL: 29/50.
- Missing or not evaluated: 0/50.

Important interpretation:

- 60/50 are evaluated FAIL under the current gate, not missing.
- 70 is partially productive but incomplete.
- 50/60/70 must continue with regeneration rounds until each bin has enough entrance-pass sequences.

## Implemented gate

Scripts:

- `projects/01-specialized-ts-aware-scorer/scripts/phase1_reset/postseq_entrance_gate.py`
- `projects/01-specialized-ts-aware-scorer/scripts/phase1_reset/run_postseq_gate_manifest.py`

Gate behavior:

- Align candidate to reference by Kabsch using all common backbone atoms.
- Evaluate global backbone RMSD, fixed-pocket backbone RMSD, catalytic heavy-atom RMSD, and protein-only key-distance deltas.
- Protein-only predictors may omit ligand at this entrance gate.
- Missing predicted PDB is `NOT_EVALUATED_MISSING_PDB`, not `FAIL`.
- Ligand/PLACER/QMMM gates are deferred and are not current-stage sequence acceptance conditions.

Thresholds:

```text
global_backbone_rmsd_max_A = 2.50
fixed_backbone_rmsd_max_A = 1.00
catalytic_heavy_rmsd_max_A = 0.75
protein_key_distance_abs_delta_max_A = 0.75
ligand_heavy_rmsd_max_A = DEFERRED
ligand_key_distance_abs_delta_max_A = DEFERRED
```

Threshold provenance:

- The Baker serine-hydrolase paper uses design-vs-experiment structural agreement as evidence that successful designs preserve the intended fold and active-site architecture.
- In this project, C-alpha/backbone RMSD belongs to the post-sequence structure gate: it asks whether the generated sequence still folds into the intended backbone/pocket before any downstream TS or QMMM calculation.
- `ligand_heavy_rmsd_max_A = 0.75` was a project-level diagnostic PLACER/full-ligand threshold, not a confirmed Baker hard cutoff, and is now deferred for the current sequence-generation stage.
- PLACER failure is not a sequence rejection reason in this stage because a TS-like conformer is not expected to be a low-energy ligand pose. The next scientific question is whether explicit TS conformer ensembles can be stabilized by the accepted protein backgrounds.

## GPU-local prediction status

GPU host:

```text
dell-PowerEdge-R760
/data/bht/project01_phase1_reset_gpu
```

ESMFold environment is now passing:

- `openfold.model.structure_module.StructureModule`: OK.
- `esm.esmfold.v1.pretrained`: OK.
- ESMFold smoke wrote `bin90_rank1_rec1_esmfold_smoke.pdb`.
- The smoke structure passed postseq entrance gate.

All-50 prediction/gate run:

```text
/data/bht/project01_phase1_reset_gpu/postseq_structure_gate/tables/all50_entrance_gate.tsv
/data/bht/project01_phase1_reset_gpu/manifests/placer_n50_entrance_pass_queue.gpu.tsv
```

The entrance-pass queue contains 21 samples.

## PLACER diagnostic status

Completed PLACER pilot:

- Inputs: one representative entrance-pass sample each from 90/80/70.
- Evaluated universe: 300 PLACER conformers from 6 sequence-run entries.
- `CROP_STRICT_PASS`: 0.
- NOT_EVALUATED: 0.

Current interpretation:

- These PLACER outputs remain useful diagnostic evidence that ligand-pose preservation is not currently reliable.
- They are not current-stage sequence acceptance evidence.
- Do not regenerate or discard sequences solely because PLACER failed the crop strict gate.
- Do not launch QMMM during the current sequence-generation stage.

Lightweight evidence files:

```text
project01/phase1_reset_20260630/two_layer_placer_crop_gate_summary.json
project01/phase1_reset_20260630/two_layer_placer_crop_gate_sequence_summary.tsv
projects/01-specialized-ts-aware-scorer/scripts/phase1_reset/two_layer_placer_crop_gate.py
```

## Regeneration policy

Active controlling plan:

- `project01/phase1_reset_20260630/regeneration_closed_loop_plan_20260630.md`

Targets:

```text
postseq_entrance_pass_sequences_per_similarity_bin = 10
defer_placer_filtering = true
defer_qmmm_calculation = true
```

Current entrance-pass shortage:

```text
90: needs 1 more pass
80: needs 2 more passes
70: needs 6 more passes
60: needs 10 more passes
50: needs 10 more passes
```

Required next action:

1. Build or refresh `regeneration_queue_round02.tsv` using postseq entrance-pass counts.
2. Generate replacement candidates for shortage bins, especially 70/60/50.
3. Predict and gate the new candidates on GPU immediately after generation.
4. Add only postseq entrance-pass rows to the sequence panel manifest.
5. Continue until every bin has 10 accepted postseq entrance-pass sequences.
6. Do not launch PLACER or QMMM as a current-stage requirement.

## GitHub sync rule

Commit only:

- markdown status files
- gate summaries
- small TSV/JSON manifests
- scripts

Do not commit:

- PDBs
- PLACER model outputs
- model checkpoints
- trajectories
- large logs
- QM/MM raw outputs
