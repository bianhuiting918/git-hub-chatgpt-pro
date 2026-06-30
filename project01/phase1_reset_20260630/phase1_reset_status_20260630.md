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

Current entrance-pass status after GPU round02/round03/round04/round05:

```text
90: complete; 16 distinct entrance-pass sequences available
80: complete; 11 distinct entrance-pass sequences available
70: complete; 38 distinct entrance-pass sequences available
60: complete; 16 distinct entrance-pass sequences available
50: complete; 11 distinct entrance-pass sequences available
```

Audited added batches:

| Batch | Evaluated | ESMFold OK | Entrance PASS | PASS bins |
| --- | ---: | ---: | ---: | --- |
| round02 controlled mutation-count candidates | 164 | 164 | 7 | 90:7 |
| round02d actual-bin LigandMPNN refilter | 66 | 66 | 18 | 70:5, 60:2, 50:11 |
| round03 actual-bin production | 152 | 152 | 33 | 80:3, 70:29, 60:1 |
| round04 empirical 60 subsets | 144 | 144 | 5 | 60:5 |
| round05 template2 60 subset | 72 | 72 | 8 | 60:8 |

Combined accepted distinct sequence-panel counts:

| Bin | Accepted distinct PASS | Target | Still needed |
| --- | ---: | ---: | ---: |
| 90 | 16 | 10 | 0 |
| 80 | 11 | 10 | 0 |
| 70 | 38 | 10 | 0 |
| 60 | 16 | 10 | 0 |
| 50 | 11 | 10 | 0 |

Previous pre-round02 shortage was:

```text
90: needed 1 more pass
80: needs 2 more passes
70: needs 6 more passes
60: needs 10 more passes
50: needs 10 more passes
```

Required next action:

1. Keep the accepted sequence panel at `/data/bht/project01_phase1_reset_gpu/manifests/postseq_entrance_pass_sequence_panel.tsv`.
2. Use the capped panel at `/data/bht/project01_phase1_reset_gpu/manifests/postseq_entrance_pass_sequence_panel_capped10.tsv` as the 50-row current-stage panel.
3. Use the dual-track scaffold strategy in `project01/phase1_reset_20260630/dual_track_scaffold_strategy_20260630.md`.
4. For natural scaffolds, use sequence alignment to fix conserved catalytic/pocket/core positions before generating new sequences. The executable protocol is now `project01/phase1_reset_20260630/natural_scaffold_msa_generation_protocol_20260630.md`.
5. For Baker-style diversity, use active-site constrained de novo scaffold generation as a separate track rather than forcing one reference backbone to carry 50/60% identity mutations.
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
