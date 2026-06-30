# Project 01 Phase1 Reset Status - 2026-06-30

## Reset target

Run the sequence-similarity series through structure-aware gates before PLACER and QMMM. Fixed-template zero RMSD is not sufficient evidence that mutant folded structures preserve the active pocket.

The active control rule is now a closed loop:

```text
generate sequences/backbones -> predict full protein structure -> postseq entrance gate
-> PLACER n=50 only for entrance-pass samples -> crop ligand RMSD/reaction gate
-> full-ligand completion only for accepted sequences -> strict QMMM manifest
-> if any bin lacks 10 accepted sequences, return to the nearest upstream stage
```

Final acceptance rule:

```text
Each similarity bin must have 10 distinct accepted sequences.
Each accepted sequence must have:
  postseq_entrance_gate == PASS
  PLACER_conformers_generated == 50
  PLACER_crop_RMSD_gate_PASS_count >= 10
  full_ligand_strict_QMMM_ready_count >= 1
```

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
- Ligand/key-distance gates become mandatory at PLACER crop and completion stages.
- Missing predicted PDB is `NOT_EVALUATED_MISSING_PDB`, not `FAIL`.

Thresholds:

```text
global_backbone_rmsd_max_A = 2.50
fixed_backbone_rmsd_max_A = 1.00
catalytic_heavy_rmsd_max_A = 0.75
protein_key_distance_abs_delta_max_A = 0.75
ligand_heavy_rmsd_max_A = 0.75
ligand_key_distance_abs_delta_max_A = 0.50
```

Threshold provenance:

- The Baker serine-hydrolase paper uses design-vs-experiment structural agreement as evidence that successful designs preserve the intended fold and active-site architecture. Reported examples include C-alpha RMSD around the sub-angstrom to low-angstrom range and active-site mean all-atom RMSD around 0.54 to 0.7 A for successful refined designs.
- In this project, C-alpha/backbone RMSD belongs mainly to the pre-PLACER post-sequence structure gate: it asks whether the generated sequence still folds into the intended backbone/pocket before ligand conformer sampling.
- `ligand_heavy_rmsd_max_A = 0.75` is our current `PROJECT_STRICT_GATE` for PLACER crop/full-ligand ligand-pose preservation. It is inspired by the sub-angstrom active-site agreement in the Baker work, but it is not currently confirmed as a direct Baker hard cutoff.
- Do not describe the 0.75 A ligand RMSD gate as `BAKER_LITERATURE_GATE` unless the exact threshold is later verified from the paper's methods or released filtering code.

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

## PLACER status

Holo input preparation:

```text
/data/bht/project01_phase1_reset_gpu/manifests/placer_n50_holo_inputs.gpu.tsv
/data/bht/project01_phase1_reset_gpu/placer_inputs_holo/
```

Policy:

```text
append_reference_bu2_refpose_before_PLACER
```

This holo input is only the initial PLACER pose. It is not final ligand evidence.

Reference bond constraint:

```text
A-128-SER-OG:X-1-bu2-C1:1.533
```

Completed PLACER pilot:

- Inputs: one representative entrance-pass sample each from 90/80/70.
- Samples: 50 PLACER conformers per representative input.
- Total conformers generated: 150.
- Output directory: `/data/bht/project01_phase1_reset_gpu/placer_runs/representative_90_80_70_n50_bonded/`.
- Runtime: 170.61 seconds.

The PLACER raw CSV/PDB outputs are not strict-pass evidence.

Crop ligand RMSD/reaction gate result for the first representative PLACER run:

- Evaluated universe: 150 PLACER conformers from 3 sequences.
- PASS: 0.
- FAIL: 150.
- NOT_EVALUATED: 0.
- Per-sequence threshold result: all 3 sequences fail `minimum_placer_crop_pass_conformers_per_sequence = 10`.
- Best observed ligand heavy-atom RMSD was still 2.12 A, above the 0.75 A gate.

Therefore the first representative PLACER run produces no QMMM-ready sequence and no full-ligand completion queue.

## Regeneration policy

Active controlling plan:

- `project01/phase1_reset_20260630/regeneration_closed_loop_plan_20260630.md`

Targets:

```text
final_accepted_sequences_per_similarity_bin = 10
placer_samples_per_postseq_pass_sequence = 50
minimum_placer_crop_pass_conformers_per_sequence = 10
minimum_qmmm_ready_strict_pass_conformers_per_accepted_sequence = 1
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

1. Treat the first 90/80/70 representative PLACER run as failed under the new per-sequence `>=10/50` crop-pass rule.
2. Use the running pose-centered PLACER batch only as a diagnostic; it must still meet `>=10/50` to count.
3. Build `regeneration_queue_round02.tsv`.
4. Generate replacement candidates for shortage bins, especially 70/60/50.
5. Predict and gate the new candidates on GPU immediately after generation.
6. Launch PLACER only for postseq entrance-pass samples.
7. Discard any sequence with fewer than 10/50 PLACER crop-pass conformers and regenerate a replacement from upstream.
8. Continue until every bin has 10 accepted sequences ready for QMMM.

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
