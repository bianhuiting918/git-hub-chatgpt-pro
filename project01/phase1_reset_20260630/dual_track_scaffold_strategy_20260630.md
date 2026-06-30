# Project 01 Phase1 Dual-Track Scaffold Strategy - 2026-06-30

## Decision

Project 01 Phase1 now uses two explicitly separated sequence/backbone-generation tracks. Route-state correction on 2026-06-30: completed existing/reference-scaffold sequence panels must not be counted as true de novo/new-backbone completion.

1. Natural-scaffold track: start from known enzyme family scaffolds. Use sequence alignment to identify conserved catalytic, pocket, and structurally conserved residues; keep these fixed or strongly constrained; generate different-similarity sequences only over tolerated background positions.
2. De novo scaffold track: follow the Baker-style serine-hydrolase logic. Fix the active-site/reactive geometry first, then generate new backbones that support that motif. Sequence identity bins in this track should compare generated sequences within a scaffold family or against the chosen design reference, but should not be interpreted as mutations on one fixed natural backbone.

This avoids mixing two different biological questions:

- Natural-scaffold question: how much background sequence variation can a known fold tolerate while preserving the active pocket?
- De novo-scaffold question: how many different folds/backgrounds can be built around the same catalytic geometry?

## Why This Change Was Needed

The single-reference mutation strategy is not robust for low-identity bins. GPU evidence from round02 showed that forcing many mutations on one reference sequence often makes ESMFold move the active pocket:

| Batch | Evaluated | ESMFold OK | Entrance PASS | Key interpretation |
| --- | ---: | ---: | ---: | --- |
| initial all50 | 50 | 50 | 21 | 90/80/70 partially productive; 60/50 failed. |
| round02 controlled mutation count | 164 | 164 | 7 | PASS only in 90%; random low-identity mutations drifted the pocket. |
| round02d actual-bin LigandMPNN refilter | 66 | 66 | 18 | Cooperative LigandMPNN designs restored some 70/60/50 passes. |

Combined accepted distinct sequence-panel counts after these evaluated batches:

| Bin | Accepted distinct PASS | Target | Still needed |
| --- | ---: | ---: | ---: |
| 90 | 16 | 10 | 0 |
| 80 | 8 | 10 | 2 |
| 70 | 9 | 10 | 1 |
| 60 | 2 | 10 | 8 |
| 50 | 11 | 10 | 0 |

That shortage has now been closed by later regeneration rounds. The current capped sequence panel has at least 10 entrance-pass sequences in every 90/80/70/60/50 bin. This historical table is retained to document why the workflow was split, but it is no longer the active shortage.

Current accepted distinct sequence-panel counts:

| Bin | Accepted distinct PASS | Target | Still needed |
| --- | ---: | ---: | ---: |
| 90 | 16 | 10 | 0 |
| 80 | 11 | 10 | 0 |
| 70 | 38 | 10 | 0 |
| 60 | 16 | 10 | 0 |
| 50 | 11 | 10 | 0 |

## Natural-Scaffold Track

Goal: generate sequence-diverse variants on natural or known-family scaffolds while preserving the active pocket.

Hard rule: do not mutate a natural scaffold by random identity-bin forcing before conservation analysis. Natural-scaffold generation must start from sequence alignment, fixed-site annotation, and a mutable-position mask.

Required steps:

1. Collect homologs for the selected enzyme family or scaffold.
2. Run MSA and annotate conserved positions.
3. Define fixed groups:
   - catalytic residues and ligand/direct-contact residues;
   - MSA-conserved structural core residues;
   - pocket shell residues required to keep catalytic geometry stable.
4. Define mutable groups:
   - nonconserved, solvent/background residues;
   - second-shell positions only when pilot gates show they are tolerated.
5. Generate candidates with LigandMPNN/ProteinMPNN using fixed conserved positions.
6. Predict full structures with ESMFold or a stronger sequence+structure predictor when available.
7. Accept only rows with `postseq_entrance_gate == PASS`.

Identity-bin definition for this track:

```text
identity_bin = sequence identity to the selected natural scaffold, family representative, or MSA consensus
fixed_positions = catalytic residues + direct ligand/TS-contact residues + MSA-conserved pocket/core residues
mutable_positions = nonconserved background residues that are not known pocket-drift drivers
```

If a requested identity bin cannot be reached without mutating fixed positions, the row is `NOT_GENERATED_CONSERVATION_CONFLICT`, not a failed biological candidate.

Detailed natural-scaffold execution protocol:

```text
project01/phase1_reset_20260630/natural_scaffold_msa_generation_protocol_20260630.md
```

The current capped sequence panel is complete for the first serine-hydrolase existing/reference-scaffold test system. It is not a true de novo/new-backbone result. New enzyme systems or new natural scaffold families must start from the MSA-constrained fixed-site mask instead of random low-identity mutation.

## De Novo Scaffold Track

Goal: reproduce the Baker-style design logic at the workflow level.

Key rule: fixed active-site/reactive geometry comes before backbone generation. The scaffold is generated around the motif; it is not a natural backbone being heavily mutated.

Required steps:

1. Define the catalytic motif and TS/reactive ligand geometry.
2. Use active-site constrained backbone generation, such as RFdiffusion-style motif scaffolding, when GPU capacity is available.
3. Sequence generated backbones with ProteinMPNN/LigandMPNN.
4. Predict or validate full structures.
5. Run the same postseq entrance gate for the protein pocket.
6. Keep de novo candidates in a separate manifest from natural-scaffold candidates until downstream comparison is explicitly defined.

This track is the correct way to obtain low sequence similarity while still preserving a designed active center. As of this update, the serine-hydrolase true new-backbone Route B has not yet started; the next step is motif extraction and a small RFdiffusion/RFdiffusionAA-style smoke batch.

Identity-bin definition for this track:

```text
identity_bin = sequence identity within a generated scaffold family or against a chosen design reference
fixed_positions = active-site motif and geometry-support residues required by the generated backbone
mutable_positions = designable positions allowed by the generated scaffold model
```

Do not directly compare natural-scaffold identity bins and de novo-scaffold identity bins as if they were the same perturbation. They answer different questions and must remain separately labeled in manifests.

## Current Acceptance Rule

For the current stage:

```text
accepted_sequence = postseq_entrance_gate == PASS
PLACER_required_for_sequence_acceptance = false
QMMM_current_stage = false
```

PLACER and QMMM remain downstream. PLACER failure is diagnostic evidence about ligand/reaction-pose handling, not a sequence-generation rejection reason in this stage.

## Lightweight Remote Evidence

Remote GPU root:

```text
/data/bht/project01_phase1_reset_gpu
```

Key files:

```text
manifests/postseq_entrance_pass_sequence_panel.tsv
manifests/postseq_entrance_pass_sequence_panel_summary.json
manifests/postseq_entrance_gate_manifest_round02.tsv
manifests/postseq_entrance_gate_manifest_round02d_actual_bins.tsv
postseq_structure_gate/tables/round02_esmfold_prediction_summary.json
postseq_structure_gate/tables/round02d_esmfold_prediction_summary.json
```

Do not upload predicted PDBs, raw PLACER outputs, trajectories, or large logs.

