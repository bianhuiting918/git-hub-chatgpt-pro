# Project 01 Phase1 Next Generation Queue - 2026-06-30

## Current Decision

Proceed with two separated generation tracks:

1. Natural-scaffold conserved-site redesign.
2. Baker-style active-site constrained de novo scaffold generation.

The current serine-hydrolase sequence panel for the local `denovo_SER_hydrolase` reference is complete at 10 accepted sequences per 90/80/70/60/50 bin. Do not keep mutating that same de novo reference as the main route for new biological diversity.

## Immediate Queue

### 1. KSI Natural-Scaffold Track

Status: `NOT_STARTED`.

Reason to start first:

- KSI is the mechanistic control system.
- It has a compact active-site electrostatic/oxyanion-hole question.
- The natural-scaffold workflow directly matches the user's requested MSA -> conserved-site fixed mask -> new sequence generation logic.

Initial candidate reference:

```text
RCSB 6UBQ: pKSI bound to 4-androstenedione
```

Reference sources:

```text
https://www.rcsb.org/structure/6UBQ
https://www.ebi.ac.uk/thornton-srv/m-csa/entry/349/
```

Initial fixed-site hypothesis:

```text
Tyr16
Asp40
Asp103
ligand/direct-contact residues
MSA conserved structural core
```

Required next files:

```text
natural_scaffold_msa_summary.tsv
natural_scaffold_fixed_positions.tsv
natural_scaffold_mutable_positions.tsv
natural_scaffold_generation_manifest.tsv
```

Do not label missing homologs, failed downloads, or unmapped residues as `FAIL`; report them as `NOT_EVALUATED`.

### 2. Serine-Hydrolase De Novo Track

Status: `SEQUENCE_PANEL_COMPLETE` for current local reference.

Use existing capped panel:

```text
/data/bht/project01_phase1_reset_gpu/manifests/postseq_entrance_pass_sequence_panel_capped10.tsv
```

Next de novo expansion should not be random mutation of this one reference. It should use active-site constrained backbone generation around the catalytic/reactive motif, then sequence and postseq-gate the generated backbones.

### 3. Natural Serine-Hydrolase Track

Status: `NOT_STARTED`.

This is separate from the Baker-style de novo serine hydrolase. Before generating sequences, select a natural hydrolase family/scaffold, build the MSA, and derive a conservation-aware fixed-site mask.

## Current Stage Acceptance

For sequence generation:

```text
accepted_sequence = postseq_entrance_gate == PASS
```

PLACER and QMMM are downstream analyses and are not sequence-stage acceptance gates.

## Active Audit Questions

Before accepting any generated sequence panel:

1. Which positions were fixed by catalytic chemistry?
2. Which positions were fixed by ligand or TS contact?
3. Which positions were fixed by MSA conservation?
4. Which positions were actually mutable?
5. Is the identity bin measured against a natural scaffold, a family reference, or a de novo design reference?
6. Does the predicted structure pass the postseq entrance gate?
