# Project 01 Phase1 Next Generation Queue - 2026-06-30

## Current Decision

Proceed with two separated generation tracks:

1. Natural-scaffold conserved-site redesign.
2. Baker-style active-site constrained de novo scaffold generation.

Route-state correction on 2026-06-30: the completed serine-hydrolase panel is an existing/reference-scaffold sequence panel, not proof that new backbones have been generated. KSI natural-scaffold sequence generation is complete for the current stage. True active-site constrained de novo backbone generation remains not started for serine hydrolase and deferred for KSI until the TS-like motif is defined.

Updated generation rule:

```text
natural_scaffold_track:
  sequence alignment first
  fixed positions = catalytic + ligand/TS-contact + pocket-shell + MSA-conserved core
  generate 90/80/70/60/50 bins only over allowed mutable positions

de_novo_scaffold_track:
  active-site/reactive geometry first
  generate new backbones around that motif
  sequence generated backbones and then apply the same postseq entrance gate
```

These tracks must remain separately labeled in manifests because their similarity bins have different meanings.

## Current Route Status

| Enzyme system | Track | Current status | Evidence | Next action |
| --- | --- | --- | --- | --- |
| serine hydrolase | existing/reference scaffold sequence panel | COMPLETE | `/data/bht/project01_phase1_reset_gpu/manifests/postseq_entrance_pass_sequence_panel_capped10.tsv` | freeze as Route A panel |
| serine hydrolase | true active-site constrained de novo/new-backbone | NOT_STARTED | no generated-backbone manifest yet | extract motif constraints and run RFdiffusion/RFdiffusionAA smoke when remote access is restored |
| serine hydrolase | natural scaffold | NOT_STARTED | no selected natural family/scaffold yet | select natural hydrolase family, build MSA/fixed mask |
| KSI | natural scaffold | SEQUENCE_PANEL_COMPLETE | `/data/bht/project01_phase1_reset_gpu/natural_scaffold/KSI/manifests/ksi_final_postseq_entrance_pass_sequence_panel_capped10.tsv` | freeze capped10 current-stage panel |
| KSI | de novo/motif scaffold | DEFERRED_TS_GEOMETRY_NOT_DEFINED | no KSI TS-like steroid/enolate motif definition yet | define ligand/protonation/motif geometry first |

## Immediate Queue

### 1. Serine-Hydrolase True De Novo/New-Backbone Track

Status: `NOT_STARTED`.

This is now the first missing route to start. The existing capped10 serine-hydrolase sequence panel should not be counted as de novo backbone completion.

Input motif/reference starting point:

```text
/data/bht/project01_phase1_reset_gpu/manifests/denovo_SER_hydrolase_full_input.pdb
```

First executable actions on the GPU host:

```text
1. Confirm the reference PDB exists and identify catalytic/motif residues plus ligand/reactive atoms.
2. Check available motif-scaffolding tools: RFdiffusion, RFdiffusionAA, or existing project scripts.
3. Generate a small smoke batch of new backbones around the fixed active-site/reactive motif.
4. Sequence generated backbones with ProteinMPNN/LigandMPNN.
5. Predict/validate full protein structures.
6. Apply the same postseq entrance gate before expanding to 90/80/70/60/50-style panels.
```

Do not use PLACER or QMMM as acceptance gates for this current sequence/backbone stage.

### 2. KSI Natural-Scaffold Track

Status: `SEQUENCE_PANEL_COMPLETE` for the current stage.

Accepted distinct PASS counts from the final postseq entrance gate:

| Bin | Accepted distinct PASS | Target |
| --- | ---: | ---: |
| 90 | 13 | 10 |
| 80 | 22 | 10 |
| 70 | 137 | 10 |
| 60 | 94 | 10 |
| 50 | 95 | 10 |

Final capped panel:

```text
/data/bht/project01_phase1_reset_gpu/natural_scaffold/KSI/manifests/ksi_final_postseq_entrance_pass_sequence_panel_capped10.tsv
```

### 3. Natural Serine-Hydrolase Track

Status: `NOT_STARTED`.

This is separate from the Baker-style de novo serine hydrolase. Before generating sequences, select a natural hydrolase family/scaffold, build the MSA, and derive a conservation-aware fixed-site mask.

## Current Stage Acceptance

For sequence generation:

```text
accepted_sequence = postseq_entrance_gate == PASS
PLACER_required_for_sequence_acceptance = false
QMMM_current_stage = false
```

PLACER and QMMM are downstream analyses and are not sequence-stage acceptance gates.

## Active Audit Questions

Before accepting any generated sequence panel:

1. What is the evaluated universe and denominator?
2. Which positions were fixed by catalytic chemistry?
3. Which positions were fixed by ligand or TS contact?
4. Which positions were fixed by MSA conservation or backbone generation constraints?
5. Which positions were actually mutable?
6. Is the identity bin measured against a natural scaffold, a family reference, a generated scaffold family, or a design reference?
7. Does the predicted structure pass the postseq entrance gate?
8. Which candidates are NOT_EVALUATED rather than FAIL?
