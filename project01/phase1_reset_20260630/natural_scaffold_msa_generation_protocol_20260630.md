# Project 01 Natural-Scaffold MSA-Constrained Generation Protocol - 2026-06-30

## Scope

This protocol is for the natural-scaffold track. It is separate from the Baker-style de novo scaffold track.

Natural-scaffold question:

```text
Given an existing enzyme fold/family, how much sequence-background variation can be introduced while preserving the catalytic pocket geometry?
```

Do not treat this as the same problem as generating a new backbone around an active-site motif.

## Inputs

Required inputs for each enzyme system:

1. Reference holo or apo scaffold with annotated catalytic residues.
2. Homologous natural sequences from the same enzyme family or mechanistically related subfamily.
3. MSA with residue numbering mapped back to the reference scaffold.
4. Pocket annotation:
   - catalytic residues;
   - ligand or TS-contact residues;
   - first-shell residues;
   - second-shell residues;
   - structurally conserved core residues.

## Fixed-Site Definition

Positions are fixed during sequence generation if any of the following is true:

1. Catalytic atom geometry would be directly affected by mutation.
2. The position contacts ligand, TS analog, catalytic water, oxyanion hole, or proton-transfer group.
3. The MSA shows strong family conservation and the position is part of the fold core or pocket support.
4. The position was already identified by post-sequence entrance-gate failures as sensitive to pocket drift.

Initial project-policy thresholds:

```text
MSA_conservation_fixed_if_major_amino_acid_frequency >= 0.80
MSA_conservation_review_if_major_amino_acid_frequency = 0.60-0.80
always_fixed = catalytic residues + direct TS/ligand contact residues
mutable_default = nonconserved solvent/background positions
```

These thresholds are project policy, not tool-native pass/fail criteria. They should be reported separately from model outputs.

## Mutable-Site Definition

Positions are mutable only when they satisfy all conditions:

1. Not catalytic.
2. Not direct ligand/TS contact.
3. Not required for oxyanion/proton-transfer geometry.
4. Not strongly conserved in the MSA.
5. Not previously associated with post-sequence pocket drift.

Second-shell residues are mutable only after a pilot round shows that their mutation does not break the post-sequence entrance gate.

## Generation Loop

For each similarity bin:

```text
50 / 60 / 70 / 80 / 90 percent identity
```

run:

```text
MSA -> fixed-site mask -> ProteinMPNN/LigandMPNN generation
-> full-structure prediction
-> post-sequence entrance gate
-> accepted sequence panel
```

Acceptance at this stage:

```text
postseq_entrance_gate == PASS
```

Current sequence-stage targets:

```text
10 accepted distinct sequences per similarity bin
50 accepted sequences total per enzyme system
```

PLACER and QMMM are downstream analyses. They are not sequence-generation acceptance gates for this stage.

## Required Evidence Files

Each natural-scaffold run should produce lightweight evidence:

```text
natural_scaffold_msa_summary.tsv
natural_scaffold_fixed_positions.tsv
natural_scaffold_mutable_positions.tsv
natural_scaffold_generation_manifest.tsv
postseq_entrance_gate_summary.json
postseq_entrance_pass_sequence_panel.tsv
```

Large files not uploaded to GitHub:

```text
raw MSA databases
predicted PDBs
PLACER outputs
trajectories
model checkpoints
QM/MM raw outputs
```

## De Novo Track Boundary

The de novo track should use active-site constrained backbone generation. Its candidate manifests must stay separate from natural-scaffold candidates until downstream comparison is intentionally defined.

De novo acceptance still uses the same post-sequence entrance gate, but the identity-bin interpretation is different:

- natural scaffold: identity is relative to a natural/reference scaffold or family representative;
- de novo scaffold: identity is relative within a generated scaffold family or chosen design reference.
