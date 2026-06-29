# 2026-06-29 serine-hydrolase LigandMPNN v2 shell-panel run

## Purpose

Execute the next GitHub plan step for the Baker-style serine-hydrolase example:

```text
fixed active pocket -> controlled sequence-background panel -> layer-wise identity audit
```

This run is a sequence-generation and mask-validation step only. It is not a PLACER/QMMM label-generation step.

## What changed

LigandMPNN sampled sequences at three temperatures while redesigning only the 5-12 A shell around the `bu2` ligand.

| Intended bin | Temperature | Sequences |
| --- | ---: | ---: |
| A1_low | 0.05 | 4 |
| B1_mid | 0.30 | 4 |
| C1_high | 0.70 | 4 |

## What stayed fixed

- Protein backbone.
- `bu2` ligand atom context.
- Catalytic/direct-contact pocket: catalytic residues plus all protein residues with any heavy atom <= 5.0 A from ligand heavy atoms.
- Distal residues outside 12 A.
- Ligand/reaction template and residue numbering.

## v2 mask decision

The original mask fixed 10 residues using <=4 A plus hard-coded catalytic contacts. Re-audit showed additional 4-5 A direct ligand contacts that should not vary in the first background-learning panel:

```text
A57 A67 A76 A99 A102 A125 A129
```

The accepted v2 mask fixes 17 residues:

```text
A57 A61 A67 A70 A71 A76 A78 A86 A95 A98 A99 A102 A125 A126 A127 A128 A129
```

The redesigned shell contains 60 residues:

```text
A6 A7 A9 A10 A11 A14 A49 A52 A53 A54 A55 A56 A58 A59 A60 A62 A63 A64 A65 A66
A68 A69 A72 A73 A74 A75 A77 A79 A80 A81 A82 A83 A84 A85 A87 A88 A89 A90 A91
A92 A93 A94 A96 A97 A100 A101 A103 A104 A105 A106 A109 A121 A123 A124 A130
A131 A132 A133 A134 A135
```

Remote evidence files:

```text
/Dell/Dell14/bianht/project01_sequence_panel_serhydrolase_20260629/fixed_residues_v2_direct5A.txt
/Dell/Dell14/bianht/project01_sequence_panel_serhydrolase_20260629/redesigned_shell_residues_v2_5to12A.txt
/Dell/Dell14/bianht/project01_sequence_panel_serhydrolase_20260629/residue_ligand_distance_layers_v2.tsv
/Dell/Dell14/bianht/project01_sequence_panel_serhydrolase_20260629/mask_summary_v2_direct5A_shell5to12A.json
```

## LigandMPNN validation

Smoke run:

- Model type: `ligand_mpnn`
- Ligand atoms parsed: 18
- Exit code: 0
- Designed sequence length: 160 aa
- Changed residues: 10
- Fixed/direct-contact changes: 0
- Changes outside redesigned shell: 0

Panel run:

```text
/Dell/Dell14/bianht/project01_sequence_panel_serhydrolase_20260629/ligandmpnn_v2_shell_panel_20260629_223527
```

Layer-wise identity table:

```text
/Dell/Dell14/bianht/project01_sequence_panel_serhydrolase_20260629/ligandmpnn_v2_shell_panel_20260629_223527/sequence_identity_layers_v2.tsv
```

Summary:

| Bin | n | Global identity range | Redesigned-shell identity range | Global changes |
| --- | ---: | ---: | ---: | ---: |
| A1_low | 4 | 0.91875-0.93750 | 0.78333-0.83333 | 10-13 |
| B1_mid | 4 | 0.91875-0.95000 | 0.78333-0.86667 | 8-13 |
| C1_high | 4 | 0.86875-0.91875 | 0.65000-0.78333 | 13-21 |

All 12 generated sequences passed hard checks:

```text
all_fixed_direct5A_unchanged = true
all_changes_within_redesigned_shell = true
```

## Interpretation

This successfully executes the first controlled shell-only sequence-panel step. The results are useful for testing PLACER and embedding extraction under fixed-pocket, varied-shell conditions.

However, the panel does not yet cover the full intended global identity bins from the research plan:

- It covers high global identity and moderate shell diversity.
- It does not produce true C1 global identity of 60-80%.
- This limitation is expected because only 60 of 160 positions are mutable in the shell-only mask.

Therefore, this panel should be treated as a B1-like shell-background pilot, not a complete A1/B1/C1 production panel.

## Label-quality status

No energy label was produced in this step. This is a technical sequence-panel result, not a QMMM/no-MM label-quality result.

## Next concrete action

1. Select 2-4 accepted sequences spanning the observed shell-identity range.
2. Run PLACER on those accepted sequence backgrounds using the same ligand/reaction template.
3. In parallel, define a broader v3 mask for true C1-style global diversity: keep catalytic/direct-contact fixed, allow selected distal/background positions to vary, and recompute global plus layer identities.
4. Do not start large QMMM scaling until PLACER reconstruction and layer-mask extraction are stable for the designed sequence backgrounds.
