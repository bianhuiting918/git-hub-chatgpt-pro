# KSI final postseq entrance-pass sequence panel - 2026-06-30

## Scope

This file records the KSI natural-scaffold sequence-generation stage after post-sequence ESMFold prediction and KSI pocket-geometry screening. PLACER was not used as a sequence gate.

## Final result

The KSI sequence panel reached the current stage target: 10 accepted sequences per similarity level.

| Bin | Distinct PASS available | Capped selected | Needed to 10 |
| --- | ---: | ---: | ---: |
| 90 | 13 | 10 | 0 |
| 80 | 22 | 10 | 0 |
| 70 | 137 | 10 | 0 |
| 60 | 94 | 10 | 0 |
| 50 | 95 | 10 | 0 |

Final capped panel: 50 sequences, all `entrance_gate_status=PASS`, with 10 unique sequences per bin.

## Gate definition

KSI-specific postseq gate used 6UBQ chain A as reference and checked:

- global backbone RMSD <= 2.5 A
- fixed pocket/shell backbone RMSD <= 1.0 A
- catalytic heavy-atom RMSD <= 0.75 A for Tyr16, Asp40, Asp103
- max active-site protein key-distance delta <= 0.75 A

## Iteration evidence

- Round01 all selected: 31 PASS / 192 evaluated. Counts: 90=9, 80=22, 70=0, 60=0, 50=0.
- Round02 expanded pocket shell: 141 PASS / 248 evaluated. Counts: 90=4, 70=137.
- Round03 low-bin expanded shell: 189 PASS / 200 evaluated. Counts: 60=94, 50=95.

The expanded shell fixed 45 extra pocket/shell residues in addition to the original 12 MSA/catalytic fixed residues. Exact 50% mutation targeting is constrained by this shell, so the accepted 50 bin is the feasible expanded-shell low-identity bin with identity 0.531746 in the capped panel.

## Remote authoritative files

- Final all-pass panel: `/data/bht/project01_phase1_reset_gpu/natural_scaffold/KSI/manifests/ksi_final_postseq_entrance_pass_sequence_panel.tsv`
- Final capped10 panel: `/data/bht/project01_phase1_reset_gpu/natural_scaffold/KSI/manifests/ksi_final_postseq_entrance_pass_sequence_panel_capped10.tsv`
- Final summary: `/data/bht/project01_phase1_reset_gpu/natural_scaffold/KSI/postseq_structure_gate/tables/ksi_final_postseq_entrance_pass_sequence_panel_summary.json`

Full FASTA, PDB structures, gate JSON directories, and long per-candidate manifests are intentionally not uploaded to GitHub.
