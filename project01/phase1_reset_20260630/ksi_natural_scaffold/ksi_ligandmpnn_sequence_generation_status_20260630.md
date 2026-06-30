# KSI natural-scaffold sequence generation status - 2026-06-30

## Scope

This status is for the reset Project 01 phase1 sequence-generation stage only. PLACER is not an entrance gate for this stage. PLACER remains downstream and optional for later conformer diagnostics.

## Input and constraints

- Scaffold: KSI natural scaffold, PDB 6UBQ chain A.
- Ligand context: ASD in the 6UBQ holo reference used by LigandMPNN ligand-aware design.
- Fixed pocket/conservation positions from PyFAMSA MSA mask: 12 positions.
- Mutable design positions: 100 positions.
- Output target: five sequence-identity bins, 90/80/70/60/50 percent, with at least 10 accepted sequences per bin.
- Acceptance at this stage: no fixed-position mutations, no mutations outside the redesignable mask, and identity in the intended bin.

## Remote run evidence

- Remote workspace: `/data/bht/project01_phase1_reset_gpu`.
- Script: `/data/bht/project01_phase1_reset_gpu/scripts/ksi_generate_ligandmpnn_bins.py`.
- Command log: `natural_scaffold/KSI/manifests/ksi_ligandmpnn_sequence_bins_chainA_rerun90.run.log`.
- Run exit code: `0`.
- Finished: `2026-06-30T15:58:54`.
- Full selected TSV on GPU: `natural_scaffold/KSI/manifests/ksi_ligandmpnn_sequence_bins_selected.tsv`.
- Full selected FASTA on GPU: `natural_scaffold/KSI/manifests/ksi_ligandmpnn_sequence_bins_selected.fasta`.

## Result

| Bin | Available pass | Capped selected | Identity | Mutation count | Fixed mutations | Outside-mask mutations |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 90 | 12 | 10 | 0.912698 | 11 | 0 | 0 |
| 80 | 30 | 10 | 0.809524 | 24 | 0 | 0 |
| 70 | 40 | 10 | 0.706349 | 37 | 0 | 0 |
| 60 | 50 | 10 | 0.611111 | 49 | 0 | 0 |
| 50 | 60 | 10 | 0.507937 | 62 | 0 | 0 |

Total capped panel: 50 KSI sequences, 10 per similarity level.

## Files synced here

- `ksi_ligandmpnn_sequence_bins_summary.json`: lightweight run summary.
- `ksi_ligandmpnn_sequence_bins_capped10_summary.json`: capped 10/bin summary.
- `ksi_ligandmpnn_sequence_bins_capped10_public_manifest.tsv`: lightweight per-sample manifest without full sequences.
- `ksi_postseq_structure_prediction_queue_capped10.tsv`: next-stage queue; `placer_required=no`.

Full sequence FASTA/TSV and model outputs are intentionally not uploaded to GitHub. They remain on the GPU server.

## Next step

Run post-sequence structure prediction and geometry screening for the 50 capped KSI sequences. This next step should verify structure/pocket preservation before any later embedding, PLACER, or QMMM decisions.
