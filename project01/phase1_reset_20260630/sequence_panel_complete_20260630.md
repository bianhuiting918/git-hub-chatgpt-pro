# Project 01 Phase1 Sequence Panel Complete - 2026-06-30

## Result

The current sequence-generation stage is complete.

Completion rule:

```text
For each similarity bin 90/80/70/60/50:
  count(distinct sequences with postseq_entrance_gate == PASS) >= 10
```

Final accepted distinct counts on the GPU working directory:

| Bin | Accepted distinct PASS | Target | Status |
| --- | ---: | ---: | --- |
| 90 | 16 | 10 | COMPLETE |
| 80 | 11 | 10 | COMPLETE |
| 70 | 38 | 10 | COMPLETE |
| 60 | 16 | 10 | COMPLETE |
| 50 | 11 | 10 | COMPLETE |

Total accepted distinct rows: 92.

A capped panel with exactly 10 rows per bin was written for downstream use:

```text
/data/bht/project01_phase1_reset_gpu/manifests/postseq_entrance_pass_sequence_panel_capped10.tsv
```

The uncapped accepted panel remains available at:

```text
/data/bht/project01_phase1_reset_gpu/manifests/postseq_entrance_pass_sequence_panel.tsv
```

## Evaluated Batches

| Batch | Evaluated | Structure prediction OK | Entrance PASS | Role |
| --- | ---: | ---: | ---: | --- |
| initial all50 | 50 | 50 | 21 | Initial panel; 90/80/70 partially successful. |
| round02 controlled mutation-count candidates | 164 | 164 | 7 | PASS only in 90; showed random low-identity mutation drifts pocket. |
| round02d actual-bin LigandMPNN refilter | 66 | 66 | 18 | Recovered 70/60/50 passes from cooperative LigandMPNN designs. |
| round03 actual-bin production | 152 | 152 | 33 | Completed 80/70 and added one 60. |
| round04 subset72 empirical 60 | 72 | 72 | 3 | Used empirical 60-pass mutation position sets. |
| round04 next72 empirical 60 | 72 | 72 | 2 | Continued best empirical 60 blocks. |
| round05 subset72 template2 60 | 72 | 72 | 8 | Completed 60 bin. |

All reported PASS/FAIL values are downstream project entrance-gate labels, not native ESMFold or LigandMPNN labels.

## Current Gate

Current accepted rows satisfy:

```text
postseq_entrance_gate == PASS
```

Entrance-gate thresholds:

```text
global_backbone_rmsd_max_A = 2.50
fixed_backbone_rmsd_max_A = 1.00
catalytic_heavy_rmsd_max_A = 0.75
protein_key_distance_abs_delta_max_A = 0.75
```

Protein-only predictors are allowed to omit ligand at this stage. Ligand/PLACER/QMMM filters remain downstream and were not used to reject current-stage sequences.

## Next Work

Do not start QMMM directly from this status update.

Recommended next step:

1. Freeze the capped10 sequence panel as the current Project 01 Phase1 protein-background panel.
2. Start the natural-scaffold MSA track separately: choose the natural enzyme scaffold/family, align homologs, fix conserved catalytic/pocket/core residues, then generate variants.
3. Start the Baker-style de novo track separately: fix active-site/reactive geometry first, generate new backbones around that motif, then sequence and gate them.
4. Only after the accepted protein-background panel and TS/reactive-state definition are explicit should PLACER/TS/QMMM preparation resume.

Remote lightweight evidence:

```text
/data/bht/project01_phase1_reset_gpu/manifests/postseq_entrance_pass_sequence_panel_summary.json
/data/bht/project01_phase1_reset_gpu/manifests/postseq_entrance_pass_sequence_panel.tsv
/data/bht/project01_phase1_reset_gpu/manifests/postseq_entrance_pass_sequence_panel_capped10.tsv
```

Do not upload predicted PDBs, PLACER outputs, trajectories, model weights, or large logs.
