# Data files

This directory should contain lightweight seed and toy evidence tables only.

Expected examples:

```text
seed_enzymes.example.tsv
catalytic_residues.example.tsv
toy_sequence_hits.tsv
toy_foldseek_hits.tsv
toy_folddisco_after_recall_hits.tsv
toy_folddisco_denovo_hits.tsv
toy_pocket_hits.tsv
```

Large sequence databases, structure databases, docking ensembles, trajectories, model checkpoints, and raw production outputs should not be committed.

For catalytic residue annotations, include a `motif_layer` column:

```text
L0_catalytic_core
L1_catalytic_environment
L2_substrate_pocket
L3_processing_interface
```

Pocket/environment residues such as T227-like positions should generally be L2/L3 unless a target-specific annotation explicitly marks them as catalytic core.
