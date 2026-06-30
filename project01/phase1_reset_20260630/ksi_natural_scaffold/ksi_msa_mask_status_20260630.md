# KSI Natural-Scaffold MSA Mask Status - 2026-06-30

## Status

KSI natural-scaffold generation has now entered the conservation-first route.

Remote GPU root:

```text
/data/bht/project01_phase1_reset_gpu/natural_scaffold/KSI
```

Reference:

```text
UniProt P07445 / RCSB 6UBQ chain A / ligand ASD
```

## Method

1. Downloaded UniProt `ec:5.3.3.1` FASTA records.
2. Filtered to KSI-like sequences by reference length, identity, and coverage against P07445.
3. Ran PyFAMSA on 200 selected homolog/reference-filtered sequences.
4. Mapped alignment columns back to P07445 residue numbering.
5. Marked fixed positions from:
   - catalytic residues: Tyr16, Asp40, Asp103;
   - ASD ligand-contact residues in 6UBQ chain A;
   - MSA-conserved positions with major amino-acid frequency >= 0.80.

## Audited Counts

```text
raw_homolog_records = 7400
selected_homolog_records_from_prefilter = 200
alignment_sequences = 200
alignment_columns = 401
fixed_positions_count = 12
review_positions_count = 19
mutable_positions_count = 100
```

Fixed P07445 positions:

```text
16, 40, 41, 43, 49, 56, 60, 86, 90, 99, 103, 120
```

Ligand-contact positions from 6UBQ chain A / ASD:

```text
16, 40, 56, 60, 86, 90, 99, 103, 120
```

## Synced Lightweight Files

```text
project01/phase1_reset_20260630/ksi_natural_scaffold/natural_scaffold_msa_pyfamsa_summary.json
project01/phase1_reset_20260630/ksi_natural_scaffold/natural_scaffold_msa_pyfamsa_fixed_positions.tsv
project01/phase1_reset_20260630/ksi_natural_scaffold/natural_scaffold_msa_pyfamsa_mutable_positions.tsv
project01/phase1_reset_20260630/ksi_natural_scaffold/natural_scaffold_msa_pyfamsa_review_positions.tsv
project01/phase1_reset_20260630/ksi_natural_scaffold/natural_scaffold_msa_pyfamsa_conservation_by_position.tsv
```

Not synced to GitHub:

```text
raw UniProt FASTA
aligned FASTA
PDB files
model outputs
large logs
```

## Next Action

Use `natural_scaffold_msa_pyfamsa_fixed_positions.tsv` as the fixed-position mask for KSI sequence generation. Generate 90/80/70/60/50 sequence-identity bins only over the mutable positions, then run structure prediction and the postseq entrance gate before accepting any sequence.
