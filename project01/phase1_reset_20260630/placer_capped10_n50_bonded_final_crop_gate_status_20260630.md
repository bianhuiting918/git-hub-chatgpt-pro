# capped10 PLACER n50 Final Crop Gate Status - 2026-06-30

## Scope

This file records the completed downstream PLACER diagnostic for the existing serine-hydrolase capped10 panel. It is not a sequence-generation gate for the current natural-scaffold or de novo scaffold sequence-design stage.

Remote root:

```text
/data/bht/project01_phase1_reset_gpu
```

Run label:

```text
capped10_n50_bonded
```

## Driver Completion Evidence

Checked at `2026-06-30 15:40:55 CST`:

```text
CSV outputs = 50
model PDB outputs = 50
driver error count = 0
last driver event = DONE 50 bin60_r05_subset72_rank18_orig18.holo_bu2_refpose
```

## Final Crop Gate Evidence

Command:

```text
scripts/run_crop_gate_capped10_final.sh
```

Exit code:

```text
0
```

Evaluated universe:

```text
total_sequences = 50
total_conformers = 2500
NOT_EVALUATED = 0
```

Project strict crop filter:

```text
FAIL_BOTH_PROTEIN_AND_LIGAND_GATES = 2500 conformers
CROP_STRICT_PASS = 0 conformers
DISCARD_SEQUENCE_REGENERATE_FROM_UPSTREAM = 50 sequences
```

Thresholds:

```text
global_backbone_rmsd_max_A = 2.5
fixed_backbone_rmsd_max_A = 1.0
catalytic_heavy_rmsd_max_A = 0.75
protein_key_distance_abs_delta_max_A = 0.75
ligand_heavy_rmsd_max_A = 0.75
ser128_og_to_bu2_c1_abs_delta_max_A = 0.5
expected_conformers_per_sequence = 50
minimum_crop_pass_conformers_per_sequence = 10
```

## Interpretation

PLACER itself produced 50 per-input output sets with no driver errors. The downstream project strict crop filter rejected every conformer because ligand/reaction-pose geometry did not pass. Therefore this capped10 PLACER batch must not proceed to full ligand/protein completion or strict QMMM manifest generation.

This result does not reject current sequence-generation work. The active sequence-generation stage continues through natural-scaffold and de novo-scaffold routes using postseq pocket/catalytic geometry gates before any later PLACER/QMMM stage.
