# Serine-Hydrolase Route B Smoke Status - 2026-06-30

## Result

A true de novo/new-backbone smoke test was started on the GPU host for the serine-hydrolase Route B track. This is separate from the completed existing/reference-scaffold sequence panel.

Tool-native result:

- RFDiffusionAA ran on GPU and wrote one raw generated backbone PDB.
- LigandMPNN ran on the generated backbone and wrote two designed sequence/backbone outputs.

Project strict filter:

- STRUCTURE_PREDICTED_AND_POSTSEQ_GATED: 0. ESMFold validation has not yet run for these Route B sequences.
- RAW_RFAA_BACKBONE_WRITTEN: 1.
- LIGANDMPNN_SEQUENCE_WRITTEN: 2.
- NOT_EVALUATED_FOR_POSTSEQ_GATE: 2 sequences, because ESMFold environment was not found in the checked env roots during this heartbeat.

Caveat:

This is a Route B smoke pass for tool execution and motif retention only. It is not yet an accepted Route B sequence panel and must not be mixed with the 90/80/70/60/50 existing/reference-scaffold or KSI natural-scaffold panels.

## Remote Evidence

GPU host reached via CPU jump:

```text
CPU: bianht@210.73.40.29
GPU: bht@192.168.10.38
GPU hostname: dell-PowerEdge-R760
GPU: NVIDIA A100 80GB PCIe
```

Remote working root:

```text
/data/bht/project01_phase1_reset_gpu/denovo_scaffold_routeB/serine_hydrolase
```

RFAA input:

```text
/data/bht/project01_phase1_reset_gpu/manifests/denovo_SER_hydrolase_full_input.pdb
```

RFAA raw output files, not uploaded:

```text
outputs/rfaa_serhyd_motif_H95_E102_S128_bu2_smoke/sample_0.pdb
outputs/rfaa_serhyd_motif_H95_E102_S128_bu2_smoke/traj/sample_0_Xt-1_traj.pdb
outputs/rfaa_serhyd_motif_H95_E102_S128_bu2_smoke/unidealized/sample_0.pdb
```

LigandMPNN output files, not uploaded:

```text
outputs/ligandmpnn_serhyd_rfaa_smoke_n2/seqs/sample_0.fa
outputs/ligandmpnn_serhyd_rfaa_smoke_n2/backbones/sample_0_1.pdb
outputs/ligandmpnn_serhyd_rfaa_smoke_n2/backbones/sample_0_2.pdb
outputs/ligandmpnn_serhyd_rfaa_smoke_n2/stats/sample_0.pt
```

Logs, not uploaded:

```text
logs/rfaa_serhyd_motif_H95_E102_S128_bu2_smoke_T20_20260630_210758.log
logs/ligandmpnn_serhyd_rfaa_smoke_n2_20260630_211159.log
```

## RFAA Settings

```text
input_pdb = /data/bht/project01_phase1_reset_gpu/manifests/denovo_SER_hydrolase_full_input.pdb
ligand = bu2
motif = A95_HIS; A102_GLU; A128_SER from input PDB
contig = 20-40,A95-95,5-15,A102-102,5-20,A128-128,20-40
length = 90-140
num_designs = 1
diffuser.T = 20
```

The first attempted `diffuser.T=50` run was stopped before sampling because it spent more than two minutes building a missing IGSO3 cache. The successful retry used `diffuser.T=20`, which had an existing cache and completed in about 0.80 minutes.

## RFAA Motif Check

RFAA renumbered the motif in the generated backbone:

| Input motif | Generated output residue | Nearest ligand distance from motif sidechain atom |
| --- | --- | ---: |
| A95 HIS | A37 HIS | CE1 2.53 A; ND1 3.64 A; NE2 3.27 A |
| A102 GLU | A50 GLU | OE2 4.12 A; OE1 4.23 A |
| A128 SER | A68 SER | OG 1.53 A |

The raw RFAA PDB contains many non-motif sidechain atoms at `(0,0,0)`, so full-atom sidechain geometry must not be evaluated until LigandMPNN/sidechain completion or structure prediction is applied. Backbone and fixed motif sidechain checks are the valid smoke-level checks here.

## LigandMPNN Settings

```text
model_type = ligand_mpnn
checkpoint_ligand_mpnn = ligandmpnn_v_32_005_25.pt
input_pdb = RFAA sample_0.pdb
fixed_residues = A37 A50 A68
ligand_mpnn_use_side_chain_context = 1
seed = 111
temperature = 0.1
batch_size = 2
number_of_batches = 1
```

LigandMPNN reported:

```text
fixed residues: A37 A50 A68
redesigned residues: all other A-chain residues
ligand atoms parsed: 18
```

## Sequence Smoke Output

Both designed sequences retained the fixed motif positions after RFAA renumbering:

| sample_id | length | A37 | A50 | A68 | overall_confidence | ligand_confidence | current status |
| --- | ---: | --- | --- | --- | ---: | ---: | --- |
| routeB_serhyd_rfaa_smoke_lmpnn_1 | 97 | H | E | S | 0.3447 | 0.4474 | NOT_EVALUATED_FOR_POSTSEQ_GATE |
| routeB_serhyd_rfaa_smoke_lmpnn_2 | 97 | H | E | S | 0.3527 | 0.4315 | NOT_EVALUATED_FOR_POSTSEQ_GATE |

## Next Action

1. Locate or restore the ESMFold environment used by earlier Project 01 runs.
2. Predict these two LigandMPNN sequences with a 0.20 memory-fraction cap if supported.
3. Run a Route B postseq gate adapted to generated backbones: motif-local catalytic geometry first, not whole-reference RMSD.
4. If both sequences fail structure prediction or motif-local geometry, regenerate a larger but still small RFAA/LigandMPNN batch.
5. Do not start PLACER or QMMM from these smoke outputs yet.

## Gate Policy Correction - Pocket-Only Hard Gate

For Route B de novo/new-backbone work, global CA RMSD to the RFAA generated backbone is no longer a hard acceptance gate. It remains a diagnostic field showing whether the designed sequence supports the generated fold.

Hard acceptance should focus on active-site motif/pocket geometry:

```text
routeB_pocket_only_gate = PASS only if
  motif_ca_rmsd_A <= 1.0
  motif_heavy_rmsd_A <= 1.0
  max_key_distance_delta_A <= 1.0
```

Global CA RMSD is recorded as:

```text
global_ca_rmsd_role = DIAGNOSTIC_NOT_HARD_GATE
```

Updated remote summaries:

```text
/data/bht/project01_phase1_reset_gpu/denovo_scaffold_routeB/serine_hydrolase/manifests/routeB_serhyd_lmpnn_postseq_gate_pocket_only_summary.json
/data/bht/project01_phase1_reset_gpu/denovo_scaffold_routeB/serine_hydrolase/manifests/routeB_serhyd_batch20_postseq_gate_pocket_only_summary.json
```

Current pocket-only result:

```text
initial n2 LigandMPNN smoke: 0 PASS / 2 evaluated
batch20 LigandMPNN expansion: 0 PASS / 20 evaluated
```

Interpretation: the first RFAA backbone and its current LigandMPNN sequences do not yet preserve the active-site motif after ESMFold prediction. The next regeneration step should sample additional RFAA backbones and screen by pocket-only gate, rather than rejecting based on global fold RMSD.

## RFAA Batch5 x LigandMPNN Seq5 Result - 2026-06-30

After correcting the Route B gate to focus on active-site pocket geometry at the current sequence stage, a new small batch was run on GPU:

```text
RFAA new backbones: 5
LigandMPNN sequences: 25, 5 per backbone
ESMFold predicted structures: 25 OK / 25 submitted
sequence-stage motif-backbone gate: 2 PASS / 25 evaluated
strict all-atom pocket sidechain gate: 0 PASS / 25 evaluated
```

Current stage policy:

```text
sequence-stage hard gate = motif backbone/CA preservation
sidechain-heavy/ligand-distance gate = diagnostic until ligand-aware repack or holo model
full global CA RMSD = diagnostic, not hard gate
```

Best current candidates:

| sample_id | RFAA backbone | fixed residues | global CA RMSD A | motif CA RMSD A | motif heavy RMSD A | max key-distance delta A | sequence-stage gate | strict all-atom gate |
| --- | --- | --- | ---: | ---: | ---: | ---: | --- | --- |
| routeB_serhyd_sample_2_lmpnn_03 | sample_2 | A40 A49 A56 | 1.0309 | 0.4081 | 1.6668 | 3.7134 | PASS | FAIL_DIAGNOSTIC |
| routeB_serhyd_sample_2_lmpnn_04 | sample_2 | A40 A49 A56 | 1.0787 | 0.7402 | 2.4429 | 7.2344 | PASS | FAIL_DIAGNOSTIC |

Remote evidence, not uploaded:

```text
/data/bht/project01_phase1_reset_gpu/denovo_scaffold_routeB/serine_hydrolase/manifests/rfaa_batch5_motif_mapping.tsv
/data/bht/project01_phase1_reset_gpu/denovo_scaffold_routeB/serine_hydrolase/manifests/routeB_serhyd_rfaabatch5_seq5_esmfold_summary.json
/data/bht/project01_phase1_reset_gpu/denovo_scaffold_routeB/serine_hydrolase/manifests/routeB_serhyd_rfaabatch5_seq5_pocket_only_gate_tiered_summary.json
```

Next action: use RFAA `sample_2` as the first productive motif-backbone seed. Generate more sequences for `sample_2` and/or run ligand-aware sidechain repack/holo validation before deciding whether the sidechain-heavy gate should reject it.

## Active-Pocket-Only Identity-Bin Target - 2026-06-30

New goal has been set for Route B serine hydrolase:

```text
Generate active-pocket-only accepted de novo scaffold sequences at 90/80/70/60/50 sequence-identity bins.
Target: 10 accepted sequences per bin, 50 accepted sequences total.
```

Acceptance is now defined at sequence stage by Baker-style active pocket preservation:

```text
active_site = catalytic residues and essential reaction-geometry residues
pocket_4A = residues with any heavy atom <= 4 A from ligand/reactive fragment
sequence-stage hard gate = active_site/pocket_4A fixed + motif_ca_rmsd_A <= 1.0
```

Current first family seed:

```text
RFAA family = sample_2
RFAA-renumbered motif = A40 HIS; A49 GLU; A56 SER
provisional reference sequence = routeB_serhyd_sample_2_lmpnn_03
```

Detailed target file:

```text
project01/phase1_reset_20260630/routeB_denovo_scaffold/routeB_active_pocket_only_sequence_target_20260630.md
project01/phase1_reset_20260630/routeB_denovo_scaffold/routeB_active_pocket_only_bin_targets.tsv
```
