# Route B Active-Pocket-Only Sequence Target - 2026-06-30

## Objective

Generate active-pocket-preserving de novo scaffold sequences for the serine-hydrolase Route B track.

Final stage target:

```text
identity_bins = 90, 80, 70, 60, 50
accepted_sequences_per_bin = 10
accepted_total = 50
```

This target is sequence-generation only. PLACER and QMMM are downstream and are not acceptance gates for this stage.

## Identity Definition

For true de novo/new-backbone Route B, identity is measured within a generated-backbone family, not against the original project-local reference scaffold.

Current implementation:

```text
identity_reference = first accepted motif-backbone sequence for the same RFAA backbone family
identity_positions = non-pocket/background positions plus fixed motif positions counted as unchanged
fixed_positions = active_site + pocket_4A residues
mutable_positions = all nonfixed positions
```

If an identity bin cannot be reached without mutating fixed active-site/pocket residues, that row is `NOT_GENERATED_FIXED_POCKET_CONFLICT`, not a biological failure.

## Baker-Style Pocket Implementation

The current project implementation follows the Baker-style logic of building/designing around a fixed catalytic/reactive pocket rather than preserving one global fold.

Operational layer definitions:

```text
active_site = catalytic residues and essential reaction-geometry residues
pocket_4A = residues with any heavy atom <= 4 A from ligand/reactive fragment
shell_diagnostic = residues with any heavy atom <= 6-8 A from ligand/reactive fragment
background = all remaining designable residues
```

For the current serine-hydrolase Route B seed:

```text
reactive ligand / fragment = bu2
input motif source = /data/bht/project01_phase1_reset_gpu/manifests/denovo_SER_hydrolase_full_input.pdb
input catalytic motif = A95 HIS; A102 GLU; A128 SER
current productive RFAA family = sample_2
RFAA-renumbered motif for sample_2 = A40 HIS; A49 GLU; A56 SER
```

The `sample_2` family is the first productive seed because it produced two ESMFold-predicted sequence variants with motif-backbone CA RMSD <= 1.0 A.

## Acceptance Gates

Sequence-stage hard gate:

```text
routeB_motif_backbone_gate = PASS if motif_ca_rmsd_A <= 1.0
```

Recommended pocket gate for this target:

```text
active_site residues unchanged = required
pocket_4A residues unchanged or explicitly fixed = required
motif_ca_rmsd_A <= 1.0 = required
```

Diagnostic fields, not hard gates at this stage:

```text
global_ca_rmsd_A = diagnostic only
motif_heavy_rmsd_A = diagnostic until ligand-aware repack/holo validation
max_key_distance_delta_A = diagnostic until ligand-aware repack/holo validation
```

Reason: ESMFold predicts apo structures and does not see ligand `bu2`, so side-chain rotamers and ligand distances cannot be used as the only sequence-stage hard rejection rule.

## Current Starting Point

Current evaluated Route B universe:

```text
RFAA new backbones evaluated = 6 total, including initial smoke sample_0 plus batch5 sample_1..sample_5
LigandMPNN/ESMFold sequence predictions evaluated = 47 total
sequence-stage motif-backbone PASS = 2
strict all-atom sidechain/ligand diagnostic PASS = 0
```

Current accepted motif-backbone candidates:

| sample_id | RFAA family | fixed motif | motif CA RMSD A | status |
| --- | --- | --- | ---: | --- |
| routeB_serhyd_sample_2_lmpnn_03 | sample_2 | A40/A49/A56 | 0.4081 | seed_candidate |
| routeB_serhyd_sample_2_lmpnn_04 | sample_2 | A40/A49/A56 | 0.7402 | seed_candidate |

## Next Executable Step

Use RFAA `sample_2` as the first productive backbone family and generate controlled identity bins:

```text
1. choose routeB_serhyd_sample_2_lmpnn_03 as the provisional family reference sequence;
2. define fixed positions from active_site + pocket_4A around bu2 in RFAA sample_2;
3. generate 90/80/70/60/50 sequence identity bins over nonfixed background positions;
4. run ESMFold with memory_fraction 0.20;
5. accept only rows with fixed pocket conserved and motif_ca_rmsd_A <= 1.0;
6. continue closed-loop generation until each bin has 10 accepted rows.
```

Do not upload PDB, trajectory, model, or raw log files. Sync only lightweight status markdown, TSV manifests, scripts, and summary JSON.
