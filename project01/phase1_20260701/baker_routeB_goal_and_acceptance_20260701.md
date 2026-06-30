# Project 01 Baker Route B Goal And Acceptance

Date: 2026-07-01

## Goal

Reset the serine-hydrolase new-backbone route to Baker-style active-site motif scaffolding:

```text
fixed object = TS/substrate/theozyme geometry
generated object = new protein backbone around that motif
current stage = sequence/scaffold generation and motif geometry screening
not current stage = PLACER acceptance, QMMM, barrier calculation
```

This supersedes using a pre-existing ligand-4A shell as the fixed object for new-backbone generation.

## Completion Criteria For This Round

A round is complete only when all items below are true:

1. Baker source inputs are recorded and staged on the GPU working filesystem:
   - `mu1.params`, `1LNS_mu1.cst`, `simple_theozyme.pdb`
   - `bn1.params`, `theozyme.cst`, `theozyme.pdb`
   - `super_af2_bu2.pdb` as substrate-bound / GS-like PLACER reference only
2. A new Route B run directory exists separately from the old `/data/bht/project01_phase1_reset_gpu/denovo_scaffold_routeB/serine_hydrolase` route.
3. The run manifest labels the route as `baker_theozyme_new_backbone`, not `fixed_reference_pocket`.
4. If GPU capacity is available, a small smoke is launched from Baker `theozyme.pdb` or `simple_theozyme.pdb`.
5. If GPU capacity is not available, the status is explicitly `BLOCKED_GPU_BUSY_OR_UNAVAILABLE`; this is not counted as biological failure.

## Pass/Fail Accounting

The evaluated universe for this reset is the set of new backbones generated from Baker theozyme inputs in the new run directory.

Do not count any of these as reset passes:

- fixed-reference sequence variants from the completed capped10 panel;
- old Route B variants derived from `denovo_SER_hydrolase_full_input.pdb`;
- PLACER crop outputs;
- QMMM preparation outputs.

Initial smoke pass condition:

```text
RFAA/RFdiffusionAA output PDB exists
fixed Baker motif residues/ligand are present and mapped
sequence design command can consume the generated backbone
lightweight manifest records seed, input, output, status
```

Production pass condition, after smoke:

```text
per identity bin 90/80/70/60/50:
  at least 10 generated sequence/scaffold candidates
  postseq structure prediction OK
  motif-local geometry gate PASS
```

Thresholds are project-policy thresholds unless directly copied from a tool-native output.

## Grill-Me Checks

- Are we comparing the same denominator across routes? No: fixed-backbone, old Route B, and Baker-theozyme Route B remain separate manifests.
- Did we turn PLACER into the current accept/reject classifier? No: PLACER is downstream for this stage.
- Are GPU/SSH failures biological failures? No: label `NOT_EVALUATED` or `BLOCKED_GPU_BUSY_OR_UNAVAILABLE`.
- Is ligand 4A the fixed object for de novo generation? No: ligand/TS/theozyme reaction geometry is fixed first; generated 4A pocket is evaluated after generation.

## Remote Working Root

Preferred GPU root:

```text
/data/bht/project01_baker_serhyd_routeB_20260701
```

Large generated PDBs, trajectories, model weights, and raw logs stay on the remote filesystem. GitHub receives only Markdown, JSON/TSV manifests, and scripts.
