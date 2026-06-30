# Public-Checkpoint CA_RFDiffusion Batch Next Step

Date: 2026-07-01 CST

## Current state

The Baker-style CA_RFDiffusion n=1 public-checkpoint smoke completed diffusion, refinement, and exploratory motif CA screening.

The exact Baker checkpoints `BFF_7.pt` and `refine_BFF_3.pt` remain unavailable from a verified source. Until those are obtained, production-scale work must be labeled:

`public CA_RFDiffusion checkpoint, Baker-parameter exploratory route`

## GPU status at last check

Host:

`dell-PowerEdge-R760`

Observed at:

`2026-07-01 04:54 CST`

GPU:

`NVIDIA A100 80GB PCIe`

Status:

- GPU utilization: `100%`
- memory used: `21390 MiB / 81920 MiB`
- active compute jobs were from other users' `pxdesign` runs

Decision:

Do not launch another CA_RFDiffusion batch while the shared GPU is fully occupied by other users. This is `BLOCKED_GPU_BUSY_OR_UNAVAILABLE`, not a biological or model failure.

## Prepared script

New lightweight script:

`project01/phase1_20260701/scripts/launch_ca_rfdiffusion_public_batch.sh`

Remote intended location:

`/data/bht/project01_baker_serhyd_routeB_20260701/scripts/launch_ca_rfdiffusion_public_batch.sh`

Default batch settings:

- `NUM_DESIGNS=5`
- `DESIGN_STARTNUM=1`
- `diffuser.T=50`
- `denoiser.noise_scale_frame=0.05`
- `denoiser.noise_scale_ca=0.0`
- `motif_only_2d=true`
- `preprocess.eye_frames=true`
- `inference.ligand=bn1`
- `inference.ij_visible=abcde`
- contig: `12,A56-60,36,A83-85,15,A113-115,73,B145-147,10`
- checkpoint: `/data/bht/design_tools/models/ca_rf_diffusion/ca_rfd_diffusion.pt`

The script refuses to launch by default if GPU compute processes are present or GPU utilization is above `MAX_GPU_UTIL=40`.

## Run command when GPU is available

```bash
cd /data/bht/project01_baker_serhyd_routeB_20260701
chmod +x scripts/launch_ca_rfdiffusion_public_batch.sh
NUM_DESIGNS=5 DESIGN_STARTNUM=1 MAX_GPU_UTIL=40 scripts/launch_ca_rfdiffusion_public_batch.sh
```

If a user explicitly authorizes shared-GPU execution:

```bash
ALLOW_SHARED_GPU=1 NUM_DESIGNS=5 DESIGN_STARTNUM=1 scripts/launch_ca_rfdiffusion_public_batch.sh
```

Do not use `FORCE=1` unless the shared-GPU risk has been explicitly accepted.

## Acceptance for this next batch

Evaluated universe:

Only CA_RFDiffusion outputs from this public-checkpoint batch.

Tool-native completion:

- output PDB/TRB exists for each generated design
- CA_RFDiffusion log reaches `Finished design`

Project exploratory motif filter:

- task-reference motif CA residues mapped
- ligand `bn1` present
- motif CA Kabsch RMSD <= `1.0 A`
- max motif pair-distance delta <= `1.0 A`

This motif filter is a project-policy exploratory gate, not a native CA_RFDiffusion pass/fail label.

## Next workflow after batch finishes

1. Run motif CA screen on each generated scaffold using Baker task-reference residue mapping.
2. Keep only motif-gate passing scaffolds for sequence design.
3. Start sequence design and bin candidates into 90/80/70/60/50 similarity levels.
4. For each level, continue generating until there are 10 accepted sequences after post-sequence structure prediction and motif-local geometry screening.

PLACER and QMMM remain downstream and are not required for this sequence-generation phase.
