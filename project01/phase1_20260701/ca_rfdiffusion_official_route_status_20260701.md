# Project 01 CA_RFDiffusion Official-Route Status

Date: 2026-07-01 CST

## Current state

The corrected Baker-style new-backbone route has been moved from the earlier RFAA smoke to the official CA_RFDiffusion code path.

Remote working root:

`/data/bht/project01_baker_serhyd_routeB_20260701`

CA_RFDiffusion source:

`/data/bht/design_tools/src/CA_RFDiffusion`

CA_RFDiffusion commit:

`11a3dfb08235c97985bc81accbc950b8848c1b06`

Installed environment:

`/data/bht/design_tools/envs/ca_rfd_release`

Import smoke passed:

`torch 2.2.0 cuda True`

## Weights

Downloaded public CA_RFDiffusion weights:

- `/data/bht/design_tools/models/ca_rf_diffusion/ca_rfd_diffusion.pt`
- `/data/bht/design_tools/models/ca_rf_diffusion/ca_rfd_refinement.pt`

Checkpoint caveat:

Baker `design_pipeline/01_diffusion/inputs/tasks.json` references `BFF_7.pt`.
That checkpoint was not present on the server after setup. The current n=1 smoke uses the public CA_RFDiffusion `ca_rfd_diffusion.pt`, and is therefore an official-parameter route smoke, not a fully checkpoint-identical Baker reproduction.

## Active CA diffusion run

Run label:

`ca_rfd_baker_theozyme_diffusion_n1_publicckpt`

Remote PID file:

`/data/bht/project01_baker_serhyd_routeB_20260701/logs/ca_rfd_baker_theozyme_diffusion_n1_publicckpt.pid`

Observed PID:

`4116728`

Remote log:

`/data/bht/project01_baker_serhyd_routeB_20260701/logs/ca_rfd_baker_theozyme_diffusion_n1_publicckpt.log`

Output prefix:

`/data/bht/project01_baker_serhyd_routeB_20260701/outputs/ca_rfd_baker_theozyme_diffusion_n1_publicckpt/sample`

Parameters:

- `--config-name=RFdiffusion_CA_inference`
- `inference.num_designs=1`
- `inference.design_startnum=0`
- `inference.input_pdb=.../design_pipeline/01_diffusion/inputs/theozyme.pdb`
- `inference.ligand=bn1`
- `inference.ij_visible=abcde`
- `inference.cautious=true`
- `inference.write_trajectory=true`
- `diffuser.T=50`
- `denoiser.noise_scale_frame=0.05`
- `denoiser.noise_scale_ca=0.0`
- `motif_only_2d=true`
- `preprocess.eye_frames=true`
- `contigmap.contigs=['12,A56-60,36,A83-85,15,A113-115,73,B145-147,10']`

GPU sharing:

This run was launched only after explicit user authorization for small shared-GPU execution. It is limited by scope (`n=1`) rather than a guaranteed hardware 10% GPU cap. `CUDA_MPS_ACTIVE_THREAD_PERCENTAGE=10` was set as a best-effort environment variable, but it only has effect if MPS is active.

## Evidence observed before commit

The second CA run passed import and configuration, loaded the public checkpoint, created first-time IGSO3 cache, and entered denoising:

- `Making design ... sample_0`
- `ckpt_path: /data/bht/design_tools/models/ca_rf_diffusion/ca_rfd_diffusion.pt`
- `Computing IGSO3`
- `ema: True`
- `making design 0 of 0:1`
- Denoising reached at least `t=tensor(44)` by the last check.

Approximate speed under shared GPU load:

About 43 seconds per CA denoising step during the first observed steps. Expected completion for n=1 is roughly half-hour scale after cache/model initialization.

## Acceptance scope

This active run can only establish:

- CA_RFDiffusion official route executes on the staged Baker theozyme input.
- The run writes a CA diffusion PDB/TRB/traj output.
- The output can then be passed into the CA_RFDiffusion refinement step.

It must not be reported as:

- A full Baker checkpoint-identical reproduction, unless `BFF_7.pt` is obtained.
- A sequence-similarity-bin completion.
- A PLACER/QMMM-ready structure.
- A motif geometry pass before post-run geometric screening.

## Next action

1. Wait for `PID 4116728` to finish.
2. Generate a lightweight completion status JSON from output paths, file sizes, log summary, and PID state.
3. If diffusion completed, run CA_RFDiffusion refinement using `ca_rfd_refinement.pt` on `sample_0`.
4. After refinement, run motif/theozyme geometry screening before any larger sequence/scaffold campaign.

Heartbeat monitor:

`project-01-ca-rfdiffusion-n1-monitor`
