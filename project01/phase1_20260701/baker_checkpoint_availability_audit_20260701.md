# Baker CA_RFDiffusion Checkpoint Availability Audit

Date: 2026-07-01 CST

## Question

Why did we not download Baker's `BFF_7.pt` checkpoint directly?

## Short answer

`BFF_7.pt` and `refine_BFF_3.pt` are referenced by the Baker serine-hydrolase task files, but I did not find a public, verifiable download source for those exact files. I therefore did not download an arbitrary same-name checkpoint from an unverified location.

The current CA_RFDiffusion smoke used the public CA_RFDiffusion release weights:

- `/data/bht/design_tools/models/ca_rf_diffusion/ca_rfd_diffusion.pt`
- `/data/bht/design_tools/models/ca_rf_diffusion/ca_rfd_refinement.pt`

This means the completed n=1 run is an official-parameter, public-checkpoint route validation, not a checkpoint-identical Baker reproduction.

## Evidence checked

Remote GPU host:

`dell-PowerEdge-R760`

Remote project root:

`/data/bht/project01_baker_serhyd_routeB_20260701`

Serine hydrolase repository staged on GPU:

`/data/bht/project01_baker_serhyd_routeB_20260701/external/serine-hydrolase-design`

Repository remote and commit:

- `https://github.com/laukoag/serine-hydrolase-design.git`
- `2442801d2c4a37994197faa8f20b459381e5dbd6`

CA_RFDiffusion source:

`/data/bht/design_tools/src/CA_RFDiffusion`

Repository remote and commit:

- `https://github.com/baker-laboratory/CA_RFDiffusion.git`
- `11a3dfb08235c97985bc81accbc950b8848c1b06`

Observed Baker task references:

- `design_pipeline/01_diffusion/inputs/tasks.json` references `../../software/ca_rf_diffusion/checkpoints/BFF_7.pt`
- `design_pipeline/02_refinement/inputs/tasks.json` references `../../software/ca_rf_diffusion/checkpoints/refine_BFF_3.pt`
- `motif_gen/02_diffusion/inputs/tasks.json` references `BFF_7.pt`
- `motif_gen/03_refinement/inputs/tasks.json` references `refine_BFF_3.pt`

Observed historical Hydra output inside the staged repository:

- `design_pipeline/01_diffusion/outputs/2024-09-04/16-59-59/run_inference.log` says it read `../../software/ca_rf_diffusion/checkpoints/BFF_7.pt`
- `.hydra/config.yaml` in that output references an internal path ending in `BFF_10.pt`

Server file search:

- strict search for `BFF_7.pt`, `BFF*.pt`, `*BFF*.ckpt`, and `*BFF*.pth` under `/data/bht` and `/data` returned no usable checkpoint file
- only permission noise observed: `/data/lost+found`

Git LFS check:

- `git lfs ls-files` in the staged `serine-hydrolase-design` checkout did not report checkpoint files

CA_RFDiffusion public README evidence:

- the installed CA_RFDiffusion README points users to public `ca_rfd_diffusion.pt`
- the inference config comments say users must point `inference.ckpt_path` to their copy of `ca_rfd_diffusion.pt`
- no README reference to public `BFF_7.pt` or `refine_BFF_3.pt` was found in the local CA_RFDiffusion checkout

## Decision

Do not download `BFF_7.pt` from an unverified source. Treat exact Baker checkpoint reproduction as blocked until one of these becomes available:

1. an author-provided release asset,
2. an official supplementary-data or Zenodo file entry,
3. a verified file already present on an authorized server path,
4. direct author-provided checkpoint access.

## Current valid route

Continue the project with the public CA_RFDiffusion route as an exploratory Baker-parameter implementation:

1. Generate CA_RFDiffusion scaffolds from the Baker theozyme/TS input using public `ca_rfd_diffusion.pt`.
2. Refine using public `ca_rfd_refinement.pt`.
3. Screen motif/theozyme geometry using task-reference motif mapping, not raw TRB numbering alone.
4. Only after motif geometry passes, enter sequence design and 90/80/70/60/50 sequence-similarity binning.

The n=1 public-checkpoint scaffold already passed the exploratory motif CA screen:

- mapped CA pairs: `14/14`
- motif CA Kabsch RMSD: `0.0976 A`
- max motif pair-distance delta: `0.1630 A`

This is a project-policy exploratory motif gate, not a tool-native pass and not final sequence-bin completion.

## Next action

Resume the main target: generate enough scaffold/sequence candidates to obtain 10 accepted sequences per similarity bin, while clearly labeling that the route uses public CA_RFDiffusion checkpoints unless `BFF_7.pt` and `refine_BFF_3.pt` are later obtained.
