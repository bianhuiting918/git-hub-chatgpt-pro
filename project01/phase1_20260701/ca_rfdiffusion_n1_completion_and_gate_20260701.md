# CA_RFDiffusion n=1 Completion and Motif Gate

Date: 2026-07-01 CST

## Run Summary

The official CA_RFDiffusion route was run on the staged Baker serine-hydrolase theozyme input.

Remote root:

`/data/bht/project01_baker_serhyd_routeB_20260701`

Input:

`/data/bht/project01_baker_serhyd_routeB_20260701/inputs/baker_serhyd/design_pipeline/01_diffusion/inputs/theozyme.pdb`

Checkpoint caveat:

Baker `tasks.json` references `BFF_7.pt`, which was not present on the server. This run used the public CA_RFDiffusion checkpoint:

`/data/bht/design_tools/models/ca_rf_diffusion/ca_rfd_diffusion.pt`

This should be treated as an official-parameter public-checkpoint smoke, not a checkpoint-identical Baker reproduction.

## Diffusion

Run label:

`ca_rfd_baker_theozyme_diffusion_n1_publicckpt`

Status:

`COMPLETED_DIFFUSION_OUTPUT_WRITTEN`

Log evidence:

`Finished design in 46.57 minutes`

Outputs written remotely:

- `outputs/ca_rfd_baker_theozyme_diffusion_n1_publicckpt/sample_0.pdb` - 56,690 bytes
- `outputs/ca_rfd_baker_theozyme_diffusion_n1_publicckpt/sample_0.trb` - 27,006 bytes
- `outputs/ca_rfd_baker_theozyme_diffusion_n1_publicckpt/unidealized/sample_0.pdb` - 56,712 bytes
- `outputs/ca_rfd_baker_theozyme_diffusion_n1_publicckpt/traj/sample_0_Xt-1_traj.pdb` - 2,835,640 bytes
- `outputs/ca_rfd_baker_theozyme_diffusion_n1_publicckpt/traj/sample_0_pX0_traj.pdb` - 2,835,640 bytes

Remote lightweight status:

`/data/bht/project01_baker_serhyd_routeB_20260701/manifests/ca_rfd_baker_theozyme_diffusion_n1_publicckpt_completion_status.json`

## Refinement

Refinement checkpoint:

`/data/bht/design_tools/models/ca_rf_diffusion/ca_rfd_refinement.pt`

Status:

`COMPLETED_REFINEMENT_OUTPUT_WRITTEN`

Log evidence:

`Finished design in 12.73 minutes`

Important output behavior:

The refinement config wrote the refined outputs beside the diffusion output, not under the separate refinement output directory supplied on the command line.

Refined outputs written remotely:

- `outputs/ca_rfd_baker_theozyme_diffusion_n1_publicckpt/sample_0_refined_0.pdb` - 56,690 bytes
- `outputs/ca_rfd_baker_theozyme_diffusion_n1_publicckpt/sample_0_refined_0.trb` - 26,394 bytes
- `outputs/ca_rfd_baker_theozyme_diffusion_n1_publicckpt/unidealized/sample_0_refined_0.pdb` - 56,712 bytes
- `outputs/ca_rfd_baker_theozyme_diffusion_n1_publicckpt/traj/sample_0_refined_0_Xt-1_traj.pdb` - 56,712 bytes
- `outputs/ca_rfd_baker_theozyme_diffusion_n1_publicckpt/traj/sample_0_refined_0_pX0_traj.pdb` - 56,712 bytes

Remote lightweight status:

`/data/bht/project01_baker_serhyd_routeB_20260701/manifests/ca_rfd_baker_theozyme_refine_n1_publicckpt_completion_status.json`

## Motif CA Gate

The first raw `.trb` mapping screen was misleading because CA_RFDiffusion internally renumbered the motif into a single chain. The valid exploratory screen uses Baker `tasks.json` original motif residues as the reference and `.trb` hal residues as the generated scaffold mapping.

Reference motif residues:

`A56-A60`, `A83-A85`, `A113-A115`, `B145-B147`

Generated mapped residues:

`A13-A17`, `A54-A56`, `A72-A74`, `A148-A150`

Project exploratory gate:

`PASS`

Metrics:

- mapped CA pairs: `14 / 14`
- motif CA Kabsch RMSD: `0.0976 A`
- max motif pair-distance delta: `0.1630 A`
- mean motif pair-distance delta: `0.0320 A`
- reference ligand HETATM count: `36`
- refined ligand HETATM count: `22`
- ligand resname retained: `bn1`

Remote gate JSON:

`/data/bht/project01_baker_serhyd_routeB_20260701/manifests/ca_rfd_baker_theozyme_refined_motif_ca_screen_tasksref.json`

## Interpretation

This proves the corrected Baker-style CA_RFDiffusion route can run end-to-end on the server for one public-checkpoint theozyme scaffold, including refinement and motif CA preservation.

This is not yet:

- a full Baker checkpoint-identical reproduction;
- a completed 90/80/70/60/50 sequence-similarity panel;
- a PLACER or QMMM-ready result;
- a full all-atom catalytic geometry pass.

Next step:

Use this validated route to launch a small batch only after deciding whether to continue with the public CA checkpoint or obtain the Baker `BFF_7.pt` checkpoint.
