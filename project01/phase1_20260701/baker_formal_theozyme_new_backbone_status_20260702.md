# Project 01 Baker Formal Theozyme New-Backbone Route

Date: 2026-07-02

## Route Decision

The main route is switched from fixed-parent 4 A pocket sequence scanning to the Baker-style new-backbone route.

Formal Baker-style constraint used here:

- input motif/theozyme: `/data/bht/project01_baker_serhyd_routeB_20260701/inputs/baker_serhyd/design_pipeline/01_diffusion/inputs/theozyme.pdb`
- ligand: `bn1`
- contig: `12,A56-60,36,A83-85,15,A113-115,73,B145-147,10`
- CA_RFDiffusion flags include `motif_only_2d=true`, `preprocess.eye_frames=true`, `denoiser.noise_scale_frame=0.05`, `denoiser.noise_scale_ca=0.0`, `diffuser.T=50`

This is the Baker theozyme/motif/stub constraint route. It is not the fixed existing-parent 4 A pocket route.

## Checkpoint Caveat

The server does not have the Baker example checkpoint path `../../software/ca_rf_diffusion/checkpoints/BFF_7.pt` under the external repo copy.

Current run uses the available public CA_RFDiffusion checkpoint:

`/data/bht/design_tools/models/ca_rf_diffusion/ca_rfd_diffusion.pt`

Therefore the constraints follow the Baker example, but the checkpoint is the public available checkpoint, not the unavailable Baker `BFF_7.pt` file.

## Fixed-Parent Route Status

The fixed-parent 4 A large ESMFold batch was stopped after the route decision changed.

Last fixed-parent evidence before stopping:

- round01 ESMFold: 200/200 OK
- round01 active gate: 90% had 3 `FINAL_QUALIFIED_ACTIVE`; 80/70/60/50 had 0
- large batch01 partial: 24 evaluated 90% candidates, 0 additional `FINAL_QUALIFIED_ACTIVE`

These are retained as exploratory evidence and are not the main route after this decision.

## New-Backbone Smoke Run

Run ID:

`ca_rfd_baker_theozyme_formal_constraints_n1_20260702`

Remote status:

`/data/bht/project01_baker_serhyd_routeB_20260701/manifests/ca_rfd_baker_theozyme_formal_constraints_n1_20260702_status.json`

Remote manifest:

`/data/bht/project01_baker_serhyd_routeB_20260701/manifests/ca_rfd_baker_theozyme_formal_constraints_n1_20260702_manifest.tsv`

Remote log:

`/data/bht/project01_baker_serhyd_routeB_20260701/logs/ca_rfd_baker_theozyme_formal_constraints_n1_20260702.log`

Remote output prefix:

`/data/bht/project01_baker_serhyd_routeB_20260701/outputs/ca_rfd_baker_theozyme_formal_constraints_n1_20260702/sample`

Launch evidence:

- status: `LAUNCHED`
- PID: `408252`
- num_designs: `1`
- design_startnum: `101`
- started on shared GPU with `FORCE=1`, because GPU utilization was 100% but memory was available
- after about 9 minutes, process was still running at the checkpoint/model-load stage, with no output PDB yet

Next action: monitor PID 408252 and the log. If output PDB/TRB appears, run a motif-preservation check against the theozyme contig/motif before launching a larger new-backbone batch.

## Batch50 Launch And Sequence-Bin Scale - 2026-07-02 Update

The n=1 smoke was stopped by user direction; the route now proceeds directly to a larger Baker-theozyme constrained backbone batch.

Stopped smoke:

- run id: `ca_rfd_baker_theozyme_formal_constraints_n1_20260702`
- stopped PID: `408252`
- reason: user requested direct large-scale prediction rather than waiting for one-design smoke completion

Started larger new-backbone batch:

- run id: `ca_rfd_baker_theozyme_formal_constraints_batch50_20260702`
- PID at launch: `413631`
- num_designs: `50`
- design_startnum: `1000`
- constraints: same Baker theozyme input, `bn1`, and contig `12,A56-60,36,A83-85,15,A113-115,73,B145-147,10`
- checkpoint caveat remains: public `/data/bht/design_tools/models/ca_rf_diffusion/ca_rfd_diffusion.pt`, not unavailable Baker `BFF_7.pt`

Remote files:

- status: `/data/bht/project01_baker_serhyd_routeB_20260701/manifests/ca_rfd_baker_theozyme_formal_constraints_batch50_20260702_status.json`
- manifest: `/data/bht/project01_baker_serhyd_routeB_20260701/manifests/ca_rfd_baker_theozyme_formal_constraints_batch50_20260702_manifest.tsv`
- log: `/data/bht/project01_baker_serhyd_routeB_20260701/logs/ca_rfd_baker_theozyme_formal_constraints_batch50_20260702.log`
- output prefix: `/data/bht/project01_baker_serhyd_routeB_20260701/outputs/ca_rfd_baker_theozyme_formal_constraints_batch50_20260702/sample`

At the last check, the batch50 process was still running at the checkpoint/model-load stage and had not yet written output PDB/TRB files.

Downstream sequence-bin scale requested by user:

| sequence identity bin | candidate sequences to generate/screen |
|---:|---:|
| 90% | 200 |
| 80% | 1000 |
| 70% | 1000 |
| 60% | 1000 |
| 50% | 1000 |

These bins apply after Baker-constrained backbones are generated/refined and sequence design begins. The CA_RFDiffusion backbone stage itself does not have sequence-identity bins.
