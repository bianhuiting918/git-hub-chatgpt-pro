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

## Batch50 Runtime Check - 2026-07-02

A later check confirmed the batch50 job is not stuck at launch. It entered denoising for the first backbone:

- PID: `413631`
- current design: `sample_1000`
- log evidence: `making design 1000 of 1000:1050`
- denoising reached at least `t=tensor(43)` / about `7/50` steps
- observed speed: about 48 seconds per denoising step on the shared GPU
- output files at check time: `0`

Interpretation: the Baker-theozyme constrained batch is running normally, but the first backbone has not completed yet. Do not start a second CA_RFDiffusion batch until this one has produced enough output or is intentionally stopped.

## Candidate Scale Rule - 2026-07-02 User-Confirmed

After a Baker-theozyme constrained new backbone passes the motif/refinement gate, sequence generation will use the following candidate counts per identity bin:

- 90% identity bin: 200 candidate sequences
- 80% identity bin: 1000 candidate sequences
- 70% identity bin: 1000 candidate sequences
- 60% identity bin: 1000 candidate sequences
- 50% identity bin: 1000 candidate sequences

These numbers are candidate-generation denominators, not final pass counts. The final target remains 10 structure/pocket-qualified sequences per identity bin.

## Batch50 Runtime Progress - 2026-07-02 12:44 CST

Latest GPU check:

- run id: `ca_rfd_baker_theozyme_formal_constraints_batch50_20260702`
- PID: `413631`
- process state: still running
- elapsed: about 38 minutes
- GPU: A100 80 GB, 100% utilization; this CA_RFDiffusion process used about 3700 MiB
- log progress: first design `sample_1000`, denoising reached `34/50`
- output files: none yet under the batch50 output directory

Next action: continue monitoring until the first PDB/TRB is written, then run a motif/theozyme CA gate before any sequence-bin generation is counted.

## First Batch50 Motif Gate Result - 2026-07-02 12:57 CST

The first CA_RFDiffusion output from the formal Baker-theozyme batch was written while the batch continued to `sample_1001`.

Remote output files:

- PDB: `/data/bht/project01_baker_serhyd_routeB_20260701/outputs/ca_rfd_baker_theozyme_formal_constraints_batch50_20260702/sample_1000.pdb`
- TRB: `/data/bht/project01_baker_serhyd_routeB_20260701/outputs/ca_rfd_baker_theozyme_formal_constraints_batch50_20260702/sample_1000.trb`

Motif/theozyme gate output:

- TSV: `/data/bht/project01_baker_serhyd_routeB_20260701/manifests/ca_rfd_baker_theozyme_formal_constraints_batch50_20260702_motif_gate.tsv`
- JSON: `/data/bht/project01_baker_serhyd_routeB_20260701/manifests/ca_rfd_baker_theozyme_formal_constraints_batch50_20260702_motif_gate_summary.json`

Evaluated universe at this check: 1 written PDB. Future samples not yet written are `NOT_EVALUATED`, not failures.

Result for `sample_1000`:

- gate: `PASS`
- motif CA Kabsch RMSD: `0.1041 A`
- max motif pair-distance delta: `0.1532 A`
- mean motif pair-distance delta: `0.0525 A`
- ligand `bn1` atom records: `22`

This is a first-layer Baker theozyme/motif pass for a generated backbone. It is not a final sequence-bin pass. The next stage is to accumulate additional motif-pass backbones from the running batch and prepare the refinement/sequence-design stage using the user-confirmed candidate scale: 90%=200; 80/70/60/50=1000 each.

## Sample1000 Public Refinement Launch - 2026-07-02 13:10 CST

Because `sample_1000` passed the Baker theozyme/motif gate, a single small public-checkpoint refinement task was launched while the main batch50 continued.

Refinement run:

- run id: `ca_rfd_baker_theozyme_refine_sample1000_publicckpt_20260702`
- PID at launch: `426780`
- input: `/data/bht/project01_baker_serhyd_routeB_20260701/outputs/ca_rfd_baker_theozyme_formal_constraints_batch50_20260702/sample_1000.pdb`
- checkpoint: `/data/bht/design_tools/models/ca_rf_diffusion/ca_rfd_refinement.pt`
- checkpoint caveat: public `ca_rfd_refinement.pt`, not unavailable Baker `refine_BFF_3.pt`
- log: `/data/bht/project01_baker_serhyd_routeB_20260701/logs/ca_rfd_baker_theozyme_refine_sample1000_publicckpt_20260702.log`
- status: `/data/bht/project01_baker_serhyd_routeB_20260701/manifests/ca_rfd_baker_theozyme_refine_sample1000_publicckpt_20260702_status.json`
- manifest: `/data/bht/project01_baker_serhyd_routeB_20260701/manifests/ca_rfd_baker_theozyme_refine_sample1000_publicckpt_20260702_manifest.tsv`

Launch policy: only one small refinement was started, not a parallel refinement batch. The main batch50 CA_RFDiffusion PID `413631` continued running.

Next gate after refinement: run the same Baker theozyme motif CA gate on the refined output before using it for sequence-bin generation.

## Sample1000 Refinement Debug And Relaunch - 2026-07-02 13:14 CST

The first refinement launch failed before model inference:

- failed run id: `ca_rfd_baker_theozyme_refine_sample1000_publicckpt_20260702`
- failure: `ModuleNotFoundError: No module named 'rf_diffusion'`
- root cause: `run_inference.py` was launched from the `rf_diffusion` package directory without the CA_RFDiffusion repository root on `PYTHONPATH`

Minimal environment test passed:

`PYTHONPATH=/data/bht/design_tools/src/CA_RFDiffusion python -c 'import rf_diffusion'`

A corrected single-sample refinement was relaunched:

- run id: `ca_rfd_baker_theozyme_refine_sample1000_publicckpt_v2_20260702`
- PID at launch: `427624`
- input: `/data/bht/project01_baker_serhyd_routeB_20260701/outputs/ca_rfd_baker_theozyme_formal_constraints_batch50_20260702/sample_1000.pdb`
- checkpoint: `/data/bht/design_tools/models/ca_rf_diffusion/ca_rfd_refinement.pt`
- environment fix: `PYTHONPATH=/data/bht/design_tools/src/CA_RFDiffusion`
- log: `/data/bht/project01_baker_serhyd_routeB_20260701/logs/ca_rfd_baker_theozyme_refine_sample1000_publicckpt_v2_20260702.log`
- status: `/data/bht/project01_baker_serhyd_routeB_20260701/manifests/ca_rfd_baker_theozyme_refine_sample1000_publicckpt_v2_20260702_status.json`
- manifest: `/data/bht/project01_baker_serhyd_routeB_20260701/manifests/ca_rfd_baker_theozyme_refine_sample1000_publicckpt_v2_20260702_manifest.tsv`

Latest v2 log showed successful import/config startup and checkpoint loading; completion still pending at this status update.

## Sample1000 Refined Motif Gate - 2026-07-02 13:30 CST

The corrected v2 refinement completed successfully:

- run id: `ca_rfd_baker_theozyme_refine_sample1000_publicckpt_v2_20260702`
- runtime: 14.57 minutes
- refined PDB: `/data/bht/project01_baker_serhyd_routeB_20260701/outputs/ca_rfd_baker_theozyme_formal_constraints_batch50_20260702/sample_1000_refined_0.pdb`
- refined TRB: `/data/bht/project01_baker_serhyd_routeB_20260701/outputs/ca_rfd_baker_theozyme_formal_constraints_batch50_20260702/sample_1000_refined_0.trb`

Updated motif/theozyme gate evaluated 2 written PDBs in the batch50 output directory:

- `sample_1000`: `PASS`, motif CA RMSD `0.1041 A`, max pair-distance delta `0.1532 A`, ligand `bn1` records `22`
- `sample_1000_refined_0`: `PASS`, motif CA RMSD `0.0865 A`, max pair-distance delta `0.1922 A`, ligand `bn1` records `22`

Remote gate files:

- TSV: `/data/bht/project01_baker_serhyd_routeB_20260701/manifests/ca_rfd_baker_theozyme_formal_constraints_batch50_20260702_motif_gate.tsv`
- JSON: `/data/bht/project01_baker_serhyd_routeB_20260701/manifests/ca_rfd_baker_theozyme_formal_constraints_batch50_20260702_motif_gate_summary.json`

This proves refinement preserved the Baker motif gate for the first backbone. It still does not complete any sequence-similarity bin; next action is sequence design/binning from `sample_1000_refined_0` while batch50 continues producing more backbones.

## Sequence-Bin Generation From Sample1000 Refined - 2026-07-02 13:53 CST

Sequence design has started from the first refined Baker-theozyme backbone:

- source backbone: `/data/bht/project01_baker_serhyd_routeB_20260701/outputs/ca_rfd_baker_theozyme_formal_constraints_batch50_20260702/sample_1000_refined_0.pdb`
- fixed residues: Baker theozyme contig motif positions only, `A13,A14,A15,A16,A17,A54,A55,A56,A72,A73,A74,A148,A149,A150`
- fixed policy: no ligand-4A shell lock for this new-backbone route
- output summary: `/data/bht/project01_baker_serhyd_routeB_20260701/manifests/baker_theozyme_sample1000_refined_ligandmpnn_bins_summary.json`
- selected manifest: `/data/bht/project01_baker_serhyd_routeB_20260701/manifests/baker_theozyme_sample1000_refined_ligandmpnn_bins_selected.tsv`

90% identity bin:

- target candidate count: 200
- status: `TARGET_MET`
- selected: `200/200`
- exact mutation count: 16 mutations over 160 residues
- fixed motif mutations: 0 in selected records
- first pass at temperature 0.10 produced only 59 unique pass records due to duplicate sequences
- supplement at temperature 0.20 with seed `991101` added enough unique pass records; the final selected count is 200

Batch50 backbone motif gate also advanced:

- evaluated PDBs: 3
- PASS: `sample_1000`, `sample_1000_refined_0`, `sample_1001`
- `sample_1001` motif CA RMSD: `0.1015 A`; max pair-distance delta: `0.1465 A`; ligand `bn1` records: `22`

80% identity bin:

- target candidate count: 1000
- launched run: `baker_theozyme_sample1000_refined_ligandmpnn_bin80_20260702`
- PID at launch: `437140`
- log: `/data/bht/project01_baker_serhyd_routeB_20260701/logs/baker_theozyme_sample1000_refined_ligandmpnn_bin80_20260702.log`

These are sequence candidates only. They are not final-qualified structures until post-sequence structure prediction and pocket/motif gates pass.

## Sequence-Bin Generation Continued - 2026-07-02 13:58 CST

Updated sequence candidate status from `sample_1000_refined_0`:

- 90% bin: `200/200` selected, target met
- 80% bin: `1000/1000` selected, target met
- selected manifest: `/data/bht/project01_baker_serhyd_routeB_20260701/manifests/baker_theozyme_sample1000_refined_ligandmpnn_bins_selected.tsv`
- summary: `/data/bht/project01_baker_serhyd_routeB_20260701/manifests/baker_theozyme_sample1000_refined_ligandmpnn_bins_summary.json`

80% bin details:

- records generated: 1251
- pass filter count: 1162
- selected count: 1000
- exact mutation count: 32 mutations over 160 residues
- fixed motif mutations: 0 in selected records by filter definition

70% identity bin was launched next:

- target candidate count: 1000
- launched run: `baker_theozyme_sample1000_refined_ligandmpnn_bin70_20260702`
- PID at launch: `438180`
- log: `/data/bht/project01_baker_serhyd_routeB_20260701/logs/baker_theozyme_sample1000_refined_ligandmpnn_bin70_20260702.log`

Reminder: these are sequence candidate manifests. Final-qualified status still requires post-sequence structure prediction and pocket/motif gates.

## Sequence-Bin Candidate Pool Complete - 2026-07-02 14:12 CST

The Baker-theozyme new-backbone sample1000 refined route now has the requested
candidate sequence pool size:

- 90% identity bin: `200/200` selected; exact mutation count `16/160`
- 80% identity bin: `1000/1000` selected; exact mutation count `32/160`
- 70% identity bin: `1000/1000` selected; exact mutation count `48/160`
- 60% identity bin: `1000/1000` selected; exact mutation count `64/160`
- 50% identity bin: `1000/1000` selected; exact mutation count `80/160`
- total selected sequence candidates: `4200`

Remote evidence:

- summary JSON: `/data/bht/project01_baker_serhyd_routeB_20260701/manifests/baker_theozyme_sample1000_refined_ligandmpnn_bins_summary.json`
- selected TSV: `/data/bht/project01_baker_serhyd_routeB_20260701/manifests/baker_theozyme_sample1000_refined_ligandmpnn_bins_selected.tsv`

Bin-level generation evidence:

- 90%: `records=401`, initial `pass_filter_count=59`; supplemented at temperature `0.20`, final selected `200`
- 80%: `records=1251`, `pass_filter_count=1162`, selected `1000`
- 70%: `records=1251`, `pass_filter_count=1250`, selected `1000`
- 60%: `records=1251`, `pass_filter_count=1250`, selected `1000`
- 50%: `records=1251`, `pass_filter_count=1250`, selected `1000`

Important audit note: these are sequence-similarity-bin candidates only. They are
not final-qualified enzyme designs until structure prediction, side-chain/ligand
rebuild as needed, and Baker motif/pocket geometry gates pass. The fixed residues
for this new-backbone route are the Baker theozyme contig motif positions only
(`A13-A17`, `A54-A56`, `A72-A74`, `A148-A150`); this route does not lock the
full ligand-4A shell.

Concurrent backbone generation status:

- CA_RFDiffusion batch50 PID `413631` was still running at 2026-07-02 14:14 CST
- current written PDBs in batch50 output: `sample_1000.pdb`, `sample_1000_refined_0.pdb`, `sample_1001.pdb`
- post-sequence structure prediction and final pocket/motif gates remain pending
