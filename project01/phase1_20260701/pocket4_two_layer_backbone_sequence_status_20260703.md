# Project 01 Serine Hydrolase Pocket4 L1/L2 Backbone And Sequence Routes

Updated: 2026-07-03 00:47 Asia/Shanghai

## Current Decision

The target is now a two-route comparison, not a single L2-only route.

Route A and Route B are both formal work products:

1. Route A, L1-direct: use pocket4 PASS L1 scaffolds directly as designable backbones, then run complete sequence design on those L1 backbones.
2. Route B, L2-extended: use a pocket4 PASS L1 scaffold as the fixed core, extend it once into a larger L2 backbone, then run complete sequence design on the L2 backbone.

For both routes, the output target is:

- 90/80/70/60/50 percent identity bins.
- At least 10 structure/pocket-gated sequences per bin if the route can support them.
- Explicit PASS/FAIL/NOT_EVALUATED accounting, with no un-evaluated rows counted as biological failures.

No L3/full third layer is planned unless both L1-direct and L2-extended routes fail to produce enough acceptable sequence/structure candidates.

## Hard Evidence So Far

L1 pocket4 PASS scaffolds:

| parent_set | run_id | output | length_aa | ligand_heavy_records | gate |
|---|---|---|---:|---:|---|
| compact | ca_rfd_baker_pocket4_layered_compact_n1_v2_20260702 | sample_7200.pdb | 76 | 22 | PASS |
| medium | ca_rfd_baker_pocket4_layered_medium_n1_20260702 | sample_7300.pdb | 82 | 22 | PASS |
| expanded | ca_rfd_baker_pocket4_layered_expanded_n1_20260702 | sample_7400.pdb | 88 | 22 | PASS |

The earlier LigandMPNN panel from these L1 scaffolds had 10620 sequence-layer candidates. Under the revised goal, this panel is the raw Route A L1-direct sequence panel, not the final accepted set. It still requires structure prediction and pocket/core gate before any per-bin success count can be claimed.

Route A existing raw panel:

| manifest | total | 90 | 80 | 70 | 60 | 50 |
|---|---:|---:|---:|---:|---:|---:|
| manifests/pocket4_layered_pass_scaffolds_ligandmpnn_20260702_2350_combined_selected.tsv | 10620 | 173 | 1667 | 2780 | 3000 | 3000 |

## Route B Active L2 Smoke

Remote root:

`/data/bht/project01_baker_serhyd_routeB_20260701`

Active run:

| field | value |
|---|---|
| run_id | ca_rfd_baker_pocket4_l2_extend_expanded_n1_20260703 |
| parent | expanded L1, sample_7400.pdb |
| parent_len | 88 aa |
| contig | `24-36,A1-88,24-36` |
| intended length | about 136-160 aa |
| PID at launch | 959968 |
| log | logs/ca_rfd_baker_pocket4_l2_extend_expanded_n1_20260703.log |
| output_prefix | outputs/ca_rfd_baker_pocket4_l2_extend_expanded_n1_20260703/sample |
| checkpoint caveat | public `ca_rfd_diffusion.pt`, not Baker `BFF_7.pt` |

This run fixes the whole L1 expanded core and extends both termini in one step. If this passes gate, it is treated as the final backbone layer for Route B sequence design.

Observed progress as of 2026-07-03 00:47 Asia/Shanghai:

- PID 959968 is still running.
- The log has entered denoising and reached at least t=48/50.
- Runtime is slow because GPU 0 is shared and already heavily occupied.
- This is still a valid active run; no second L2 job should be launched in parallel.

## Gate For Route B L2

Because the N-terminal extension length is sampled from `24-36`, the parent core offset in the L2 output is not fixed. The L2 gate therefore scans possible offsets and selects the best parent-core alignment.

Script:

`scripts/gate_pocket4_l2_parent_core.py`

Gate thresholds:

| metric | cutoff |
|---|---:|
| parent-core CA RMSD | <= 1.0 A |
| parent-core pair max delta | <= 1.0 A |
| ligand records | required |

This is a downstream project gate, not a native CA_RFDiffusion binary criterion.

## Next Actions

1. Route A: use the existing L1-direct raw LigandMPNN panel as input for structure prediction and pocket/core gate; report accepted counts per bin.
2. Route B: monitor PID 959968 until `sample_7500.pdb`/`.trb` appear or the log reports failure.
3. Route B: run `gate_pocket4_l2_parent_core.py` against the L2 output directory.
4. Route B: if L2 gate PASS, run LigandMPNN on the L2 backbone and generate 90/80/70/60/50 sequence bins.
5. Compare Route A versus Route B by per-bin accepted sequence counts and pocket geometry stability.

Do not upload PDB/TRB/model/trajectory/log files to GitHub.

## Live Progress Update 2026-07-03 00:55 Asia/Shanghai

Route A L1-direct ESMFold structure-screening smoke has been launched:

| field | value |
|---|---|
| run_id | routeA_l1_direct_esmfold_lpb10_20260703_0052 |
| PID | 966399 |
| input panel | manifests/pocket4_layered_pass_scaffolds_ligandmpnn_20260702_2350_combined_selected.tsv |
| intended evaluated universe | 50 rows, 10 per bin from 90/80/70/60/50 |
| current output count | 0 PDB at first check |
| current status | running/loading; no PASS/FAIL assigned yet |

Route B L2 extension remains active:

| field | value |
|---|---|
| PID | 959968 |
| current log evidence | denoising reached at least t=39/50 |
| current status | running; no L2 output PDB yet |

Audit note: unevaluated Route A rows and unfinished Route B output are NOT_EVALUATED, not FAIL.

## Live Progress Update 2026-07-03 01:42 Asia/Shanghai

Route A L1-direct ESMFold smoke completed 50/50 structures and strict pocket gate was run.

Project strict filter result for Route A run `routeA_l1_direct_esmfold_lpb10_20260703_0052`:

| bin | evaluated | PASS | FAIL | NOT_EVALUATED |
|---:|---:|---:|---:|---:|
| 50 | 10 | 0 | 10 | 2990 |
| 60 | 10 | 0 | 10 | 2990 |
| 70 | 10 | 0 | 10 | 2770 |
| 80 | 10 | 0 | 10 | 1657 |
| 90 | 10 | 0 | 10 | 163 |

Route B L2 extension completed and passed parent-core gate:

| sample | gate | parent_core_CA_RMSD_A | parent_core_pair_max_delta_A | ligand_records | best_n_flank |
|---|---|---:|---:|---:|---:|
| sample_7500 | PASS | 0.1212 | 0.3864 | 22 | 27 |

Route B L2 LigandMPNN raw panel completed:

| selected_tsv | total | 90 | 80 | 70 | 60 | 50 |
|---|---:|---:|---:|---:|---:|---:|
| manifests/pocket4_l2_ligandmpnn_from_sample7500_20260703_0135_selected.tsv | 1701 | 102 | 399 | 400 | 400 | 400 |

Route B L2 ESMFold structure screening was launched:

| field | value |
|---|---|
| run_id | routeB_l2_esmfold_lpb10_20260703_0140 |
| PID | 983495 |
| intended evaluated universe | 50 rows, 10 per bin from 90/80/70/60/50 |
| gate script | scripts/gate_pocket4_l2_esmfold_core.py |
| current status | running/loading at first check; no PDB yet |

Audit note: Route A FAIL counts are for evaluated ESMFold outputs only. Unpredicted rows remain NOT_EVALUATED. Route B raw sequence counts are not final accepted sequence counts until ESMFold/core gate is run.

## Live Progress Update 2026-07-03 01:51 Asia/Shanghai

Route B L2 ESMFold run `routeB_l2_esmfold_lpb10_20260703_0140` has started producing structures.

Current partial status:

| metric | value |
|---|---:|
| ESMFold status rows | 3 |
| PDB outputs | 3 |
| evaluated bin | 50 only so far |
| bin50 PASS | 0 |
| bin50 FAIL | 3 |

The first three Route B bin50 failures are due to strict core geometry and pLDDT filters, with core CA RMSD about 13-14 A, pair max delta about 35-40 A, and core mean pLDDT about 42.

Audit note: this does not evaluate bins 60/70/80/90 yet, and all unpredicted rows remain NOT_EVALUATED.
