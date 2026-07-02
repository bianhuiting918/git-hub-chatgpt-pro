# Project 01 Serine Hydrolase Pocket4 Two-Layer Backbone Route

Updated: 2026-07-03 00:35 Asia/Shanghai

## Current Decision

Do not continue the premature branch that designs Ala-rich sequences directly on the L1 mini-scaffolds and sends them to ESMFold.

The route is now capped at two backbone layers:

1. L1: generate and gate pocket-core scaffolds around the Baker serine hydrolase theozyme.
2. L2: use a PASS L1 scaffold as the fixed core and directly extend it into the final designable backbone.
3. Sequence design: only after an L2 backbone passes parent-core/pocket geometry gate, run LigandMPNN to generate complete sequences and bin them at 90/80/70/60/50 percent identity.

No L3/full third layer is planned unless L2 repeatedly fails.

## Hard Evidence So Far

L1 pocket4 PASS scaffolds:

| parent_set | run_id | output | length_aa | ligand_heavy_records | gate |
|---|---|---|---:|---:|---|
| compact | ca_rfd_baker_pocket4_layered_compact_n1_v2_20260702 | sample_7200.pdb | 76 | 22 | PASS |
| medium | ca_rfd_baker_pocket4_layered_medium_n1_20260702 | sample_7300.pdb | 82 | 22 | PASS |
| expanded | ca_rfd_baker_pocket4_layered_expanded_n1_20260702 | sample_7400.pdb | 88 | 22 | PASS |

The earlier LigandMPNN panel from these L1 scaffolds had 10620 sequence-layer candidates, but these are **not final accepted sequences** because the parent scaffolds are only 76/82/88 aa and are not the complete backbone target.

## Active L2 Smoke

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

This run fixes the whole L1 expanded core and extends both termini in one step. If this passes gate, it is treated as the final backbone layer for sequence design.

## Gate For L2

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

1. Monitor PID 959968 until `sample_7500.pdb`/`.trb` appear or the log reports failure.
2. Run `gate_pocket4_l2_parent_core.py` against the L2 output directory.
3. If PASS, run LigandMPNN on the L2 backbone and generate 90/80/70/60/50 sequence bins.
4. Only after L2 sequence design, evaluate predicted structures/pocket geometry and count accepted final sequences.

Do not upload PDB/TRB/model/trajectory/log files to GitHub.
