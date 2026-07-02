# attempt037 deacylation path probe status (2026-07-02)

Purpose: low-cost restrained QM/MM path probes from blind HOH272 preorganized seeds. This is a path-probe step only; no structure here is claimed as a transition state without committor/path validation.

Directory:
`blind_work/04_qmmm_exploration/attempt037_deacylation_path_probe`

Inputs:
- chainA_v01_h1: attempt036 v01_h1_mid seed, selected water H1 to His.
- chainB_v05_h2: attempt036 v05_h2_tight seed, selected water H2 to His.
- QM region includes the acylation QM region plus HOH272 atoms 3863/3864/3865 (Amber 1-based).

Run status:
- chainA_v01_h1: w00-w04 completed.
- chainB_v05_h2: w00-w04 completed.
- No serious Amber/QMMM error strings were found in the window logs by the final grep scan.
- No residual `run_chain.sh` or `sander -O -i w..in` processes remained after completion.

Key geometry summary (Angstrom):

| chain/window | OW-C | selected H-His | C-SerO | C-Oleave | triage |
|---|---:|---:|---:|---:|---|
| chainA w00 | 3.117 | 2.610 | 1.394 | 2.485 | path_probe_intermediate |
| chainA w01 | 2.492 | 1.255 | 2.104 | 2.273 | path_probe_intermediate |
| chainA w02 | 2.489 | 1.919 | 2.032 | 2.062 | path_probe_intermediate |
| chainA w03 | 2.118 | 1.363 | 2.840 | 1.648 | deacylation_ts_like_probe |
| chainA w04 | 1.904 | 1.198 | 2.530 | 1.852 | deacylation_ts_like_probe |
| chainB w00 | 3.171 | 3.343 | 1.402 | 2.448 | path_probe_intermediate |
| chainB w01 | 2.928 | 3.107 | 1.393 | 2.384 | path_probe_intermediate |
| chainB w02 | 2.650 | 2.210 | 1.454 | 2.428 | path_probe_intermediate |
| chainB w03 | 2.283 | 1.308 | 1.390 | 2.801 | path_probe_intermediate |
| chainB w04 | 1.911 | 1.447 | 1.449 | 3.307 | water_pt_advanced_acyl_bond_intact |

Triage interpretation:
- chainA w03/w04 are the current best deacylation path candidates because water O approaches the acyl carbon while the selected water H is close to His and C-SerO is lengthened.
- chainB w04 should not be used directly as a deacylation TS candidate because C-SerO remains short (~1.45 A), even though OW-C and H-His are advanced.
- The next scientific step is not more blind pulling; it is short unbiased/weakly restrained committor-style tests seeded from chainA w03/w04, with early-product/reactant classification based on C-SerO, OW-C, and H-His evolution.

Machine-readable table:
`blind_work/04_qmmm_exploration/attempt037_deacylation_path_probe/attempt037_deacylation_path_probe_summary.csv`
