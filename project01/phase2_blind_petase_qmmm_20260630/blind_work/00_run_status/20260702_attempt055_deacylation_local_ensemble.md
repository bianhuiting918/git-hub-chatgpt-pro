# Deacylation attempt055: local ensemble around bridge_cser240_h145 (2026-07-02)

Purpose: generate a local deacylation TS-core ensemble around the blind best candidate `bridge_cser240_h145`, without using literature coordinates.

CPU policy: user limit is total CPU <=64. Minimizations were serial with OMP_NUM_THREADS=2. The committor smoke batch used 24 Amber/sander jobs with OMP_NUM_THREADS=2, total process threads = 48 by ps NLWP check.

## Minimized local candidates

| candidate | OW-C A | H2-His A | C-SerO A | C-Oleave A | SerH-His A |
|---|---:|---:|---:|---:|---:|
| n00_center | 2.134 | 1.555 | 2.236 | 3.379 | 1.730 |
| n01_more_attack | 2.116 | 1.516 | 2.262 | 3.347 | 1.704 |
| n02_less_attack | 2.142 | 1.571 | 2.256 | 3.383 | 1.741 |
| n03_more_break | 2.146 | 1.557 | 2.272 | 3.382 | 1.742 |
| n04_less_break | 2.136 | 1.563 | 2.180 | 3.402 | 1.736 |
| n05_proton_late | 2.134 | 1.559 | 2.236 | 3.379 | 1.731 |

## 3000-step committor smoke results

All 24/24 trajectories completed; post-run check found 0 live `sander` processes. Fatal/SCF/NaN/segmentation scan returned no failure matches.

| candidate | product | reactant | undecided | interpretation |
|---|---:|---:|---:|---|
| n00_center | 1 | 1 | 2 | balanced but under-sampled; expand |
| n01_more_attack | 2 | 1 | 1 | product-leaning but still dual-commitment; expand |
| n02_less_attack | 0 | 2 | 2 | lower/reactant-side bracket |
| n03_more_break | 0 | 3 | 1 | reactant-side; do not expand first |
| n04_less_break | 0 | 1 | 3 | mostly undecided/reactant-side |
| n05_proton_late | 1 | 1 | 2 | balanced but under-sampled; expand |

Next action: extend n00/n01/n05 to 8 replicas each and compare with the 16-rep center `bridge_cser240_h145` result.

## Attempt056 extension of n00/n01/n05

Purpose: extend the dual-commitment local neighbors n00/n01/n05 from 4 to 8 replicas each. CPU use: 12 concurrent Amber/sander jobs with OMP_NUM_THREADS=2, total process threads = 24 by ps NLWP check.

All 12 added trajectories completed; post-run check found 0 live `sander` processes. Fatal/SCF/NaN/segmentation scan returned no failure matches.

Combined 8-rep results for expanded local candidates:

| candidate | product | reactant | undecided | product fraction excluding undecided | interpretation |
|---|---:|---:|---:|---:|---|
| n00_center | 1 | 4 | 3 | 0.20 | lower/reactant-side bracket after expansion |
| n01_more_attack | 3 | 2 | 3 | 0.60 | accepted deacylation local TS-core ensemble member |
| n05_proton_late | 3 | 2 | 3 | 0.60 | accepted deacylation local TS-core ensemble member |

Together with `bridge_cser240_h145` (6 product / 4 reactant / 6 undecided over 16 replicas), n01 and n05 provide a small blind deacylation TS-core ensemble.

## Selected local deacylation TS-core ensemble artifact

Directory: `blind_work/06_ts_ensemble/deacylation/cser240_local_ts_core_20260702`

Contents:
- 3 selected `.rst7` representatives.
- 3 corresponding `.pdb` files for visual inspection/alignment.
- `deacylation_cser240_local_ts_core_manifest.tsv`.
- `deacylation_ts_ensemble_snapshot.tsv`.
