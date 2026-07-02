# Attempt020/021 PETase acylation TS-candidate update

## Context

The earlier product-side AS points had C-Oleave around 2.45 A and were judged too late for the current blind mechanism objective. The branch was redirected to find earlier reaction-state / TS-like structures before full ester C-O cleavage.

## Attempt020: parallel expansion from the best early-state seed

Directory:
`blind_work/05_atesa_acylation/rank022_manual_shooting_attempt008_direct_seed/attempt020_parallel_frame039_expansion`

Design:
- Primary seed: `a025_early_rep2_2000step_frame039.rst7`, expanded into 8 independent 0.2 ps DFTB3/MM replicas.
- Controls: frame041 and frame002, 2 replicas each.
- Total: 12 parallel `sander` jobs, all completed normally with 12 `.rst7` and 12 `.nc` files.

Classification outputs:
- `attempt020_all_frame_classification.csv`
- `attempt020_trajectory_summary.csv`

Key result:
- `frame039_primary` final states: 3 reactant-like, 4 other, 1 tetrahedral-intermediate-like.
- `frame041_boundary_ctrl`: 2/2 reactant-like.
- `frame002_balanced_ctrl`: 1 reactant-like, 1 other.

Interpretation:
`frame039` is the most useful dividing-region seed because independent replicas can return to reactant or proceed toward tetrahedral-intermediate-like geometry. The controls are more reactant-biased.

Selected attempt021 seed frames:
- `ts_seed_rep7_f003`: SerO-C 2.017 A, C-Oleave 2.043 A.
- `ts_seed_rep2_f009`: SerO-C 1.999 A, C-Oleave 2.087 A.
- `ts_seed_rep3_f004`: SerO-C 2.038 A, C-Oleave 2.052 A.
- `ts_seed_rep1_f010`: SerO-C 2.047 A, C-Oleave 2.055 A.

## Attempt021: aimless shooting toward tetrahedral intermediate, not late product

Directories:
- `blind_work/05_atesa_acylation/rank022_as_attempt021_tetint_seed_rep1_f010`
- `blind_work/05_atesa_acylation/rank022_as_attempt021_tetint_seed_rep2_f009`
- `blind_work/05_atesa_acylation/rank022_as_attempt021_tetint_seed_rep3_f004`
- `blind_work/05_atesa_acylation/rank022_as_attempt021_tetint_seed_rep7_f003`

Basin definition changed from the old late/product-side leaving-O-protonated condition:

Forward basin, tetrahedral-intermediate-like:
- SerO-C < 1.8 A
- SerH-SerO > 1.4 A
- SerH-HisNE2 < 1.2 A
- C-Oleave > 1.8 A

Backward basin, Michaelis/reactant-like:
- SerO-C > 2.4 A
- SerH-SerO < 1.2 A
- SerH-HisNE2 > 1.4 A
- C-Oleave < 1.8 A

This follows the mechanistic objective of finding the early acylation TS / tetrahedral intermediate crossing, without forcing the later leaving-O protonation/product-side state.

## Current accepted TS-candidate ensemble

Best branch: `rank022_as_attempt021_tetint_seed_rep1_f010`.

Current ATESA history at extraction time:
- total moves: 3
- accepted moves: 3
- results: `[['fwd','bwd'], ['bwd','fwd'], ['bwd','fwd']]`

Extracted candidates:
`blind_work/05_atesa_acylation/rank022_as_attempt021_tetint_seed_rep1_f010/accepted_ts_candidates/`

Geometry CSV:
`accepted_ts_candidate_geometries.csv`

Accepted candidate geometries:

| move | result | SerO-C | C-Oleave | SerH-SerO | SerH-His | SerH-Oleave | His-Oleave |
|---:|---|---:|---:|---:|---:|---:|---:|
| 0 | fwd/bwd | 2.046 | 2.055 | 1.045 | 1.707 | 2.525 | 3.529 |
| 1 | bwd/fwd | 1.989 | 1.976 | 1.048 | 1.588 | 2.763 | 3.527 |
| 2 | bwd/fwd | 2.013 | 1.933 | 0.994 | 1.576 | 2.736 | 3.497 |

Interpretation:
These are earlier than the old product-side accepted points (`C-Oleave ~2.45 A`) and cluster around simultaneous SerO attack / ester weakening before full leaving-group protonation. They are TS-candidates from shooting, not final validated TS ensemble yet.

## Other attempt021 seeds

At the last checked history:
- `seed_rep2_f009`: biased bwd; no accepted point yet.
- `seed_rep3_f004`: biased bwd; no accepted point yet.
- `seed_rep7_f003`: three bwd/bwd moves; de-prioritize for expansion.

## Next actions

1. Let attempt021 finish all currently running moves.
2. Re-read `restart.pkl` in each attempt021 directory and update accepted counts.
3. Extract any new accepted points into the same candidate table.
4. Run geometry clustering and short committor validation on accepted candidates.
5. Only after the attempt021 branch is coherent, sync this milestone to GitHub.
