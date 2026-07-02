# Attempt021 update: rep1 branch completed with 3 accepted TS candidates

Check time: 2026-07-02 12:42-12:44 Asia/Shanghai.

## Completed branch

Directory:
`blind_work/05_atesa_acylation/rank022_as_attempt021_tetint_seed_rep1_f010`

ATESA history:
- total moves: 4
- accepted moves: 3
- results: `[['fwd','bwd'], ['bwd','fwd'], ['bwd','fwd'], ['fwd','fwd']]`

The fourth move is not accepted and is not included in the TS-candidate table.

Accepted candidates were extracted to:
`blind_work/05_atesa_acylation/rank022_as_attempt021_tetint_seed_rep1_f010/accepted_ts_candidates/`

Geometry table:
`accepted_ts_candidate_geometries.csv`

| move | result | SerO-C | C-Oleave | SerH-SerO | SerH-His | SerH-Oleave | His-Oleave |
|---:|---|---:|---:|---:|---:|---:|---:|
| 0 | fwd/bwd | 2.046 | 2.055 | 1.045 | 1.707 | 2.525 | 3.529 |
| 1 | bwd/fwd | 1.989 | 1.976 | 1.048 | 1.581 | 2.763 | 3.527 |
| 2 | bwd/fwd | 2.013 | 1.933 | 0.994 | 1.575 | 2.736 | 3.497 |

## Still running at check time

At the initial check, 5 `sander` processes were still running and ATESA workers remained active for the other seeds.

Other seed histories at check time:
- `seed_rep2_f009`: moves 2, accepted 0, results `[['bwd','bwd'], ['', 'bwd']]`.
- `seed_rep3_f004`: moves 3, accepted 0, results `[['bwd','bwd'], ['bwd','bwd'], ['bwd','bwd']]`.
- `seed_rep7_f003`: moves 3, accepted 0, results `[['bwd','bwd'], ['bwd','bwd'], ['bwd','bwd']]`.

## Interpretation

The rep1 branch is currently the only productive attempt021 branch. Its accepted candidates cluster around early acylation crossing geometries with SerO-C around 1.99-2.05 A and C-Oleave around 1.93-2.05 A, earlier than the old product-side points near C-Oleave ~2.45 A.

This is a TS-candidate ensemble from aimless shooting, not yet a final validated TS ensemble. Next required validation is committor-style short shooting from the extracted candidates plus geometry clustering.
