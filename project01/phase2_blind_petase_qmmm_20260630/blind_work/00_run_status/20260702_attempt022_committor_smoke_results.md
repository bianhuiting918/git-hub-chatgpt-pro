# Attempt022 committor-smoke results

Directory:
`blind_work/05_atesa_acylation/attempt022_committor_smoke_from_attempt021`

Input:
`blind_work/05_atesa_acylation/attempt021_all_accepted_ts_candidates.csv`

Run status:
- 10/10 DFTB3/MM trajectories completed.
- No real error signatures were found.
- Output classification table: `attempt022_committor_smoke_classification.csv`.

First-hit classification over 0.15 ps:
- reactant basin: 4
- tetrahedral-intermediate basin: 4
- no basin hit: 2

Per-candidate result:
- `cand00_rep1_f010_move0`: 1 reactant, 1 none.
- `cand01_rep1_f010_move1`: 1 reactant, 1 none.
- `cand02_rep1_f010_move2`: 2 tetrahedral-intermediate; likely product/tetint biased for this basin definition.
- `cand03_rep3_f004_move3`: 1 reactant, 1 tetrahedral-intermediate; best TS-like split.
- `cand04_rep7_f003_move3`: 1 reactant, 1 tetrahedral-intermediate; best TS-like split.

Interpretation:
The full attempt021 accepted set is not homogeneous. `cand03` and `cand04` are the best current TS-candidate structures because even with only two random-velocity replicas each, they split to opposite basins. `cand02` appears downstream/tetrahedral-biased; `cand00` and `cand01` are reactant/undecided biased.

Next action:
Run attempt023 focused committor expansion on `cand03` and `cand04` only.
