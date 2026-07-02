# Attempt023 focused committor results

Directory:
`blind_work/05_atesa_acylation/attempt023_focused_committor_cand03_cand04`

Run status:
- 12/12 focused DFTB3/MM trajectories completed.
- No `sander` process remained after completion.
- No real error signatures were found.
- Classification output: `attempt023_focused_committor_classification.csv`.
- Combined attempt022+attempt023 output: `attempt023_combined_attempt022_attempt023_cand03_cand04.csv`.

Attempt023 first-hit counts:
- reactant basin: 7
- tetrahedral-intermediate basin: 0
- no basin hit: 5

Combined with attempt022:
- `cand03_rep3_f004_move3`: 8 total replicas = 6 reactant, 1 tetint, 1 none. Committed tetint fraction = 1/7 = 0.143.
- `cand04_rep7_f003_move3`: 8 total replicas = 3 reactant, 1 tetint, 4 none. Committed tetint fraction = 1/4 = 0.25.

Interpretation:
The two initially promising candidates are still reactant-side biased under this short committor test. `cand02_rep1_f010_move2` from attempt022 is tetint-biased (2/2 tetint), while cand03/cand04 are reactant-biased. Therefore the next productive search should not simply add more replicas to cand03/cand04; it should extract/rank pre-commit boundary frames between these regions and launch a new focused test from those frames.
