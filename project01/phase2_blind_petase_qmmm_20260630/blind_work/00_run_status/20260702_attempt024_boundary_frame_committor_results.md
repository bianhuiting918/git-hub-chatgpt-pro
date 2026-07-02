# Attempt024 boundary-frame committor results

Directory:
`blind_work/05_atesa_acylation/attempt024_boundary_frame_committor`

Run status:
- 12/12 DFTB3/MM trajectories completed.
- No `sander` process remained after completion.
- No real error signatures were found.
- Classification output: `attempt024_boundary_committor_classification.csv`.

Overall first-hit counts:
- tetrahedral-intermediate basin: 6
- reactant basin: 6

Per selected boundary seed:
- `boundary00_cand02_rep1_f010_move2_rep1_1500step_f006`: 2 tetint / 0 reactant, tetint-biased.
- `boundary01_cand03_rep3_f004_move3_rep2_1500step_f012`: 2 tetint / 0 reactant, tetint-biased.
- `boundary02_cand01_rep1_f010_move1_rep1_1500step_f019`: 0 tetint / 2 reactant, reactant-biased.
- `boundary03_cand00_rep1_f010_move0_rep2_1500step_f007`: 1 tetint / 1 reactant, best split.
- `boundary04_cand02_rep1_f010_move2_rep2_1500step_f002`: 1 tetint / 1 reactant, best split.
- `boundary05_cand03_rep3_f004_move3_rep7_1500step_f005`: 0 tetint / 2 reactant, reactant-biased.

Interpretation:
The pre-commit frame extraction improved the TS search. Boundary03 and boundary04 are better current TS candidates than the original accepted structures, because they split evenly in the first short two-replica test. They should be expanded with more random-velocity replicas before being reported as a TS ensemble.
