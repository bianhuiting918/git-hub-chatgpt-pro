# Attempt025 focused boundary committor results

Directory:
`blind_work/05_atesa_acylation/attempt025_focused_boundary03_boundary04`

Run status:
- 12/12 DFTB3/MM trajectories completed.
- No `sander` process remained after completion.
- No real error signatures were found.
- Classification output: `attempt025_focused_boundary_committor_classification.csv`.
- Combined output with attempt024: `attempt025_combined_attempt024_attempt025_boundary03_boundary04.csv`.

Attempt025 first-hit counts:
- tetint basin: 7
- reactant basin: 4
- none: 1

Combined attempt024 + attempt025, 8 replicas per candidate:

| seed | reactant | tetint | none | committed tetint fraction | interpretation |
|---|---:|---:|---:|---:|---|
| boundary03_cand00_rep1_f010_move0_rep2_1500step_f007 | 5 | 3 | 0 | 0.375 | best current TS ensemble core |
| boundary04_cand02_rep1_f010_move2_rep2_1500step_f002 | 1 | 6 | 1 | 0.857 | tetint/product-side biased |

Interpretation:
Boundary03 is now the best current acylation TS-candidate ensemble member. Its 5/3 reactant/tetint split over 8 short DFTB3/MM replicas is not exactly 50/50, but it is the closest validated candidate obtained so far. Boundary04 is too tetint-biased and should not be treated as the central TS structure.

Next recommended step:
- Extract frames around boundary03 trajectories immediately before first basin commitment and run a final smaller refinement/committor expansion, or use boundary03 as the current TS ensemble core for a milestone GitHub sync.
