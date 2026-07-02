# Attempt023 focused committor expansion launched

Rationale: Attempt022 showed that `cand03_rep3_f004_move3` and `cand04_rep7_f003_move3` each split 1 reactant / 1 tetrahedral-intermediate in two short random-velocity replicas, making them the best current TS-like candidates.

Directory:
`blind_work/05_atesa_acylation/attempt023_focused_committor_cand03_cand04`

Design:
- Candidates: `cand03_rep3_f004_move3`, `cand04_rep7_f003_move3`.
- Additional replicas: 6 per candidate, reps 3-8.
- Total new trajectories: 12.
- Each trajectory: DFTB3/MM, 1500 steps, dt = 0.0001 ps, 0.15 ps.
- Classification will use the same reactant/tetrahedral-intermediate basin definitions as attempt021/022.

Expected use:
Combine attempt022 and attempt023 counts to estimate whether each candidate is near a 50/50 committor region or biased to one basin.
