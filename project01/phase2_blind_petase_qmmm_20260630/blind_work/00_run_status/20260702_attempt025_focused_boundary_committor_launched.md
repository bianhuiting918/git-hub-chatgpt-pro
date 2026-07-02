# Attempt025 focused boundary committor launched

Rationale: Attempt024 identified two boundary frames with 1 reactant / 1 tetrahedral-intermediate split in the first two replicas.

Directory:
`blind_work/05_atesa_acylation/attempt025_focused_boundary03_boundary04`

Candidates:
- `boundary03_cand00_rep1_f010_move0_rep2_1500step_f007`
- `boundary04_cand02_rep1_f010_move2_rep2_1500step_f002`

Design:
- Additional replicas: 6 per candidate, reps 3-8.
- Total new trajectories: 12.
- DFTB3/MM, 1500 steps, dt = 0.0001 ps.

Goal: combine attempt024 and attempt025 to estimate whether either boundary frame remains near a 50/50 committor after 8 replicas total.
