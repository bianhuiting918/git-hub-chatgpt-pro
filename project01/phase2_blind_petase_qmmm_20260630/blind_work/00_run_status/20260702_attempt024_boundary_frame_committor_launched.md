# Attempt024 boundary-frame committor launched

Rationale: Attempt023 showed that cand03/cand04 accepted structures were still reactant-side biased, while cand02 was tetrahedral-intermediate biased. Therefore, all attempt022/023 trajectories were scanned for pre-commit early-region frames with balanced SerO-C/C-Oleave and partial Ser-H/His proton transfer.

Inputs:
- Ranking: `blind_work/05_atesa_acylation/attempt024_boundary_frame_candidates/attempt024_boundary_frame_ranking.csv`
- Selected frames: `blind_work/05_atesa_acylation/attempt024_boundary_frame_candidates/attempt024_selected_boundary_frames.csv`

Run directory:
`blind_work/05_atesa_acylation/attempt024_boundary_frame_committor`

Design:
- 6 selected boundary frames.
- 2 random-velocity DFTB3/MM replicas per frame.
- 12 total trajectories.
- 1500 steps each, dt = 0.0001 ps.

Goal: identify whether any extracted pre-commit frame gives a better reactant/tetrahedral split than the original accepted structures.
