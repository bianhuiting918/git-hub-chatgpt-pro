# PETase acylation TS search milestone through attempt025

This milestone follows the blind mechanism-search strategy: the literature was used for method direction only, while seeds and TS candidates were selected from our own structure/trajectory evidence.

## Main outcome

Current best TS-candidate core:
`boundary03_cand00_rep1_f010_move0_rep2_1500step_f007`

Evidence:
- It was extracted as a pre-commit boundary frame from short committor trajectories.
- Combined attempt024 + attempt025 short DFTB3/MM committor test gives 8 replicas:
  - reactant basin: 5
  - tetrahedral-intermediate basin: 3
  - none: 0
  - committed tetint fraction: 0.375

This is the closest validated candidate so far to a reactant/tetrahedral-intermediate dividing region.

## Important negative result

`boundary04_cand02_rep1_f010_move2_rep2_1500step_f002` is tetrahedral-intermediate biased:
- reactant basin: 1
- tetrahedral-intermediate basin: 6
- none: 1
- committed tetint fraction: 0.857

It should not be treated as the central TS candidate.

## Attempt sequence

- attempt020: expanded the best early-state seed (`frame039`) and confirmed mixed reactant/tetrahedral behavior.
- attempt021: ATESA aimless shooting from four independent early-region seeds; obtained 5 accepted TS-candidates.
- attempt022: light committor smoke over all 5 accepted candidates; selected cand03/cand04 for focused testing.
- attempt023: focused test showed cand03/cand04 were still reactant-biased or undecided, not final 50/50 points.
- attempt024: extracted pre-commit boundary frames from attempt022/023 trajectories; found boundary03 and boundary04 as 1/1 split candidates.
- attempt025: expanded boundary03/boundary04 to 8 replicas each; boundary03 is the current best TS-candidate core.

## Key files

- `blind_work/05_atesa_acylation/attempt021_all_accepted_ts_candidates.csv`
- `blind_work/05_atesa_acylation/attempt022_committor_smoke_from_attempt021/attempt022_committor_smoke_classification.csv`
- `blind_work/05_atesa_acylation/attempt023_focused_committor_cand03_cand04/attempt023_combined_attempt022_attempt023_cand03_cand04.csv`
- `blind_work/05_atesa_acylation/attempt024_boundary_frame_candidates/attempt024_selected_boundary_frames.csv`
- `blind_work/05_atesa_acylation/attempt024_boundary_frame_committor/attempt024_boundary_committor_classification.csv`
- `blind_work/05_atesa_acylation/attempt025_focused_boundary03_boundary04/attempt025_combined_attempt024_attempt025_boundary03_boundary04.csv`

## Next step

Use boundary03 as the current TS ensemble core for a higher-replica committor refinement and eventual comparison to literature TS geometry. Do not continue from the earlier product-side C-Oleave ~2.45 A points as central TS candidates.
