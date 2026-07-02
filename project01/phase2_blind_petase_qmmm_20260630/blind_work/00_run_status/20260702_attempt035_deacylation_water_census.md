# Attempt035 deacylation water census - 2026-07-02

Purpose: identify candidate hydrolytic water for the next PETase deacylation/hydrolysis stage from our own acylation/tetrahedral-intermediate trajectories, not from literature coordinates.

Input trajectories:

- `blind_work/05_ts_refinement/attempt033_micro_relax_committor/`
- `blind_work/05_atesa_acylation/attempt028_refined02_extended_committor/`

Method:

- Inspect tetint-committing frames and final tetint frames.
- Rank waters by proximity to acyl carbon, His NE2, Ocarb, and SerO.
- Repeat with H-aware scoring that includes closest water-H to His NE2.

Result:

- Both O-only and H-aware censuses identify the same recurring active-site water: `HOH272`.
- Amber 1-based atoms: O `3863`, H1 `3864`, H2 `3865`.
- Best O-only frame: source `attempt028_refined02`, water O-C `3.26` A, water O-His `3.99` A, water O-Ocarb `3.09` A.
- Best H-aware frame: source `attempt028_refined02`, water O-C `3.69` A, closest water-H to His `3.44` A.

Interpretation:

`HOH272` is the first-principles candidate hydrolytic water because it is the only water repeatedly ranked near the reactive acyl region in our own tetint/acylation trajectories. Its orientation is not yet fully reactive: the closest water-H to His remains >3.4 A, so the next calculation should not claim a deacylation TS. The correct next step is to build a deacylation reactant/preorganization branch with HOH272 included in the QM region and sample/relax water orientation before TS search.

Generated files:

- `blind_work/04_qmmm_exploration/deacylation_water_census_attempt035/attempt035_tetint_water_census_top10_per_frame.csv`
- `blind_work/04_qmmm_exploration/deacylation_water_census_attempt035/attempt035_water_candidate_aggregate.csv`
- `blind_work/04_qmmm_exploration/deacylation_water_census_attempt035/attempt035_haware_water_census.csv`
- `blind_work/04_qmmm_exploration/deacylation_water_census_attempt035/attempt035_haware_water_candidate_aggregate.csv`
- `blind_work/04_qmmm_exploration/deacylation_water_census_attempt035/attempt035_deacylation_qm_region_with_h2o272.txt`
