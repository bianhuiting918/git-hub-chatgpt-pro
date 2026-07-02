# Attempt034 acylation TS-core ensemble expansion - 2026-07-02

Purpose: expand the single refined02 acylation TS-core into a small blind TS-core ensemble using only our own micro-relaxed committor trajectories.

Input evidence: attempt033 micro-relaxed refined02 committor gave tetint 3 / reactant 4 / none 1, so the micro-relaxed core remains close to the dividing surface.

Neighbor frame extraction:

- Extracted pre-commit early-region frames from attempt033 trajectories.
- Selected 5 diverse neighboring frames: 2 from tetint-first trajectories, 2 from reactant-first trajectories, and 1 from the none/undecided trajectory.
- Selection did not use literature TS geometries.

Micro-relaxation:

- Ran 5-cycle DFTB3/MM constrained micro-relax for all 5 selected neighbors.
- Initial long filenames failed because Amber truncated `DISANG=` paths; reran with short n00-n04 file names.
- Short-name rerun completed 5/5 with no serious errors.
- All 5 retained TS-core geometry after micro-relax.

Kept ensemble members: 5

Evidence files:

- `blind_work/06_ts_ensemble/acylation/refined02_committor_ts_core_20260702/attempt034_selected_neighbor_boundary_frames.csv`
- `blind_work/06_ts_ensemble/acylation/refined02_committor_ts_core_20260702/attempt034_neighbor_micro_relax_geometry_summary.csv`
- `blind_work/06_ts_ensemble/acylation/refined02_committor_ts_core_20260702/attempt034_relaxed_neighbor_members/attempt034_acylation_ts_core_members.csv`

Interpretation: acylation now has one central micro-relaxed TS-core plus five neighboring micro-relaxed TS-core ensemble members. The current Amber environment still lacks direct MPI-NEB/TS optimizer support, so these are constrained/dynamical TS-core ensemble structures, not one-imaginary-frequency optimized saddle points.
