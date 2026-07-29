# NylC A1 Step1 Adaptive Inherited Gap-Fill Design

## Objective

Continue only the eight incomplete intermediate routes from job 62267272:

- two endpoint seeds: seed26723 and seed26737;
- two intermediate anchors per seed: I1 and I2;
- two directions per anchor: REACTANT and PRODUCT.

The calculation fills the first missing interval in each inherited chain. It does not recompute R/P controls, replace the frozen Step1 Hamiltonian, or claim a TS, PMF, barrier, or mechanism.

## Frozen authority

- Source job: 62267272, output root `a1_activated_nac_20260726/qmmm/a1_step1_four_anchor_bidirectional_gapfill`.
- Relevant source task indices: 2, 3, 4, 5, 10, 11, 12, 13.
- Every route has 4 technically complete inherited stages and stopped at the next scientific guard.
- Step1 Hamiltonian remains 146 QM atoms, q0, 510 electrons including 6 link H, 6 link atoms, no QM water, DFTB3/3OB-3-1, Amber18 sander.MPI.
- Frozen prmtop SHA256: `a61d15bf0bf78675be93275d45f274e808ed6ae450fc1ca21a8e14aee8c12ca0`.
- Gate loop remains residues 261-266; Thr267 is excluded.
- Coordinates from job 62021985 are forbidden.

## Selected approach

Use adaptive midpoint continuation. One Slurm array task represents one of the eight routes. Within a task, windows run serially and each accepted restart is the sole input to the next window.

For each route:

1. Resolve the final persisted technically complete and guard-safe restart from the source manifest, verify its recorded SHA256, and bind the original next planned joint-coordinate target.
2. Define the interval between the last accepted joint target and the first rejected/planned target.
3. Attempt the midpoint target while moving the same coupled coordinate set used by job 62267272: attack C12-OG1, cleavage C12-N3, carbonyl C12-O2/OOP response, and proton-transfer coordinates.
4. If the midpoint is technical PASS and all inherited chemistry/direction guards pass, accept its restart and advance the lower bound toward the destination.
5. If it is technical PASS but a scientific guard fails, keep the last accepted restart and move the upper bound back to the failed midpoint.
6. Repeat until the destination gate is reached, or the bracket width is at most 1/16 of the original failed interval, or eight new windows have been attempted.
7. A failed or partial window never becomes the parent of a later window.

The midpoint rule is deterministic and recorded in every stage manifest. There is no manual target substitution and no cross-route inheritance.

## Parallelism and persistence

- Slurm array: 0-7 without throttle.
- Each task: one node, 8 MPI ranks.
- The eight routes may run concurrently; windows inside each route are serial.
- Full engine scratch remains temporary.
- Persist only the route manifest, accepted restart chain, compact per-window geometry/guards, RESULT/PASS/NOT_EVALUATED, SHA256, RUNBOOK, and run_history.
- New output roots are immutable and keyed by the new array job ID; job 62267272 is read-only.

## Gates

Technical completion and scientific interpretation remain separate.

- `PASS_TECHNICAL_ADAPTIVE_GAPFILL`: all attempted windows have complete Amber output, valid restart/SHA inheritance, the frozen engine banner, finite geometry, and no hard/SCC/vlimit/overflow errors.
- `PASS_GUIDED_ROUTE_REACHED_DESTINATION`: the inherited route reaches its predefined basin with every chemistry and direction guard satisfied.
- `PARTIAL_BRACKET_REFINED`: technical PASS and the original failed interval was reduced to at most 1/16, but the destination was not reached.
- `NOT_EVALUATED_TECHNICAL`: missing source restart/SHA authority or runtime failure.

A restrained endpoint or refined bracket is not a TS candidate. Only after a continuous guard-safe bridge exists may central frames be subjected to a separate fully unrestrained release. Paired shooting is not launched in this stage.

## Alternatives rejected

- Fixed dense windows: simpler but spends equal work in easy and failing regions.
- Immediate unrestrained release from current I1/I2 frames: likely collapses to one basin and does not fill the known missing interval.
