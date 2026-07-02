# Deacylation attempt053: bridge committor from attempt041_w04 toward w239 (2026-07-02)

Purpose: close the blind deacylation TS-core bracket between the early/reactant-side attempt041_w04 and the product-emergent w239 candidate. This attempt does not use literature geometry; it starts from our own attempt041_w04 seed and applies small, explicit restraint shifts toward the observed boundary.

CPU policy: user limit is total CPU <=64. The launched committor batch uses 24 Amber/sander jobs with OMP_NUM_THREADS=2, total process threads = 48 by ps NLWP check.

## Starting candidates after restrained DFTB3 minimization

| candidate | OW-C A | H2-His A | C-SerO A | C-Oleave A | SerH-His A |
|---|---:|---:|---:|---:|---:|
| bridge_cser240_h145 | 2.096 | 1.477 | 2.339 | 3.287 | 1.669 |
| bridge_cser242_h148 | 2.063 | 1.499 | 2.366 | 3.323 | 1.750 |
| bridge_cser244_h150 | 2.157 | 1.552 | 2.364 | 3.343 | 1.719 |

## Committor setup

- 3 candidates x 8 random-velocity replicas = 24 trajectories.
- Each trajectory: Amber sander DFTB3 QM/MM, 3000 steps, dt=0.0001 ps, 310 K Langevin.
- Launched directory: `blind_work/05_ts_refinement/attempt053_deacylation_w04_w239_bridge_committor`.
- Initial launch verification: 24 live `sander` jobs, 48 total NLWP threads.

## Pending

Classify product/reactant/undecided after all 24 trajectories finish, then compare to attempt052 lower bracket and attempt050 w239 product-biased probe.

## Results

All 24/24 trajectories completed. Post-run check found 0 live `sander` processes. Fatal/SCF/NaN/segmentation scan returned no failure matches; generic `ERROR` only matched Amber's standard Ewald relative-error banner.

| candidate | product | reactant | undecided | interpretation |
|---|---:|---:|---:|---|
| bridge_cser240_h145 | 2 | 3 | 3 | best current deacylation TS-core candidate; product fraction excluding undecided = 0.40 |
| bridge_cser242_h148 | 5 | 2 | 1 | product-shifted upper bracket; product fraction excluding undecided = 0.71 |
| bridge_cser244_h150 | 2 | 4 | 2 | early/reactant-shifted bracket; product fraction excluding undecided = 0.33 |

Interpretation: attempt053 successfully brackets the deacylation dividing surface from our own blind seed. `bridge_cser240_h145` is the closest balanced candidate so far and should be used as the next seed for micro-relax/expanded committor ensemble.

## Attempt054 extension of bridge_cser240_h145

Purpose: extend the best attempt053 candidate from 8 to 16 random-velocity committor trajectories. CPU use: 8 concurrent Amber/sander jobs with OMP_NUM_THREADS=2, total process threads = 16 by ps NLWP check.

All 8 added trajectories completed; post-run check found 0 live `sander` processes. Fatal/SCF/NaN/segmentation scan returned no failure matches.

Combined cser240 result across 16 x 3000-step trajectories:

| candidate | product | reactant | undecided | product fraction excluding undecided | interpretation |
|---|---:|---:|---:|---:|---|
| bridge_cser240_h145 | 6 | 4 | 6 | 0.60 | current best deacylation TS-core candidate; close enough for ensemble expansion/micro-relax rather than further broad bracketing |

Next action: generate a small local ensemble around bridge_cser240_h145, then compare with literature only after our blind candidate set is internally consistent.
