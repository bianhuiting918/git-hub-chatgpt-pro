# Attempt029-033 acylation TS-core refinement status - 2026-07-02

Purpose: advance the refined02 acylation committor TS-core toward Stage 5/6 refinement evidence without using literature TS coordinates.

## Direct TS/NEB capability check

Attempt029 tested `sander` `imin=5` on the refined02 seed. It did not provide a usable serial TS optimizer in the current environment: the run entered trajectory-energy behavior and stopped with `File "inptraj" is missing or unreadable`. Binary inspection also showed NEB code but indicated NEB requires MPI, while the active environment has serial `sander` and no `sander.MPI`.

Decision: do not keep guessing `imin=5` parameters. Use the protocol-allowed constrained saddle-like refinement plus dynamical committor evidence.

## Constrained relaxation tests

Attempt030: weak 100-cycle constrained relaxation ran, but slid product-side: C-Oleave moved from 1.847 A to 2.249 A.

Attempt031: stronger 50-cycle constrained relaxation also slid product-side: SerO-C 2.332 A and C-Oleave 2.164 A. This shows longer minimization is downhill away from the dividing surface.

Attempt032: 5-cycle strong micro-relax preserved the TS-core geometry:

- SerO-C: 2.055 A
- C-Oleave: 1.848 A
- SerH-SerO: 1.119 A
- SerH-His: 1.451 A
- restraint penalty: 0.022

## Micro-relaxed committor validation

Attempt033 launched 8 unbiased DFTB3/MM short trajectories from the attempt032 micro-relaxed structure.

First-hit counts:

- tetint_basin: 3
- reactant_basin: 4
- none: 1
- valid directional replicas: 7
- tetint fraction excluding none: 0.429

Interpretation: the micro-relaxed refined02 structure retains bidirectional committor behavior and is the best current acylation TS-core refinement artifact. It is not a fully optimized one-imaginary-frequency TS because the available serial Amber setup does not support direct TS/NEB optimization. It is a defensible constrained/dynamical TS-core structure for the current reproduction stage and should seed later higher-level TS optimization if MPI-NEB/geomeTRIC/ChemShell/CP2K tools are enabled.

Key directories:

- `blind_work/05_ts_refinement/attempt029_refined02_imin5_probe/`
- `blind_work/05_ts_refinement/attempt032_refined02_micro_relax/`
- `blind_work/05_ts_refinement/attempt033_micro_relax_committor/`
- `blind_work/06_ts_ensemble/acylation/refined02_committor_ts_core_20260702/`
