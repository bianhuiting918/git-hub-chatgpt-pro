# 2026-07-01 HID237 + WAT273 segmented QMMM protocol check

Purpose: determine how far the clean HID237 + WAT273 DFTB3/MM minimization can be advanced before SCC or reaction-geometry failure.

Starting point:

- Topology: `hid237_topology_prep/TOPOLOGY_AMBER_TOPO_HID237_RANK022/complex.prmtop`
- Starting coordinates for QMMM: HID237 strong300 MM geometry
- Clean checkpoint: `hid237_qmmm_wat273_smoke/qmmm_hid_ser_his_lig_wat273_10step.rst7`
- QM region: Ser160 sidechain + HID237 sidechain + full LIG + WAT273
- QM charge/spin: `qmcharge=0`, `spin=1`

## Segment-length test from the clean 10-step checkpoint

| run | SCC warning | SerOG-C005 A | attack angle deg | O004-HID NE2 A | O006-MetN A | WAT O-HID NE2 A | WAT O-C005 A | WAT O-O006 A | WAT O-O004 A |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 10-step | no | 2.898 | 91.645 | 3.074 | 2.915 | 2.932 | 3.690 | 3.540 | 2.907 |
| 20-from-10 | no | 3.004 | 81.488 | 3.097 | 2.891 | 2.909 | 3.737 | 3.658 | 2.798 |
| 30-from-10 | no | 2.874 | 81.807 | 3.337 | 2.850 | 2.889 | 3.677 | 3.824 | 2.683 |
| 40-from-10 | no | 2.891 | 65.738 | 4.123 | 2.893 | 2.846 | 3.550 | 3.784 | 3.126 |
| 50-from-10 | yes | 2.848 | 69.598 | 3.976 | 2.865 | 2.848 | 3.548 | 3.904 | 3.017 |

## Interpretation

- 20- and 30-step segments from the clean 10-step checkpoint retain acceptable Michaelis-like geometry and have no SCC warning.
- 40-step remains electronically converged but loses the His/leaving-group and attack-angle gate.
- 50-step both drifts geometrically and produces SCC-DFTB non-convergence warnings.

Operational rule for the next stage:

- Use 10-step QMMM as the clean checkpoint.
- Advance QMMM minimization in <=20-30 step segments.
- After each segment, verify both SCC warning absence and reaction geometry before continuing.
- Do not run long unmonitored 50+ step QMMM minimizations from this state.

This protocol is still reactant-preparation only. It is not a TS search and does not use paper coordinates or paper numeric results.
