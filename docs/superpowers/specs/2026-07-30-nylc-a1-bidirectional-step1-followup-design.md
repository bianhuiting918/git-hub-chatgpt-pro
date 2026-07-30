# NylC A1 Bidirectional Step1 Follow-up Design

Status: user-approved on 2026-07-30.

## Objective

Run two calculations in parallel without changing the frozen Step1 Hamiltonian:

1. Continue the successful product-to-tetra response from reverse-force tasks 6 and 8 through inherited windows.
2. Restore the missing reactant-to-product direction with matched addition-first and fully concerted first-window calibration.

## Frozen authority

- Amber topology SHA256: `a61d15bf0bf78675be93275d45f274e808ed6ae450fc1ca21a8e14aee8c12ca0`.
- Step1 QM contract: 146 QM atoms, charge 0, 510 electrons including six link H, six links, no QM water, DFTB3/3OB-3-1, Amber18.
- Computation runs only on SCNet.
- A technical or response PASS is not a TS, PMF, barrier, or mechanism PASS.

## Reverse continuation

Two independent chains, both from seed26737:

- LOW_BIAS source: task6 `stage.rst7`, SHA256 `54ecabec3cbb993a81b36ce86630355e2143797481d71726ca83edd4d65862df`, force scale 1.
- TARGET_CLOSE source: task8 `stage.rst7`, SHA256 `2f30ffd0746cc7c53f36ac102b446256ea575e24a2d5bfa9da15a1455db96a7f`, force scale 4.

Each array task runs at most eight serial windows. Every accepted window supplies the exact restart for the next window. A failed or guard-rejected restart is never inherited. Each window continues C12-N3 shortening, reverse qPT, and C12-O2 lengthening while OG1-C12 remains unrestrained. No release or shooting is automatic.

## Forward calibration

Twelve independent first-window tasks: two reactant seeds x two mechanisms x force scales 1, 2, and 4.

Frozen reactant sources:

- seed26723 SHA256 `9d93b60aee9e8d6fa97493396d757f3bf4d0a47595591968e3debaebb2640e86`.
- seed26737 SHA256 `a12c018e195c32b34239004aea36ebc80de8275c9f9332b44084ad8d88d4db0e`.

Mechanisms:

- `ADDITION_FIRST_FORWARD`: shorten OG1-C12 and lengthen C12-O2. qPT and C12-N3 are unrestrained and monitored, directly testing whether proton transfer responds spontaneously.
- `FULLY_CONCERTED_FORWARD`: shorten OG1-C12, lengthen C12-O2 and C12-N3, and drive qPT from Nalpha toward N3.

The forward job stops after the first window. Cross-audit selects the weakest passing scale per seed and mechanism; it does not launch inherited windows.

## Outputs and gates

Each task persists only source/target/actual geometry, stage restart, compact engine tail, manifest, SHA table, and RESULT/PASS or NOT_EVALUATED. Reverse chains report every window's SHA inheritance. Forward tasks report active versus monitored coordinates. Chemical guards remain mandatory. No old output is overwritten or reclassified.
