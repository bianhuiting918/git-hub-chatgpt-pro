# 2026-07-01 HID237 + WAT273 QMMM gate status

Purpose: continue blind PETase reproduction from structure-derived rank022 candidate, correcting catalytic His protonation and identifying the important active-site water from our own geometry rather than using literature water identifiers.

## Solvation and active-site water

The Amber system is explicitly solvated. Sander previously reported 13,706 three-point waters. In rank022 weakback geometry, only one water oxygen simultaneously sits near catalytic His and the ligand reactive region:

- WAT273, oxygen atom serial 3863.

In the corrected HID237 topology, WAT273 persists as the only near-bridge candidate in strong300, weakback200, and release100 geometries:

| stage | WAT O-HID NE2 A | WAT O-C005 A | WAT O-O006 A | WAT O-O004 A | WAT O-SerOG A |
|---|---:|---:|---:|---:|---:|
| HID strong300 | 2.930 | 3.678 | 3.531 | 2.919 | 4.137 |
| HID weakback200 | 2.921 | 3.664 | 3.614 | 2.829 | 4.086 |
| HID release100 | 2.859 | 3.648 | 3.314 | 3.397 | 4.102 |

Interpretation: WAT273 is our structure-derived important-water candidate. It is not imported from the paper.

## Catalytic His correction

The original rank022 topology had His237 converted by tleap to HIE209, with HE2 protonated. For the acylation/reactant-base geometry, His should accept the Ser proton through NE2, so a HID237 topology was rebuilt:

- source PDB: `blind_work/04_qmmm_rank022_reactant_prep/hid237_rebuild/complex_for_amber_hid237.pdb`
- topology prep: `blind_work/04_qmmm_rank022_reactant_prep/hid237_topology_prep/TOPOLOGY_AMBER_TOPO_HID237_RANK022/`
- HID209 contains HD1 and neutral NE2.

## HID237 MM gates

HID topology reaction geometry from MM gates:

| stage | SerOG-C005 A | attack angle deg | O004-HID NE2 A | O006-MetN A | pass |
|---|---:|---:|---:|---:|---|
| weak100 | 2.903 | 60.550 | 4.335 | 2.919 | no |
| strong300 | 2.869 | 93.701 | 3.073 | 2.921 | yes |
| weakback200 | 2.958 | 72.746 | 3.152 | 2.927 | no |
| release100 | 3.361 | 52.596 | 4.399 | 3.324 | no |

Interpretation: with chemically correct HID237, rank022 can be brought into Michaelis geometry under stronger generic restraints, but does not stay in the gate under weakback/release. This is a useful reactant-prep candidate, not a final unrestrained MC ensemble.

## QMMM smoke with WAT273

A neutral 39-atom QMMM region was tested from HID237 strong300 coordinates:

- QM atoms: Ser160 sidechain, HID237 sidechain, full LIG, WAT273.
- MM atoms: Asp206 remains MM for this smoke to avoid the odd-electron/charge problem seen when including deprotonated Asp in the QM region.
- `qmcharge=0`, `spin=1`, `qm_theory='DFTB3'`.
- work dir: `blind_work/04_qmmm_rank022_reactant_prep/hid237_qmmm_wat273_smoke/`

Results:

| QMMM gate | exit | SCC warning | SerOG-C005 A | attack angle deg | O004-HID NE2 A | WAT O-HID NE2 A | WAT O-O004 A |
|---|---:|---|---:|---:|---:|---:|---:|
| 1-step | 0 | no | not tabulated | not tabulated | not tabulated | not tabulated | not tabulated |
| 10-step | 0 | no | 2.898 | 91.645 | 3.074 | 2.932 | 2.907 |
| 50-step unrestrained | 0 | yes | 2.848 | 69.598 | 3.976 | 2.848 | 3.017 |
| 50-step weak restrained | 0 | yes | 2.613 | 81.082 | 3.465 | 2.841 | 2.843 |

Interpretation: the 10-step HID237 + WAT273 QMMM smoke is currently the cleanest QMMM entry point: no SCC warning, WAT273 retained, and Michaelis geometry retained. Longer 50-step minimization begins to show SCC convergence warnings, so production QMMM needs smaller staged minimization/equilibration settings or improved QM-region/charge treatment before TS search.

## Next action

Do not proceed to TS search from the 50-step outputs. Use the 10-step clean QMMM output as a diagnostic checkpoint, then develop a staged QMMM protocol: shorter segments, possible geometry restraints, and a separate treatment of Asp206 charge/protonation if it must enter the QM region.
