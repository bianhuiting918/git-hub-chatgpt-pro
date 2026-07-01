# 2026-07-01 reference initial-geometry calibration and rank022 stage3 reactant gate

## Purpose

The paper's public acylation dataset was downloaded only for geometry calibration and comparison. The coordinates are not used as inputs for the blind reproduction calculations.

Zenodo record: `10.5281/zenodo.10854763`

Downloaded file: `1.acylation.zip`

Zenodo MD5: `0ac86bc4b2f4148ae7480aa3bb5c98bc`

Verified local MD5: `0AC86BC4B2F4148AE7480AA3BB5C98BC`

## Paper acylation files identified

From the dataset README:

- `vmd-md-b.inpcrd`: initial reactant state coordinates in Amber format.
- `minimization.rst7`: output Amber coordinate file for QM/MM minimization of the reactant state.
- `minimization.rst7_0_ts_guess_960.0.rst7`: putative transition-state coordinate produced by ATESA during `find_ts`.

From `find_ts.config`:

- Ser O: atom 1911, SER132 OG.
- Ser H: atom 1912, SER132 HG1.
- His N: atom 2991, HSD209 NE2.
- PET reactive carbon: atom 3862, UNK266 C30.
- PET reactive/leaving oxygen: atom 3848, UNK266 O16.
- Closest carbonyl oxygen candidate to C3862: atom 3846, UNK266 O14.

## Paper reference geometry

| structure | SerO-C | SerO-C-Oleaving | SerO-C-Ocarbonyl | C-Oleaving | SerH-HisN | HisN-Oleaving | SerH-Oleaving |
|---|---:|---:|---:|---:|---:|---:|---:|
| paper initial `vmd-md-b.inpcrd` | 3.676 | 142.507 | 41.783 | 1.485 | 2.256 | 4.042 | 4.298 |
| paper QM/MM minimized reactant `minimization.rst7` | 3.739 | 142.651 | 36.832 | 1.331 | 2.929 | 4.154 | 3.889 |
| paper ATESA TS guess | 1.939 | 75.553 | 107.128 | 1.900 | 2.382 | 3.186 | 1.195 |

## Gate correction

The measurement shows that the acylation reactant in the paper does not require His237 NE2 to be close to the PET leaving oxygen. `HisN-Oleaving` is about 4.0 A in the paper initial/minimized reactant. Therefore, `His-leavingO` should be monitored, not used as a strong Michaelis-complex gate.

The useful reactant-gate constraints are:

1. Keep SerHG-HisNE2 preorganized for Ser proton abstraction.
2. Keep WAT273 out of the HisNE2 acylation proton-acceptor site.
3. Keep SerO-C in a nonbonded reactant-like range, not over-compressed toward TS.
4. Keep oxyanion-hole contact as a loose stabilizing criterion.
5. Leave His-leavingO approach to the TS/path-search stage.

## rank022 HID237 stage3 test

Work directory:

`blind_work/04_qmmm_rank022_reactant_prep/hid237_qmmm_acylation_reactant_gate/`

Stage3 inputs:

- `acylation_reactant_gate_stage3_reactant_like.RST`
- `mm_acylation_reactant_gate_stage3_300.in`
- `qmmm_acylation_reactant_gate_stage3_10step.in`

Stage3 result:

- The MM stage completed with `FINAL RESULTS`.
- The QMMM 10-step smoke completed with `FINAL RESULTS` and no `QMMM SCC-DFTB` warning.

| structure | SerO-C | SerO-C-Ocarbonyl | SerO-C-Oleaving | HisNE2-Oleaving | Ocarbonyl-MetN | SerH-HisNE2 | WAT273O-HisNE2 | WAT273O-C | SerH-Oleaving |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| stage2 QMMM input | 2.871 | 83.346 | 97.323 | 3.966 | 3.133 | 1.937 | 3.767 | 3.829 | 3.241 |
| stage3 MM | 3.142 | 73.314 | 116.889 | 4.530 | 3.289 | 2.054 | 3.735 | 3.612 | 3.849 |
| stage3 QMMM 10-step | 3.142 | 73.392 | 116.773 | 4.530 | 3.339 | 2.045 | 3.703 | 3.656 | 3.845 |

## Interpretation

The rank022 stage3 structure is now a cleaner reactant-like Michaelis candidate than the over-compressed stage2 structure:

- SerO-C is nonbonded and no longer TS-like.
- SerH-HisNE2 remains preorganized.
- WAT273 remains displaced from HisNE2 and should stay MM for acylation.
- HisNE2-Oleaving is far, which is consistent with the paper initial/minimized reactant being far.
- The geometry is electronically stable in a short Amber/Sander DFTB3 QM/MM smoke test.

This is not a TS ensemble. It is a calibrated acylation reactant candidate for subsequent TS seed/path generation.

## Next calculation step

Use the stage3 QMMM structure as the current best rank022 acylation reactant candidate:

`blind_work/04_qmmm_rank022_reactant_prep/hid237_qmmm_acylation_reactant_gate/qmmm_acylation_reactant_gate_stage3_10step.rst7`

Next, generate controlled acylation TS seeds by scanning/steering along first-principles reaction coordinates:

- SerO-C formation.
- SerH transfer from SerO to HisNE2.
- C-Oleaving elongation and later His-to-leavingO proton relay.

Do not force HisNE2-Oleaving in the reactant gate.
