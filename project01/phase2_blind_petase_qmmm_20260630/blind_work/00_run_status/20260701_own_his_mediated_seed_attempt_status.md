# 2026-07-01 own His-mediated acylation seed attempt

## Context

User clarified the workflow: do not follow the paper TS/commitment coordinates as the main route. The correct route is:

1. Try our own mechanistic hypothesis from the structure.
2. Compare with paper only afterward.
3. If different, identify the reason and retry.

Therefore, the prior paper-style net proton-to-leavingO seed is treated only as a comparison/control branch, not as the main reproduction route.

## Own hypothesis tested

Acylation begins with His237 NE2 accepting Ser160 HG while Ser160 OG approaches the PET carbonyl carbon. The leaving oxygen is monitored but not forced.

Starting coordinate:

`blind_work/04_qmmm_rank022_reactant_prep/hid237_qmmm_acylation_reactant_gate/qmmm_acylation_reactant_gate_stage3_10step.rst7`

## Runs

Work directory:

`blind_work/04_qmmm_rank022_reactant_prep/hid237_qmmm_acylation_reactant_gate/`

| step | coordinate | SCC warning | SerO-C | C-Oleaving | SerH-SerO | SerH-HisNE2 | SerH-Oleaving | HisNE2-Oleaving | interpretation |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| stage3 reactant | `qmmm_acylation_reactant_gate_stage3_10step.rst7` | no | 3.142 | 1.371 | 0.973 | 2.045 | 3.845 | 4.530 | reactant-like candidate |
| own His step1 | `qmmm_own_his_mediated_ts_seed_step1_30step.rst7` | no | 2.751 | 1.376 | 1.012 | 1.636 | 3.650 | 4.422 | His approaches H, but SerO-H does not break |
| own His step2 | `qmmm_own_his_mediated_ts_seed_step2_30step.rst7` | no | 2.874 | 1.366 | 1.000 | 1.686 | 3.418 | 4.388 | stronger SerO-H lengthening still fails; SerO-C relaxes away |

## Interpretation

The own His-mediated seed did not reach a productive proton-abstraction/attack geometry. It remains electronically stable, but the proton stays bonded to SerO and the SerO-C attack coordinate does not continue forward under the current restraints.

Likely causes to test next:

1. The His/Asp electrostatic environment is not represented strongly enough because Asp206 is outside the QM region.
2. The His237 side-chain orientation is still not suitable for simultaneous SerH abstraction and SerO-C approach.
3. Rank022 may be a reactant-like pose but not a good acylation TS precursor.
4. The restraint scheme tries to move distances but does not enforce the angular alignment of O-H-N proton transfer.

## Next retry options

1. Add an O-SerH-HisNE2 angle gate for the His-mediated branch, still without paper TS coordinates.
2. Test a QM region including Asp206 with a carefully defined charge/spin setup.
3. Re-rank poses using the corrected own gate: SerH-His preorganization, WAT-His exclusion, non-overcompressed SerO-C, and proton-transfer angle.
4. Only after our own branch is understood should we compare its geometry against the paper TS/AS pathway.
