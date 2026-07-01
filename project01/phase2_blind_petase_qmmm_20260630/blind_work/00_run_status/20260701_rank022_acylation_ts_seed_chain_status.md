# 2026-07-01 rank022 acylation TS-seed chain status

## Scope

This status records short Amber/Sander DFTB3 QM/MM seed-generation tests from the calibrated rank022 HID237 stage3 reactant candidate. These are seed-generation tests, not an aimless-shooting ensemble and not a validated TS collection.

Starting reactant candidate:

`blind_work/04_qmmm_rank022_reactant_prep/hid237_qmmm_acylation_reactant_gate/qmmm_acylation_reactant_gate_stage3_10step.rst7`

QM region: Ser160 sidechain + HID237 sidechain + PET dimer fragment used in the rank022 model. WAT273 remains MM and is restrained away from HID NE2 for acylation.

## Seed chain summary

| step | coordinate file | SCC warning | SerO-C | C-Oleaving | SerH-SerO | SerH-Oleaving | SerH-HisNE2 | HisNE2-Oleaving | WAT273O-HisNE2 | status |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| stage3 reactant | `qmmm_acylation_reactant_gate_stage3_10step.rst7` | no | 3.142 | 1.371 | 0.973 | 3.845 | 2.045 | 4.530 | 3.703 | current reactant candidate |
| early seed | `qmmm_acylation_ts_seed_early_20step.rst7` | no | 2.593 | 1.376 | 0.989 | 3.674 | 1.955 | 4.450 | 3.743 | clean early attack seed |
| mid seed | `qmmm_acylation_ts_seed_mid_20step.rst7` | no | 2.619 | 1.410 | 0.986 | 3.249 | 2.070 | 4.267 | 3.800 | clean but limited progress |
| net seed | `qmmm_acylation_ts_seed_net_30step.rst7` | no | 2.225 | 1.484 | 1.004 | 1.894 | 2.295 | 3.622 | 3.796 | best current clean TS-seed precursor |
| late net seed | `qmmm_acylation_ts_seed_late_net_30step.rst7` | yes | 2.194 | 1.564 | 0.991 | 1.817 | 2.393 | 3.607 | 3.791 | rejected: SCC convergence warnings |

## Interpretation

The stage3 reactant can be driven along the acylation direction without SCC warnings up to the net seed:

- SerO-C formation advances from 3.142 A to 2.225 A.
- C-Oleaving elongates from 1.371 A to 1.484 A.
- Net H transfer toward the leaving oxygen begins: SerH-Oleaving decreases from 3.845 A to 1.894 A.
- WAT273 remains outside the HisNE2 acylation proton-acceptor site.

The next stronger late net seed produced SCC-DFTB convergence warnings and should not be used as a trusted TS seed. The latest reliable coordinate is therefore:

`qmmm_acylation_ts_seed_net_30step.rst7`

## Next actions

1. Continue from `qmmm_acylation_ts_seed_net_30step.rst7` using smaller increments, not the rejected late-net jump.
2. Consider increasing DFTB SCC robustness only if needed, but do not hide SCC failures.
3. Generate multiple seed variants around the net seed by varying SerO-C, C-Oleaving, and proton-transfer restraints in small increments.
4. Only after a clean set of near-TS seeds exists should aimless-shooting / committor tests be launched.
