# 2026-07-01 rank021 Asp-QM and angle diagnostics

## Context

This continues the blind PETase acylation mechanism workflow. The purpose was to test why rank021, despite being a better screened Michaelis-complex-like pose than rank022, still does not transfer Ser132 HG to HID209 NE2 under short His-mediated QMMM TS-seed restraints.

No paper TS coordinates were used. These are own-mechanism diagnostics based on serine hydrolase first principles.

## Starting point

Rank021 QMMM reactant gate:

`blind_work/04_qmmm_rank021_reactant_prep/hid237_qmmm_acylation_reactant_gate/qmmm_rank021_reactant_gate_stage1_10step.rst7`

Prior rank021 no-Asp tests showed:

| state | SerO-C | attack O006 | attack O004 | O006-MetN | SerH-His | SerH-SerO | O-H-N angle | His-O004 | C-O004 | nearWAT-His |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| reactant_qmmm | 3.221 | 67.5 | 106.7 | 3.332 | 1.872 | 1.012 | 166.9 | 4.946 | 1.372 | 3.589/res273 |
| own_step1 | 2.925 | 58.1 | 106.3 | 3.298 | 1.718 | 1.009 | 160.2 | 4.762 | 1.365 | 3.676/res273 |
| own_step2 | 3.155 | 66.0 | 102.4 | 3.125 | 1.551 | 1.015 | 160.8 | 4.719 | 1.443 | 3.847/res273 |

Interpretation before this diagnostic: rank021 improves pose quality over rank022 but SerO-H remains locked near 1.01 A, suggesting a missing electronic/environmental degree of freedom or an inadequate reaction coordinate.

## New tests run

### 1. Asp206 sidechain included in QM region

Input:

`qmmm_rank021_own_his_step2_with_asp_30step.in`

QM mask added Asp178/Asp206 sidechain atoms `2533-2538` to Ser132 sidechain, HID209 sidechain, and ligand. `qmcharge=-1`.

Output:

`qmmm_rank021_own_his_step2_with_asp_30step.rst7`

Status:

- `FINAL RESULTS` present.
- No `QMMM SCC-DFTB: SCC failed` string found by grep.
- Some DFTBESCF discontinuities/spikes were seen in the output, so this should be treated as a diagnostic result, not a production-quality TS seed.

### 2. Asp-QM plus O-H-N linearity

Input:

`qmmm_rank021_own_his_step3_angle_with_asp_30step.in`

Added SerO-SerH-HisNE2 angle restraint toward linear geometry.

Output:

`qmmm_rank021_own_his_step3_angle_with_asp_30step.rst7`

Status:

- `FINAL RESULTS` present.
- No `QMMM SCC-DFTB: SCC failed` string found by grep.

### 3. Asp-QM plus O-H-N linearity plus carbonyl attack angle

Input:

`qmmm_rank021_own_his_step4_attack_angle_with_asp_30step.in`

Added generic nucleophilic carbonyl attack angle restraint, still not a paper TS coordinate.

Output:

`qmmm_rank021_own_his_step4_attack_angle_with_asp_30step.rst7`

Status:

- `FINAL RESULTS` present.
- No `QMMM SCC-DFTB: SCC failed` string found by grep.

## Geometry comparison

| state | SerO-C | attack O006 | attack O004 | O006-MetN | SerH-His | SerH-SerO | O-H-N angle | His-O004 | C-O004 | HisND1-AspOD1 | HisND1-AspOD2 | nearWAT-His |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| reactant_qmmm | 3.221 | 67.5 | 106.7 | 3.332 | 1.872 | 1.012 | 166.9 | 4.946 | 1.372 | 3.080 | 2.720 | 3.589/res273 |
| step2_AspQM | 3.000 | 58.4 | 104.6 | 3.287 | 1.663 | 1.031 | 162.9 | 4.745 | 1.372 | 3.065 | 2.720 | 3.692/res273 |
| step3_angle_AspQM | 3.043 | 68.1 | 98.0 | 3.145 | 1.476 | 1.045 | 167.5 | 4.678 | 1.350 | 2.955 | 2.696 | 3.959/res273 |
| step4_attack_AspQM | 2.768 | 63.6 | 97.6 | 3.147 | 1.548 | 1.055 | 164.8 | 4.628 | 1.353 | 3.018 | 2.687 | 3.931/res273 |

## Interpretation

- Adding Asp206 to QM does not solve the Ser HG transfer problem. SerH-SerO only increases from about 1.015 A to 1.031 A.
- Adding O-H-N angle alignment improves SerH-His contact and slightly lengthens SerO-H, but not enough for proton transfer.
- Adding the carbonyl attack angle shortens SerO-C to 2.768 A, but the O006 attack angle remains poor at 63.6 degrees and SerO-H remains only 1.055 A.
- WAT273 remains outside the His acceptor site, so the present failure is not water competition.
- rank021 therefore remains a useful MC-like candidate but these short restrained minimizations have not produced a credible TS seed.

## Decision

Stop escalating simple distance/angle restraints on rank021. The repeated failure across rank022 and rank021 indicates that the current TS-seed strategy is underparameterized or using an inadequate path coordinate.

Next scientifically appropriate move:

1. Build a small 2D relaxed QMMM scan from rank021 MC over `(SerO-C)` and proton-transfer coordinate `(SerH-HisNE2 - SerH-SerO)`, while monitoring attack angle and oxyanion contact.
2. Alternatively set up a string/NEB-style path with multiple images rather than single-endpoint restrained minimizations.
3. Only after a continuous path gives a saddle-like region should umbrella sampling or TS ensemble sampling be started.

Current stage: TS seed/path probing, not umbrella sampling and not final TS ensemble.
