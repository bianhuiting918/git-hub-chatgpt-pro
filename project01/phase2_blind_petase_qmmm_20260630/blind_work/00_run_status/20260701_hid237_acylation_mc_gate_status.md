# 2026-07-01 HID237 acylation Michaelis-complex gate update

## Scope

This status records an own-structure/own-gate acylation reactant preparation attempt for rank022 after correcting His237 to HID. It does not use paper transition-state coordinates, paper barriers, paper umbrella windows, or paper product coordinates.

## Starting issue

The earlier HID237/WAT273 QMMM smoke tests kept WAT273 in the QM region. That geometry is useful for tracking an active-site water candidate, but it is not an acylation reactant model: WAT273 occupies the His237 NE2 region and competes with Ser160 HG as the proton donor to His237.

Measured from the HID237 MM/QMMM structures:

| structure | SerOG-C005 | attack angle | HisNE2-O004 | O006-MetN | SerHG-HisNE2 | WAT273O-HisNE2 |
|---|---:|---:|---:|---:|---:|---:|
| HID strong300 MM | 2.869 | 93.701 | 3.074 | 2.921 | 4.316 | 2.930 |
| WAT273-QM 10-step | 2.898 | 91.645 | 3.074 | 2.915 | 4.334 | 2.932 |

Therefore, the previous gate was missing an essential acylation MC criterion: Ser160 HG must be positioned for His237 NE2 proton acceptance.

## New acylation MC gate strategy

1. Treat WAT273 as MM, not QM, for acylation reactant preparation.
2. Add a Ser160 HG - His237 NE2 gate.
3. Add a loose WAT273 O - His237 NE2 lower-bound restraint to prevent WAT273 from occupying the acylation proton-acceptor site.
4. Keep SerOG-C005, attack-angle, and oxyanion-hole gates.
5. Do not force a paper-derived TS geometry.

## Runs completed

Work directory:

`blind_work/04_qmmm_rank022_reactant_prep/hid237_qmmm_acylation_reactant_gate/`

### Direct QMMM pull from HID strong300

Input: `qmmm_acylation_reactant_gate_20step.in`

Result: failed as a clean reactant gate. It produced SCC-DFTB convergence warning and only reduced SerHG-HisNE2 from 4.316 to 3.676 A.

### Classical MM preorganization

Input: `mm_acylation_reactant_gate_500.in`

Result: formed the SerHG-HisNE2 contact and displaced WAT273 from His237 NE2, but attack angle degraded.

| structure | SerOG-C005 | attack angle | HisNE2-O004 | O006-MetN | SerHG-HisNE2 | WAT273O-HisNE2 |
|---|---:|---:|---:|---:|---:|---:|
| before MM gate | 2.869 | 93.701 | 3.074 | 2.921 | 4.316 | 2.930 |
| after MM gate | 2.901 | 77.600 | 3.478 | 3.076 | 2.013 | 3.764 |

### QMMM smoke after MM preorganization

Input: `qmmm_acylation_reactant_gate_from_mm_10step.in`

Result: no SCC-DFTB warning; SerHG-HisNE2 contact is stable, but attack angle remains low.

| structure | SerOG-C005 | attack angle | HisNE2-O004 | O006-MetN | SerHG-HisNE2 | WAT273O-HisNE2 |
|---|---:|---:|---:|---:|---:|---:|
| after QMMM 10 | 2.958 | 74.190 | 3.441 | 3.033 | 2.029 | 3.752 |

### Stage2 angle refinement

Inputs: `mm_acylation_reactant_gate_stage2_300.in`, `qmmm_acylation_reactant_gate_stage2_10step.in`

Result: no SCC-DFTB warning in the QMMM stage2 smoke. This is the best current acylation MC candidate, but it is not yet a full-pass reactant gate.

| structure | SerOG-C005 | attack angle | HisNE2-O004 | O006-MetN | SerHG-HisNE2 | WAT273O-HisNE2 |
|---|---:|---:|---:|---:|---:|---:|
| after MM stage2 | 2.923 | 79.683 | 3.959 | 3.108 | 1.936 | 3.784 |
| after QMMM stage2 10 | 2.871 | 83.346 | 3.966 | 3.133 | 1.937 | 3.767 |

## Interpretation

Rank022 can be converted into an acylation-like Michaelis candidate with the correct His237 tautomer and SerHG-HisNE2 proton-acceptor contact, but satisfying Ser attack geometry and His/leaving oxygen geometry simultaneously remains strained.

Current candidate is suitable for further restrained refinement or for testing as an acylation TS-seed precursor, but it is not yet sufficient to claim a validated reactant state or any TS ensemble.

## Next actions

1. Do not use WAT273 as a QM atom for acylation reactant preparation; keep it as a monitored MM water unless a later deacylation AEI water-selection gate identifies it independently.
2. Add SerHG-HisNE2 and WAT-His exclusion to the formal reactant-pose gate.
3. Either continue staged refinement from `qmmm_acylation_reactant_gate_stage2_10step.rst7`, or rebuild/rerank poses using the expanded gate.
4. Only start acylation TS-seed generation after a reactant candidate passes: SerOG-C005 near 2.7-3.1 A, attack angle at least near the accepted gate, SerHG-HisNE2 around 1.7-2.2 A, oxyanion-hole contact retained, and no SCC warnings in short QMMM.
