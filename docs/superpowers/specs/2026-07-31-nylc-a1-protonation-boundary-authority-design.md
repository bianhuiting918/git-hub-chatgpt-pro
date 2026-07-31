# NylC A1 protonation and QM-boundary authority design

## Objective

Separate two explanations for the reaction-coordinate-free Thr267 covalent-integrity failure:

1. the preformed neutral-overall zwitterion Thr267 NalphaH3+/OG1- is unstable under the current DFTB3/MM Hamiltonian;
2. the current QM/MM cut across the Thr267 C--Thr268 N peptide bond destabilizes Thr267.

## Controlled comparison

Reuse the same two PBC-reconstructed heavy-atom source frames already frozen for seeds 26723 and 26737.

- Direction A, NEUTRAL_THR267: use NalphaH2/OG1H, retain the established 146-QM-atom core, qmcharge 0, 510 electrons, and six link atoms.
- Direction B, ACTIVATED_EXPANDED_THR268: retain NalphaH3+/OG1-, add complete Thr268 (global residue 623, atoms 8964-8977) to the QM region. The expected contract is 160 QM atoms, qmcharge 0, 564 electrons, and six link atoms. The peptide-boundary link moves from Thr267 C8962--Thr268 N8964 to Thr268 C8976--Ile269 N8978.

This is one four-task array: two seeds times two directions. The previous activated-state/original-boundary result is the unchanged baseline and is not rerun.

## Calculation

Each task performs four strictly inherited, reaction-coordinate-free and position-restraint-free minimization blocks of 250 cycles. No attack, carbonyl, proton-transfer, or angle restraint is present. A failed block restart is never inherited.

For the neutral direction, HG1 is deterministically placed on OG1 while preserving the source heavy-atom coordinates. The topology bond graph is changed only from Nalpha-HG1 to OG1-HG1. Because all Thr267 atoms are QM, the electronic state remains governed by qmcharge 0 and the QM geometry; the topology patch supplies correct fragment reconstruction and audit connectivity.

## Gates

Technical interpretation requires:

- exact source and topology SHA lineage;
- no PBC-split covalent bond;
- expected engine banner for QM atoms, charge, electrons, and link atoms;
- finite numerical output without Amber/DFTB overflow;
- intact Thr267 Nalpha-CA, CA-CB, and CB-OG1 heavy-atom skeleton;
- intact Thr267--Thr268 peptide bond; for the expanded direction, intact Thr268 skeleton as well;
- strict block-to-block restart SHA inheritance.

HG1 is classified by its final distances to OG1, Nalpha, and substrate N3. A chemically plausible proton transfer is reported as a state change, not mislabeled as technical failure.

No reaction window, release, shooting, PMF, NEB, string, or mechanism conclusion is started automatically.
