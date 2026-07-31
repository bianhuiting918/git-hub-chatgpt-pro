# NylC-GYAQ A1 addition-first / free-proton pilot design

Date: 2026-07-31  
Branch: `codex/nylc-l4-nac-to-l2-rebalance`

## Objective

Test whether a persistent Michaelis/NAC geometry can undergo carbonyl addition before proton transfer. The pilot must induce tetrahedralization while leaving Thr267 HG1 chemically free to choose among all acceptors represented in the QM region.

This is a restrained response experiment. A passing result is not a transition state, barrier, minimum-free-energy path, or mechanism proof.

## Starting structures and provenance

Use coordinates extracted directly from the existing unbiased MM trajectories, not any restart produced by jobs 62471131, 62477118, or 62500360.

- seed26723: evt25 trajectory frame 240, time 480 ps; OG1-C12 2.972 A; O2-C12-OG1 107.16 degrees; 8 of 11 neighboring frames satisfy strict NAC.
- seed26737: evt25 trajectory frame 103, time 206 ps; OG1-C12 2.947 A; O2-C12-OG1 106.32 degrees; 7 of 11 neighboring frames satisfy strict NAC.

Before calculation, extract each frame once and record TPR/XTC provenance, frame/time, coordinate SHA256, topology SHA256, atom mapping, and periodic box. Reject the frame if the extracted Thr267 covalent graph is not intact.

## Protonation and Hamiltonian

- Start from neutral Thr267 OG1-HG1.
- Do not pre-form OG1-/NalphaH3+.
- Retain the established 146-QM, total charge 0, 510 electrons, six-link-atom DFTB3 Hamiltonian, subject to actual Amber engine-banner verification.
- HG1 and every acceptor used for proton-transfer interpretation must be in the QM region.
- Water-mediated proton transfer cannot be claimed unless the relevant water is explicitly included in the QM region. This pilot does not enlarge the QM region.

## Alternatives considered

1. Recommended: two-seed inherited addition-first pilot with one conservative force schedule. Lowest cost and directly tests the requested ordering.
2. Per-seed force grid before continuation. Safer against force mismatch but repeats the engineering-heavy calibration that previously delayed progress.
3. Fully concerted attack plus qPT restraints. Directly writes the mechanism into the bias and therefore cannot answer whether proton transfer follows tetrahedralization.

Use approach 1. If the first window lacks healthy same-direction response, stop that seed and diagnose; do not automatically increase all forces.

## Reaction-coordinate policy

Biased coordinates:

- decrease OG1-C12 progressively;
- increase C12-O2 progressively to promote carbonyl tetrahedralization.

Unbiased/monitored coordinates:

- OG1-HG1;
- HG1 to substrate N3;
- HG1 to Thr Nalpha;
- HG1 to carbonyl O2;
- HG1 to every other QM heteroatom within the defined monitoring cutoff;
- C12-N3;
- O2-C12-OG1 angle, carbonyl pyramidalization, and local OOP metrics.

No restraint may be applied to qPT, HG1, C12-N3, or a selected proton acceptor.

## Window inheritance

Run the two seeds independently. Within each seed:

1. Start window 0 from the SHA-fixed extracted MM frame.
2. Move the two biased targets by a small conservative increment.
3. Accept a window only when the engine is numerically healthy, Thr267/substrate covalent integrity passes, the biased coordinates respond in the required directions, and all chemical guards pass.
4. Only an accepted restart may seed the next window.
5. On failure, stop that chain; do not inherit the failed restart and do not silently increase force.
6. Persist source and output restart SHA256 for every window.

The initial implementation should be bounded to a short pilot chain sufficient to show response and early pyramidalization. Extending to a full tetrahedral endpoint requires a separately recorded decision based on the pilot response.

## Proton destination audit

At every accepted window, report:

- OG1-HG1 distance;
- distance and donor-H-acceptor angle for every QM acceptor;
- nearest acceptor identity;
- whether proton ownership remains OG1, becomes ambiguous, or transfers;
- the window at which any change first occurs;
- simultaneous OG1-C12, C12-O2, C12-N3, attack-angle, and OOP values.

Classify the outcome as one of:

- TETRAHEDRALIZATION_WITH_H_ON_OG1;
- TETRAHEDRALIZATION_WITH_H_TO_N3;
- TETRAHEDRALIZATION_WITH_H_TO_NALPHA;
- TETRAHEDRALIZATION_WITH_H_TO_OTHER_QM_ACCEPTOR;
- NO_TETRAHEDRAL_RESPONSE;
- NOT_EVALUATED_TECHNICAL.

## Success boundary

A seed qualifies only if it shows numerically healthy, chemically intact, inherited movement toward a tetrahedral geometry. Proton transfer is an observed response, not a required gate for the addition stage.

No result from this pilot automatically authorizes release, shooting, PMF, NEB, string, Step2 product construction, or a TS/barrier/mechanism claim.
