# NylC-GYAQ A1 addition-first / free-proton pilot design

Date: 2026-07-31  
Branch: `codex/nylc-l4-nac-to-l2-rebalance`

## Objective

Test whether a persistent Michaelis/NAC geometry can undergo carbonyl addition before proton transfer. The pilot induces early tetrahedralization while leaving Thr267 HG1 chemically free to choose among all acceptors represented in the QM region.

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

- OG1-C12: harmonic flat-bottom center moves by -0.04 A per accepted window;
- C12-O2: harmonic flat-bottom center moves by +0.03 A per accepted window.

Use the established base restraint constants without scale multiplication:

- OG1-C12: 24.0 kcal mol-1 A-2;
- C12-O2: 30.0 kcal mol-1 A-2.

Do not apply attack-angle or proton-position protection restraints in this pilot. The stable MM starting frames supply the initial geometry; angle loss is an observed failure mode.

Unbiased/monitored coordinates:

- OG1-HG1;
- HG1 to substrate N3;
- HG1 to Thr Nalpha;
- HG1 to carbonyl O2;
- HG1 to every QM N/O/S acceptor, retaining the ranked list for acceptors within 4.0 A;
- C12-N3;
- O2-C12-OG1 angle, carbonyl pyramidalization, and local OOP metrics.

No restraint may be applied to qPT, HG1, C12-N3, a proton acceptor, or the attack angle.

## Pilot size and window inheritance

Run two independent array tasks, one per seed, with 8 MPI ranks per task. Each task is bounded to six sequential minimization windows. The six target increments are defined relative to the SHA-fixed starting geometry: window n uses n times -0.04 A for OG1-C12 and n times +0.03 A for C12-O2, for n = 1 through 6.

Within each seed:

1. Start window 1 from the SHA-fixed extracted MM frame.
2. Accept a window only when the Amber engine is numerically healthy, the 146-QM/q0/510e/6link banner is present, Thr267 and substrate covalent integrity pass, both biased coordinates respond in the required directions, and all chemical guards pass.
3. Only an accepted restart may seed the next window.
4. On failure, stop that chain; do not inherit the failed restart and do not silently change force or target step.
5. Persist source and output restart SHA256, input restraints, engine tail, RMS/GMAX/energy, and compact geometry audit for every attempted window.

The six-window pilot tests response and early pyramidalization. It does not claim a completed tetrahedral intermediate. Extending the same accepted chain toward a full tetrahedral endpoint requires review of these six windows and a separately recorded continuation decision.

## Proton destination audit

At the starting frame and every attempted window, report:

- OG1-HG1 distance;
- distance and donor-H-acceptor angle for every QM N/O/S acceptor;
- ranked acceptors within 4.0 A and the nearest acceptor identity;
- whether proton ownership remains OG1, becomes ambiguous, or transfers;
- the first window at which any ownership change occurs;
- simultaneous OG1-C12, C12-O2, C12-N3, O2-C12-OG1 angle, pyramidalization, and OOP values.

Classify the outcome as one of:

- EARLY_TETRAHEDRALIZATION_WITH_H_ON_OG1;
- EARLY_TETRAHEDRALIZATION_WITH_H_TO_N3;
- EARLY_TETRAHEDRALIZATION_WITH_H_TO_NALPHA;
- EARLY_TETRAHEDRALIZATION_WITH_H_TO_OTHER_QM_ACCEPTOR;
- NO_TETRAHEDRAL_RESPONSE;
- NOT_EVALUATED_TECHNICAL.

## Success boundary

A seed qualifies only if it shows numerically healthy, chemically intact, inherited movement toward tetrahedral geometry. Proton transfer is an observed response, not a required gate for the addition stage. Cross-seed status must report the denominator explicitly; one passing seed cannot be labeled a two-seed result.

No result from this pilot automatically authorizes release, shooting, PMF, NEB, string, Step2 product construction, or a TS/barrier/mechanism claim.
