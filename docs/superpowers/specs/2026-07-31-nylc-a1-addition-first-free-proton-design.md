# NylC-GYAQ A1 activated-Thr addition-first pilot design

Date: 2026-07-31  
Branch: `codex/nylc-l4-nac-to-l2-rebalance`

## Objective

Test whether the already proton-transferred A1 state can undergo productive carbonyl addition while preserving a chemically intact Thr267 backbone and the correct attack geometry.

The starting chemistry is explicitly `Nalpha-H1/H2/HG1+` and `OG1-`. HG1 is already bonded to Nalpha, so this calculation does not test where an OG1 proton transfers. It tests only whether preactivated OG1 can approach C12 and induce early tetrahedralization.

This is a restrained response experiment. A passing result is not a transition state, barrier, minimum-free-energy path, or mechanism proof.

## Starting structures and provenance

Use coordinates directly from the existing activated-A1 unbiased MM trajectories. Never inherit a restart produced by jobs 62471131, 62477118, or 62500360.

- seed26723: evt25 frame 240, time 480 ps; OG1-C12 2.972 A; O2-C12-OG1 107.16 degrees; 8 of 11 neighboring frames satisfy strict NAC. TPR: `/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/a1_activated_nac_20260726/equilibration/nac_evt25_time1462ps/seed26723/attempt_61970146_4_61970151/npt300free/run.tpr`, SHA256 `c60078a92c2ace51facde4ef64e453f690177fc88b4d6363427f935944fa2e43`; XTC: sibling `run.xtc`, SHA256 `1a54f1b5b9f139b746c22d9e0f7e9a4a94eb8154bf2b881888986eedca933d89`.
- seed26737: evt25 frame 103, time 206 ps; OG1-C12 2.947 A; O2-C12-OG1 106.32 degrees; 7 of 11 neighboring frames satisfy strict NAC. TPR: `/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/a1_activated_nac_20260726/equilibration/nac_evt25_time1462ps/seed26737/attempt_61970146_5_61970146/npt300free/run.tpr`, SHA256 `dbd19a399547319d10630430ed494d33f6271bab0f30cb0de5a466c6af13ba20`; XTC: sibling `run.xtc`, SHA256 `fcba14da98b331368061dcd990f2467628ad77b9b7a9c4ce88090f92e0831b05`.

Before calculation, extract each frame once and record TPR/XTC provenance, frame/time, coordinate SHA256, topology SHA256, atom mapping, and periodic box. The extracted topology must prove:

- Nalpha bonded to H1, H2, HG1, and CA;
- HG1 bonded to Nalpha and not to OG1;
- OG1 bonded to CB and not to a hydrogen;
- Nalpha-CA, CA-CB, CB-OG1, and the substrate covalent graph intact.

A source failing this graph is `NOT_EVALUATED_TECHNICAL_SOURCE_INTEGRITY` and is not minimized or inherited.

## Hamiltonian

Retain the established activated-A1 146-QM, total charge 0, 510-electron, six-link-atom DFTB3 Hamiltonian, subject to actual Amber engine-banner verification.

The QM region and protonation state are immutable during this pilot. No water is added to the QM region. No conclusion about water-mediated proton transfer is permitted.

## Alternatives considered

1. **Recommended: authority minimization followed by addition-first windows in the same two-seed job.** Each raw MM frame first undergoes a fully unrestrained short, strict-convergence QM/MM minimization. Only a numerically healthy and chemically intact authority restart may enter the restrained windows. This distinguishes a bad starting Hamiltonian/topology from a failure to respond to the addition bias.
2. **Direct restrained windows from the raw frames.** Faster, but a covalent or numerical failure cannot be separated cleanly from the reaction-coordinate bias.
3. **New activated-A1 MM equilibration before QM/MM.** Provides more sampling, but changes the approved source ensemble and costs substantially more before checking the immediate question.

Use approach 1.

## Phase 0: unrestrained authority minimization

For each seed, run one short strict-convergence QM/MM minimization with no attack, carbonyl, proton, or angle restraint.

Persist the full engine banner, NSTEP, RMS, GMAX, energies, input and output restart SHA256, and the complete Thr267/substrate covalent-integrity audit.

The source qualifies for Phase 1 only if:

- Amber exits normally and no DFTB/EAMBER overflow or hard error occurs;
- the 146-QM/q0/510e/6link/DFTB3 contract is present;
- Thr267 and substrate covalent graphs remain intact;
- OG1-C12 remains within 3.5 A;
- O2-C12-OG1 remains within 95-115 degrees;
- no output restart SHA is duplicated spuriously.

A failure is technical or source-state evidence, not a mechanism failure. It stops that seed without retry or force adjustment.

## Phase 1: reaction-coordinate policy

Biased coordinates:

- OG1-C12: target moves by -0.04 A per accepted window;
- C12-O2: target moves by +0.03 A per accepted window.

Use the established base force constants without scale multiplication:

- OG1-C12: 24.0 kcal mol-1 A-2;
- C12-O2: 30.0 kcal mol-1 A-2.

No restraint is applied to qPT, HG1, C12-N3, proton acceptors, or the attack angle. Because HG1 is already on Nalpha, qPT is monitored only as a state/integrity coordinate, not as a free-proton reaction coordinate.

Monitored coordinates:

- OG1-C12;
- C12-O2;
- C12-N3;
- Nalpha-H1/H2/HG1 and OG1-H1/H2/HG1 distances;
- O2-C12-OG1 angle;
- carbonyl pyramidalization, angle sum, and local OOP;
- Nalpha-CA, CA-CB, and CB-OG1 covalent bonds.

## Pilot size and strict inheritance

Run two independent array tasks, one per seed, with 8 MPI ranks per task. Each task is bounded to six sequential minimization windows.

Within each seed:

1. Phase 0 starts from the SHA-fixed extracted MM frame.
2. Window 1 starts only from the accepted Phase-0 restart.
3. Window n uses targets relative to the Phase-0 geometry: `attack_0 - 0.04*n` and `carbonyl_0 + 0.03*n`, for n=1..6.
4. A window is accepted only when the engine is numerically healthy, the Hamiltonian banner passes, all covalent-integrity gates pass, both biased coordinates respond in the required directions, and the attack-angle guard remains 95-115 degrees.
5. Only an accepted output restart may seed the next window.
6. On failure, stop that chain. Never inherit the failed restart and never silently change force, target, topology, or protonation state.
7. Persist source and output restart SHA256, physical restraint lines, engine tail, RMS/GMAX/energy, and compact geometry audit for every attempted stage.

## Outcome classification

Per seed, classify exactly one terminal state:

- `ACTIVATED_A1_EARLY_TETRAHEDRAL_RESPONSE`;
- `NO_TETRAHEDRAL_RESPONSE_UNDER_TESTED_BIAS`;
- `NOT_EVALUATED_TECHNICAL_SOURCE_INTEGRITY`;
- `NOT_EVALUATED_TECHNICAL_NUMERICAL`;
- `NOT_EVALUATED_TECHNICAL_CHEMICAL_INTEGRITY`.

Report the two-seed denominator explicitly. A one-seed response is not a two-seed result.

## Success boundary

A seed qualifies only if it passes Phase 0 and then shows numerically healthy, chemically intact, strictly inherited movement toward tetrahedral geometry. This establishes only response eligibility for a later unconstrained release test.

No result from this pilot automatically authorizes release, shooting, PMF, NEB, string, Step2 product construction, or a TS/barrier/mechanism claim.
