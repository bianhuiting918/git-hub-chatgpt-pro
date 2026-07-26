# NylC activated-NAC rescreen design

Date: 2026-07-26  
Branch: `codex/nylc-l4-nac-to-l2-rebalance`

## Objective

Starting from recurrent, fully unrestrained M1 NAC conformations, construct and validate a post-activation Thr267 microstate before any Step1 QM/MM reaction-path calculation.

This workflow deliberately does not calculate how the Thr267 side-chain proton reaches the N-terminal alpha amino group. The calculated Step1 barrier will therefore be conditional on an already activated reactant state and must not be reported as the complete barrier from the resting enzyme.

## Authoritative parent ensemble

SCNet final audit job `61841413` completed with:

- technical status: `PASS`
- scientific status: `PASS_UNRESTRAINED_M1_ENSEMBLE_NAC`
- Stage B array: `61841378`, 18/18 replicas completed
- no final-audit failure-ledger entries

The three nominated classical-MM NAC medoids are:

1. `nac_evt18_time1206ps`, seed 26711, 258 ps
   - coordinate SHA256: `c8b9c8bf8988b7c7b1b470579894eb9b7501cc71bda6bb856e1c291d792e23dc`
   - Thr Oγ–substrate C distance: 0.3117068011 nm
   - O–C–Oγ angle: 109.6471457 degrees
2. `nac_evt25_time1462ps`, seed 26711, 606 ps
   - coordinate SHA256: `e1c292b7da52d68bdc3c41213d934f8fabcdc5f3dd953282274d793598690fae`
   - distance: 0.3460361793 nm
   - angle: 114.7746798 degrees
3. `nac_evt08_time1086ps`, seed 26711, 100 ps
   - coordinate SHA256: `6cddb322d02425c040afaa86dd9b2939e8dd933b9651b793b0a2ffa907611124`
   - distance: 0.3134983074 nm
   - angle: 102.0559099 degrees

These are classical-MM medoids, not QM/MM-optimized ground states.

## Activated microstate

Parent M1 state:

`Thr267 NalphaH2 / OgammaH / Asp306H / Asp308-`

Activated A1 state:

`Thr267 NalphaH3+ / Ogamma- / Asp306H / Asp308-`

The existing Thr267 HG1 atom is transferred from Oγ to Nalpha. Atom count and total system charge must remain unchanged. Asp306 and Asp308 protonation must not be changed during this transformation.

## Considered approaches

### A. Three medoids, three seeds each, 1 ns fixed-topology A1 screening — selected

Advantages: tests conformational and velocity-seed reproducibility at manageable cost and matches the user's requested workflow.

Limitation: fixed-topology MM cannot permit proton return or alternative proton transfer.

### B. All six Stage B candidates, three seeds each

Advantages: broader conformational coverage.

Limitation: doubles the A1 parameterization and simulation cost without strong evidence that lower-ranked M1 parents add useful activated-state diversity.

### C. Direct QM/MM optimization/MD from the three medoids

Advantages: avoids a nonstandard fixed-topology alkoxide patch.

Limitation: substantially more expensive and provides less classical ensemble screening before QM/MM path calculations.

Approach A is selected. Approach C is the fallback if a defensible force-field-native or independently audited A1 patch cannot be produced.

## Build and topology requirements

The A1 builder must operate on parsed topology and bonded graphs, not textual atom-name replacement.

It must:

- identify Thr267 by audited residue instance and atom connectivity;
- remove the Oγ–HG1 bond and add the Nalpha–HG1 bond;
- assign force-field-consistent bonded and nonbonded parameters for NalphaH3+/Ogamma-;
- preserve Asp306H, Asp308-, substrate, waters, ions, and all coordinates except the transferred proton;
- place the transferred proton using a tetrahedral Nalpha geometry;
- emit atom mapping, bond-difference, atom-count, total-charge, SHA256, and minimum-distance audits;
- refuse to run if the force field lacks a defensible parameter source for this combined terminal-threoninium/alkoxide state.

No trajectory, large topology, checkpoint, or secret is committed to GitHub.

## A1 simulation protocol

For each of the three parent medoids:

1. topology and charge audit;
2. `grompp` preflight;
3. energy minimization;
4. short low-temperature heavy-atom restrained relaxation;
5. gradual heating and restraint release;
6. three independent velocity seeds;
7. at least 1 ns fully unrestrained NPT per seed.

The restrained stages are technical preparation only and cannot satisfy the scientific gate.

## Scientific gate

Only the fully unrestrained A1 window is evaluated. A parent passes only if:

- the substrate remains bound in all required replicas;
- at least two of three replicas contain NAC frames after burn-in;
- a recurrent NAC conformational cluster is reproduced across at least two replicas;
- distance is at most 0.35 nm and angle is 95–115 degrees;
- temperature and volume are stable;
- no LINCS, SETTLE, NaN, or FATAL condition occurs;
- no severe nonbonded clash occurs.

MM potential energy is used only as a same-Hamiltonian screening descriptor, not as a free energy or activation barrier.

## QM/MM handoff

For each passing A1 parent, select a recurrent-cluster medoid. Retain two or three independent medoids for QM/MM preflight.

QM/MM preflight must:

- use explicit water;
- include Thr267 Nalpha/Ogamma and the transferred proton in the QM region;
- include the reactive PA66-L2 carbonyl/amide atoms and scientifically justified nearby proton-network residues or waters;
- perform an optimization without a proton-transfer coordinate constraint;
- reject an A1 medoid if the activated microstate or NAC collapses during unconstrained optimization.

A surviving structure is labeled `GS_act`, not the resting-state GS. Step1 energies are conditional on prior Thr activation. Step1 TS/PMF and Step2/deacylation are separate later gates.

## Reproducibility

All SCNet tasks must write rerunnable scripts, update `RUNBOOK.md`, and append `run_history.tsv` and `run_history.jsonl`. Failures remain explicit and do not stop unrelated candidates.
