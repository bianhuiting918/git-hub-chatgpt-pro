# NylC A1 Step2 water-reorganization sampling design

Date: 2026-07-29  
Branch: `codex/nylc-l4-nac-to-l2-rebalance`

## Objective

Determine whether an attack-competent water can enter naturally around the stable Step1 acyl enzyme before constructing a Step2 pathway. The prior fixed-water calculation showed 6/6 technically valid release legs but 0/6 attack-water residence; it did not test water-identity exchange.

## Alternatives considered

1. Pure MM water sampling: cheapest, but the frozen topology does not encode the Step1 acyl connectivity and proton state. Rejected unless a separate acyl-state classical topology is built.
2. 146-QM/no-QM-water DFTB3 QM/MM sampling: preserves the acyl chemistry while all 40,990 waters remain exchangeable MM waters. Selected.
3. Product-side reverse construction: retained as the next route if natural water-entry sampling produces no reproducible candidates.

## Frozen inputs and Hamiltonian

- Sources: the two SHA-fixed Step1 acyl endpoint restarts from attempts `62216380_0/1`.
- Frozen prmtop SHA256: `a61d15bf0bf78675be93275d45f274e808ed6ae450fc1ca21a8e14aee8c12ca0`.
- QM region: existing Step1 146-atom core, q0, 510 electrons, 6 link atoms, no QM water.
- Amber18 `sander.MPI`, DFTB3/3OB-3-1.
- Engine banner must reproduce the complete contract before production dynamics.

## Pilot calculation

Run eight independent trajectories in one unthrottled array:

- 2 endpoint seeds x 4 independent velocity seeds.
- 8 MPI ranks per task.
- 4,000 steps at 0.5 fs = 2 ps per trajectory.
- No reactive-coordinate or water-position restraints.
- Preserve only final restart, compact water-hit records, manifest, hashes and logs; full trajectory remains scratch-only.

## Water-entry analysis

At every saved frame, examine complete waters near C12 and record both strict and near-hit geometry. A strict candidate requires:

- C12-OW 2.7-3.5 A;
- O2-C12-OW angle 95-125 degrees;
- donor-H to Nalpha <=2.5 A;
- OW-H-Nalpha angle >=130 degrees;
- Step1 acyl-connectivity and proton-site guards passing;
- at least three consecutive saved frames.

Candidate identity is not fixed in advance. Keep at most two temporally independent strict candidates per endpoint seed. Record the complete local water universe and denominators even if no candidate passes.

## Downstream gate

Only a strict natural-entry candidate may be promoted as one complete QM water (149 QM atoms, q0, expected 518 electrons, 6 links) for two independent fully unrestrained A2 validation legs. If either seed has no strict candidate, do not force a replacement water; advance instead to the separately versioned product-side reverse coupled path.

## Minimal implementation and verification

Add one driver, one runner and one Slurm wrapper, plus a focused contract test. Reuse existing authority, geometry, triclinic-PBC, engine-banner and compact-audit helpers. Validate input hashes, the 146/q0/510e/6-link/no-QM-water contract, array mapping, absence of reactive restraints, and non-overwriting output roots before the single submission.

Technical completion, water entry, or A2 persistence is not evidence of a transition state, PMF, barrier or mechanism.
