# Blind PETase Stage 1 Protonation Setup Notes

Date: 2026-06-30

## Boundary

This protonation setup uses only cleaned PETase coordinates, pH-7 residue chemistry, and generic serine-hydrolase catalytic logic. It does not use the PETase paper's TS structures, reaction-coordinate terms, selected CVs, trajectory data, barriers, rate-limiting assignment, or mechanistic conclusions.

## Current Input

The scan uses `*_chainA_initial_clean_v2.pdb` files generated from RCSB structures. These structures are not production-ready: hydrogens are not added, protonation is not assigned, missing residues are not rebuilt, and local minimization has not been run.

## Default pH-7 States

Use these as the initial nonreactive defaults before automated pKa prediction:

| residue | default |
| --- | --- |
| Asp | deprotonated, charge -1 |
| Glu | deprotonated, charge -1 |
| Lys | protonated, charge +1 |
| Arg | protonated, charge +1 |
| Cys | neutral thiol unless disulfide-linked |
| Tyr | neutral phenol |
| His | neutral tautomer must be assigned per local H-bond network |

## Primary Template Hypotheses

For the primary template `6EQE`:

- Catalytic `ASP206`: deprotonated is the primary state. A neutral-Asp sensitivity branch is required only if the His-Asp hydrogen-bond network or pKa prediction supports it.
- Catalytic `HIS237`: neutral tautomer is not fixed yet. Test HID/HIE tautomer alternatives before committing to QM/MM reaction scans. HIP should be included only if pKa prediction or local electrostatics supports a cationic histidine.
- `CYS203` and `CYS239`: near the triad but part of a geometric disulfide pair; keep the disulfide assignment unless preparation software contradicts it.
- Remote histidines such as `HIS104` and `HIS293`: require tautomer naming for topology generation, but they are not automatically reaction-mechanism sensitivity branches because they are far from the catalytic triad.

## Required External Check

Before production MD or QM/MM:

1. Run PROPKA, H++, pdb2pqr, Amber `reduce`, or an equivalent pKa/protonation workflow on the initial cleaned structure.
2. Record exact software, version, command, input SHA256, output SHA256, and changed residue states.
3. Compare automated predictions to `protonation_hypothesis_manifest.tsv`.
4. If the catalytic His/Asp assignment disagrees with the automated tool, run explicit sensitivity branches rather than choosing silently.

## Grill Gates

1. Can we defend the catalytic His tautomer without already assuming the reaction mechanism?
2. Does neutral Asp remain chemically plausible after seeing the His-Asp H-bond geometry and pKa estimate?
3. Are disulfide cysteines locked correctly before hydrogen placement?
4. Are remote His residues assigned for topology generation without inflating the reaction sensitivity tree?
5. Is every protonation branch linked to a concrete coordinate/topology file and not just a note?
