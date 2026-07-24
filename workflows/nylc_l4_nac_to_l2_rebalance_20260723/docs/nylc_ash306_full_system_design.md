# NylC-C18 ASH306 full-system preflight design

## Scope

Build one neutral full-system microstate from the immutable 434 ps fully-unrestrained PA66-L2 NAC frame. This is a topology and numerical-entry preflight only. It is not a ground-state proof, transition state, reaction coordinate, PMF, or barrier.

## Immutable inputs

- Source GRO SHA256: `182f68a41c72c490e0046323c9b8a0f8514ccf7af7dc4d578dcb87342c1c8e2d`.
- Source atom count: 133589.
- Protein chain H spans source global atoms 8949-10272 and contains 1324 atoms.
- PA66-L2 spans source global atoms 10273-10351 and contains 79 atoms.
- Audited ASH306 probe: SCNet job 61717760.
- Probe chain GRO SHA256: `a3998b6b5f5df5a40c78861dfd4b65983dff5d96bba1cd7ed238c7d57528ea44`.
- Probe inline topology SHA256: `630b45eb20713b1f604ae11a193ded778edcfb52f44b70c9377d6858bd7690f9`.

## Considered approaches

1. Recommended: replace complete chain H with the audited 1325-atom ASH306 chain, extract its molecule topology, and remove one distant Na+ so the periodic system remains neutral.
2. Rejected: patch one hydrogen, charges, bonds, angles, and dihedrals manually in the old topology. This is unnecessarily error-prone.
3. Rejected: retain the +1 periodic system and rely on uniform neutralizing background. That would make later energetic comparisons and PMF baselines less clean.

## Build algorithm

1. Verify every immutable SHA and atom range before writing output.
2. Extract the chain-H molecule topology from the first `[ moleculetype ]` through the end of its bonded sections, excluding force-field, solvent, ion, system, and molecules blocks.
3. Splice the audited chain-H GRO atoms into the source GRO.
4. Treat terminal Lys355 `OC1/OC2` as chemically interchangeable only for coordinate comparison; all other heavy atoms must retain coordinates within 0.0011 nm.
5. Select one Na+ with the largest minimum periodic distance from protein plus L2 using the triclinic box. Remove that exact ion and record its source global atom, residue, coordinate, and distance.
6. Reduce the `NA` molecule count from 144 to 143. The inserted HD2 and removed Na+ keep the final atom count at 133589.
7. Preserve the PA66-L2 reactive C/O/N coordinates exactly except for the unavoidable GRO text precision already present in the immutable source.

## Required audits

- Exactly one ASH306 residue with 13 atoms including HD2; Asp308 remains ASP.
- Thr267 retains H1/H2/H3/OG1/HG1.
- Chain-H topology charge changes from -4 to -3; full-system ParmEd charge is 0 within 1e-4.
- Molecule counts are protein A/E/D/H each 1, PA66_L2 1, SOL 40990, NA 143, CL 124.
- No missing or extra bonded heavy-atom keys; terminal oxygen symmetry is explicitly reported.
- NylC NAC remains distance <=0.35 nm and angle 95-115 degrees.
- Minimum ligand-protein heavy-atom distance is reported and must be >=0.20 nm.
- GROMACS grompp succeeds with no atom-count mismatch or fatal topology error.
- ParmEd loads the topology plus coordinates and saves Amber prmtop/rst7 with 133589 atoms and total charge 0.
- Failure writes a compact FAIL audit and run-history record; it does not delete or overwrite prior jobs.

## Output boundary

Large GRO, TPR, prmtop, and rst7 files remain only on SCNet. GitHub receives scripts, tests, RUNBOOK updates, compact JSON audits, and small run-history excerpts. No credentials or trajectories are stored in GitHub.
