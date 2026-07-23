# Five-candidate PA66-L4 NAC to PA66-L2 rebalance design

## Scope

Build and compare five independent PA66-L2 systems from fully unrestrained PA66-L4 NAC frames:

1. NylC-GYAQ C18 at 11854 ps.
2. NylC-GYAQ C23 at 29684 ps.
3. Nyl50 C18 at 70792 ps.
4. Nyl12 J1 at 37848 ps.
5. Nyl12 J2 at 87418 ps.

All five remain in the formal candidate universe. A build or simulation failure in one candidate does not stop later candidates and is recorded as NOT_EVALUATED or FAIL with an exact reason. No docking pose is used as a fallback.

## Source and graph contract

Each source trajectory and topology is immutable. The representative frame is extracted from its named FREE_MONITOR segment at the local XTC time recorded in the candidate manifest.

The source L4 ITP is parsed into an element-labelled bonded graph. The audited PA66_L2_GMX.itp is the only target chemistry. The builder maps the 32 inherited L2 heavy atoms (all heavy atoms except the new terminal carboxylate O5) onto a connected source subgraph while enforcing:

- source reactive carbonyl C, carbonyl O and bonded amide N map together to one internal L2 amide;
- all non-terminal mapped heavy atoms have exact graph degree and element agreement;
- extra source heavy-atom edges are allowed only at L2 terminal amine N1 and terminal carboxyl C24;
- both endpoints must be genuine cuts for the formal mapping;
- source C/O/N coordinates are copied exactly;
- all other matched hydrogens are assigned through bonded-parent equivalence classes, never atom-name order;
- terminal ammonium hydrogens and the second carboxylate oxygen are rebuilt with audited L2 bond lengths;
- every source-to-L2 atom mapping and each cut edge are written to JSON.

The historical NylC-C23 fixed mapping is a regression oracle only. It is not applied to NylC-C18, Nyl50 or Nyl12.

## System replacement and audit

The builder replaces exactly one source L4 molecule in the GRO, changes the topology include and molecule entry to PA66_L2, and preserves the rest of the system in order. Candidate outputs include input hashes, atom counts, total L2 charge, cut edges, atom mapping, reactive geometry, minimum nonbonded distance, and topology preprocessing status.

A candidate may enter equilibration only after:

- target L2 has 79 atoms, 33 heavy atoms and the audited integer charge;
- the system atom-count delta equals source-L4 atoms minus 79;
- the target C/O/N coordinates match source coordinates within GRO precision;
- topology and coordinate molecule counts agree;
- grompp succeeds without undocumented warnings;
- energy minimization terminates without FATAL, NaN, LINCS or SETTLE errors;
- the post-EM minimum nonbonded distance and reactive geometry are reported.

## Re-equilibration

Each passing build follows an independent bounded schedule:

1. restrained energy minimization;
2. 100 ps NVT at 50 K with 1000 kJ mol-1 nm-2 L2-heavy restraint;
3. 100 ps NVT at 150 K with 500 restraint;
4. 200 ps NVT at 300 K with 100 restraint;
5. 200 ps NPT at 300 K with 100 restraint;
6. 200 ps NPT at 300 K with 10 restraint;
7. at least 1 ns fully unrestrained NPT at 300 K and 1 bar.

No restrained stage counts toward scientific NAC success.

## Scientific comparison

Only the final fully unrestrained NPT window is compared. For each candidate report the same PBC-aware metrics:

- reactive C to catalytic Thr OG1 distance <= 0.35 nm;
- C-centered attack angle 95-115 degrees;
- joint-frame NAC occupancy, never separate marginal occupancies;
- validated branch-specific gate opening (NylC gate is residues 261-266 and excludes Thr267);
- temperature and pressure stability;
- LINCS, SETTLE, NaN and FATAL scans.

A candidate is eligible for DFTB3/3OB-3-1 QM/MM preflight only from this unconstrained evidence. Technical incompleteness remains NOT_EVALUATED and geometry failure remains FAIL; neither is silently converted to biological inactivity.
