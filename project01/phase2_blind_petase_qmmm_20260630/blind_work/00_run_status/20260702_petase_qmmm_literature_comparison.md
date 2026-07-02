# PETase QM/MM literature-comparison checkpoint (2026-07-02)

Purpose: compare the current blind QM/MM TS-core ensemble against the PETase paper at the level that is currently justified: workflow, reaction-coordinate definitions, committor evidence, and selected TS-core geometry. This is not yet a full reproduction of the paper free-energy profiles because umbrella sampling/MBAR barriers have not been run for our selected RCs.

## Literature workflow extracted from the paper

- Starting point: the paper began from a hydroxyethyl-capped PET dimer Michaelis complex from prior docking, followed by QM/MM equilibration with DFTB3.
- TS discovery: the paper used ATESA `find_ts` / aimless shooting to obtain putative TS structures without assuming a single hand-built coordinate.
- RC selection: linear combinations of CVs were optimized and validated by committor probability histograms centered near 0.5.
- Free energies: after RC discovery and validation, umbrella sampling along the discovered RC was used for acylation and deacylation profiles.
- Mechanistic motif: both steps use a moving His237 proton shuttle; no His flip or double-proton-transfer mechanism is required in the paper interpretation.

## Literature CV definitions used for comparison

Acylation, from the Fig. 2 caption:
- CV1 = d(Ser160 H, PET carboxyl O) - d(Ser160 O, Ser160 H).
- CV2 = d(Ser160 O, PET carboxyl C) - d(PET carboxyl C, PET carboxyl O).
- CV3 = angle Ser160 H - PET ether O - BHET C, in radians.

Deacylation, from the Fig. 4 caption:
- CV1 = d(Ser160 O, water H) - d(water H, water O).
- CV2 = d(Ser160 O, MHET C) - d(MHET C, water O).
- CV3 = d(MHET C, water O) - d(water O, water H).

Our acylation CV3 is currently a proxy because the reduced atom map tracks the scissile ether O and carboxyl C, but not the exact adjacent BHET C used in the paper caption. Deacylation CVs use the H2 proton that was identified in our blind deacylation path as the transferring water proton.

## Selected TS-core geometry

| member | step | SerO-C A | C-Oleave A | SerH-HisNE2 A | WatO-C A | WatH2-HisNE2 A | acyl attack deg | deacyl attack deg | paper-like acyl CV2 | paper deacyl CV1 | paper deacyl CV2 | paper deacyl CV3 |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| acyl_refined02_core | acylation | 2.057 | 1.847 | 1.450 | 4.099 | 3.622 | 99.6 | 64.6 | 0.881 | 2.838 | -2.042 | 3.142 |
| deacyl_bridge_cser240_h145 | deacylation | 2.339 | 3.287 | 1.669 | 2.096 | 1.477 | 99.6 | 101.6 | 1.090 | 1.349 | 0.243 | 1.034 |
| deacyl_n01_more_attack | deacylation | 2.262 | 3.347 | 1.704 | 2.116 | 1.516 | 102.0 | 102.2 | 1.095 | 1.369 | 0.145 | 1.076 |
| deacyl_n05_proton_late | deacylation | 2.236 | 3.379 | 1.731 | 2.134 | 1.559 | 101.8 | 101.2 | 1.071 | 1.400 | 0.103 | 1.110 |

## Comparison and interpretation

Acylation: our accepted acylation representative is chemically TS-like by the same qualitative criteria as the paper: SerO-C is partially formed at 2.057 A, C-Oleave is elongated to 1.847 A, and the Ser proton is transferred toward His237 with SerH-HisNE2 = 1.450 A. The attack angle is 99.6 degrees, which is compatible with a tetrahedralizing carbon but not identical to the user-quoted literature initial/pose angle of 79.7 degrees; that difference is expected because this row is a TS-core member, not the literature initial Michaelis complex.

Deacylation: all three selected deacylation representatives reproduce the paper qualitative TS motif: catalytic water is near the acyl carbon (WatO-C = 2.096-2.134 A), the water proton is partly transferred to His237 (WatH2-HisNE2 = 1.477-1.559 A), and the SerO-C acyl bond is elongated (2.236-2.339 A). The deacylation CV2 values are positive but small (0.103-0.243 A), consistent with a dividing-surface region where water attack and Ser-acyl bond cleavage compete.

Committor agreement: the paper validates RCs by committor histograms near 0.5. Our current evidence is a lower-throughput analog: acylation selected member has near-balanced tetint/reactant commitment after excluding none outcomes, and deacylation selected members have product fractions of 0.60 after excluding undecided outcomes. This supports TS-core selection, but does not replace the paper-level ATESA histogram plus umbrella sampling validation.

Main discrepancy: our workflow has not yet generated paper-equivalent free-energy profiles, transmission coefficients, or literature-coordinate overlays. It has produced a blind TS-core ensemble that is ready for those next stages.

## Next actions required for full paper reproduction

1. Recover or reconstruct the paper reference structures/coordinates for Michaelis complex, acylation TS, AEI, deacylation TS, and product states.
2. Align our selected PDBs to the reference active site and compute RMSD plus the same key CVs.
3. Use ATESA or an equivalent linear-CV search around our TS-core members to generate formal pB histograms.
4. Run umbrella sampling along the discovered RCs and compute barriers with MBAR before claiming full reproduction of the paper energy profiles.
