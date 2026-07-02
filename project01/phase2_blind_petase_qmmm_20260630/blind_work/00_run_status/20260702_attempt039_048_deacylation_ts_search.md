# Deacylation TS search status: attempts 039-048 (2026-07-02)

Purpose: continue blind PETase deacylation mechanism reconstruction from self-identified HOH272, without using literature geometries as seeds.

## Key evidence

- attempt039 tested attempt037 chainA w02. Result: 8/8 complete short trajectories stayed undecided; water moved away (final OW-C 3.5-4.8 A) while C-SerO stayed long (2.9-3.4 A). This branch is not a clean deacylation TS path.
- attempt040 tested attempt037 chainB w03/w04 using H2 as the transferring water proton. Result: 8/8 returned to acyl-enzyme/reactant basin. This established a clean reactant-side bracket.
- attempt041 restrained C-SerO breaking from chainB w04 while keeping OW-C and H2-His preorganized. It produced TS-core probe geometries w01-w04; w04 was the most advanced clean early/undecided seed.
- attempt042 smoke committor from attempt041 w01-w04: w01-w03 were 6/6 reactant-side; w04 was 2/2 undecided/advanced. This moved the boundary later than w04.
- attempt043 extended later from w04. w05/w06 preserved late deacylation geometry; w07 was overdriven.
- attempt044 smoke committor from w05/w06: w05 was product 3/4, undecided 1/4; w06 was product 4/4. This established a product-side bracket.
- attempt045 generated fine windows between w04 and w05: w045, w050, w055.
- attempt046 smoke committor from w045/w050/w055: w045 was product 0/3, reactant 1/3, undecided 2/3; w050 was product 1/3, reactant 1/3, advanced/undecided 1/3; w055 was product 3/3.
- attempt047 extended w050 to 12 total replicas by adding 9 more. Combined w050 result: product 9/12, reactant 1/12, undecided 2/12. Therefore w050 is product-biased, not final balanced TS.
- attempt048 generated a finer bracket between w045 and w050: w239, w241, w243.

## Current bracket

Reactant/early side:
- attempt045 w045: OW-C 2.074 A, H2-His 1.466 A, C-SerO 2.376 A; committor product 0/3, reactant 1/3, undecided 2/3.

Product-biased side:
- attempt045 w050: OW-C 2.057 A, H2-His 1.458 A, C-SerO 2.439 A; combined committor product 9/12, reactant 1/12, undecided 2/12.

New finer candidates pending committor:
- attempt048 w239: OW-C 2.067 A, H2-His 1.461 A, C-SerO 2.341 A, C-Oleave 3.331 A.
- attempt048 w241: OW-C 2.006 A, H2-His 1.282 A, C-SerO 2.395 A, C-Oleave 3.801 A.
- attempt048 w243: OW-C 1.746 A, H2-His 1.060 A, C-SerO 2.412 A, C-Oleave 3.690 A.

## Interpretation

The current deacylation TS-core search has narrowed to the w045-w050 interval, with the next committor tests focused on attempt048 w239/w241/w243. The strongest current seed for a balanced TS is not w050 itself because it is product-biased after 12 replicas; the next likely center is one of the finer attempt048 windows.

## Next action

Run committor smoke tests from attempt048 w239/w241/w243, then expand whichever gives mixed product/reactant commitment into the deacylation TS-core ensemble.
