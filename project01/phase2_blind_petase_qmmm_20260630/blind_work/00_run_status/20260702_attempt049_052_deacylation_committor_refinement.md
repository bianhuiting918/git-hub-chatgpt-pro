# Deacylation committor refinement: attempts 049-052 (2026-07-02)

Purpose: refine the blind deacylation TS-core bracket after attempt048 generated finer candidates between w045 and w050.

CPU policy: after user instruction, keep total QMMM CPU threads <= 64. The runs here used 8-12 concurrent jobs with OMP_NUM_THREADS=2, i.e. <=24 threads during the largest batch.

## Attempt049: finer candidates w239/w241/w243, 1500-step

All 12/12 trajectories completed, no fatal/QMMM failure strings.

- w239: product 1/4, reactant 0/4, undecided 3/4.
- w241: product 4/4, reactant 0/4, undecided 0/4.
- w243: product 4/4, reactant 0/4, undecided 0/4.

Interpretation: w241/w243 are product-side. w239 is closest to the boundary but mostly undecided on the 1500-step timescale.

## Attempt050: w239 extended to 3000-step

All 8/8 trajectories completed, no fatal/QMMM failure strings.

- attempt050 alone: product 5/8, reactant 1/8, undecided 2/8.
- combined w239 with attempt049: product 6/12, reactant 1/12, undecided 5/12.

Interpretation: w239 has both product and reactant commitment, but excluding undecided it is product-biased. It remains the best current TS-core probe, not yet a final balanced TS ensemble member.

## Attempt051: w045 extended to 3000-step

All 8/8 trajectories completed, no fatal/QMMM failure strings.

- product 4/8, reactant 0/8, undecided 4/8.

Interpretation: the earlier w045 also shows product commitment on the longer timescale and no reactant commitment; it is not a clean lower bracket at 3000-step.

## Attempt052: attempt041 w04 extended to 3000-step

All 8/8 trajectories completed, no fatal/QMMM failure strings.

- product 0/8, reactant 3/8, undecided 5/8.

Interpretation: attempt041 w04 is the best current early/reactant-side bracket under the same 3000-step criterion.

## Current deacylation bracket

Lower/early bracket:
- attempt041 w04: OW-C 2.033 A, H2-His 1.499 A, C-SerO 2.382 A; 3000-step product 0/8, reactant 3/8, undecided 5/8.

Upper/product-emergent bracket:
- w239: OW-C 2.067 A, H2-His 1.461 A, C-SerO 2.341 A; 3000-step product 5/8, reactant 1/8, undecided 2/8; combined product 6/12, reactant 1/12, undecided 5/12.

This indicates the TS-core dividing surface is close to this geometry range, but more refinement is needed because w239 is product-biased and attempt041 w04 is mostly undecided/reactant-side.

## Next action

Generate and test an intermediate relaxed candidate between attempt041 w04 and w239, or vary H2-His/OW-C slightly while keeping C-SerO near 2.34-2.38 A. Use 3000-step committor for comparability.
