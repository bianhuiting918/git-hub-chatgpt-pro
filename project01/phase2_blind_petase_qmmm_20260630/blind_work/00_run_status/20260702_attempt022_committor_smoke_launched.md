# Attempt022 committor smoke launched from attempt021 accepted candidates

Purpose: perform a light committor-style validation from the attempt021 accepted TS-candidate ensemble before claiming a final TS set.

Input candidate table:
`blind_work/05_atesa_acylation/attempt021_all_accepted_ts_candidates.csv`

Attempt022 directory:
`blind_work/05_atesa_acylation/attempt022_committor_smoke_from_attempt021`

Design:
- 5 accepted attempt021 TS-candidates.
- 2 random-velocity DFTB3/MM trajectories per candidate.
- 10 total `sander` jobs launched in parallel.
- Each trajectory: 1500 steps, dt = 0.0001 ps, 0.15 ps total.
- Same QM region and DFTB3/MM setup as attempt021.

This is a smoke validation only. It should classify whether candidates trend toward the reactant basin, tetrahedral-intermediate basin, late/product side, or remain undecided over the short window.

Basin definitions for classification should match attempt021:
- Reactant/Michaelis-like: SerO-C > 2.4 A, SerH-SerO < 1.2 A, SerH-His > 1.4 A, C-Oleave < 1.8 A.
- Tetrahedral-intermediate-like: SerO-C < 1.8 A, SerH-SerO > 1.4 A, SerH-His < 1.2 A, C-Oleave > 1.8 A.

The attempt022 results are not a final committor calculation; they are intended to decide which candidates deserve more replicas.
