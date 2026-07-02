# PETase blind QM/MM TS-core milestone (2026-07-02)

Objective: reproduce the PETase mechanism workflow using the literature methodology as a guide, while deriving reactant geometry, water involvement, and TS-core candidates from first-principles/blind structure-based reasoning rather than copying literature coordinates.

## Methodological scope

- System: PETase reactive complex prepared from structural/reactive-pose workflow.
- QM/MM engine used here: Amber `sander` DFTB3 QM/MM.
- TS search logic: generate chemically motivated reaction-coordinate brackets, run short unbiased/random-velocity committor trajectories, then locally refine and expand dual-commitment candidates.
- CPU policy: after user instruction, all QMMM batches kept total process threads <=64. The largest post-instruction batches used 24 jobs x OMP_NUM_THREADS=2 = 48 threads.
- Literature use: methodology only. Candidate geometries and water participation were selected from blind workflow observations and then organized for later comparison with the paper.

## Acylation TS-core status

Selected acylation center: `refined02_boundary03_cand00_rep1_f010_move0_rep2_1500step_f007_rep9_1500step_f004`.

Evidence:
- 16-rep committor: tetint 6, reactant 7, none 3; excluding none tetint fraction 0.462.
- Micro-relax attempt033: tetint 3, reactant 4, none 1; excluding none tetint fraction 0.429.
- attempt034 generated five local relaxed neighbor members preserving TS-core geometry.

Artifact directory:
- `blind_work/06_ts_ensemble/acylation/refined02_committor_ts_core_20260702`

## Deacylation TS-core status

Blind hydrolytic water candidate: HOH272.

Main deacylation seed found by bracketing:
- `bridge_cser240_h145`: product 6/16, reactant 4/16, undecided 6/16; excluding undecided product fraction 0.60.

Local ensemble expansion:
- `n01_more_attack`: product 3/8, reactant 2/8, undecided 3/8; excluding undecided product fraction 0.60.
- `n05_proton_late`: product 3/8, reactant 2/8, undecided 3/8; excluding undecided product fraction 0.60.
- `n00_center` became lower/reactant-side after expansion: product 1/8, reactant 4/8, undecided 3/8.

Artifact directory:
- `blind_work/06_ts_ensemble/deacylation/cser240_local_ts_core_20260702`

## Current interpretation

This milestone establishes a two-step PETase TS-core set from the blind workflow:

1. Acylation: a refined TS-core candidate with near-balanced commitment and local relaxed neighbors.
2. Deacylation: a hydrolytic-water-assisted TS-core cluster with three selected representatives and lower/product-side brackets.

The next scientific step is not another broad conformer search. It is a controlled literature comparison: align selected PDB/rst7 representatives to the reported initial/intermediate/TS geometries, quantify key distances/angles, and explain any deviations from differences in starting Michaelis complex, water preorganization, QM region, or committor criteria.

## Verification performed

- Server process check after calculation: no active `sander`, `pmemd`, `cp2k`, `gmx`, or `sqm` processes.
- `deacylation_ts_ensemble.tsv`: all rows verified as 10 tab-separated columns.
- `deacylation_cser240_local_ts_core_manifest.tsv`: all rows verified as 8 tab-separated columns.
- Raw attempt055 working directory is approximately 5.6G and intentionally excluded from Git milestone upload except for selected summaries/artifacts.
