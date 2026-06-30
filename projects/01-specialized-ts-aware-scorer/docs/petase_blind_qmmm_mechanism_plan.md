# PETase Blind QM/MM Mechanism Plan

Date: 2026-06-30

## Boundary

This plan treats the PETase paper as methodological inspiration only. The blind workflow must not use the paper's concrete mechanistic results as inputs.

Allowed:

- PETase structures from public structural databases.
- General serine hydrolase chemistry.
- Generic QM/MM, TS-search, PMF, committor, and sampling methods.
- The paper's broad research style: unbiased sampling, likelihood/RC validation, committor tests, and free-energy calculation.

Not allowed before final validation:

- Paper TS coordinates.
- Paper reaction-coordinate formulas.
- Paper selected CVs.
- Paper umbrella window setup.
- Paper aimless-shooting trajectories.
- Paper reported barriers, rate constants, or rate-limiting assignment.
- Paper's conclusion about moving/flipping histidine, Trp185 role, Asp206 role, or tetrahedral TS/intermediate status.

The final comparison against the paper happens only after our blind QM/MM workflow has produced its own TS ensembles and mechanistic conclusions.

## Goal

Starting only from PETase structure and substrate chemistry, determine:

- the acylation mechanism and TS ensemble;
- the deacylation mechanism and TS ensemble;
- whether the tetrahedral configuration is a TS or metastable intermediate;
- how His/Asp proton transfer behaves;
- the relative barriers of the two steps;
- whether active-site conformational features such as Trp185 affect catalysis.

## Phase 1 - Structure-Only System Definition

Inputs:

- PETase structure candidates, including apo and ligand/analog-bound structures if available.
- PET model fragments such as BHET, MHET, PET dimer, or PET tetramer.
- Standard residue chemistry and pH around 7.

Tasks:

1. Select 2-4 PETase structural templates by resolution, active-site completeness, catalytic triad geometry, and substrate/analog relevance.
2. Repair missing atoms/residues and assign disulfides.
3. Assign protonation states using PROPKA/H++ plus manual catalytic-site inspection.
4. Generate PET fragment parameters for MM pre-equilibration.
5. Build multiple substrate-bound Michaelis candidates:
   - docking;
   - PLACER-like local conformer generation if useful;
   - manually restrained poses based on generic serine-hydrolase attack geometry.
6. Reject poses where Ser O-gamma cannot plausibly attack the ester carbonyl or the oxyanion hole cannot stabilize the carbonyl oxygen.

Deliverable:

- `blind_work/01_system_setup/gs_pose_manifest.tsv`
- 10-50 chemically plausible GS pose candidates.

## Phase 2 - Classical Ensemble Preparation

Tasks:

1. Solvate and neutralize each selected GS pose.
2. Run restrained minimization and staged equilibration.
3. Run short classical MD replicates for each pose.
4. Cluster active-site conformations by:
   - Ser O-gamma to ester carbonyl C;
   - attack angle;
   - His position relative to Ser and leaving oxygen;
   - oxyanion-hole H-bonds;
   - water positions for deacylation;
   - Trp185 rotamer and substrate aromatic pose.
5. Select representative productive conformers and explicitly record rejected nonproductive basins.

Deliverable:

- `blind_work/02_classical_md/productive_conformer_manifest.tsv`
- representative equilibrated GS structures for QM/MM exploration.

## Phase 3 - Mechanism Hypothesis Tree

Hypotheses to test without using the paper's conclusions:

### Acylation

- concerted Ser attack plus proton transfer;
- stepwise proton transfer before attack;
- tetrahedral intermediate followed by leaving-group protonation;
- direct vs indirect His-mediated leaving-group proton transfer.

### Deacylation

- water attack assisted by His;
- water deprotonation before carbonyl attack;
- concerted water attack plus Ser-acyl bond cleavage;
- stable tetrahedral intermediate vs single TS.

### Proton Relay

- His shuttles proton without ring flip;
- His flips before donating proton;
- Asp participates only by H-bond stabilization;
- Asp participates in double proton transfer.

Deliverable:

- `blind_work/03_mechanism_tree/mechanism_hypotheses.yaml`
- explicit CV candidates for each hypothesis.

## Phase 4 - Low-Cost QM/MM Exploration

Purpose: find reactive paths, not final barriers.

Initial QM region:

- Ser side chain;
- His side chain;
- Asp side chain;
- reactive PET ester fragment;
- leaving group atoms;
- acylation/deacylation reactive water as needed;
- link atoms at C-alpha/C-beta boundaries.

Exploration methods:

- DFTB3/MM, xTB/MM, or low-cost DFT/MM;
- 1D and 2D relaxed scans;
- restrained short QM/MM MD;
- string/NEB if available.

Candidate acylation coordinates:

- Ser O-gamma - carbonyl C forming distance;
- carbonyl C - leaving O breaking distance;
- Ser H - His acceptor distance;
- His H - leaving O distance;
- nucleophilic attack angle;
- carbonyl oxygen H-bonds.

Candidate deacylation coordinates:

- water O - acyl carbonyl C forming distance;
- water H - His acceptor distance;
- Ser O-gamma - acyl carbonyl C breaking distance;
- His H - Ser O-gamma distance;
- tetrahedralization geometry.

Deliverable:

- `blind_work/04_qmmm_exploration/path_screening_table.tsv`
- candidate TS-like structures for each viable pathway.

## Phase 5 - TS Search and Refinement

For each viable pathway:

1. Extract TS guesses from scan maxima, string images, or constrained path points.
2. Refine using QM/MM TS optimization if available.
3. If full TS optimization is unavailable, refine by constrained optimization at the dividing surface plus vibrational/IRC-like checks where possible.
4. Expand or adjust QM region and repeat for top candidates.
5. Reject candidates that:
   - relax to reactant or product without a saddle region;
   - show chemically wrong protonation;
   - require impossible active-site geometry;
   - are unstable to modest QM-region expansion.

Deliverable:

- `blind_work/05_ts_refinement/acylation_ts_candidates/`
- `blind_work/05_ts_refinement/deacylation_ts_candidates/`
- `blind_work/05_ts_refinement/ts_refinement_manifest.tsv`

## Phase 6 - TS Ensemble Construction

Tasks:

1. Repeat TS search across multiple productive classical conformer clusters.
2. Cluster TS candidates by reaction geometry and active-site conformation.
3. Preserve ensemble diversity instead of keeping only the lowest-energy structure.
4. For each TS cluster, select representative structures for committor testing.

Deliverable:

- `blind_work/06_ts_ensemble/acylation_ts_ensemble.tsv`
- `blind_work/06_ts_ensemble/deacylation_ts_ensemble.tsv`

## Phase 7 - Dynamical Validation

Required checks:

- short forward/backward trajectories from TS candidates;
- committor estimates for representative candidates;
- product/reactant endpoint classification;
- IRC-like verification where feasible;
- no reliance on energy maximum alone.

Acceptance criteria:

- representative TS candidates commit to reactants and products with roughly balanced probability;
- trajectories reach chemically correct acylation/deacylation endpoints;
- key reaction atoms move along the expected bond-making, bond-breaking, and proton-transfer directions.

Deliverable:

- `blind_work/07_committor/committor_summary.tsv`
- accepted and rejected TS candidate labels.

## Phase 8 - Free-Energy Barriers

Only after TS candidates pass dynamical checks:

1. Define reaction coordinates from our own successful trajectories and CV analysis.
2. Run umbrella sampling/string-window sampling/metadynamics as appropriate.
3. Analyze PMFs with MBAR/WHAM.
4. Estimate barriers and uncertainty for acylation and deacylation.
5. Check sensitivity to QM region, protonation, starting conformer, and substrate pose.

Deliverable:

- `blind_work/08_free_energy/acylation_pmf.tsv`
- `blind_work/08_free_energy/deacylation_pmf.tsv`
- `blind_work/08_free_energy/barrier_summary.md`

## Phase 9 - Final Paper Validation

Only here may we compare against the paper.

Compare:

- mechanism sequence;
- acylation TS geometry;
- deacylation TS geometry;
- tetrahedral TS vs intermediate;
- His relay behavior;
- Asp role;
- oxyanion hole role;
- Trp185/conformational effects;
- barrier ordering and rate-limiting step.

If our result differs, audit:

- starting structure;
- substrate model;
- protonation state;
- QM region;
- sampling depth;
- reaction coordinate bias;
- solvent and water placement;
- force-field and QM method.

## Grill Gates

Gate 1: Which structure/substrate scope?

Recommended: WT PETase with a PET dimer or BHET/MHET-like ester fragment, acylation and deacylation both included.

Gate 2: Which compute tier first?

Recommended: DFTB3/MM or xTB/MM exploratory scans first, then DFT/MM refinement only for surviving pathways.

Gate 3: What counts as success in the first round?

Recommended: chemically validated TS ensembles and relative barrier ordering for WT, not mutant ranking.

Gate 4: When to open the paper result?

Recommended: only after Phase 7 has accepted TS candidates and Phase 8 has preliminary PMFs.
