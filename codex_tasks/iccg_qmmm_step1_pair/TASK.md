# Project1 continuation: ICCG Step1 LG1/LG2 QM/MM pair gate

Work only inside `codex_tasks/iccg_qmmm_step1_pair` on branch
`codex/iccg-qmmm-step1-pair-gate`. Read-only use of
`codex_tasks/iccg_lg1_dimer/inputs` and its geometry audit is allowed.
Do not open or merge a pull request.

## Final scientific context

Project1 phase 1 is an auditable four-state fixed-geometry QM/MM
single-point-energy dataset, not a free-energy-barrier calculation. The four
states are Step1 ground/NAC (LG1), Step1 TS (LG2), Step2 initial (LG3), and
Step2 TS (LG4). This task advances the next dependency only: make and gate one
ICCG Step1 LG1/LG2 paired pilot before any production array.

The prior PHL7 Step2 pose is invalid and must never be reused as a label:
`FAIL_GEOMETRY_CLASH_NOT_LABEL`. The verified ICCG LG1 rigid transfer has
heavy-atom geometry PASS (minimum contact 2.669 A, maximum vdW overlap 0.551 A)
but is only a viewer/heavy-atom artifact, not a protonated topology or energy
input.

## Inputs and provenance

The scripts must take explicit CLI paths rather than embed local paths.

Required production inputs on the remote scientific server:

- active ICCG protein:
  `/work/home/acshdt1dks/iccg_sequence_design_20260723/inputs/ICCG_active_chainA.pdb`
- literature LG1:
  `/work/home/acshdt1dks/petase_qmmm_pilot8_20260721/inputs/reference/IsPETase_WT_LG1.pdb`
- literature LG2:
  `/work/home/acshdt1dks/petase_qmmm_pilot8_20260721/inputs/reference/IsPETase_WT_LG2.pdb`
- released literature Amber topology/parameter sources already audited under
  `/work/home/acshdt1dks/petase_qmmm_pilot8_20260721/inputs/reference`
  and the corresponding literature-source project paths recorded there.
- Amber18 runtime and the existing SHA-256-pinned AmberTools20 ff14SB data
  snapshot from the pilot project.

The GitHub sibling directory supplies the ICCG heavy-atom and LG1 reference
fixtures for unit tests. Do not copy generated local debug artifacts into the
repository.

## Pair construction contract

1. Parse and validate the 258-residue active ICCG chain A and Ser165 OG.
2. Independently align each released literature complex to ICCG using mapped
   active-pocket N/CA/C atoms, with a proper rigid transform. Apply exactly one
   rigid transform to each complete 54-atom ligand. Do not move ICCG heavy atoms
   and do not minimize ligand or protein heavy coordinates.
3. Build LG1 and LG2 as a strict pair: same protein atom order, same ligand atom
   names/order, same total atom count, same QM atom selection, and only
   state-specific ligand coordinates.
4. Protonate and parameterize the pair with one shared Amber topology. Retain
   catalytic Ser165-H in both Step1 states. Use ff14SB for protein and the
   released literature ligand parameters; no invented generic ligand charges.
5. Use mbondi3 radii and Amber18 DFTB3/3OB-3-1 with GBN2
   (`igb=8`, `qmgb=2`, `gbsa=1`, `saltcon=0.10`) for the fast paired
   pilot. This is the high-throughput pilot method, not a PBE-D3 substitute.
6. The formal QM region is the complete 54-atom ligand plus complete sidechains
   of the five paper-equivalent pocket residues, cut at CA-CB with link H:
   ICCG Trp164, Ser165, Asp210, His242, and Ile243. Ile243 is the three-dimensional
   structural equivalent of paper Ser238: transformed CA separation 0.939 A.
   Record this non-identical structural mapping explicitly. Never silently omit
   the fifth sidechain or fall back to QM65.
7. Set QM charge to -1 and singlet multiplicity. Audit the actual QM atom count,
   link boundaries, residue identities, atom names, and total charge rather than
   assuming the WT paper count of 104 atoms, because Ser238->Ile243 changes the
   sidechain composition.
8. Generate frozen-geometry single-point inputs only. No optimization, dynamics,
   umbrella sampling, PMF, or barrier claim.

## Hard gates before submission

For each state and for the pair, produce machine-readable audits and stop before
Slurm submission unless all gates pass:

- exact input provenance and SHA-256 hashes;
- 258 protein residues and Ser165 OG/HG present;
- 54 ligand atoms and 32 ligand heavy atoms;
- two six-membered carbon rings;
- ligand atom names/order identical between LG1/LG2;
- total topology atom count and ordering identical between states;
- all five mapped QM sidechains present and all ligand atoms in QM;
- QM charge -1, singlet, mbondi3, GBN2, DFTB3/3OB-3-1;
- no NaN/Inf or duplicate atoms;
- no nonbonded heavy-atom distance below 1.20 A;
- maximum vdW overlap <= 0.80 A, excluding true topology bonds and declared
  link-boundary pairs;
- catalytic geometry fields non-empty, including Ser165 OG to reactive ligand
  carbon distance and the Step1 RC components used by the literature;
- paired heavy-coordinate invariance for the protein;
- `sander.MPI -O` input smoke validation succeeds.

A gate failure must be recorded as `NOT_SUBMITTED_<REASON>` or
`FAIL_GEOMETRY_CLASH_NOT_LABEL`. Exit 0 alone is not scientific PASS.

## Post-run scientific PASS

The two-state Slurm pilot may be submitted only after the preflight audit PASS.
Each state requires normal Amber termination, converged DFTB SCC, a finite
`DFTBESCF`, finite `EGB`, finite total energy, and a post-run geometry audit
that still passes. Only then write `PASS.json` for that state. Write one pair
summary only if both states pass; label the energy difference as a fixed-geometry
electronic/implicit-solvent energy difference, never a free-energy barrier.

## Required files

Create exactly these source-controlled files:

- `build_iccg_step1_pair.py`
- `audit_iccg_step1_pair.py`
- `test_iccg_step1_pair.py`
- `run_iccg_step1_pair.sbatch`
- `collect_iccg_step1_pair.py`
- `deploy_and_run_remote.sh`
- `RUNBOOK.md`
- `RUN_HISTORY.tsv`

Generated topologies, coordinates, logs, energies, and JSON run artifacts belong
only in the remote project directory
`/work/home/acshdt1dks/iccg_qmmm_step1_pair_20260723`, not in GitHub.

## TDD and verification

Write failing tests first and record the RED command/status in RUN_HISTORY.tsv.
Implement the smallest code satisfying the contract. Tests must cover at least:

- structural residue mapping fixes paper Ser238 to ICCG Ile243 and rejects the
  sequence-gap null mapping;
- strict LG1/LG2 atom-order and topology pairing;
- five-residue/full-ligand QM selection;
- bonded exclusions versus nonbonded geometry;
- prevention of submission on any hard-gate failure;
- PASS.json cannot be written from mere process exit 0.

Run `python test_iccg_step1_pair.py` in Codex Web. Commit the completed files to
the current branch. Report the commit and test status, but do not claim that the
remote scientific calculation ran; the desktop continuation will deploy and
submit it separately after reviewing the generated scripts.
