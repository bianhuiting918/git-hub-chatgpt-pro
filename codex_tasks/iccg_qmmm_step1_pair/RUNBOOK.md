# ICCG Step1 LG1/LG2 QM/MM Pair Gate Runbook

This directory contains only source-controlled builders, auditors, tests, and
run instructions. Generated Amber topologies, coordinates, logs, energies, and
JSON run artifacts must stay in `/work/home/acshdt1dks/iccg_qmmm_step1_pair_20260723`.

## Scope

- Build one frozen-geometry Step1 LG1/LG2 paired pilot.
- Use one shared Amber topology and state-specific ligand coordinates.
- Use ff14SB, released literature ligand parameters, mbondi3 radii, and
  Amber18 DFTB3/3OB-3-1 with GBN2 (`igb=8`, `qmgb=2`, `gbsa=1`, `saltcon=0.10`).
- Do not claim a free-energy barrier; any reported difference is a fixed-geometry
  electronic/implicit-solvent energy difference.

## Hard stops

Any failed geometry, topology, QM-region, provenance, or smoke-validation gate
must write `NOT_SUBMITTED_<REASON>` or `FAIL_GEOMETRY_CLASH_NOT_LABEL` and must
not submit Slurm. A process exit code of 0 is not a scientific PASS.

## Remote deployment

Review the scripts, then from this directory on the scientific server run:

```bash
./deploy_and_run_remote.sh
```

The deployment script stages source files to the remote project directory,
builds explicit-provenance manifests from CLI paths, runs preflight submission
checks, and submits `run_iccg_step1_pair.sbatch` only when preflight status is
`PASS`.

## PASS criteria

Each state requires normal Amber termination, converged DFTB SCC, finite
`DFTBESCF`, finite `EGB`, finite total energy, and a passing post-run geometry
audit before a state `PASS.json` can be written. The pair summary is permitted
only after both state PASS files exist.

## Remote CYX Stage-B notes

The validated remote CYX topology has two disulfides (Cys238-Cys283 and
Cys275-Cys292), 258 protein residues plus one ligand residue, 3848 protein
atoms, 3902 total atoms, ligand54, QM atom count 100, system charge about +6,
qmcharge -1, His242 as HID with HD1/no HE2, and Ser165 HG present. Reported
GMAX values are diagnostics only, not low-gradient hard gates. The hard geometry
gate remains protein-ligand heavy-atom max vdW overlap <= 0.80 A, with verbose
SCC convergence and finite DFTBESCF/EGB/ENERGY required for PASS.json.
