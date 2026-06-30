# Stage 1 Environment Enablement Status

Date: 2026-06-30

## Boundary

This status file concerns environment readiness for blind PETase Stage 1 only. It does not use paper TS coordinates, reaction-coordinate terms, selected CVs, paper trajectories, barriers, or mechanistic conclusions.

## Current Status

- Local bundled Python lacks RDKit, Open Babel, and Biopython.
- Local ligand 3D construction is therefore blocked to avoid producing unreliable coordinates.
- Non-interactive SSH probe to `bianht@210.73.40.29` failed with `Permission denied (publickey,password)`.
- No password was written into commands, scripts, GitHub files, logs, or terminal history.
- A remote-safe probe script has been added: `project01/phase2_blind_petase_qmmm_20260630/scripts/probe_stage1_compute_environment.sh`.
- Remote execution instructions have been added: `stage1_remote_execution_instructions.md`.

## Next Required Action

Run an interactive SSH session or configure key-based login, then execute the environment probe on the compute server.

The minimum usable Stage 1 environment needs:

- RDKit or Open Babel for ligand 3D generation and atom-label preservation;
- PROPKA, pdb2pqr, Amber reduce, H++, or equivalent for protonation/pKa handling;
- AmberTools or equivalent for ligand/protein parameter preparation;
- later, Amber/Sander QM/MM tooling for DFTB3/MM execution.

## Grill Gate

No docking, MD, or QM/MM input should be accepted until the environment probe output is committed or otherwise recorded with exact tool paths and versions.

