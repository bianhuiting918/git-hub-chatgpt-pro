# 2026-07-02 GitHub Upload Scope

Purpose: upload the PETase blind QM/MM reproduction task in a usable, auditable form without pushing large simulation artifacts.

## Include

- README and protocol/status Markdown files.
- Small TSV manifests that define gates, queues, and screening state.
- Python and shell scripts used to generate queues, reactive pose seeds, AMBER/Gromacs-CP2K inputs, and status audits.
- Unit tests for generated manifests and scripts.
- `.gitignore` rules that keep heavy Amber/ATESA outputs on the CPU server.

## Exclude

- Amber/ATESA trajectory and restart outputs: `.nc`, `.rst7`, `.mdcrd`, `.ncrst`.
- Topology/coordinate binaries: `.prmtop`, `.inpcrd`, `.parm7`.
- Large prepared system directories and live work directories under `blind_work/04_qmmm_*` and `blind_work/05_atesa_acylation`.
- Passwords, private keys, and server-local authentication data.

## Scientific State to Upload

The GitHub state should make clear that:

- The project is following a blind first-principles PETase mechanism workflow.
- rank022 is the active acylation route.
- `lam26.55625` and `lam26.6` are TS-near but product-biased seeds from independent QM/MM exploration.
- ATESA now runs numerically with fresh velocities, but no validated TS ensemble exists yet.
- `attempt011` centered on `lam26.6` also gave one-sided product-biased commitment; the next decision point is generation of a cleaner reactant-side TS-near seed before further ATESA scaling.
