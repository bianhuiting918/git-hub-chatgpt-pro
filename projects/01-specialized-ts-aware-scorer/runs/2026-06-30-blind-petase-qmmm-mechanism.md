# 2026-06-30 Blind PETase QM/MM Mechanism Task

## Status

A blind first-principles PETase QM/MM mechanism task was created and uploaded. Stage 1 has progressed from structure-only setup to executable environment, ligand, and protonation gates. The actual ligand 3D structures, protonated production coordinates, Michaelis complexes, QM/MM paths, and TS ensembles are not yet complete.

## Scope correction

The task is not to reproduce the paper by consuming the paper's concrete trajectory/TS/RC results. The task is to reproduce the scientific discovery process from PETase structure and substrate chemistry, using only the paper's broad methodology as inspiration.

## Uploaded documents

- `projects/01-specialized-ts-aware-scorer/docs/petase_blind_qmmm_mechanism_plan.md`
- `project01/phase2_blind_petase_qmmm_20260630/README.md`
- `project01/phase2_blind_petase_qmmm_20260630/blind_work/01_system_setup/structure_selection.tsv`
- `project01/phase2_blind_petase_qmmm_20260630/blind_work/01_system_setup/stage1_system_setup_protocol.md`
- `project01/phase2_blind_petase_qmmm_20260630/blind_work/01_system_setup/stage1_ligand_and_protonation_execution_protocol.md`
- `project01/phase2_blind_petase_qmmm_20260630/blind_work/01_system_setup/gs_pose_manifest.tsv`
- `project01/phase2_blind_petase_qmmm_20260630/blind_work/01_system_setup/rejected_pose_manifest.tsv`
- `project01/phase2_blind_petase_qmmm_20260630/blind_work/01_system_setup/ligand_model_manifest.tsv`
- `project01/phase2_blind_petase_qmmm_20260630/blind_work/01_system_setup/ligand_model_definitions.md`
- `project01/phase2_blind_petase_qmmm_20260630/blind_work/01_system_setup/ligand_smiles.tsv`
- `project01/phase2_blind_petase_qmmm_20260630/blind_work/01_system_setup/ligand_build_environment_report.md`
- `project01/phase2_blind_petase_qmmm_20260630/blind_work/01_system_setup/stage1_environment_enablement_status.md`
- `project01/phase2_blind_petase_qmmm_20260630/blind_work/01_system_setup/stage1_remote_execution_instructions.md`
- `project01/phase2_blind_petase_qmmm_20260630/blind_work/01_system_setup/structure_download_manifest.csv`
- `project01/phase2_blind_petase_qmmm_20260630/blind_work/01_system_setup/ser_his_asp_triad_candidates.tsv`
- `project01/phase2_blind_petase_qmmm_20260630/blind_work/01_system_setup/pdb_preparation_log.md`
- `project01/phase2_blind_petase_qmmm_20260630/blind_work/01_system_setup/pdb_preparation_audit.tsv`
- `project01/phase2_blind_petase_qmmm_20260630/blind_work/01_system_setup/template_chain_decisions.tsv`
- `project01/phase2_blind_petase_qmmm_20260630/blind_work/01_system_setup/retained_water_candidates.tsv`
- `project01/phase2_blind_petase_qmmm_20260630/blind_work/01_system_setup/prepared_structure_manifest.tsv`
- `project01/phase2_blind_petase_qmmm_20260630/blind_work/01_system_setup/altloc_resolution_decisions.tsv`
- `project01/phase2_blind_petase_qmmm_20260630/blind_work/01_system_setup/protonation_site_scan.tsv`
- `project01/phase2_blind_petase_qmmm_20260630/blind_work/01_system_setup/protonation_hypothesis_manifest.tsv`
- `project01/phase2_blind_petase_qmmm_20260630/blind_work/01_system_setup/protonation_setup_notes.md`
- `project01/phase2_blind_petase_qmmm_20260630/scripts/download_stage1_rcsb_structures.ps1`
- `project01/phase2_blind_petase_qmmm_20260630/scripts/identify_ser_his_asp_triads.py`
- `project01/phase2_blind_petase_qmmm_20260630/scripts/audit_stage1_pdb_preparation.py`
- `project01/phase2_blind_petase_qmmm_20260630/scripts/prepare_stage1_initial_structures_v2.py`
- `project01/phase2_blind_petase_qmmm_20260630/scripts/write_ligand_build_environment_report.py`
- `project01/phase2_blind_petase_qmmm_20260630/scripts/scan_stage1_protonation_sites.py`
- `project01/phase2_blind_petase_qmmm_20260630/scripts/probe_stage1_compute_environment.sh`
- `project01/phase2_blind_petase_qmmm_20260630/scripts/build_stage1_ligands_rdkit.py`
- `project01/phase2_blind_petase_qmmm_20260630/scripts/run_stage1_protonation_gate.sh`

## Boundary

Before final validation, do not use:

- paper TS coordinates;
- paper reaction-coordinate formulas;
- paper selected CVs;
- paper umbrella windows;
- paper aimless-shooting trajectories;
- paper barrier values, rate constants, or rate-limiting-step assignment;
- paper conclusions about His motion, Asp role, Trp185 role, oxyanion-hole role, or tetrahedral TS/intermediate status.

Use only:

- PETase structures;
- substrate chemistry;
- generic serine-hydrolase chemistry;
- QM/MM, TS-search, PMF, committor, and trajectory-validation methods.

## Compute note

A previous Zenodo paper-trajectory download on the CPU server was stopped after the boundary correction, because paper shooting trajectories are concrete results and should not enter the blind workflow. The partial archive fragments were not deleted.

Non-interactive SSH to the known compute server was tested without embedding a password and failed with `Permission denied (publickey,password)`. A remote-safe environment probe script and execution instructions were added; running them requires interactive login or key-based authentication.

## Stage 1 progress

RCSB metadata was queried on 2026-06-30 for PETase structure candidates. The initial blind production template set is:

- primary: `6EQE`;
- secondary sensitivity templates: `5XJH`, `5YFE`, `6ILW`;
- backup WT-like templates: `6EQD`, `6EQF`, `6EQG`, `6EQH`, `6QGC`.

Mutant, non-PETase, and failed pre-query IDs were explicitly marked as excluded in `structure_selection.tsv`.

The selected WT-like/backup RCSB coordinate files were downloaded locally. A coordinate-only Ser-His-Asp geometric scan found the same active-site triad geometry across the WT-like templates, with `6EQE`, `5XJH`, and `6ILW` using chain A `SER160-HIS237-ASP206`, and `5YFE` using a numbering-shifted chain A `SER134-HIS211-ASP180`.

A PDB preparation audit was completed. It records selected chain decisions, missing residues/atoms from PDB headers, alternate conformers, non-water heterogens, geometric disulfide candidates, water counts, and crystallographic waters within 4 A of catalytic hetero atoms.

Initial cleaned chain-A PDB files were generated locally from RCSB inputs. The manifest records nine generated files, retained water counts, altloc decision counts, SHA256 hashes, and explicit `not_assigned` protonation status. A validation check found zero nonblank altloc indicators in the v2 cleaned PDB files.

Protonation setup was started from cleaned coordinates. The scan records all titratable residues near the catalytic triad and all histidines requiring tautomer assignment. For primary template `6EQE`, catalytic `ASP206` is deprotonated in the primary hypothesis but has a conditional neutral-Asp sensitivity gate; catalytic `HIS237` requires HID/HIE tautomer sensitivity before QM/MM reaction scans.

Blind substrate model definitions were added for:

- `PET_dimer_capped` as the primary PET-like acylation substrate;
- `BHET_like` as a small neutral docking/pose control;
- `MHET_like` as a product-side/deacylation reference fragment with explicit pH-7 carboxylate state;
- `MHET_like_acyl_enzyme_precursor` as a protein-covalent Ser160 acyl-enzyme precursor for deacylation setup.

Executable Stage 1 gates were added:

- `build_stage1_ligands_rdkit.py` reads `ligand_smiles.tsv`, builds RDKit conformers when RDKit is present, writes SDF/PDB outputs, and enumerates ester atom-label candidates without choosing a paper-derived scissile ester.
- `run_stage1_protonation_gate.sh` records input hashes and pKa/protonation/hydrogen-tool outputs for manual catalytic Asp/His review.
- `stage1_ligand_and_protonation_execution_protocol.md` defines pass/fail grill gates for ligand, protonation, and pose acceptance.
- `rejected_pose_manifest.tsv` is now present so failed Michaelis candidates must be recorded rather than silently discarded.

Local verification completed:

- `build_stage1_ligands_rdkit.py` passed Python compile check with bundled Python.
- Running it in the local missing-RDKit environment wrote a `blocked_missing_rdkit` manifest instead of generating fake ligand coordinates.
- Bash syntax validation was not possible locally because `bash` is not installed in the Windows environment; the shell script must be checked on the Linux compute server.

Current structure-prep decisions:

- use chain A for the initial blind setup across selected templates;
- keep geometric disulfide candidates unless preparation software contradicts connectivity;
- retain listed catalytic-site water candidates for first local relaxation, then test sensitivity;
- use v2 altloc decisions before substrate placement;
- preserve scissile ester atom labels through ligand 3D generation and topology conversion;
- run external pKa/protonation tooling before production topology generation and compare against `protonation_hypothesis_manifest.tsv`;
- do not accept docking/MD/QM-MM inputs until the compute environment probe, ligand build gate, and protonation gate record exact tool paths, versions, and hashes.

This is structure/substrate-derived preparation evidence and is not a paper-derived mechanism result.

## Next action

Continue Stage 1 by running `probe_stage1_compute_environment.sh` on the compute server through interactive/key-based login. If RDKit and pKa/protonation tools are available, run `build_stage1_ligands_rdkit.py` and `run_stage1_protonation_gate.sh`, then generate repaired/protonated coordinate files, ligand conformers/parameters, and accepted/rejected Michaelis-complex manifests.
