# 2026-06-30 Blind PETase QM/MM Mechanism Task

## Status

A blind first-principles PETase QM/MM mechanism task was created and uploaded. Stage 1 has structure, ligand, protonation, and pose-gate scaffolding. Stage 3/4 now has a blind mechanism hypothesis tree, candidate CV table, path screening manifest, and TS-like guess manifest. The actual ligand 3D structures, protonated production coordinates, Michaelis complexes, low-cost QM/MM scans, refined TS structures, committors, and PMFs are not yet complete.

## Scope correction

The task is not to reproduce the paper by consuming the paper's concrete trajectory/TS/RC results. The task is to reproduce the scientific discovery process from PETase structure and substrate chemistry, using only the paper's broad methodology as inspiration.

Before final validation, do not use paper TS coordinates, paper reaction-coordinate formulas, paper selected CVs, paper umbrella windows, paper aimless-shooting trajectories, paper barrier values, rate constants, rate-limiting-step assignment, or paper conclusions about His/Asp/Trp/oxyanion-hole roles.

Allowed inputs are PETase structures, substrate chemistry, generic serine-hydrolase chemistry, QM/MM methods, TS search, PMF, committor, and trajectory-validation ideas.

## Uploaded documents and scripts

Core entry points:

- `projects/01-specialized-ts-aware-scorer/docs/petase_blind_qmmm_mechanism_plan.md`
- `project01/phase2_blind_petase_qmmm_20260630/README.md`

Stage 1 setup and gates:

- `blind_work/01_system_setup/structure_selection.tsv`
- `blind_work/01_system_setup/stage1_system_setup_protocol.md`
- `blind_work/01_system_setup/stage1_ligand_and_protonation_execution_protocol.md`
- `blind_work/01_system_setup/stage1_pose_geometry_filter_protocol.md`
- `blind_work/01_system_setup/gs_pose_manifest.tsv`
- `blind_work/01_system_setup/rejected_pose_manifest.tsv`
- `blind_work/01_system_setup/ligand_model_manifest.tsv`
- `blind_work/01_system_setup/ligand_model_definitions.md`
- `blind_work/01_system_setup/ligand_smiles.tsv`
- `blind_work/01_system_setup/ligand_build_environment_report.md`
- `blind_work/01_system_setup/stage1_environment_enablement_status.md`
- `blind_work/01_system_setup/stage1_remote_execution_instructions.md`
- `blind_work/01_system_setup/structure_download_manifest.csv`
- `blind_work/01_system_setup/ser_his_asp_triad_candidates.tsv`
- `blind_work/01_system_setup/pdb_preparation_log.md`
- `blind_work/01_system_setup/pdb_preparation_audit.tsv`
- `blind_work/01_system_setup/template_chain_decisions.tsv`
- `blind_work/01_system_setup/retained_water_candidates.tsv`
- `blind_work/01_system_setup/prepared_structure_manifest.tsv`
- `blind_work/01_system_setup/altloc_resolution_decisions.tsv`
- `blind_work/01_system_setup/protonation_site_scan.tsv`
- `blind_work/01_system_setup/protonation_hypothesis_manifest.tsv`
- `blind_work/01_system_setup/protonation_setup_notes.md`

Stage 3/4 mechanism and low-cost QM/MM scaffolds:

- `blind_work/03_mechanism_tree/mechanism_hypotheses.yaml`
- `blind_work/03_mechanism_tree/candidate_cv_sets.tsv`
- `blind_work/03_mechanism_tree/stage3_stage4_mechanism_qmmm_protocol.md`
- `blind_work/04_qmmm_exploration/path_screening_table.tsv`
- `blind_work/04_qmmm_exploration/ts_like_guess_manifest.tsv`

Scripts and tests:

- `scripts/download_stage1_rcsb_structures.ps1`
- `scripts/identify_ser_his_asp_triads.py`
- `scripts/audit_stage1_pdb_preparation.py`
- `scripts/prepare_stage1_initial_structures_v2.py`
- `scripts/write_ligand_build_environment_report.py`
- `scripts/scan_stage1_protonation_sites.py`
- `scripts/probe_stage1_compute_environment.sh`
- `scripts/build_stage1_ligands_rdkit.py`
- `scripts/run_stage1_protonation_gate.sh`
- `scripts/score_stage1_pose_geometry.py`
- `scripts/generate_stage3_mechanism_tree.py`
- `tests/test_generate_stage3_mechanism_tree.py`

## Stage 1 progress

RCSB metadata was queried on 2026-06-30 for PETase structure candidates. The initial blind production template set is primary `6EQE`, secondary sensitivity templates `5XJH`, `5YFE`, `6ILW`, and backups `6EQD`, `6EQF`, `6EQG`, `6EQH`, `6QGC`.

The selected WT-like/backup RCSB coordinate files were downloaded locally. A coordinate-only Ser-His-Asp geometric scan found the same active-site triad geometry across the WT-like templates, with `6EQE`, `5XJH`, and `6ILW` using chain A `SER160-HIS237-ASP206`, and `5YFE` using a numbering-shifted chain A `SER134-HIS211-ASP180`.

Initial cleaned chain-A PDB files were generated locally from RCSB inputs. The manifest records nine generated files, retained water counts, altloc decision counts, SHA256 hashes, and explicit `not_assigned` protonation status. A validation check found zero nonblank altloc indicators in the v2 cleaned PDB files.

Protonation setup was started from cleaned coordinates. For primary template `6EQE`, catalytic `ASP206` is deprotonated in the primary hypothesis but has a conditional neutral-Asp sensitivity gate; catalytic `HIS237` requires HID/HIE tautomer sensitivity before QM/MM reaction scans.

Blind substrate model definitions were added for `PET_dimer_capped`, `BHET_like`, `MHET_like`, and `MHET_like_acyl_enzyme_precursor`.

Executable Stage 1 gates were added for ligand conformer generation, protonation review, and pose geometry scoring. Local verification showed the RDKit ligand script compiles and writes `blocked_missing_rdkit` in the local missing-RDKit environment instead of fabricating coordinates. The pose geometry scorer compiles locally. Bash validation for `run_stage1_protonation_gate.sh` remains server-side because local Windows lacks `bash`.

## Stage 3/4 progress

A TDD check was used for `generate_stage3_mechanism_tree.py`:

- Red: `test_generate_stage3_mechanism_tree_unittest.py` failed because the generator did not exist.
- Green: after implementation, the test passed with bundled Python: `Ran 1 test in 0.465s OK`.

The generator produced:

- `mechanism_hypotheses.yaml`: blind acylation/deacylation path hypotheses and explicit forbidden inputs.
- `candidate_cv_sets.tsv`: generic bond-formation, bond-breaking, proton-transfer, and protonation-sensitivity CV candidates sourced only from generic serine-hydrolase chemistry.
- `path_screening_table.tsv`: Stage 4 low-cost QM/MM screening rows, all `not_started` until accepted poses exist.
- `ts_like_guess_manifest.tsv`: placeholder manifest that must not be populated from paper TS coordinates or paper trajectories.

The Stage 3/4 protocol defines grill gates: paths must map to accepted Stage 1/2 structures or documented sensitivity branches; protonation branches need protonation-gate evidence; water-network branches need classical ensemble evidence; TS-like guesses can only come from scan maxima, string images, or saddle-like structures generated by this blind workflow.

## Compute note

A previous Zenodo paper-trajectory download on the CPU server was stopped after the boundary correction, because paper shooting trajectories are concrete results and should not enter the blind workflow. The partial archive fragments were not deleted.

Non-interactive SSH to the known compute server was tested without embedding a password and failed with `Permission denied (publickey,password)`. Running the environment and chemistry gates requires interactive login or key-based authentication.

## Next action

Continue by running `probe_stage1_compute_environment.sh` on the compute server. If RDKit and pKa/protonation tools are available, run `build_stage1_ligands_rdkit.py`, `run_stage1_protonation_gate.sh`, build accepted/rejected Michaelis or acyl-enzyme poses, score them with `score_stage1_pose_geometry.py`, then activate supported Stage 3 paths in `path_screening_table.tsv` for low-cost QM/MM scans.
