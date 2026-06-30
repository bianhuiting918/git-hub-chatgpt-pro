# 2026-06-30 Blind PETase QM/MM Mechanism Task

## Status

A blind first-principles PETase QM/MM mechanism task was created and uploaded. Stages 1 through 9 now have blind workflow scaffolds, manifests, protocols, generator tests, a Stage 1/2 gate runner for compute-server execution, and a Stage 1 pose-generation queue generator. The actual ligand 3D structures, protonated production coordinates, Michaelis complexes, classical MD trajectories, productive conformer clusters, low-cost QM/MM scans, refined TS structures, committors, PMFs, barriers, and final paper comparison are not yet complete.

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
- `blind_work/01_system_setup/pose_generation_queue.tsv`
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

Stage 2 classical-MD ensemble gate:

- `blind_work/02_classical_md/md_replicate_queue.tsv`
- `blind_work/02_classical_md/productive_conformer_manifest.tsv`
- `blind_work/02_classical_md/rejected_pose_manifest.tsv`
- `blind_work/02_classical_md/stage2_classical_md_protocol.md`

Stage 3/4 mechanism and low-cost QM/MM scaffolds:

- `blind_work/03_mechanism_tree/mechanism_hypotheses.yaml`
- `blind_work/03_mechanism_tree/candidate_cv_sets.tsv`
- `blind_work/03_mechanism_tree/stage3_stage4_mechanism_qmmm_protocol.md`
- `blind_work/04_qmmm_exploration/path_screening_table.tsv`
- `blind_work/04_qmmm_exploration/ts_like_guess_manifest.tsv`

Stage 5/6/7 TS and committor scaffolds:

- `blind_work/05_ts_refinement/ts_refinement_manifest.tsv`
- `blind_work/05_ts_refinement/stage5_stage6_ts_refinement_protocol.md`
- `blind_work/06_ts_ensemble/acylation_ts_ensemble.tsv`
- `blind_work/06_ts_ensemble/deacylation_ts_ensemble.tsv`
- `blind_work/07_committor/committor_queue.tsv`

Stage 8/9 free-energy and final-validation scaffolds:

- `blind_work/08_free_energy/acylation_pmf.tsv`
- `blind_work/08_free_energy/deacylation_pmf.tsv`
- `blind_work/08_free_energy/barrier_summary.md`
- `blind_work/09_paper_validation/blind_vs_paper_comparison.md`
- `blind_work/09_paper_validation/discrepancy_audit.md`

Scripts and tests:

- `scripts/download_stage1_rcsb_structures.ps1`
- `scripts/identify_ser_his_asp_triads.py`
- `scripts/audit_stage1_pdb_preparation.py`
- `scripts/prepare_stage1_initial_structures_v2.py`
- `scripts/write_ligand_build_environment_report.py`
- `scripts/scan_stage1_protonation_sites.py`
- `scripts/probe_stage1_compute_environment.sh`
- `scripts/build_stage1_ligands_rdkit.py`
- `scripts/generate_stage1_pose_generation_queue.py`
- `scripts/run_stage1_protonation_gate.sh`
- `scripts/score_stage1_pose_geometry.py`
- `scripts/run_blind_stage1_stage2_gates.py`
- `scripts/generate_stage2_classical_md_manifests.py`
- `scripts/generate_stage3_mechanism_tree.py`
- `scripts/generate_stage5_ts_manifests.py`
- `scripts/generate_stage8_stage9_manifests.py`
- `tests/test_generate_stage1_pose_generation_queue.py`
- `tests/test_run_blind_stage1_stage2_gates.py`
- `tests/test_generate_stage2_classical_md_manifests.py`
- `tests/test_generate_stage3_mechanism_tree.py`
- `tests/test_generate_stage5_ts_manifests.py`
- `tests/test_generate_stage8_stage9_manifests.py`

## Stage 1 progress

RCSB metadata was queried on 2026-06-30 for PETase structure candidates. The initial blind production template set is primary `6EQE`, secondary sensitivity templates `5XJH`, `5YFE`, `6ILW`, and backups `6EQD`, `6EQF`, `6EQG`, `6EQH`, `6QGC`.

The selected WT-like/backup RCSB coordinate files were downloaded locally. A coordinate-only Ser-His-Asp geometric scan found the same active-site triad geometry across the WT-like templates, with `6EQE`, `5XJH`, and `6ILW` using chain A `SER160-HIS237-ASP206`, and `5YFE` using a numbering-shifted chain A `SER134-HIS211-ASP180`.

Initial cleaned chain-A PDB files were generated locally from RCSB inputs. The manifest records nine generated files, retained water counts, altloc decision counts, SHA256 hashes, and explicit `not_assigned` protonation status. A validation check found zero nonblank altloc indicators in the v2 cleaned PDB files.

Protonation setup was started from cleaned coordinates. For primary template `6EQE`, catalytic `ASP206` is deprotonated in the primary hypothesis but has a conditional neutral-Asp sensitivity gate; catalytic `HIS237` requires HID/HIE tautomer sensitivity before QM/MM reaction scans.

Blind substrate model definitions were added for `PET_dimer_capped`, `BHET_like`, `MHET_like`, and `MHET_like_acyl_enzyme_precursor`.

Executable Stage 1 gates were added for ligand conformer generation, protonation review, pose generation, and pose geometry scoring. Local verification showed the RDKit ligand script compiles and writes `blocked_missing_rdkit` in the local missing-RDKit environment instead of fabricating coordinates. The pose geometry scorer compiles locally. Bash validation for `run_stage1_protonation_gate.sh` remains server-side because local Windows lacks `bash`.

A Stage 1/2 compute-server runner was added. It records gate status to `blind_work/00_run_status/stage1_stage2_gate_status.tsv` and next actions to `blind_work/00_run_status/stage1_stage2_next_actions.md`; it refuses to queue Stage 2 MD without an accepted GS pose.

A TDD check was used for `run_blind_stage1_stage2_gates.py`:

- Red: `test_run_blind_stage1_stage2_gates.py` failed because the runner did not exist and no status file was produced.
- Green: after implementation and skip-probe status ordering correction, the test passed with bundled Python: `Ran 1 test in 0.275s OK`.
- Syntax check: `python -m py_compile work/run_blind_stage1_stage2_gates.py` exited 0.

A TDD check was used for `generate_stage1_pose_generation_queue.py`:

- Red: `test_generate_stage1_pose_generation_queue.py` failed because the pose queue generator did not exist.
- Green: after implementation, the test passed with bundled Python: `Ran 1 test in 0.508s OK`.
- The generator builds docking queue rows and Vina-style config files from prepared structures, Ser-His-Asp triad geometry, and ligand atom-label manifests. The docking box center is computed from catalytic Ser OG coordinates in the prepared structure, not from paper coordinates or paper trajectories.
- With the current local missing-RDKit ligand manifest, the generated `pose_generation_queue.tsv` remains `not_ready` with `no_ready_ligand_build_manifest`; no pose or docking result was fabricated.

## Stage 2 progress

A TDD check was used for `generate_stage2_classical_md_manifests.py`:

- Red: `test_generate_stage2_classical_md_manifests.py` failed because the generator did not exist.
- Green: after implementation and one boundary-text correction, the test passed with bundled Python: `Ran 1 test in 0.605s OK`.

The generator accepts only Stage 1 GS pose rows with `pass_fail=pass` and a non-pending relaxed structure path. It writes an MD replicate queue, productive conformer manifest, rejected pose manifest, and Stage 2 classical-MD protocol. Because Stage 1 has no accepted GS pose yet, the generated Stage 2 queue remains `not_ready_no_accepted_gs_pose` and no productive conformers are claimed.

## Stage 3/4 progress

A TDD check was used for `generate_stage3_mechanism_tree.py`:

- Red: `test_generate_stage3_mechanism_tree_unittest.py` failed because the generator did not exist.
- Green: after implementation, the test passed with bundled Python: `Ran 1 test in 0.465s OK`.

The generator produced blind mechanism hypotheses, generic CV candidates, Stage 4 path screening rows, and a placeholder TS-like guess manifest that must not be populated from paper TS coordinates or paper trajectories.

## Stage 5/6/7 progress

A TDD check was used for `generate_stage5_ts_manifests.py`:

- Red: `test_generate_stage5_ts_manifests.py` failed because the generator did not exist.
- Green: after implementation, the test passed with bundled Python: `Ran 1 test in 0.386s OK`.

The generator accepts only Stage 4 TS-like guesses with `validation_status=endpoint_checked`, real structure paths, and non-placeholder IDs. Because the current Stage 4 manifest is still a placeholder, the generated Stage 5/6/7 manifests remain empty or `pending/not_generated`. No TS coordinates or TS claims were fabricated.

## Stage 8/9 progress

A TDD check was used for `generate_stage8_stage9_manifests.py`:

- Red: `test_generate_stage8_stage9_manifests.py` failed because the generator did not exist.
- Green: after implementation and one boundary-text correction, the test passed with bundled Python: `Ran 1 test in 0.292s OK`.

The generator accepts only accepted TS entries from Stage 7. It writes acylation/deacylation PMF manifests, a barrier summary scaffold, and final paper-validation documents. Because there is no accepted TS manifest yet, the generated PMF and validation outputs correctly remain `not_ready` and contain no barrier values or paper-derived comparison results.

## Compute note

A previous Zenodo paper-trajectory download on the CPU server was stopped after the boundary correction, because paper shooting trajectories are concrete results and should not enter the blind workflow. The partial archive fragments were not deleted.

Non-interactive SSH to the known compute server was tested without embedding a password and failed with `Permission denied (publickey,password)`. Running the environment and chemistry gates requires interactive login or key-based authentication.

## Next action

Run this from the repository root on the compute server:

```text
python project01/phase2_blind_petase_qmmm_20260630/scripts/run_blind_stage1_stage2_gates.py
```

If the runner reports a blocked ligand/protonation/pose/GS-pose gate, fix that upstream input or toolchain issue first. Only after an accepted GS pose exists should the runner queue Stage 2 classical MD; only after productive conformers exist should Stage 4 low-cost QM/MM scans be activated.