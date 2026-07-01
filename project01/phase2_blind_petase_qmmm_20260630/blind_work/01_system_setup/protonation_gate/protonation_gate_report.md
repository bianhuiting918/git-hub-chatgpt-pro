# Stage 1 Protonation Gate Report

Generated: 2026-07-01T04:47:07+08:00

Input PDB: project01/phase2_blind_petase_qmmm_20260630/blind_work/01_system_setup/prepared_initial_pdbs/6EQE_chainA_initial_clean_v2.pdb
Input SHA256: 3ec12a00836d5b5525ae78f7db1321fb19811d619c9a887be0726f5e4e04c543

Boundary: this gate uses the cleaned PETase structure and pH/protonation tools only. It does not use paper TS structures, paper reaction coordinates, selected CVs, trajectories, barriers, or mechanism conclusions.

## Required Manual Review

- Catalytic Asp state.
- Catalytic His tautomer/protonation.
- Nearby Cys/disulfide consistency.
- Remote His tautomer naming for topology.
- Any branch that differs from the primary hypothesis manifest.

## Probe Table

```tsv
category	tool	status	path_or_version	command	output_path	output_sha256	note
input	pdb	available	project01/phase2_blind_petase_qmmm_20260630/blind_work/01_system_setup/prepared_initial_pdbs/6EQE_chainA_initial_clean_v2.pdb	provided_input	project01/phase2_blind_petase_qmmm_20260630/blind_work/01_system_setup/prepared_initial_pdbs/6EQE_chainA_initial_clean_v2.pdb	3ec12a00836d5b5525ae78f7db1321fb19811d619c9a887be0726f5e4e04c543	blind_stage1_clean_structure
pka	propka3	ran_check_output	propka3 3.5.1	propka3 project01/phase2_blind_petase_qmmm_20260630/blind_work/01_system_setup/prepared_initial_pdbs/6EQE_chainA_initial_clean_v2.pdb	project01/phase2_blind_petase_qmmm_20260630/blind_work/01_system_setup/protonation_gate/propka_stdout_stderr.log	b5380ebdc18dc7c0b82454c29e67dea183217f3177862dab41e3b14fc02238d1	inspect pka file for catalytic Asp/His states
protonation	pdb2pqr	produced_output	pdb2pqr 3.6.1	pdb2pqr --ff=AMBER --with-ph=7.0 project01/phase2_blind_petase_qmmm_20260630/blind_work/01_system_setup/prepared_initial_pdbs/6EQE_chainA_initial_clean_v2.pdb project01/phase2_blind_petase_qmmm_20260630/blind_work/01_system_setup/protonation_gate/protonated_pdb2pqr_ph7.pqr	project01/phase2_blind_petase_qmmm_20260630/blind_work/01_system_setup/protonation_gate/protonated_pdb2pqr_ph7.pqr	474eb4383eab57ca754d34a0d3c91a216551954496ddb826db321fe55b9298ef	review residue states before topology generation
hydrogen	reduce	produced_output	reduce.4.10.230211	reduce -BUILD project01/phase2_blind_petase_qmmm_20260630/blind_work/01_system_setup/prepared_initial_pdbs/6EQE_chainA_initial_clean_v2.pdb > project01/phase2_blind_petase_qmmm_20260630/blind_work/01_system_setup/protonation_gate/reduce_build_h.pdb	project01/phase2_blind_petase_qmmm_20260630/blind_work/01_system_setup/protonation_gate/reduce_build_h.pdb	8feebd4c2e9efc2a47305164f3e0b2ae59a35049297551cb2d521d69bf3d01c4	compare His tautomer choices with pKa output
```
