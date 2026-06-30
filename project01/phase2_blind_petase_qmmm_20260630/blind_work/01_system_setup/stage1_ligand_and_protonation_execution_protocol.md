# Blind PETase Stage 1 Ligand And Protonation Execution Protocol

Date: 2026-06-30

## Boundary

This protocol advances the blind workflow from structure selection to executable ligand/protonation gates. It does not use the PETase paper's TS coordinates, reaction-coordinate formulas, selected CVs, umbrella windows, aimless-shooting trajectories, barrier values, rate constants, or residue/mechanism conclusions.

## Inputs

- `ligand_smiles.tsv`: structure-derived substrate model definitions and SMILES seeds.
- `prepared_structure_manifest.tsv`: cleaned structure hashes from Stage 1.
- `protonation_hypothesis_manifest.tsv`: pre-tool hypothesis table for catalytic and nearby titratable residues.
- `ligand_model_manifest.tsv`: required atom-label and validation targets.

## Ligand Build Command

Run only after the environment probe confirms RDKit is available:

```bash
python3 project01/phase2_blind_petase_qmmm_20260630/scripts/build_stage1_ligands_rdkit.py \
  --smiles-table project01/phase2_blind_petase_qmmm_20260630/blind_work/01_system_setup/ligand_smiles.tsv \
  --out-dir work/blind_work/01_system_setup/ligand_build \
  --conformers 20 \
  --seed 20260630
```

Expected outputs:

```text
work/blind_work/01_system_setup/ligand_build/ligand_build_manifest.tsv
work/blind_work/01_system_setup/ligand_build/qm_atom_labels/*_atoms.tsv
work/blind_work/01_system_setup/ligand_build/*/*_conformers.sdf
work/blind_work/01_system_setup/ligand_build/*/*_conf001_named_atoms.pdb
```

## Ligand Grill Gate

The ligand build passes only if:

1. `PET_dimer_capped`, `BHET_like`, and `MHET_like` produce 3D conformers.
2. The formal charge in `ligand_smiles.tsv` matches the generated molecule.
3. Ester candidates are enumerated from chemistry and not chosen from paper TS/RC data.
4. Atom names and RDKit indices are recorded in `qm_atom_labels/*_atoms.tsv`.
5. The selected scissile ester candidate is justified later by accepted Michaelis pose geometry, not by paper coordinates.

## Protonation Gate Command

Run this after recreating the v2 cleaned PDB files on the compute server:

```bash
bash project01/phase2_blind_petase_qmmm_20260630/scripts/run_stage1_protonation_gate.sh \
  work/blind_work/01_system_setup/prepared_initial_pdbs/6EQE_chainA_initial_clean_v2.pdb \
  work/blind_work/01_system_setup/protonation_gate/6EQE
```

Expected outputs:

```text
work/blind_work/01_system_setup/protonation_gate/6EQE/protonation_gate_report.tsv
work/blind_work/01_system_setup/protonation_gate/6EQE/protonation_gate_report.md
```

## Protonation Grill Gate

The protonation stage passes only if:

1. The report records the exact input PDB SHA256.
2. At least one pKa/protonation tool runs successfully, or the failure is recorded with tool/version evidence.
3. Catalytic Asp and His assignments are explicitly reviewed against `protonation_hypothesis_manifest.tsv`.
4. Any catalytic Asp/His branch disagreement is carried forward as a sensitivity branch.
5. No topology or docking step starts before this gate is resolved.

## Pose Gate

After ligand and protonation gates pass, generate candidate Michaelis complexes by at least two independent routes. Record accepted candidates in `gs_pose_manifest.tsv` and failures in `rejected_pose_manifest.tsv`.

Stage 2 remains blocked until the Stage 1 system setup protocol's accepted-pose count or written justification gate is satisfied.
