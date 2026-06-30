# Blind PETase Stage 1 Pose Geometry Filter Protocol

Date: 2026-06-30

## Boundary

This filter scores Michaelis-complex candidates using generic serine-hydrolase geometry only. It does not use paper TS coordinates, paper reaction coordinates, selected CVs, shooting trajectories, umbrella windows, barriers, or mechanistic conclusions.

## Required Inputs

- A relaxed or restrained-relaxed complex PDB containing protein and named ligand atoms.
- A ligand label table produced by `build_stage1_ligands_rdkit.py`, such as `qm_atom_labels/PET_dimer_capped_atoms.tsv`.
- A selected ester `candidate_id` from the label table. The candidate must be selected from pose chemistry and not from paper data.

## Command

```bash
python3 project01/phase2_blind_petase_qmmm_20260630/scripts/score_stage1_pose_geometry.py \
  --complex-pdb work/blind_work/01_system_setup/poses/AC_GS_0001_relaxed.pdb \
  --label-tsv work/blind_work/01_system_setup/ligand_build/qm_atom_labels/PET_dimer_capped_atoms.tsv \
  --model-id PET_dimer_capped \
  --candidate-id E01 \
  --pose-id AC_GS_0001 \
  --template-pdb 6EQE \
  --generation-method constrained_docking \
  --out-tsv work/blind_work/01_system_setup/gs_pose_manifest.tsv
```

For numbering-shifted templates such as `5YFE`, pass explicit residue IDs:

```bash
--ser-og A:134:OG --his A:211 --trp A:159
```

## Filters

The script writes the same fields as `gs_pose_manifest.tsv` and applies these broad Stage 1 checks:

- Ser O-gamma to ester carbonyl C <= 3.6 A.
- O-gamma - carbonyl C - carbonyl O attack angle between 80 and 125 degrees.
- Carbonyl O has at least one protein N/O donor within 3.2 A.
- Catalytic His ring N is within 3.4 A of Ser O-gamma.
- Leaving ester O is within 4.0 A of catalytic His ring N.

## Grill Gate

A pose can move to classical MD only if it passes the automated filters and the operator records:

1. Which ligand ester candidate was selected and why.
2. Whether the pose depends on a single retained crystal water.
3. Whether the same substrate construction is still plausible on at least one secondary template.
4. Whether a rejected pose was written to `rejected_pose_manifest.tsv` instead of being silently discarded.

The filter is a triage tool, not a mechanistic result.
