# PETase Blind QM/MM Stage 1 System Setup Protocol

Date: 2026-06-30

## Boundary

This stage must be performed without using the PETase paper's concrete mechanistic results. Do not use its TS structures, RC equations, selected CVs, umbrella windows, shooting trajectories, barrier values, rate-limiting assignment, or residue-role conclusions.

Allowed inputs are public PETase structures, substrate chemistry, generic serine-hydrolase chemistry, and standard QM/MM setup methods.

## Stage 1 Goal

Build chemically plausible WT PETase Michaelis-complex candidates for blind acylation and deacylation mechanism discovery.

The output of this stage is not a mechanism. It is a defensible ensemble of starting structures whose active-site geometry is good enough to enter classical MD and later QM/MM exploration.

## Selected Structural Templates

Primary template:

- `6EQE`: high-resolution WT-like PETase structure, 0.92 A.

Secondary WT-like templates for sensitivity checks:

- `5XJH`: WT PETase, 1.54 A.
- `5YFE`: WT-like PETase, 1.39 A.
- `6ILW`: WT PETase, 1.575 A.

Backup WT-like structures:

- `6EQD`, `6EQF`, `6EQG`, `6EQH`, `6QGC`.

Excluded from production WT setup:

- Mutants or analog-bound mutant controls: `5XH3`, `5XH2`, `5XFZ`, `5XFY`, `5YNS`, `6ILX`, `7SH6`.
- Non-PETase or failed pre-query IDs: `7SH7`, `7SH8`, `7SH9`, `7SHA`, `8JY2`, `8JY3`.

The full selection table is `structure_selection.tsv`.

## Deliverable Layout

Use this layout in the repository:

```text
project01/phase2_blind_petase_qmmm_20260630/blind_work/01_system_setup/
  structure_selection.tsv
  stage1_system_setup_protocol.md
  gs_pose_manifest.tsv
  rejected_pose_manifest.tsv
  pdb_preparation_log.md
  ligand_model_manifest.tsv
```

## Structure Preparation

For each included WT-like template:

1. Download mmCIF/PDB from RCSB.
2. Keep the biological monomer unless assembly inspection shows a required crystallographic interface near the active site.
3. Inspect missing residues and alternate conformers.
4. Retain crystallographic waters only if they hydrogen-bond to catalytic residues, oxyanion-hole donors, or the bound substrate region; otherwise mark them as removable candidates.
5. Assign disulfides from structure connectivity, then verify geometrically.
6. Assign protonation states at pH 7 using an automated predictor, followed by manual review of catalytic triad and nearby residues.

Record every decision in `pdb_preparation_log.md`.

## Substrate Models

Build at least three neutral or chemically justified substrate fragments before choosing one for production:

| model_id | purpose | required atoms near scissile ester |
| --- | --- | --- |
| `PET_dimer_capped` | Minimal PET chain mimic for acylation | ester carbonyl C/O, leaving alkoxy O, aromatic ring on both sides if retained |
| `BHET_like` | Small soluble PET-like ester control | same scissile ester, hydroxyethyl caps |
| `MHET_like` | Product-side and deacylation reference | ester/carboxylate state explicitly defined |

Each ligand model must include:

- 2D structure source or construction script.
- Protonation/charge state.
- MM parameter source.
- QM atom labels for the future QM/MM region.
- Atom names for scissile carbonyl C, carbonyl O, leaving O, and nucleophile-side O.

Record these in `ligand_model_manifest.tsv`.

## Michaelis Pose Generation

Generate candidate poses by at least two independent routes:

1. Docking or constrained docking into the active-site groove.
2. Local conformer generation around the catalytic triad.
3. Manual restrained construction from generic serine-hydrolase geometry only.

Do not use the PETase paper's pose, atom distances, RC terms, or TS snapshots.

## Geometry Filters

A pose can enter classical MD only if it passes all filters below after restrained relaxation:

| filter | accept criterion | reason |
| --- | --- | --- |
| nucleophile distance | Ser O-gamma to ester carbonyl C is less than or equal to 3.6 A | allows plausible acyl substitution approach |
| attack angle | O-gamma - C=O vector is within a broad near-Burgi-Dunitz window, 80-125 deg | screens impossible attack geometries without imposing a paper RC |
| oxyanion stabilization | carbonyl O has at least one backbone/side-chain H-bond donor within 3.2 A | required for generic serine hydrolase chemistry |
| His accessibility | catalytic His has one ring N within 3.4 A of Ser O-gamma or a candidate relay water | permits proton abstraction alternatives |
| leaving group path | leaving O has no severe steric clash and can reach His/water network within 4.0 A after local relaxation | permits product formation |
| substrate stability | substrate aromatic/ester fragment remains in the binding groove during short restrained relaxation | removes poses that are only docking artifacts |

These thresholds are deliberately broad. They are filters for chemical plausibility, not assumed reaction coordinates.

## Manifest Schemas

`gs_pose_manifest.tsv`:

```text
pose_id	template_pdb	substrate_model	generation_method	relaxed_structure_path	ser_og_to_c_A	attack_angle_deg	oxyanion_hbond_count	his_acceptor_distance_A	leaving_group_relay_distance_A	trp_rotamer_label	pass_fail	rejection_reason	next_step
```

`rejected_pose_manifest.tsv`:

```text
pose_id	template_pdb	substrate_model	generation_method	failed_filter	measured_value	threshold	rejection_reason
```

`ligand_model_manifest.tsv`:

```text
model_id	formal_charge	protonation_state	parameter_source	qm_atom_label_file	scissile_carbonyl_c	scissile_carbonyl_o	leaving_o	nucleophile_target_o	build_command	validation_status
```

## Grill Gates

Before Stage 2 starts, answer these questions with evidence:

1. If the primary template is `6EQE`, what changes when the same substrate construction is repeated on `5XJH` and `5YFE`?
2. Are any accepted poses dependent on a single crystal-water placement?
3. Does each accepted pose satisfy generic serine-hydrolase geometry without using paper-specific TS or RC information?
4. Are rejected poses recorded clearly enough that we can defend why they were not taken into QM/MM?
5. Is the chosen substrate model chemically close enough to PET to answer the mechanism question, but small enough for QM/MM sampling?

Stage 2 may begin only after at least 10 accepted acylation Michaelis candidates and at least 3 deacylation precursor candidates are recorded, or after a written justification explains why fewer candidates are chemically sufficient.
