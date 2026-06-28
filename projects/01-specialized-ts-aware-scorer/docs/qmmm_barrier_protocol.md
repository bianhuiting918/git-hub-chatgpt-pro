# QM/MM and DFT barrier-label protocol — Project 01

Created: 2026-06-28

This document defines the metadata and boundary conditions for external QM/MM or DFT labels consumed by Project 01.

Project 01 does **not** run heavy QM/MM, DFT, PLACER, or TS-search jobs inside this repository. It consumes externally generated labels and stores only lightweight manifests, summaries, schemas, and synthetic examples.

## 1. Purpose

Project 01 trains a true-TS teacher model:

```text
protein + GS + true/refined TS + QM/MM or DFT label
→ ΔG‡ / ΔΔG‡ / ranking / mechanism components
```

Every computed label must be linked to a reproducible protocol. Barrier labels without protocol metadata should not be used for training except as toy or synthetic examples.

## 2. State definitions

Each quantum label should identify the reaction state:

```yaml
state_type: GS | TS | intermediate | product
reaction_template_id: serine_hydrolase_acylation_TS1
reaction_step: acylation
```

Recommended first-stage reaction template:

```text
serine_hydrolase_acylation_TS1
```

Core geometry:

```yaml
forming_bond:
  - Ser_Ogamma
  - substrate_carbonyl_C
breaking_or_weakening_bond:
  - carbonyl_C
  - leaving_group_O_or_N
proton_transfer:
  - Ser_Ogamma_H
  - His_N
oxyanion_stabilization:
  - substrate_carbonyl_O
  - oxyanion_hole_donors
```

## 3. True TS definition

A true/refined TS in Project 01 may be one of:

```text
QM/MM transition-state optimization
QM/MM reaction-coordinate scan maximum or constrained TS-like point
QM/MM free-energy barrier state
DFT cluster TS optimization
DFT cluster constrained TS-like state
validated TS template or TS analog geometry refined by a consistent protocol
```

A raw PLACER pose is not considered a true TS. PLACER can seed or screen conformers, but the label source must be quantum or protocol-validated.

## 4. Role of PLACER

Recommended Project 01 workflow:

```text
true/refined TS template + GS complex
→ PLACER-generated or PLACER-screened conformer ensemble
→ geometry filtering
→ representative conformer selection
→ restrained relaxation / QM/MM / DFT refinement
→ barrier label
```

Required PLACER metadata, if applicable:

```yaml
placer_run_id:
placer_model_version:
placer_input_structure:
placer_input_ligand_or_ts_template:
placer_score:
conformer_id:
cluster_id:
geometry_filter_passed:
```

## 5. QM/MM or DFT protocol metadata

Every external label should include:

```yaml
protocol_id: qmmm_protocol_v001
label_tier: QM_MM_scan | QM_MM_free_energy | DFT_cluster_optimized | DFT_cluster_single_point | xTB | proxy
software:
  program:
  version:
method:
  functional:
  basis_set:
  dispersion:
  solvent_model:
  embedding:
  charge:
  spin_multiplicity:
qm_region:
  atom_selection:
  catalytic_residues:
  substrate_atoms:
  waters:
  cofactors_or_metals:
  capping_scheme:
mm_environment:
  force_field:
  fixed_atoms:
  boundary_treatment:
protonation:
  histidine_states:
  catalytic_residue_states:
  ligand_charge_state:
optimization:
  optimized_or_single_point:
  constraints:
  reaction_coordinate:
  imaginary_frequency_checked:
  scan_points:
sampling:
  source_conformer_id:
  source_trajectory:
  representative_cluster:
```

## 6. Energy labels

Required or recommended fields:

```yaml
E_GS:
E_TS:
G_GS:
G_TS:
delta_E_TS_GS:
delta_G_dagger:
reference_variant_id:
delta_delta_E_vs_reference:
delta_delta_G_vs_reference:
units: kcal_per_mol
```

If only electronic energies are available, use `delta_E_TS_GS` and set `delta_G_dagger` to null.

## 7. Recommended label hierarchy

Use label tiers rather than mixing all labels as one truth source:

```text
Tier 0: geometry / electrostatic proxy only
Tier 1: xTB or semiempirical single-point
Tier 2: DFT cluster single-point
Tier 3: DFT cluster restrained optimization
Tier 4: QM/MM reaction-coordinate scan
Tier 5: QM/MM free-energy barrier
```

The first model can train on multiple tiers, but should preserve `label_tier` and `protocol_id` so later models can do calibration or multi-fidelity learning.

## 8. Preferred targets

Use these targets in order of robustness:

```text
1. ΔΔG‡ or ΔΔE‡ vs WT/reference within the same reaction template
2. pairwise ranking within the same enzyme/family/template
3. absolute ΔG‡ or ΔE_TS-GS under a single protocol
```

Do not directly mix absolute barriers across unrelated catalytic classes unless protocol and reaction-template effects are explicitly modeled.

## 9. Minimal external label package

A minimal external calculation should return:

```text
protocol.yaml
energies.json
optimized_GS.pdb or relaxed_GS.pdb
optimized_TS.pdb or refined_TS.pdb
geometry_features.json
charges_GS.csv, optional
charges_TS.csv, optional
delta_q.csv, optional
electrostatic_features.json, optional
```

Large trajectories and raw wavefunction files should stay outside this repository.
