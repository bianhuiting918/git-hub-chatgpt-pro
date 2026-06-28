# External compute interface — Project 01

Created: 2026-06-28

This document defines how external PLACER, QM/MM, and DFT workflows should pass lightweight results into Project 01.

The repository should not contain raw trajectories, large PLACER ensembles, wavefunction files, checkpoint files, or large quantum-chemistry output directories.

## 1. Expected directory layout outside git

External compute results can live in a non-versioned workspace such as:

```text
/external/enzyme_barrier_runs/
  project01/
    PETase/
      S238F/
        acylation_TS1/
          conf_017/
            protocol.yaml
            manifest.json
            energies.json
            relaxed_GS.pdb
            refined_TS.pdb
            geometry_features.json
            electrostatic_features.json
            charges_GS.csv
            charges_TS.csv
            delta_q.csv
```

Only small manifest examples should be committed to this repository.

## 2. Manifest format

Each externally computed conformer/state pair should provide a lightweight manifest:

```yaml
sample_id: PETase_S238F_conf_017_acylation_TS1
enzyme_id: PETase
variant_id: S238F
catalytic_class: serine_hydrolase
reaction_template_id: serine_hydrolase_acylation_TS1
reaction_step: acylation
conformer_id: conf_017
conformer_source: PLACER_screened_true_TS
protocol_id: qmmm_protocol_v001
label_tier: QM_MM_scan

paths:
  protocol: protocol.yaml
  energies: energies.json
  GS_structure: relaxed_GS.pdb
  TS_structure: refined_TS.pdb
  geometry_features: geometry_features.json
  electrostatic_features: electrostatic_features.json
  charges_GS: charges_GS.csv
  charges_TS: charges_TS.csv
  delta_q: delta_q.csv
```

## 3. `energies.json`

Recommended fields:

```json
{
  "units": "kcal_per_mol",
  "E_GS": null,
  "E_TS": null,
  "G_GS": null,
  "G_TS": null,
  "delta_E_TS_GS": null,
  "delta_G_dagger": null,
  "reference_variant_id": "WT",
  "delta_delta_E_vs_reference": null,
  "delta_delta_G_vs_reference": null,
  "protocol_id": "qmmm_protocol_v001",
  "label_tier": "QM_MM_scan"
}
```

Null values are allowed when a label is unavailable.

## 4. `geometry_features.json`

Recommended fields:

```json
{
  "forming_bond_length": null,
  "breaking_bond_length": null,
  "nucleophilic_attack_angle": null,
  "proton_transfer_distance": null,
  "oxyanion_hole_distances": [],
  "catalytic_triad_geometry": {},
  "leaving_group_alignment": null,
  "geometry_filter_passed": true
}
```

## 5. Charge and electrostatic files

`charges_GS.csv` and `charges_TS.csv` should use a simple format:

```text
atom_id,atom_name,residue_id,residue_name,x,y,z,charge
```

`delta_q.csv` should use:

```text
atom_id,atom_name,residue_id,residue_name,delta_q
```

`electrostatic_features.json` may include:

```json
{
  "reaction_projected_field": null,
  "reaction_projected_potential": null,
  "probe_points": [],
  "dielectric": null,
  "method": "classical_point_charge_or_external"
}
```

## 6. Validation expectations

Codex should implement validators that check:

- required manifest fields exist.
- `protocol_id` is present.
- `label_tier` is present.
- units are defined.
- GS and TS structures are both available when a barrier label is present.
- `delta_delta` labels define a reference variant.

## 7. Export to Project 02

After training Project 01, export a teacher table containing:

```text
sample_id
enzyme_id
variant_id
catalytic_class
reaction_template_id
reaction_step
protein_embedding_path
GS_embedding_path
true_TS_embedding_path
project01_delta_G_pred
project01_delta_delta_G_pred
project01_field_score
project01_geometry_score
project01_ensemble_score
project01_uncertainty
qmmm_delta_G_dagger
qmmm_delta_delta_G_dagger_vs_reference
protocol_id
```

This export is the main interface consumed by Project 02.
