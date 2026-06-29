# Project 01 Phase1 gate diagnosis after strict-pass resampling

## Evaluated universe

This diagnosis combines all currently evaluated fixed-pocket-grafted strict geometry outputs:

- pilot rank1 representatives: 15 structures, 15 evaluated with output
- bonded Ser128OG-C1 resampling: 5 structures, 5 evaluated with output
- fixedlig target-res smoke: 1 structure, 1 evaluated with output
- fixedlig no-target run: 1 row, failed before evaluation due to PLACER crop-center error

Total rows: 22. Evaluated with output: 21. Failed before evaluation: 1.

## Gate decomposition

Strict gate was decomposed into four thresholds:

- catalytic heavy RMSD <= 0.5 A
- fixed-pocket heavy RMSD <= 0.75 A
- max key-distance delta <= 0.35 A
- full ligand heavy RMSD <= 0.5 A

`reaction_geometry_gate` is defined as catalytic heavy RMSD pass + fixed-pocket heavy RMSD pass + key-distance pass. It does not include full ligand RMSD.

`strict_all_gate` is reaction geometry pass + full ligand heavy RMSD pass.

## Result

- strict all-gate pass: 0/21 evaluated structures
- reaction-geometry pass: 1/21 evaluated structures
- full ligand RMSD pass: 0/21 evaluated structures

The only reaction-geometry-pass row is:

- universe: `resample_bonded_ser128_c1_n5`
- sample: `bin80_rank1_rec2_model_005`
- catalytic heavy RMSD: 0.1892 A
- fixed-pocket heavy RMSD: 0.1879 A
- max key-distance delta: 0.2455 A
- full ligand heavy RMSD: 8.5047 A
- current QMMM launch gate: `BLOCKED_STRICT_GEOMETRY_FAIL`

## Interpretation

The current PLACER bonded strategy can satisfy the local reaction-geometry criteria for at least one structure, but the full ligand heavy RMSD threshold blocks all candidates. This means the current strict all-gate is dominated by global/full-ligand pose preservation, not by catalytic-distance recovery.

## Recommendation

Do not launch production QMMM under the current strict all-gate. The next project decision should be explicit:

1. Keep full ligand RMSD as a hard gate. Then upstream input construction must be changed to preserve the full ligand pose more tightly before/through PLACER.
2. Reclassify full ligand RMSD as an ensemble-diversity descriptor rather than a launch blocker. Then `bin80_rank1_rec2_model_005` can enter a clearly labeled exploratory +/-MM comparison under a reaction-geometry-only gate.
3. Replace full ligand RMSD with a local ligand/contact RMSD around the reactive atoms and oxyanion/contact region, then rerank all evaluated candidates by that local metric.

No PDB, trajectory, model, or large log files were uploaded. Remote evidence is under `/Dell/Dell14/bianht/project01_similarity_phase1_20260630/gate_diagnosis_20260630_0343`.
