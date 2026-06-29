# Project 01 Phase1 strict-pass resampling outcome

## Evaluated universe

This resampling round intentionally used only the current closest background:

- sequence/background: bin80 rank1 record2
- original input PDB: `/Dell/Dell14/bianht/project01_similarity_phase1_20260630/bin80_t002_seed8102_redesign100/backbones/denovo_SER_hydrolase_full_input_2.pdb`
- run root: `/Dell/Dell14/bianht/project01_similarity_phase1_20260630/placer_resample_strictgate_bin80_rank1_20260630_0312`

## Variants

1. `fixedlig_n5`
   - status: `FAILED_BEFORE_EVALUATION`
   - reason: PLACER raised `IndexError: list index out of range` because the only ligand was fixed and no predicted ligand index remained for automatic crop/corruption center selection.
   - action: do not repeat this exact command.

2. `bonded_ser128_c1_n5`
   - status: PLACER output generated, 5/5 split, mapped to full-bu2, sidechain completed, fixed-pocket sidechain grafted, strict geometry evaluated.
   - result: 0/5 strict pass.
   - best row: `bin80_rank1_rec2_model_005`, key-distance max delta 0.2455 A, but ligand heavy RMSD 8.5047 A, so still `BLOCKED_STRICT_GEOMETRY_FAIL`.

3. `fixedlig_targetres_smoke_n1`
   - command idea: fixed ligand, zero ligand noise, explicit target residue A128.
   - status: PLACER output generated and fully evaluated.
   - result: 0/1 strict pass.
   - reason: after full-length mapping and fixed-pocket graft, ligand heavy RMSD 2.0484 A and key-distance max delta 1.8271 A.

## Current gate result

- total evaluated-with-output in this resampling round: 6
- failed-before-evaluation: 1
- strict pass: 0
- QMMM-ready: 0

The bonded strategy improved the reactive-distance criterion for one sample: `bonded_ser128_c1_n5/bin80_rank1_rec2_model_005` passes the current key-distance threshold but fails the ligand RMSD threshold. This suggests the strict ligand-heavy RMSD gate is now the main blocker after reactive-distance targeting.

## Recommendation

Do not start production QMMM yet under the current strict gate. The next decision should be one of:

1. Treat ligand-heavy RMSD as an ensemble-diversity metric rather than a hard launch gate, and define a reaction-geometry-only exploratory gate for +/-MM comparison.
2. Keep the current strict gate and change the upstream input construction so the full ligand reference pose is preserved more tightly before PLACER sampling.
3. Add a stronger post-PLACER selection criterion focused on catalytic distances and local ligand atoms instead of full-ligand RMSD.

All current structures remain on the CPU server only; no PDB, model, trajectory, or large log files were uploaded.
