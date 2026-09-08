# Multiscale material shape matching v1 — implementation plan
Approved scope: PET/PA6/PA66; all original dense exterior sites and deduplicated sites accessible to CURRENT 101-core panel. Never call this the exhaustive historical library.
Radii 10,15,20 A. Same cohort/coordinate frame at each radius. Input normals uncertain remain NOT_EVALUATED in original denominator.
## Definition
Actual matching criterion: symmetric area-weighted surface distance D95 (max of two directional 95th percentiles); sampled surface quadrature, also record worst sampled distance. Report required A tolerance for ceil(N*0.9/0.7/0.5), plus achieved rates at 1/2/3 A. Not enzyme collision or chemical energy. No manufactured curvature or artificial flattening.
Choose candidate surfaces at R20 from fixed finite library of real representatives and median SES fields; same chosen 3D field across smaller radii. Record searched library and empirical in-sample status, not global optimum. Chemical annotation post hoc with existing Crippen/HBA/HBD.
## Tasks / gates
- [ ] Verify metric against identical surface, 2 A planar shift, asymmetric extra patch, empty mesh, ceil coverage with missing sites.
- [ ] Recompute original exterior-water SES fields on padded >=43^3 grid for R20, CPU existing environment. Preserve full atom sources and XY periodicity.
- [ ] PET 8-site pilot: generate real connected patches at R10/15/20, check finite mesh, radius, topology, metric refinement, and dense-body-vs-tail scope. STOP publication if dense external face requirement is not met.
- [ ] On pilot success run full fixed cohorts and retain 90/70/50 achieved tolerances, same membership across radii; build chemistry views.
- [ ] Verify data/figures, publish only verified results and exact scope on existing GitHub branch; append remote RUN_LOG/RUNBOOK.
Inputs: panel_surface_scan_v1/sites.json, coverage.npz, polymer_exposure_full_v2/*/sites.json, inputs_v2/*_parent_atoms.npz, existing chemical_contact_surface_v3.py.
Outputs: multiscale_material_match_v1.py, test_multiscale_material_match_v1.py, new multiscale_material_match_pilot_v1; never overwrite old outputs.
Minimal resources: sequential CPU, original m.ChemEnv dependencies, no GPU, no new MD or donor sampling, no local project writes.
