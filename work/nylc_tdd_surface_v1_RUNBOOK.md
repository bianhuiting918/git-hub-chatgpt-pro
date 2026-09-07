# NylC TDD / PA6 and PA66 surface screening v1

## Approved objective
Use the user-selected v10 ge90max uncapped PA6 and PA66 structures. Analyze top/bottom dense-body surfaces, intact amide targets and connected 15 A patches. Keep NylC 3AXG A Thr267/Asp306/Asp308 rigid; sample NH fragments and jointly feasible donor pairs. Rank fixed whole motifs by distinct site coverage and existing-point 1 A neighborhood density conditional on identical TDD pose. Produce aligned PyMOL core/surface views. Geometry is not activity validation.

## Input stage, 2026-09-08
Remote root: /data/bht2/polymer_material_reference_20260827/simulation_slabs/rectangular_400chain_v1/nylc_tdd_surface_v1
Local source root: D:/Codex/polymer_material_reference_20260827/simulation_slabs/rectangular_400chain_v1/pymol/axis_density_sixface_anchor_filtered_v10
Source file in each material subdirectory: MATERIAL_axis_density_ge90max_sixface_all_external_amide_trim_bound_terminal_trim_anchor_filtered_uncapped.pdb
Remote copies: inputs/PA6_v10_ge90_uncapped.pdb and inputs/PA66_v10_ge90_uncapped.pdb.
PA6 SHA256: 1e93b32c62f442c34531ca65767cffbeffd26375b9c0f786db8a48e4b1a630bf
PA66 SHA256: 125963a9af3f6d1cb3b62e567cb6247bf7adde330f6489933bbd2a5c12087fe4
Both material anchor_filtered_short_molecules_audit.json source records copied alongside PDBs.

## Reproduce input staging
Use native PowerShell Get-FileHash -LiteralPath on the two source files. Verify against the hashes above. Use scp with literal absolute source arguments into the remote inputs directory. Never replace a differing existing file; create a new version instead. On the remote server, from the remote root, run `sha256sum inputs/*.pdb` and require exact matches above.

## Template
/data/bht2/new_seed_structure_sequence_folddisco_20260825_v3/structure_collection_20260826_v1/structures/folddisco/PDB/PDB_3AXG_A.pdb
Coordinate audit found all 23 TDD heavy atoms; chain A residue 266 absent. Mature-terminal state and protonation require verification before sampling. Do not assign Thr267 as an ordinary internal peptide nitrogen.

## Method gates
1. Use CONECT plus original untrimmed parent mapping to mark missing cut bonds. Do not treat artificial exposed cuts as native surface accessibility.
2. ge90max is axis-density cropping, not the later 90% inward-support criterion. Keep denominators separate.
3. Coarse density boundary is a shape summary; use atomic geometry for clash checks. Color chemistry as a proxy, not a calibrated hydrophobic potential.
4. PET reaction thresholds are not automatically validated for NylC amide chemistry. Verify mechanism before formal scan.
5. Preserve fixed TDD plus both donor relative coordinates in transfer rankings. Existing-neighbor density must not pool unrelated TDD poses.
6. Run diagnostic subset before full screening, save seeds and failures. Keep native donor arrangements distinct from independently sampled NH fragments.

## Environment and state
CPU Python: /data/bht2/polymer_material_reference_20260827/simulation_slabs/unbiased_100chain_v1/env/bin/python
No GPU/CUDA needed for input audit. No existing jobs altered.
Input transfer: exit 0; local and remote hashes match. Scientific scan: NOT_STARTED. All surface/core results: NOT_EVALUATED.
