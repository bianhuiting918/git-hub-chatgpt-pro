# Exploratory QMMM region summary

This package is prepared under `reaction_geometry_only_exploratory`, not under the strict all-gate.

## Candidate

- candidate: `bin80_rank1_rec2_model_005`
- universe: `resample_bonded_ser128_c1_n5`
- remote structure: `/Dell/Dell14/bianht/project01_similarity_phase1_20260630/placer_resample_strictgate_bin80_rank1_20260630_0312/fixed_pocket_sidechain_grafted_no_min/bin80_rank1_rec2_model_005/full_length_full_bu2_fixed_pocket_sidechains_no_min.pdb`
- strict all-gate: `FAIL`
- reaction-geometry gate: `PASS`
- QMMM label: `EXPLORATORY_NOT_STRICT_PASS`

## Gate Metrics

- catalytic heavy RMSD: 0.1892 A
- fixed-pocket heavy RMSD: 0.1879 A
- max key-distance delta: 0.2455 A
- full ligand heavy RMSD: 8.5047 A
- strict fail reason: `ligand_heavy_rmsd`
- missing atoms: none in catalytic, fixed, ligand, or key-distance sets

## QM/MM Definition

QM atom count: 42

- full `bu2` ligand: 32 atoms
- A95 HIS sidechain/imidazole fragment: 6 atoms
- A126 SER sidechain fragment: 2 atoms
- A128 SER sidechain fragment: 2 atoms

MM atom count: 1201

- plus-MM run: full enzyme environment outside the QM atoms
- no-MM control: same QM atoms only, no protein MM environment

Link atom hints for later engine-specific input:

- A95: CA-CB link atom for HIS95 sidechain fragment
- A126: CA-CB link atom for SER126 sidechain fragment
- A128: CA-CB link atom for SER128 sidechain fragment

Remote generated files:

- selection TSV: `/Dell/Dell14/bianht/project01_similarity_phase1_20260630/exploratory_qmmm_reaction_geometry_20260630_0353/selection/qm_atom_selection.tsv`
- index NDX: `/Dell/Dell14/bianht/project01_similarity_phase1_20260630/exploratory_qmmm_reaction_geometry_20260630_0353/selection/qmmm_index.ndx`
- QM XYZ without link atoms: `/Dell/Dell14/bianht/project01_similarity_phase1_20260630/exploratory_qmmm_reaction_geometry_20260630_0353/selection/qm_region.xyz`
- QM PDB without link atoms: `/Dell/Dell14/bianht/project01_similarity_phase1_20260630/exploratory_qmmm_reaction_geometry_20260630_0353/selection/qm_region_no_link_atoms.pdb`
- machine-readable summary: `/Dell/Dell14/bianht/project01_similarity_phase1_20260630/exploratory_qmmm_reaction_geometry_20260630_0353/selection/qmmm_region_summary.json`

No production QMMM has been launched from this package.
