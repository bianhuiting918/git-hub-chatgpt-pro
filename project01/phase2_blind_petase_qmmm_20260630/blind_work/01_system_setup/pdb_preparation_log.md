# Blind PETase Stage 1 PDB Preparation Log

Date: 2026-06-30

## Boundary

This log records structure-only preparation steps for the blind PETase QM/MM mechanism workflow. No paper TS coordinates, RCs, selected CVs, umbrella windows, shooting trajectories, barrier values, or mechanistic conclusions were used.

## Completed

1. Queried RCSB metadata for candidate PETase structures.
2. Selected `6EQE` as the primary WT-like template, with `5XJH`, `5YFE`, and `6ILW` as secondary WT-like sensitivity templates.
3. Downloaded nine WT-like or backup PDB coordinate files from RCSB:
   - `6EQE`, `5XJH`, `5YFE`, `6ILW`, `6EQD`, `6EQF`, `6EQG`, `6EQH`, `6QGC`.
4. Ran a coordinate-only Ser-His-Asp geometric scan with these thresholds:
   - Ser OG to His ND1/NE2 <= 4.0 A;
   - His ND1/NE2 to Asp OD1/OD2 <= 4.0 A.
5. Ran PDB preparation audit for chains, missing residues, missing atoms, alternate conformers, non-water heterogens, crystallographic water, and geometric disulfide candidates.
6. Wrote explicit template-chain decisions and retained-water candidate tables.
7. Defined blind substrate model inputs: `PET_dimer_capped`, `BHET_like`, `MHET_like`, and `MHET_like_acyl_enzyme_precursor`.
8. Generated initial cleaned chain-A PDB files locally with retained catalytic-site waters, non-water heterogens removed, disulfide records written, and altlocs resolved by highest altloc-atom occupancy.
9. Verified the generated v2 PDB files have zero nonblank altloc indicators after cleaning.
10. Checked the local ligand-build environment and confirmed RDKit/Open Babel/Biopython are not available in the bundled Python runtime.
11. Scanned v2 cleaned structures for titratable residues near the catalytic triad and all histidines requiring tautomer assignment.
12. Wrote the primary-template protonation hypothesis manifest for `6EQE`.

## Triad Scan Result

The scan identifies a consistent catalytic triad geometry across the WT-like templates:

- `6EQE`: chain A `SER160-HIS237-ASP206`.
- `5XJH`: chain A `SER160-HIS237-ASP206`.
- `6ILW`: chain A `SER160-HIS237-ASP206`.
- `5YFE`: chain A `SER134-HIS211-ASP180`, likely a residue-numbering offset relative to the other templates.
- Multi-chain crystal forms report the same local geometry per chain.

This is a structure-derived active-site assignment, not a paper-derived mechanism result.

## PDB Preparation Audit Result

- Initial chain decision: use chain A for the primary blind setup because chain A contains the consistent Ser-His-Asp triad in every selected template.
- `6EQE` has a complete triad, 33 missing residues reported in the PDB header, no missing atoms reported, 26 alternate-conformer residues, two geometric Cys-Cys pairs, and two crystallographic water candidates within 4 A of catalytic hetero atoms.
- `5XJH` has a complete triad, 37 missing residues reported, no missing atoms or alternate-conformer residues, two geometric Cys-Cys pairs, and three retained-water candidates.
- `5YFE` uses shifted residue numbering, has four missing C-terminal His residues reported, glycerol and sulfate heterogens, two geometric Cys-Cys pairs, and no retained-water candidates within the 4 A triad cutoff.
- Backup multi-chain structures contain repeated chain copies and should remain sensitivity/backup templates rather than primary production templates.

## Initial Cleaned Structure Result

- `prepared_structure_manifest.tsv` records nine generated initial-cleaned PDB files and SHA256 checksums.
- These are not final simulation-ready coordinates: protonation is not assigned, missing residues are not rebuilt, and no local minimization has been run.
- These files are intended as inputs to the next preparation layer: protonation, hydrogen placement, residue completion where needed, and local relaxation.

## Protonation Setup Result

- Default pH-7 assumptions are explicit: Asp/Glu deprotonated, Lys/Arg protonated, Cys/Tyr neutral, His neutral tautomer to assign.
- For primary template `6EQE`, catalytic `ASP206` is deprotonated in the primary hypothesis but requires a neutral-Asp sensitivity branch if pKa or the His-Asp network supports it.
- Catalytic `HIS237` requires HID/HIE tautomer sensitivity before committing to QM/MM reaction scans; HIP is only a branch if pKa/local electrostatics supports cationic histidine.
- Remote histidines such as `HIS104` and `HIS293` require tautomer naming for topology generation but are not automatically reaction-mechanism sensitivity branches.
- `CYS203` and `CYS239` are near the triad but treated as disulfide-linked unless preparation software contradicts the geometric disulfide assignment.

## Ligand Model Definition Result

- `PET_dimer_capped` is the primary acylation substrate model because it preserves PET-like ester/aromatic repetition while remaining tractable.
- `BHET_like` is a small neutral control for docking and pose-generation sanity checks.
- `MHET_like` is the default product-side/deacylation reference fragment, with pH-7 carboxylate state explicit.
- `MHET_like_acyl_enzyme_precursor` is explicitly a protein-covalent Ser160 acyl-enzyme model, not a standalone free ligand.
- Stable atom-label files are required before topology conversion so the scissile ester atoms can be traced into QM/MM setup.
- Ligand 3D conformer generation is blocked in the current local bundled Python because RDKit/Open Babel are unavailable; see `ligand_build_environment_report.md`.

## Current Decisions

1. Primary production template remains `6EQE`, chain A.
2. Secondary sensitivity templates remain `5XJH`, `5YFE`, and `6ILW`, chain A.
3. Geometric disulfide candidates should be kept unless later preparation software contradicts the connectivity.
4. Water molecules listed in `retained_water_candidates.tsv` should be retained for the first local relaxation pass, then tested by sensitivity runs.
5. Alternate conformers in `6EQE` are resolved in the v2 initial-cleaned PDB by highest altloc-atom occupancy; this remains subject to active-site geometry inspection after ligand placement.
6. Ligand construction must preserve the named scissile ester atom labels from `ligand_model_manifest.tsv`.
7. Protonation production setup must run an external pKa/protonation tool and compare its output with `protonation_hypothesis_manifest.tsv` before topology generation.

## Not Yet Done

- Production repaired/protonated PDB/mmCIF generation.
- Missing-residue rebuilding where needed.
- Hydrogen placement and protonation-state assignment with PROPKA/H++/pdb2pqr/Amber reduce or equivalent.
- 3D substrate conformer generation and parameterization in an environment with RDKit/Open Babel or an equivalent chemistry stack.
- Michaelis-complex generation and relaxation.

## Evidence Files

- `structure_selection.tsv`
- `structure_download_manifest.csv`
- `ser_his_asp_triad_candidates.tsv`
- `pdb_preparation_audit.tsv`
- `template_chain_decisions.tsv`
- `retained_water_candidates.tsv`
- `prepared_structure_manifest.tsv`
- `altloc_resolution_decisions.tsv`
- `protonation_site_scan.tsv`
- `protonation_hypothesis_manifest.tsv`
- `protonation_setup_notes.md`
- `ligand_model_manifest.tsv`
- `ligand_model_definitions.md`
- `ligand_build_environment_report.md`
- `stage1_system_setup_protocol.md`
