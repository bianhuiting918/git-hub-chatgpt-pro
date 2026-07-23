# Codex Web task: ICCG with literature LG1 PET dimer

Work only inside this directory on branch `codex/iccg-lg1-dimer-viewer`.
Do not use or recreate the 16-atom HEMT/856 monomer visualization.

## Goal

Generate a verified display PDB containing repaired active ICCG as a protein cartoon,
the complete literature Step1 ground-state LG1 PET dimer as sticks, and the
corrected provisional LG1 5 A fixed-site mask in red.

## Supplied inputs

- `inputs/ICCG_active_chainA_heavy.pdb`: 258-residue repaired 6THT chain A,
  heavy atoms only; residue 165 is SER and contains OG.
- `inputs/IsPETase_WT_LG1_reference.pdb`: audited literature reference complex;
  exactly 54 LG1 atoms, including 32 heavy atoms and two aromatic carbon rings.
- `inputs/6THT_secondary_structure_records.pdb`: 21 HELIX/SHEET records.
- `inputs/ICCG_legacy_mask_components.tsv`: mapped ICCG flags from the previous
  audit. Reuse only `high_conservation`, `native_cysteine`, and catalytic
  identities. Ignore the old `hemt_5A` and old `fixed` columns.
- `inputs/phl7_wt_final_fixed_positions.tsv`: provenance-only conservation source.

## Scientific placement contract

1. Parse protein residues from the literature reference and ICCG.
2. Build a global sequence alignment and residue correspondence.
3. In the literature reference, select protein residues with any heavy atom
   within 5.0 A of an LG1 heavy atom. Union those residues with source anchors
   159, 160, 206, 237, and 238.
4. Fit mapped N/CA/C backbone atoms from the literature protein onto ICCG with a
   proper Kabsch transform.
5. Apply exactly that one rigid transform to all 54 LG1 atoms. Do not optimize,
   minimize, or change any protein or ligand internal coordinate.
6. Infer intra-LG1 bonds with covalent radii plus 0.45 A tolerance. Verify the
   carbon bond graph contains exactly two distinct six-membered rings.
7. Recompute ICCG residues within 5.0 A of transformed LG1 heavy atoms.
8. Viewer fixed set:
   `catalytic {165,210,242} | new_LG1_5A_pocket | high_conservation | native_cysteines`.
   Treat this as a provisional visualization mask if the pose fails geometry.

## Geometry gate

Calculate and record:

- minimum nonbonded protein-LG1 heavy-atom distance;
- maximum van der Waals overlap;
- worst protein/ligand atom pair;
- alignment anchor residue and backbone-atom counts;
- backbone fit RMSD.

Use a maximum overlap threshold of 0.8 A. If exceeded:

- still generate the requested viewer;
- mark it `PROVISIONAL_DISPLAY_ONLY_LG1_CLASH`;
- do not call the pose computation-ready;
- do not claim the provisional ligand-contact mask replaces any production
  sequence-generation mask.

## Required outputs

Create only:

- `build_iccg_lg1_viewer.py`
- `test_iccg_lg1_viewer.py`
- `ICCG_active_cartoon_fixed_sites_with_LG1.pdb`
- `ICCG_active_fixed_sites.tsv`
- `ICCG_active_fixed_sites_audit.json`
- `RUNBOOK.md`
- `RUN_HISTORY.tsv`

PDB requirements:

- include all 21 supplied HELIX/SHEET records;
- protein uses ATOM records, chain A;
- LG1 uses HETATM records, residue name LG1, chain L, residue 301;
- all 54 ligand atoms are present;
- explicit symmetric intra-ligand CONECT records are present;
- protein fixed residues have B=0.00, variable residues B=99.00, LG1 B=50.00;
- add a clear display-only REMARK if the geometry gate fails.

## TDD and verification

Write the test first and run it once to demonstrate the expected RED state.
Then implement the smallest builder that satisfies the contract.

The final test must verify:

- 258 protein residues;
- Ser165 OG;
- 54 total LG1 atoms and 32 LG1 heavy atoms;
- two six-carbon rings;
- 21 HELIX/SHEET records;
- non-empty symmetric ligand connectivity;
- fixed mask equals the exact required union;
- audit counts match PDB/TSV;
- geometry status and worst-contact fields are non-empty;
- no HEMT or residue 856 appears in the primary viewer.

Run:

`python test_iccg_lg1_viewer.py`

Record commands, timestamps, inputs, outputs, exit status, counts, and geometry
summary in RUN_HISTORY.tsv. Commit all completed files to the current branch.
Do not open or merge a PR automatically; report the final commit and scientific
geometry status.
