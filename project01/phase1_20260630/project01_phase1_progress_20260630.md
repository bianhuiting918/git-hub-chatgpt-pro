# Project 01 Phase 1 Progress - 2026-06-30

## Current status

Phase 1 sequence-similarity panel is complete for the serine hydrolase design target.

- Target bins completed: 90%, 80%, 70%, 60%, 50%
- Selected designs per bin: 10
- Total selected designs: 50
- Active/direct-contact pocket residues fixed: 17 residues
- PLACER-ready lightweight manifest: `placer_manifest_phase1_top50.tsv` and `placer_manifest_phase1_top50.json`

## Hard evidence files on CPU server

- Phase root: `/Dell/Dell14/bianht/project01_similarity_phase1_20260630`
- ProTrek smoke: `/Dell/Dell14/bianht/project01_protrek_smoke_latest.json`
- Similarity/bin spec: `/Dell/Dell14/bianht/project01_similarity_phase1_20260630/phase1_similarity_bins.json`
- Top50 PLACER manifest:
  - `/Dell/Dell14/bianht/project01_similarity_phase1_20260630/placer_manifest_phase1_top50.tsv`
  - `/Dell/Dell14/bianht/project01_similarity_phase1_20260630/placer_manifest_phase1_top50.json`

## Completed checks

1. ProTrek_650M sequence+structure smoke passed.
   - Sequence embedding shape: `[1, 1024]`
   - Structure embedding shape: `[1, 1024]`
   - Sequence-structure cosine on reference: `0.4478713870`

2. LigandMPNN bins were generated and filtered.
   - 90%: 10 selected, identity 0.89375-0.925
   - 80%: 10 selected, identity 0.78125-0.825
   - 70%: 10 selected, identity 0.69375-0.725
   - 60%: 10 selected, identity 0.60000
   - 50%: 10 selected, identity 0.50000

3. Geometry checks passed at current LigandMPNN stage.
   - LigandMPNN `backbones/` PDBs omit protein side chains.
   - Therefore the stage-valid check is catalytic/fixed backbone RMSD plus ligand heavy-atom RMSD.
   - Full side-chain key-distance checks are deferred to packed/full/PLACER structures.

4. ProTrek global/pocket/shell cosines were computed for all selected designs.
   - Pocket 5A sequence cosine is 1.0 for all selected records because pocket/direct-contact residues are fixed.
   - Global sequence cosine decreases across bins as intended.

## Files uploaded in this lightweight report bundle

- `project01_phase1_progress_20260630.md`
- `phase1_similarity_bins.json`
- `placer_manifest_phase1_top50.tsv`
- `placer_manifest_phase1_top50.json`
- `compute_protrek_layer_cosines.py`
- `filter_ligandmpnn_sequences.py`
- `pocket_backbone_ligand_filter.py`
- `pocket_geometry_filter.py`

No model weights, PDB batches, `.pt` stats, or large outputs are included.

## Next action

Use `placer_manifest_phase1_top50.tsv` to start small PLACER conformer generation, beginning with the 90% bin, then apply full side-chain/key-distance geometry filtering before QMMM preparation.
