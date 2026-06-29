# Project 01 Phase1 QMMM preparation notes

Evaluated universe: 15 PLACER conformers from five sequence-similarity bins (90/80/70/60/50), three conformers per selected rank-1 sequence.

Current strict launch gate: no structure is `READY_STRICT_PASS`. All 15 were evaluated with output after sidechain completion and fixed-pocket sidechain grafting; all fail only `ligand_heavy_rmsd` and `key_distance_delta`.

QM definition for future inputs: full `bu2` ligand, chain X, residue 1, 32 atoms, plus catalytic residue fragments A95 HIS, A126 SER, and A128 SER. Link atoms should be placed at CA-CB cuts for sidechain fragments unless a larger backbone-inclusive QM region is chosen.

MM definition: all remaining protein atoms in the full enzyme structure. The no-MM control uses the same QM atoms without the protein MM point-charge/mechanical environment.

Do not launch production QMMM from this manifest until the strict geometry gate is relaxed by project decision or new PLACER/conformer sampling produces strict-pass structures. The top5 exploratory rows are ranked only by closest reactive-distance deviation, not by strict pass.
