# Project 01 Phase1 strict-pass PLACER resampling pilot

## Why this pilot was started

The previous PLACER pilot produced 15 evaluated conformers across 90/80/70/60/50 sequence-similarity bins, but strict-ready QMMM structures remained 0/15 after full ligand mapping, sidechain completion, and fixed-pocket sidechain grafting.

Best current candidate before this resampling round:

- sample: `bin80_rank1_rec2_model_002`
- max key-distance delta: 0.5598 A
- strict threshold: 0.35 A
- ligand heavy RMSD: 8.0213 A
- strict status: `BLOCKED_STRICT_GEOMETRY_FAIL`

## New resampling scope

Evaluated universe for this run is intentionally small:

- one input sequence/background: bin80 rank1 record2
- input PDB: `/Dell/Dell14/bianht/project01_similarity_phase1_20260630/bin80_t002_seed8102_redesign100/backbones/denovo_SER_hydrolase_full_input_2.pdb`
- two PLACER strategies, five samples each

Strategies:

1. `fixedlig_n5`: PLACER with `--fixed_ligand X-bu2-1 --fixed_ligand_noise 0`. This is a reactive-geometry sanity control, not a ligand-ensemble discovery run.
2. `bonded_ser128_c1_n5`: PLACER with `--bonds A-128-75I-OG:X-1-bu2-C1:1.533`. This tests whether PLACER's bond graph can pull the catalytic Ser128OG-to-ligand-C1 distance toward the reference geometry.

Both variants keep the previous `--mutate 128A:75I`, custom `75I.json`, and `--rerank prmsd` settings.

## Remote run state

- run root: `/Dell/Dell14/bianht/project01_similarity_phase1_20260630/placer_resample_strictgate_bin80_rank1_20260630_0312`
- driver script: `scripts/run_bin80_rank1_resample.sh`
- driver PID at launch: `2492441`
- run manifest: `resample_run_manifest.tsv`
- driver log: `logs/driver.nohup.log`

## Pass condition

A new output is not QMMM-ready just because PLACER runs. It must go through the same downstream chain:

1. PLACER output exists.
2. Split conformers.
3. Map back to full-length full-bu2 structure without minimization.
4. Complete sidechains.
5. Fixed-pocket sidechain graft.
6. Strict all-atom pocket geometry pass.

Only rows with `READY_STRICT_PASS` should be used for production QMMM. Otherwise they remain `BLOCKED_STRICT_GEOMETRY_FAIL` or tool-state failures with explicit reason.
