# Project 01 Phase1 PLACER Pilot Status - 2026-06-30

## Evaluated universe

This pilot evaluates 5 representative rank-1 structures from the Phase1 Top50 manifest, one per sequence identity bin.

- Source manifest: `/Dell/Dell14/bianht/project01_similarity_phase1_20260630/placer_manifest_phase1_top50.tsv`
- Pilot root: `/Dell/Dell14/bianht/project01_similarity_phase1_20260630/placer_pilot_rank1_20260630_0148`
- PLACER input list: `/Dell/Dell14/bianht/project01_similarity_phase1_20260630/placer_pilot_rank1_20260630_0148/placer_input_files.txt`
- PLACER command log: `/Dell/Dell14/bianht/project01_similarity_phase1_20260630/placer_pilot_rank1_20260630_0148/logs/placer_pilot_rank1_n3_retry1.log`

## Pilot inputs

| bin | rank | record_index | identity | mutation_count | pilot_input_pdb |
|---:|---:|---:|---:|---:|---|
| 90 | 1 | 1 | 0.9125 | 14 | `/Dell/Dell14/bianht/project01_similarity_phase1_20260630/placer_pilot_rank1_20260630_0148/inputs/bin90_rank1_rec1.pdb` |
| 80 | 1 | 2 | 0.79375 | 33 | `/Dell/Dell14/bianht/project01_similarity_phase1_20260630/placer_pilot_rank1_20260630_0148/inputs/bin80_rank1_rec2.pdb` |
| 70 | 1 | 1 | 0.725 | 44 | `/Dell/Dell14/bianht/project01_similarity_phase1_20260630/placer_pilot_rank1_20260630_0148/inputs/bin70_rank1_rec1.pdb` |
| 60 | 1 | 1 | 0.6000 | 64 | `/Dell/Dell14/bianht/project01_similarity_phase1_20260630/placer_pilot_rank1_20260630_0148/inputs/bin60_rank1_rec1.pdb` |
| 50 | 1 | 1 | 0.5000 | 80 | `/Dell/Dell14/bianht/project01_similarity_phase1_20260630/placer_pilot_rank1_20260630_0148/inputs/bin50_rank1_rec1.pdb` |

## Current PLACER run

- Tool: `/data/bht/design_tools/src/PLACER/run_PLACER.py`
- Environment: `/data/bht/design_tools/envs/placer_env/bin/python`
- Samples: `-n 3` per input
- Mutation setting: `--mutate 128A:75I`
- Custom residue JSON: `/Dell/Dell14/bianht/project01_conformer_queue/placer_denovo_ser_hydrolase_n5_20260629T011253/inputs/75I.json`
- Retry PID: `1635100`
- Status at last check: running on CPU, no output CSV/PDB yet flushed to disk.

## Failed attempt recorded

The first launch failed before structure evaluation because the mutation argument was given as `A128:75I`. `run_PLACER.py` expects the token format used in its help string, e.g. `128A:75I`. This is recorded as a command-format failure, not as a PLACER structural failure.

Failure record: `/Dell/Dell14/bianht/project01_similarity_phase1_20260630/placer_pilot_rank1_20260630_0148/logs/failed_attempt_001.txt`

## Pass condition for this pilot

A bin is counted as evaluated only when PLACER writes both score CSV and multi-model PDB for that bin input. Missing output while the process is still running is `IN_PROGRESS`; command/setup failure is `FAILED_BEFORE_EVALUATION`; input/model failure after PLACER starts is recorded with the log error text.
