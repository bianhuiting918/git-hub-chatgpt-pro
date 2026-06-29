# Project 01 Phase1 PLACER Pilot Status - 2026-06-30 02:01 CST

## Current state

PLACER pilot is running on the CPU server. The initial all-input process produced the bin90 output, then failed before evaluating bin80 because the custom residue `75I` was already registered inside the same PLACER process. This is recorded as a tool-state issue, not a structural failure for bins 80/70/60/50.

A serial per-input recovery run is now active for bins 80, 70, 60, and 50.

## Hard evidence

- Pilot root: `/Dell/Dell14/bianht/project01_similarity_phase1_20260630/placer_pilot_rank1_20260630_0148`
- Failed attempt 001: `/Dell/Dell14/bianht/project01_similarity_phase1_20260630/placer_pilot_rank1_20260630_0148/logs/failed_attempt_001.txt`
- Failed attempt 002: `/Dell/Dell14/bianht/project01_similarity_phase1_20260630/placer_pilot_rank1_20260630_0148/logs/failed_attempt_002.txt`
- Serial recovery script: `/Dell/Dell14/bianht/project01_similarity_phase1_20260630/placer_pilot_rank1_20260630_0148/run_remaining_bins_serial.sh`
- Serial recovery PID: `1695461`
- Current active child: `run_PLACER.py` PID `1695471`, input `bin80_rank1_rec2.pdb`

## Bin-level status

| bin | input | status | evidence |
|---:|---|---|---|
| 90 | `bin90_rank1_rec1.pdb` | `EVALUATED_WITH_OUTPUT` | CSV and multi-model PDB written |
| 80 | `bin80_rank1_rec2.pdb` | `IN_PROGRESS` | child `run_PLACER.py` PID `1695471` |
| 70 | `bin70_rank1_rec1.pdb` | `NOT_EVALUATED_YET` | queued in serial recovery script |
| 60 | `bin60_rank1_rec1.pdb` | `NOT_EVALUATED_YET` | queued in serial recovery script |
| 50 | `bin50_rank1_rec1.pdb` | `NOT_EVALUATED_YET` | queued in serial recovery script |

## Bin90 output

- CSV: `/Dell/Dell14/bianht/project01_similarity_phase1_20260630/placer_pilot_rank1_20260630_0148/outputs/bin90_rank1_rec1_phase1pilot_n3.csv`
- Multi-model PDB: `/Dell/Dell14/bianht/project01_similarity_phase1_20260630/placer_pilot_rank1_20260630_0148/outputs/bin90_rank1_rec1_phase1pilot_n3_model.pdb`
- Split manifest: `/Dell/Dell14/bianht/project01_similarity_phase1_20260630/placer_pilot_rank1_20260630_0148/placer_pilot_split_manifest.tsv`
- Split model count so far: 3

## Audit note

The evaluated universe for this pilot is 5 rank-1 representatives, one per similarity bin. At this checkpoint only bin90 has completed PLACER output. The remaining bins are not failed; they are still in progress or queued. The bin90 three-model order is the PLACER `--rerank prmsd` output order.
