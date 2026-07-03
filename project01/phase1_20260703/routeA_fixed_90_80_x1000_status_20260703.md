# Project 01 RouteA fixed-scaffold 90/80 x1000 status

Updated: 2026-07-03 12:50 Asia/Shanghai

## Scope

This run is the fixed/reference-scaffold route only. It does not use the Baker
new-backbone `bn1` route, sample21/refine route, PLACER/QMMM exploratory route,
or any old strict-ligand RMSD gate outputs.

Remote root:

`/data/bht/project01_phase1_reset_gpu`

Run id:

`routeA_90_80_x1000_20260703_123422`

## Denominator audit

LigandMPNN generation requested:

- 90% identity bin: 1000 requested, 1001 raw records written.
- 80% identity bin: 1000 requested, 1001 raw records written.

Sequence-level filter output:

- 90% bin: 73 sequence-filter pass.
- 80% bin: 227 sequence-filter pass.
- Total entering ESMFold pocket/structure evaluation: 300.

Rows excluded at the sequence-filter stage are not pocket-geometry FAIL and are
not PLACER FAIL. They are excluded before structure evaluation because they did
not satisfy the fixed-residue, target-identity/mutation-count, duplicate, or
related sequence-level filters.

## Current live step

ESMFold structure prediction is running on the 300 sequence-filter-pass rows.

Evidence files:

- Input manifest:
  `/data/bht/project01_phase1_reset_gpu/routeA_fixed_90_80_x1000/routeA_90_80_x1000_20260703_123422/manifests/routeA_90_80_x1000_sequence_filter_pass_manifest.tsv`
- ESMFold status TSV:
  `/data/bht/project01_phase1_reset_gpu/routeA_fixed_90_80_x1000/routeA_90_80_x1000_20260703_123422/manifests/routeA_90_80_x1000_esmfold_status.tsv`
- ESMFold log:
  `/data/bht/project01_phase1_reset_gpu/logs/routeA_90_80_x1000_20260703_123422.esmfold.log`

Snapshot at 2026-07-03 12:49 Asia/Shanghai:

- ESMFold PID: 1317361.
- Status TSV lines: 22 including header.
- Completed predictions: 21/300 OK.
- No ESMFold failures observed in the snapshot.

## Pass criteria for next steps

1. ESMFold complete means 300 evaluated rows are represented in the status TSV
   and the summary JSON is written.
2. Pocket/entrance geometry gate will be run with the existing fixed-scaffold
   Project 01 post-sequence gate scripts.
3. Only pocket-gate PASS rows will be converted into holo PLACER inputs.
4. PLACER screening will be reported with explicit denominators:
   sequence-filter pass, ESMFold OK/FAIL/NOT_EVALUATED, pocket-gate PASS/FAIL,
   and PLACER conformer crop-gate PASS/FAIL/NOT_EVALUATED.

Large PDB, model, trajectory, and raw log files are intentionally not uploaded
to GitHub.
