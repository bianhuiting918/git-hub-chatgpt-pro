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

## 2026-07-03 14:50 Asia/Shanghai Update

ESMFold final status:

- Evaluated universe: 300 sequence-filter-pass rows.
- OK: 300.
- FAIL: 0.
- NOT_EVALUATED: 0.
- Summary JSON: `/data/bht/project01_phase1_reset_gpu/routeA_fixed_90_80_x1000/routeA_90_80_x1000_20260703_123422/manifests/routeA_90_80_x1000_esmfold_summary.json`

Post-sequence fixed-pocket/entrance gate:

- Evaluated universe: 300 ESMFold OK rows.
- PASS: 15.
- FAIL: 285.
- NOT_EVALUATED: 0.
- PASS by bin: 90% = 15, 80% = 0.
- Gate manifest: `/data/bht/project01_phase1_reset_gpu/routeA_fixed_90_80_x1000/routeA_90_80_x1000_20260703_123422/manifests/routeA_90_80_x1000_postseq_gate.tsv`

PLACER n50:

- Holo PLACER inputs were prepared only from the 15 postseq PASS rows.
- The 285 postseq gate FAIL rows are not PLACER FAIL.
- PLACER ifile: `/data/bht/project01_phase1_reset_gpu/routeA_fixed_90_80_x1000/routeA_90_80_x1000_20260703_123422/manifests/routeA_90_80_x1000_placer.ifile`
- PLACER queue manifest: `/data/bht/project01_phase1_reset_gpu/routeA_fixed_90_80_x1000/routeA_90_80_x1000_20260703_123422/manifests/routeA_90_80_x1000_placer_queue.tsv`
- PLACER serial driver PID: 1415493.
- PLACER status TSV: `/data/bht/project01_phase1_reset_gpu/routeA_fixed_90_80_x1000/routeA_90_80_x1000_20260703_123422/manifests/routeA_90_80_x1000_placer_n50_status.tsv`
- Expected maximum PLACER conformer denominator: 15 inputs x 50 = 750.

## 2026-07-03 17:55 Asia/Shanghai Final PLACER Crop-Gate Update

PLACER n50 driver final status:

- PLACER holo inputs: 15.
- Expected conformer denominator: 15 inputs x 50 = 750 conformers.
- PLACER input status: OK = 15, FAIL = 0, NOT_EVALUATED = 0.
- PLACER output files: 15 CSV score files and 15 model PDB files on the GPU host.
- PLACER status TSV: `/data/bht/project01_phase1_reset_gpu/routeA_fixed_90_80_x1000/routeA_90_80_x1000_20260703_123422/manifests/routeA_90_80_x1000_placer_n50_status.tsv`

Crop ligand geometry screen:

- Evaluated conformers: 750/750.
- Conformer strict-pass: 0/750.
- Conformer status: `FAIL_BOTH_PROTEIN_AND_LIGAND_GATES` = 750.
- Evaluated sequences: 15/15.
- Sequence strict-pass: 0/15.
- Sequence status: `DISCARD_SEQUENCE_REGENERATE_FROM_UPSTREAM` = 15.
- Summary JSON: `/data/bht/project01_phase1_reset_gpu/routeA_fixed_90_80_x1000/routeA_90_80_x1000_20260703_123422/placer_crop_gate/two_layer_placer_crop_gate_summary.json`
- Sequence summary TSV: `/data/bht/project01_phase1_reset_gpu/routeA_fixed_90_80_x1000/routeA_90_80_x1000_20260703_123422/placer_crop_gate/two_layer_placer_crop_gate_sequence_summary.tsv`
- Conformer TSV: `/data/bht/project01_phase1_reset_gpu/routeA_fixed_90_80_x1000/routeA_90_80_x1000_20260703_123422/placer_crop_gate/two_layer_placer_crop_gate_routeA_90_80_x1000.tsv`

Important denominator notes:

- The 285 postseq gate FAIL rows are not PLACER FAIL.
- The 80% bin has no postseq PASS rows in this run, so it has no PLACER-evaluated conformers.
- The PLACER crop-gate result applies to the 15 postseq PASS rows from the 90% bin only.

Applied crop thresholds:

- global backbone RMSD <= 2.5 A.
- fixed backbone RMSD <= 1.0 A.
- catalytic heavy RMSD <= 0.75 A.
- protein key-distance max absolute delta <= 0.75 A.
- ligand heavy RMSD <= 0.75 A.
- Ser128-OG to bu2-C1 delta <= 0.5 A.
- sequence-level pass requires at least 10 crop-pass conformers out of 50.

Interpretation:

- This run did not produce any sequence that should proceed downstream from PLACER.
- The correct next action is upstream regeneration or altered constraints, not QMMM.
