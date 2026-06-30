# Project 01 Phase1 Reset Status - 2026-06-30

## Reset target

Run the sequence-similarity series through structure-aware gates before PLACER and QMMM. Fixed-template zero RMSD is not sufficient evidence that mutant folded structures preserve the active pocket.

## Current audited denominator

- 50 selected sequences: 90/80/70/60/50 identity bins x 10 sequences per bin.
- Post-sequence structure gate: NOT_EVALUATED for all 50 because predicted PDB files are not yet available.
- PLACER n=50 queue: 0 entrance-pass samples, so no PLACER launch is allowed yet.

## Implemented gate

New scripts:

- `projects/01-specialized-ts-aware-scorer/scripts/phase1_reset/postseq_entrance_gate.py`
- `projects/01-specialized-ts-aware-scorer/scripts/phase1_reset/run_postseq_gate_manifest.py`

Gate behavior:

- Align candidate to reference by Kabsch using all common backbone atoms.
- Evaluate global backbone RMSD, fixed-pocket backbone RMSD, catalytic heavy-atom RMSD, and protein-only key-distance deltas.
- Protein-only predictors may omit ligand at this entrance gate.
- Ligand/key-distance gates become mandatory at PLACER crop and completion stages.
- Missing predicted PDB is NOT_EVALUATED_MISSING_PDB, not FAIL.

## Verification

- CPU reference-self smoke with `/Dell/Dell14/bianht/miniconda3/envs/colabfold/bin/python` returned PASS.
- GPU-local reference-self smoke under `/data/bht/project01_phase1_reset_gpu` returned PASS.
- Current `postseq_entrance_gate_manifest.tsv` status counts: NOT_EVALUATED = 50.
- Current `placer_n50_entrance_pass_queue.tsv` has header only; `entrance_pass_count = 0`.

## GPU-local immediate screening

Copied to GPU host `dell-PowerEdge-R760` under `/data/bht/project01_phase1_reset_gpu/`:

- scripts
- 50-sequence FASTA and manifest
- reference PDB
- fixed active/direct-contact residues
- GPU-local target manifest with predicted PDB paths rewritten to GPU-local paths

GPU policy for structure prediction remains: single process, `CUDA_VISIBLE_DEVICES=0`, PyTorch per-process memory fraction 0.50, and ESMFold chunk size 64 after the ESMFold/ColabFold environment smoke passes.

## Current blocker to prediction

GPU A100 is reachable, but ESMFold smoke is not yet passing because `protrek_gpu` lacks a complete `openfold` package. The bundled LigandMPNN `openfold` is partial and lacks `openfold.model`. PyPI `openfold==0.0.1` failed in Python 3.12 metadata generation. A clean Python 3.9/3.10 ESMFold or ColabFold environment is the next setup step.
