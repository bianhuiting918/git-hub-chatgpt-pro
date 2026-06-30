# PLACER capped10 n50 bonded run status - 2026-06-30

## Evaluated Universe

Input universe for this run:

```text
50 postseq entrance-pass capped10 sequences
10 each from 90/80/70/60/50 similarity bins
```

Remote input manifest:

```text
/data/bht/project01_phase1_reset_gpu/manifests/placer_n50_holo_inputs_capped10.gpu.tsv
```

Remote holo input directory:

```text
/data/bht/project01_phase1_reset_gpu/placer_inputs_holo_capped10
```

Generation result before PLACER:

```text
50/50 holo inputs written
32 bu2 HETATM appended per input
missing predicted PDB errors: 0
```

## Smoke Test

Smoke sample:

```text
bin60_r2d_r2b_from80_rank9_rec12
```

Command policy:

```text
--mutate 128A:75I
--predict_ligand X-bu2-1
--bonds A-128-75I-OG:X-1-bu2-C1:1.533
--residue_json /data/bht/design_tools/src/PLACER/examples/ligands/75I.json
--rerank prmsd
```

Smoke result:

```text
PASS: PLACER wrote one model PDB and one CSV under
/data/bht/project01_phase1_reset_gpu/placer_runs/smoke_capped10_bin60_n1_bonded
```

This only validates command/input compatibility. It is not a structural pass label.

## n50 Batch

Run label:

```text
capped10_n50_bonded
```

Requested sampling:

```text
50 inputs x 50 conformers = 2500 requested PLACER conformers
```

Remote output directory:

```text
/data/bht/project01_phase1_reset_gpu/placer_runs/capped10_n50_bonded
```

Observed active process at launch check:

```text
run_PLACER.py --ifile /data/bht/project01_phase1_reset_gpu/manifests/placer_n50_holo_inputs_capped10.ifile --odir /data/bht/project01_phase1_reset_gpu/placer_runs/capped10_n50_bonded --suffix bonded_bu2_n50 --nsamples 50 --mutate 128A:75I --predict_ligand X-bu2-1 --bonds A-128-75I-OG:X-1-bu2-C1:1.533 --residue_json /data/bht/design_tools/src/PLACER/examples/ligands/75I.json --rerank prmsd --cautious
```

Status at 2026-06-30 14:21 CST:

```text
SINGLE_PROCESS_RUN_FAILED_AFTER_FIRST_INPUT
CSV outputs present: 1
model PDB outputs present: 1
```

Failure mode:

```text
AssertionError: Residue 75I already in database, please choose a different name3.
```

Interpretation:

```text
This is a PLACER process-state/custom-residue registration issue after the first input completed.
It is not a structural failure for the remaining 49 inputs.
```

Recovery action:

```text
Started per-input driver:
/data/bht/project01_phase1_reset_gpu/scripts/run_placer_capped10_n50_bonded_per_input.sh

Driver behavior:
- skip existing CSV/PDB outputs;
- run one PLACER Python process per input;
- keep the same bonded 128A:75I command policy;
- write per-input logs under /data/bht/project01_phase1_reset_gpu/placer_runs/capped10_n50_bonded/per_input_logs.
```

Status at 2026-06-30 14:23 CST:

```text
PER_INPUT_DRIVER_RUNNING
driver PID: 3692495
first completed input skipped
second input started: bin90_rank2_rec3.holo_bu2_refpose
```

Status at 2026-06-30 14:24 CST:

```text
PER_INPUT_DRIVER_VERIFIED
CSV outputs present: 2
second input completed: bin90_rank2_rec3.holo_bu2_refpose
third input started: bin90_rank3_rec4.holo_bu2_refpose
```

This verifies the per-input recovery path across at least one input boundary after the single-process `75I` registration failure.

Status at 2026-06-30 14:31 CST:

```text
PER_INPUT_DRIVER_RUNNING
CSV outputs present: 6
model PDB outputs present: 6
current input started: 7 / 50
```

Partial crop gate smoke:

```text
remote out dir: /data/bht/project01_phase1_reset_gpu/placer_crop_gate_capped10_partial
evaluated universe: 6 completed sequences, 300 conformers
CROP_STRICT_PASS: 0 / 300
```

This partial result validates the crop gate interface only. It is not the final 50-sequence crop gate result. Missing PLACER outputs remain `NOT_EVALUATED`.

Final crop gate script prepared:

```text
project01/phase1_reset_20260630/scripts/run_crop_gate_capped10_final.sh
```

The final script refuses to run until 50 CSV and 50 model PDB files exist.

## Gate Caveat

These PLACER outputs are not complete until the downstream crop ligand geometry screen has run. The project strict filter remains:

```text
screen PLACER crop ligand geometry before full ligand/protein completion
only passing crops proceed to full completion and strict QMMM manifest preparation
```

Do not upload PLACER PDB/CSV outputs to GitHub.
