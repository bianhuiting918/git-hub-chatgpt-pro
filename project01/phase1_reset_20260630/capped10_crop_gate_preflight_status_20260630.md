# capped10 crop gate preflight status - 2026-06-30

## Purpose

Prepare the PLACER crop ligand geometry screen for the active `capped10_n50_bonded` run without falsely treating incomplete PLACER output as a final pass/fail result.

## Current PLACER Sampling State

Remote run:

```text
/data/bht/project01_phase1_reset_gpu/placer_runs/capped10_n50_bonded
```

At 2026-06-30 14:31 CST:

```text
CSV outputs: 6
model PDB outputs: 6
driver PID: 3692495
current input: 7 / 50
```

At the final preflight check:

```text
NOT_READY csv=7 pdb=7 expected=50
```

This means PLACER sampling is still running. Missing outputs are `NOT_EVALUATED`, not structural failures.

## Partial Gate Smoke

A partial crop gate smoke was run on the outputs available at that moment:

```text
remote out dir: /data/bht/project01_phase1_reset_gpu/placer_crop_gate_capped10_partial
evaluated sequences: 6
evaluated conformers: 300
```

Partial strict result:

```text
CROP_STRICT_PASS: 0 / 300
status_counts: FAIL_BOTH_PROTEIN_AND_LIGAND_GATES = 300
```

Interpretation:

```text
This validates that the crop gate script runs on the capped10 output directory.
It is not the final 50-sequence result.
The final crop gate must wait for 50 CSV and 50 model PDB outputs.
```

## Final Gate Script

Prepared remote script:

```text
/data/bht/project01_phase1_reset_gpu/scripts/run_crop_gate_capped10_final.sh
```

GitHub copy:

```text
project01/phase1_reset_20260630/scripts/run_crop_gate_capped10_final.sh
```

Protection:

```text
The script exits with NOT_READY unless csv_count >= 50 and pdb_count >= 50.
```

Only the final crop gate output can determine which conformers proceed to full ligand/protein completion and strict QMMM manifest preparation.
