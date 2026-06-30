# 2026-07-01 Blind PETase Stage 1 Compute Launcher

## Scope

This log continues the blind PETase QM/MM reproduction task. The paper may be used only for broad methodological organization. It is still not allowed to use paper TS coordinates, reaction coordinates, selected CVs, umbrella windows, aimless-shooting trajectories, barriers, rate assignments, or mechanism conclusions before final validation.

## Change

Added a compute-server launcher:

```text
project01/phase2_blind_petase_qmmm_20260630/scripts/launch_blind_stage1_stage2_compute.py
```

The launcher runs the existing Stage 1/2 gate runner and preserves execution evidence under:

```text
project01/phase2_blind_petase_qmmm_20260630/blind_work/00_run_status/
```

Expected evidence files:

```text
compute_launch_summary.md
stage1_stage2_runner.stdout.log
stage1_stage2_runner.stderr.log
stage1_stage2_runner_command.txt
stage1_stage2_gate_status.tsv
stage1_stage2_next_actions.md
```

The launcher does not contain, request, or record passwords. It records the runner exit code and the blind boundary, then returns the runner exit code unchanged.

## Test Artifact

Added:

```text
project01/phase2_blind_petase_qmmm_20260630/tests/test_launch_blind_stage1_stage2_compute.py
```

The test defines the required behavior with a fake runner: preserve exit code 2, write stdout/stderr logs, write the command log, write a launch summary, include the blind boundary, and avoid the word `password` in the summary.

## Verification Status

Local command execution was unstable in this turn: attempts to launch Python intermittently failed at process startup with `helper_unknown_error: apply deny-read ACLs`. Therefore the new launcher test was written but not successfully executed locally in this turn. The GitHub upload and readback verification must be used as the current evidence for file presence and content; the compute server should run the test before using the launcher for production Stage 1/2 execution.

## README Update

The phase README now recommends the launcher as the first compute-server action instead of invoking the runner directly. The runner remains the authoritative gate engine; the launcher only captures evidence and preserves the exit code.

## Next Compute Action

From the repository root on the compute server:

```text
python project01/phase2_blind_petase_qmmm_20260630/tests/test_launch_blind_stage1_stage2_compute.py
python project01/phase2_blind_petase_qmmm_20260630/scripts/launch_blind_stage1_stage2_compute.py
```

If the launcher returns 2, inspect `blind_work/00_run_status/stage1_stage2_gate_status.tsv` and fix the first upstream `blocked` gate before moving downstream to docking, MD, or QM/MM.