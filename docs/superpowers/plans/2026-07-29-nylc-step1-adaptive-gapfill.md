# NylC A1 Step1 Adaptive Inherited Gap-Fill Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Continue the eight incomplete I1/I2-to-basin routes from job 62267272 using deterministic midpoint refinement and strict restart inheritance.

**Architecture:** A small Python driver imports the existing four-anchor gap-fill authority, resolves each source route's last accepted restart and original target sequence, and generates one midpoint window at a time. A Bash runner keeps each route serial while an unthrottled Slurm array runs all eight routes concurrently.

**Tech Stack:** Python 3, ParmEd, Amber18 `sander.MPI`, Bash, Slurm, JSON manifests, SHA256.

## Global Constraints

- SCNet computation only; local workspace remains read-only.
- Array0-7 without throttle; 8 MPI ranks per task.
- Source task mapping: `(2,3,4,5,10,11,12,13)`.
- Frozen Step1 contract: 146 QM atoms, q0, 510 electrons including 6 link H, 6 link atoms, no QM water.
- Frozen prmtop SHA256: `a61d15bf0bf78675be93275d45f274e808ed6ae450fc1ca21a8e14aee8c12ca0`.
- Source job 62267272 and all old results are read-only.
- Every accepted window must inherit the preceding accepted restart SHA.
- At most eight new windows; stop at destination or bracket width <=1/16 of the original failed interval.
- No unrestrained release, shooting, PMF, NEB, string, TS, barrier, or mechanism claim in this stage.

---

### Task 1: Contract test and RED

**Files:**
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_nylc_a1_step1_adaptive_gapfill.py`
- Test: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_nylc_a1_step1_adaptive_gapfill.py`

**Interfaces:**
- Consumes: no new production code.
- Produces: contract for `task_spec(index)`, `midpoint(low, high)`, `next_bracket(low, high, passed)`, and `describe()`.

- [ ] **Step 1: Write the failing test**

The test must assert:

```python
SOURCE_TASKS = (2, 3, 4, 5, 10, 11, 12, 13)
assert describe()["array_task_count"] == 8
assert describe()["mpi_ranks_per_task"] == 8
assert describe()["maximum_new_windows"] == 8
assert describe()["minimum_bracket_fraction"] == 1 / 16
assert task_spec(0)["source_task_index"] == 2
assert task_spec(7)["source_task_index"] == 13
assert midpoint({"attack": 2.0, "cn": 1.5}, {"attack": 1.8, "cn": 1.7}) == {
    "attack": 1.9, "cn": 1.6
}
assert next_bracket(0.0, 1.0, True) == (0.5, 1.0)
assert next_bracket(0.0, 1.0, False) == (0.0, 0.5)
```

It must also assert that the runner uses `attempt_${ARRAY_JOB}_${INDEX}`, verifies source and prmtop SHA, serially inherits `accepted.rst7`, uses `mpirun --bind-to none -np 8 sander.MPI`, and that the sbatch file contains `#SBATCH --array=0-7` without a throttle.

- [ ] **Step 2: Run the test on a fresh SCNet snapshot**

Run:

```bash
python tests/test_nylc_a1_step1_adaptive_gapfill.py
```

Expected: RED only because the driver, runner, and sbatch files do not exist.

- [ ] **Step 3: Commit the RED test**

Commit only the test to branch `codex/nylc-l4-nac-to-l2-rebalance`.

### Task 2: Adaptive driver

**Files:**
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/prepare_audit_nylc_a1_step1_adaptive_gapfill.py`
- Reuse: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/prepare_audit_nylc_a1_step1_four_anchor_bidirectional_gapfill.py`
- Test: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_nylc_a1_step1_adaptive_gapfill.py`

**Interfaces:**
- Consumes: source `MANIFEST.json`, `RESULT.json`, persisted accepted restart, original route targets and guard functions.
- Produces:
  - `task_spec(task_index: int) -> dict`
  - `midpoint(low: Mapping[str,float], high: Mapping[str,float]) -> dict[str,float]`
  - `next_bracket(low: float, high: float, passed: bool) -> tuple[float,float]`
  - CLI modes `describe`, `initialize`, `prepare-window`, `audit-window`, `finalize`, `merge-if-ready`.

- [ ] **Step 1: Resolve source authority**

For each new array index, map to one old task, require old technical PASS, exactly four completed stages, `SCIENTIFIC_GUARD_STOP`, source restart existence, and recorded SHA match. A missing persisted restart yields `NOT_EVALUATED_SOURCE_AUTHORITY`; it must not fall back to a different coordinate.

- [ ] **Step 2: Implement deterministic midpoint brackets**

Represent the original failed interval as scalar progress `[0.0, 1.0]`. Interpolate every original joint target field with:

```python
def interpolate(low, high, fraction):
    return {key: low[key] + fraction * (high[key] - low[key]) for key in low}
```

After an accepted guard-safe midpoint set `low=midpoint`; after a technically complete guard failure set `high=midpoint`. Never inherit the failed restart.

- [ ] **Step 3: Reuse original window generation and audit**

Call the existing authority's restraint/MDIN, geometry, engine-banner, and guard helpers. Do not redefine thresholds. Record per window:

```json
{
  "window": 0,
  "fraction": 0.5,
  "target": {},
  "input_restart_sha256": "",
  "output_restart_sha256": "",
  "technical_pass": true,
  "guard_pass": true,
  "accepted_for_inheritance": true
}
```

- [ ] **Step 4: Implement terminal states**

Use exactly:

- `PASS_GUIDED_ROUTE_REACHED_DESTINATION`
- `PARTIAL_BRACKET_REFINED`
- `NOT_EVALUATED_TECHNICAL_ADAPTIVE_GAPFILL`

Persist `MANIFEST.json`, `RESULT.json`, accepted restart chain, `PASS.json` or `NOT_EVALUATED.json`, and `SHA256.tsv`.

- [ ] **Step 5: Run GREEN locally on SCNet**

Run:

```bash
python tests/test_nylc_a1_step1_adaptive_gapfill.py
python -m py_compile scripts/prepare_audit_nylc_a1_step1_adaptive_gapfill.py
```

Expected: all tests PASS with exit 0.

- [ ] **Step 6: Commit driver implementation**

Commit the driver after GREEN.

### Task 3: Runner and Slurm wrapper

**Files:**
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/run_nylc_a1_step1_adaptive_gapfill.sh`
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/slurm/run_nylc_a1_step1_adaptive_gapfill.sbatch`
- Test: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_nylc_a1_step1_adaptive_gapfill.py`

**Interfaces:**
- Consumes: driver CLI and immutable snapshot environment variables.
- Produces: eight immutable attempt directories and one cross-task compact audit.

- [ ] **Step 1: Implement serial inheritance loop**

Use:

```bash
for WINDOW in 0 1 2 3 4 5 6 7; do
    prepare-window --input-rst7 "$ACCEPTED_RST7"
    mpirun --bind-to none -np 8 sander.MPI ...
    audit-window
    if accepted; then
        ACCEPTED_RST7="$WINDOW_DIR/stage.rst7"
    fi
    if terminal; then break; fi
done
```

The audit result, not engine exit alone, controls inheritance.

- [ ] **Step 2: Enforce immutable paths and compact persistence**

Output:
`a1_activated_nac_20260726/qmmm/a1_step1_adaptive_gapfill/attempt_${ARRAY_JOB}_${INDEX}`.

Scratch:
`${SLURM_TMPDIR:-/tmp}/nylc_a1_step1_adaptive_gapfill_${ARRAY_JOB}_${INDEX}`.

Reject pre-existing outputs, verify source restart SHA and prmtop SHA, preserve first failure, and append `run_history.tsv/jsonl` under `flock`.

- [ ] **Step 3: Add sbatch contract**

Use one node, 8 MPI ranks, array0-7 without throttle, CPU-only, and a bounded wall time. Verify snapshot commit/SHA and run the targeted contract test before executing the runner.

- [ ] **Step 4: Run static verification**

Run:

```bash
bash -n scripts/run_nylc_a1_step1_adaptive_gapfill.sh
bash -n slurm/run_nylc_a1_step1_adaptive_gapfill.sbatch
python tests/test_nylc_a1_step1_adaptive_gapfill.py
```

Expected: exit 0.

- [ ] **Step 5: Commit runner and sbatch**

Commit only the two production wrappers.

### Task 4: Immutable deployment and one submission

**Files:**
- Create remotely: new immutable SCNet snapshot under `code_snapshots/`.
- Modify remotely: `RUNBOOK.md`, `run_history.tsv`, `run_history.jsonl`.

**Interfaces:**
- Consumes: final GitHub commit and frozen inputs.
- Produces: one Slurm job ID for array0-7.

- [ ] **Step 1: Deploy a fresh immutable snapshot**

Use old-Git-compatible commands with `cwd`, verify exact HEAD, write `GITHUB_COMMIT` and `SNAPSHOT_SHA256.tsv`.

- [ ] **Step 2: Verify before submission**

Require targeted test GREEN, Python compile, Bash syntax, snapshot SHA, prmtop SHA, all eight source-restart SHA bindings, array0-7, 8 ranks, no throttle, and no output collisions.

- [ ] **Step 3: Check one unique job name and submit exactly once**

Submit only after one read-only duplicate check. Record exact job ID, commit, snapshot, inputs, resources, and output root. Do not perform a post-submit polling loop.

- [ ] **Step 4: Scientific handoff**

Report technical denominator, destinations reached, brackets refined, and remaining guard failures separately. Do not launch unrestrained release or shooting automatically.
