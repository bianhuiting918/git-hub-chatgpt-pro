# NylC A1 Addition-First Free-Proton Pilot Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build and submit a two-seed, six-window inherited QM/MM pilot that biases only OG1-C12 and C12-O2 while auditing the unconstrained HG1 destination.

**Architecture:** Add one focused Python driver that extracts the two approved MM frames, creates SHA-bound restarts, prepares six serial Amber minimization windows, audits numerical/chemical integrity and all-QM proton acceptors, and writes compact immutable results. Add one runner and one Slurm wrapper following the existing preorganized-calibration launch pattern.

**Tech Stack:** Python 3, ParmEd, MDAnalysis, Amber18 sander.MPI, Bash, Slurm, unittest.

## Global Constraints

- Start from seed26723 evt25 frame 240/480 ps and seed26737 evt25 frame 103/206 ps.
- Start from neutral Thr267 OG1-HG1; do not pre-form OG1-/NalphaH3+.
- Preserve the 146-QM/q0/510e/6link/DFTB3 Hamiltonian and verify the actual engine banner.
- Bias only OG1-C12 and C12-O2 with 24.0 and 30.0 kcal mol-1 A-2.
- Use six serial windows with target increments -0.04 A and +0.03 A per window.
- Do not restrain HG1, qPT, C12-N3, proton acceptors, or attack angle.
- A failed or chemically invalid restart is never inherited.
- Do not start release, shooting, PMF, NEB, string, or Step2 product work.

---

### Task 1: Contract test for the approved sources and free-proton protocol

**Files:**
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_nylc_a1_step1_mmframe_addition_freeproton.py`
- Test: same file

**Interfaces:**
- Consumes: new driver constants `SOURCES`, `ARRAY_TASKS`, `WINDOWS`, `ATTACK_DELTA_A`, `CARBONYL_DELTA_A`, `FORCE_ATTACK`, `FORCE_CARBONYL`.
- Produces: executable contract that locks the two source trajectories, frames, schedule, restraint count, and no-proton-bias rule.

- [ ] **Step 1: Write the failing test**

```python
def test_two_approved_mm_sources():
    assert mod.ARRAY_TASKS == 2
    assert mod.SOURCES[26723]["frame"] == 240
    assert mod.SOURCES[26723]["time_ps"] == 480.0
    assert mod.SOURCES[26737]["frame"] == 103
    assert mod.SOURCES[26737]["time_ps"] == 206.0

def test_six_window_schedule_and_base_forces():
    assert mod.WINDOWS == 6
    assert mod.ATTACK_DELTA_A == -0.04
    assert mod.CARBONYL_DELTA_A == 0.03
    assert mod.FORCE_ATTACK == 24.0
    assert mod.FORCE_CARBONYL == 30.0

def test_restraints_bias_only_attack_and_carbonyl():
    stage = mod.stage_spec(1, {"attack_A": 2.972, "c12_o2_A": 1.220})
    text = mod.restraints(stage)
    assert text.count("&rst") == 2
    assert str(mod.REACTIVE["hg1"]) not in text
    assert str(mod.REACTIVE["n3"]) not in text
```

- [ ] **Step 2: Run test to verify RED**

Run:
```bash
python tests/test_nylc_a1_step1_mmframe_addition_freeproton.py
```

Expected: FAIL because the production driver does not exist.

- [ ] **Step 3: Commit the RED test**

Commit message: `test: define A1 MM-frame free-proton pilot`.

### Task 2: Implement frame extraction, serial windows, and proton-destination audit

**Files:**
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/prepare_audit_nylc_a1_step1_mmframe_addition_freeproton.py`
- Reuse: existing helpers from `prepare_audit_nylc_a1_step1_raw_nac_preorganized_calibration.py`, `prepare_audit_nylc_a1_step1_raw_nac_forward_calibration.py`, and triclinic PBC helper.

**Interfaces:**
- Produces:
  - `describe() -> dict`
  - `task_spec(task_index: int) -> dict`
  - `extract_source(task_index: int, root: pathlib.Path, commit: str) -> dict`
  - `stage_spec(window_index: int, source_geometry: Mapping[str, float]) -> dict`
  - `restraints(stage: Mapping[str, Any]) -> str`
  - `prepare_window(root, scratch, window_index, input_restart) -> dict`
  - `audit_window(root, scratch, window_index) -> dict`
  - `audit_proton_destination(restart: pathlib.Path) -> dict`
  - `finalize(root: pathlib.Path) -> dict`

- [ ] **Step 1: Implement source constants and deterministic extraction**

Use exact TPR/XTC paths from the design, load the requested frame with MDAnalysis, transplant coordinates and triclinic box into the frozen Amber topology, save `source.rst7`, and record SHA256. Validate atom count and Thr267 covalent integrity before any minimization.

- [ ] **Step 2: Implement exactly two restraint lines**

```python
def restraints(stage):
    return "".join([
        tight_distance(REACTIVE["og1"], REACTIVE["c12"],
                       stage["attack_target_A"], FORCE_ATTACK),
        tight_distance(REACTIVE["c12"], REACTIVE["o2"],
                       stage["carbonyl_target_A"], FORCE_CARBONYL),
    ])
```

Reject output unless it contains exactly two physical `&rst` lines and no HG1/N3 atom number.

- [ ] **Step 3: Implement six-window schedule**

For window `n=1..6`:
```python
attack_target = source_attack + n * ATTACK_DELTA_A
carbonyl_target = source_carbonyl + n * CARBONYL_DELTA_A
```

The accepted restart from window `n-1` is the only permitted input to window `n`.

- [ ] **Step 4: Implement multi-acceptor HG1 audit**

Load the restart with ParmEd, enumerate QM-region N/O/S atoms excluding covalently bound OG1, calculate triclinic minimum-image HG1-acceptor distances and OG1-HG1-acceptor angles, sort by distance, and persist every acceptor within 4.0 A. Report nearest identity and classify proton ownership without making proton transfer a window-pass gate.

- [ ] **Step 5: Implement acceptance and stop conditions**

A window passes only if engine banner, numerical health, Thr/substrate integrity, and same-direction attack/carbonyl response pass. Otherwise write `NOT_EVALUATED_TECHNICAL` or `NO_TETRAHEDRAL_RESPONSE`, stop the chain, and do not copy the failed restart forward.

- [ ] **Step 6: Run targeted test to verify GREEN**

Run:
```bash
python tests/test_nylc_a1_step1_mmframe_addition_freeproton.py
python -m py_compile scripts/prepare_audit_nylc_a1_step1_mmframe_addition_freeproton.py
```

Expected: all tests PASS and compile exit 0.

- [ ] **Step 7: Commit production driver**

Commit message: `feat: add A1 MM-frame free-proton pilot`.

### Task 3: Add the minimal runner and Slurm wrapper

**Files:**
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/run_nylc_a1_step1_mmframe_addition_freeproton.sh`
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/slurm/run_nylc_a1_step1_mmframe_addition_freeproton.sbatch`

**Interfaces:**
- Runner consumes `A1_MMFRAME_FREEPT_CODE_ROOT`, `A1_MMFRAME_FREEPT_GITHUB_COMMIT`, and `SLURM_ARRAY_TASK_ID=0|1`.
- Outputs unique `attempt_<job>_<task>` directories and compact manifests under `a1_step1_mmframe_addition_freeproton`.

- [ ] **Step 1: Implement runner**

Initialize once, then loop windows 1 through 6. Run `mpirun --bind-to none -np 8 sander.MPI`. After each engine call, audit immediately and stop unless the window is accepted. Persist run-history STARTED and TERMINAL rows.

- [ ] **Step 2: Implement Slurm wrapper**

Use `#SBATCH --array=0-1`, `-n 8`, 24-hour wall time, unique job name/output pattern, isolated `PYTHONPYCACHEPREFIX`, snapshot SHA verification, and `bash -n` before exec.

- [ ] **Step 3: Verify shell and describe contracts**

Run:
```bash
bash -n scripts/run_nylc_a1_step1_mmframe_addition_freeproton.sh
bash -n slurm/run_nylc_a1_step1_mmframe_addition_freeproton.sbatch
python scripts/prepare_audit_nylc_a1_step1_mmframe_addition_freeproton.py --mode describe
```

Expected: two shell checks exit 0; describe reports two tasks, six windows, two restraints, and no proton bias.

- [ ] **Step 4: Commit launch files**

Commit message: `ops: launch A1 MM-frame free-proton pilot`.

### Task 4: Deploy immutable SCNet snapshot and submit once

**Files:**
- Update on SCNet only: workflow `GITHUB_COMMIT`, `SNAPSHOT_SHA256.tsv`, `RUNBOOK.md`, and append-only run history.

- [ ] **Step 1: Create immutable snapshot from the final GitHub commit**

Use a new path containing the short commit and date. Do not modify existing snapshots.

- [ ] **Step 2: Verify before submission**

Run targeted tests, `py_compile`, two `bash -n` checks, snapshot SHA verification, frozen topology SHA, source TPR/XTC SHA, extracted coordinate SHA, and `describe`.

- [ ] **Step 3: Check duplicate job name once**

Use one `squeue` query for the unique new job name. If an equivalent commit/snapshot job exists, reuse it; if a conflicting job exists, stop.

- [ ] **Step 4: Submit exactly once**

Submit array `0-1` with 8 MPI ranks per task and export the immutable code root plus full GitHub commit. Record the exact job ID. Do not query the queue again in the submission turn.

- [ ] **Step 5: Update heartbeat**

Track the exact new job ID. Audit each task immediately when it leaves the queue; do not wait for the other task.
