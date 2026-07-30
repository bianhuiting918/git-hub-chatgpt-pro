# NylC A1 Bidirectional Step1 Follow-up Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Submit a two-chain inherited reverse continuation and a twelve-task matched forward first-window calibration on SCNet.

**Architecture:** Two independent route implementations reuse the frozen Step1 geometry, restraint, QM-authority, technical-audit, and hashing functions already used by `prepare_audit_nylc_a1_step1_reverse_force_calibration.py`. Each route has one Python driver, one runner, and one sbatch wrapper. Outputs use new job-ID-scoped roots and never alter job 62441145.

**Tech Stack:** Python 3, ParmEd, Amber18 `sander.MPI`, Bash, Slurm, GitHub immutable snapshots.

## Global Constraints

- Frozen prmtop SHA256: `a61d15bf0bf78675be93275d45f274e808ed6ae450fc1ca21a8e14aee8c12ca0`.
- Frozen Step1 contract: 146 QM atoms, q0, 510 electrons, six link H, no QM water.
- Reverse sources: task6 SHA `54ecabec3cbb993a81b36ce86630355e2143797481d71726ca83edd4d65862df`; task8 SHA `2f30ffd0746cc7c53f36ac102b446256ea575e24a2d5bfa9da15a1455db96a7f`.
- Forward sources: seed26723 SHA `9d93b60aee9e8d6fa97493396d757f3bf4d0a47595591968e3debaebb2640e86`; seed26737 SHA `a12c018e195c32b34239004aea36ebc80de8275c9f9332b44084ad8d88d4db0e`.
- Failed or guard-rejected restarts are never inherited.
- No automatic release, shooting, PMF, NEB, string, or Step2 action.

---

### Task 1: Contract tests

**Files:**
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_nylc_a1_step1_bidirectional_followup_contract.py`

**Interfaces:**
- Consumes both new drivers as importable modules.
- Produces assertions for task maps, SHA-pinned sources, active coordinate sets, force scales, and inheritance policy.

- [ ] Write a test that asserts reverse tasks are exactly LOW_BIAS/task6/scale1 and TARGET_CLOSE/task8/scale4, each with `max_windows == 8` and strict serial inheritance.
- [ ] Write a test that asserts forward tasks are `2 x 2 x 3 == 12`; `ADDITION_FIRST_FORWARD` activates only attack and carbonyl while qPT/C-N are monitored; `FULLY_CONCERTED_FORWARD` activates attack, carbonyl, C-N, and qPT.
- [ ] Run on the current branch snapshot:
  `python tests/test_nylc_a1_step1_bidirectional_followup_contract.py`
  Expected: import/file failure because production drivers do not exist.
- [ ] Commit the RED test.

### Task 2: Reverse inherited continuation

**Files:**
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/prepare_audit_nylc_a1_step1_reverse_inherited_followup.py`
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/run_nylc_a1_step1_reverse_inherited_followup.sh`
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/slurm/run_nylc_a1_step1_reverse_inherited_followup.sbatch`

**Interfaces:**
- Driver modes: `describe`, `initialize`, `prepare-window`, `audit-window`, `finalize`, `merge-if-ready`.
- Runner maps Slurm tasks 0/1 to task6/task8 sources and loops windows 1..8 serially.
- Sbatch uses array `0-1`, eight MPI ranks per task, no throttle.

- [ ] Implement exact source SHA validation and copy the accepted calibration metadata into a new chain manifest.
- [ ] Generate each target relative to the last accepted restart: C12-N3 -0.08 A, C12-O2 +0.03 A, Nalpha-H -0.05 A, H-N3 +0.05 A; leave OG1-C12 unrestrained.
- [ ] Audit technical completion, actual direction, chemical guards, restart SHA, and stop without inheriting when any required gate fails.
- [ ] Persist per-window compact files plus CHAIN_MANIFEST/RESULT/SHA256.
- [ ] Run the contract test, Python AST parse, and `bash -n` on runner/sbatch. Expected: PASS.
- [ ] Commit the reverse route.

### Task 3: Forward first-window matrix

**Files:**
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/prepare_audit_nylc_a1_step1_forward_force_calibration.py`
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/run_nylc_a1_step1_forward_force_calibration.sh`
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/slurm/run_nylc_a1_step1_forward_force_calibration.sbatch`

**Interfaces:**
- Driver modes: `describe`, `initialize`, `prepare`, `audit`, `merge-if-ready`.
- Runner maps array tasks 0..11 to two seeds, two mechanisms, and scales 1/2/4.
- Sbatch uses array `0-11`, eight MPI ranks per task, no throttle.

- [ ] Implement source SHA and frozen QM/prmtop checks.
- [ ] For `ADDITION_FIRST_FORWARD`, target OG1-C12 -0.04 A and C12-O2 +0.03 A only; record unrestrained qPT and C12-N3 changes without requiring their direction.
- [ ] For `FULLY_CONCERTED_FORWARD`, additionally target C12-N3 +0.08 A, Nalpha-H +0.05 A, and H-N3 -0.05 A.
- [ ] Apply force scales 1/2/4 to the frozen baseline constants and require direction plus chemical guards only for active coordinates.
- [ ] Persist compact results and select the weakest passing scale per seed/mechanism without launching continuation.
- [ ] Run the contract test, Python AST parse, and `bash -n`. Expected: PASS.
- [ ] Commit the forward route.

### Task 4: SCNet deployment and submission

**Files:**
- Create remotely in the immutable snapshot: `GITHUB_COMMIT`, `GITHUB_BLOBS.tsv`, `SNAPSHOT_SHA256.tsv`.
- Append remotely: `RUNBOOK.md`, `run_history.tsv`, `run_history.jsonl`.

**Interfaces:**
- Consumes the final GitHub commit and exact input SHA values.
- Produces two unique Slurm job IDs.

- [ ] Deploy a new immutable snapshot and verify GitHub blobs, full snapshot SHA, test GREEN, AST, four `bash -n` checks, two reverse input SHA, two forward input SHA, and prmtop SHA.
- [ ] Query the user's queue once for the two new unique job names; reuse an equivalent job or block on a conflicting one.
- [ ] Submit reverse array 0-1 and forward array 0-11 exactly once.
- [ ] Record job IDs and update the existing heartbeat to monitor both jobs.
- [ ] Do not poll again in the submission turn.

## Self-review

The plan covers both approved routes, exact source identities, coordinate definitions, failure inheritance, testing, deployment, and scientific boundaries. It contains no placeholders and introduces no PMF, TS, release, or Step2 scope.
