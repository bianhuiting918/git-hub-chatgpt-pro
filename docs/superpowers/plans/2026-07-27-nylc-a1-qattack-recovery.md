# NylC A1 q_attack Recovery Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Run three controlled recovery scans that correct the Amber outer-tail restraint problem and identify a technically clean Step1 attack bracket seed.

**Architecture:** Add recovery-specific prepare, audit, runner and Slurm array files so completed job 62012919 remains immutable. All array elements consume the same q03 restart and differ only in distance force constant.

**Tech Stack:** Python 3, unittest, Bash, Slurm, Amber18 sander.MPI, ParmEd.

## Global Constraints

- Compute only on SCNet.
- Keep the frozen 146-atom unified QM core, charge 0, 510 electrons and six link H.
- Do not overwrite job 62012919.
- Do not call a constrained scan a TS, PMF, barrier or mechanism.
- GitHub stores scripts, tests, RUNBOOK and compact audits only.

---

### Task 1: Lock the recovery contract

**Files:**
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_nylc_a1_unified_step1_qattack_recovery_contract.py`

- [ ] Assert the source is job 62012919 q03, targets are 2.65 through 1.65 A, steps are 150, r4 is 4.50 A, and force constants are 50/100/200.
- [ ] Assert the frozen QM count, charge, electron/link counts, Step1 zero-QM-water rule and hard-error gates.
- [ ] Run the new test and confirm it fails because recovery scripts do not exist.

### Task 2: Implement prepare and audit scripts

**Files:**
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/prepare_nylc_a1_unified_step1_qattack_recovery.py`
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/audit_nylc_a1_unified_step1_qattack_recovery.py`

- [ ] Validate the immutable PASS and q03 restart from job 62012919.
- [ ] Generate six sequential windows for one force constant, with r4=4.50 A.
- [ ] Audit geometry, completion and all hard numerical patterns.
- [ ] Emit technical status separately from the attack-bracket seed gate.

### Task 3: Implement immutable runner and Slurm array

**Files:**
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/run_nylc_a1_unified_step1_qattack_recovery.sh`
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/slurm/run_nylc_a1_unified_step1_qattack_recovery.sbatch`

- [ ] Map array indices 0-2 to force constants 50/100/200.
- [ ] Use eight MPI ranks and distinct output directories.
- [ ] Append START and terminal records to TSV/JSONL run history.
- [ ] Hash the immutable snapshot and all compact authority files.

### Task 4: Verify, deploy and submit

**Files:**
- Modify: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/RUNBOOK.md`
- Modify: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/audit/nylc_a1_unified_core_definition_20260727.json`

- [ ] Run unit tests, Python compilation and Bash syntax checks on SCNet.
- [ ] Freeze the GitHub commit into a verified remote code snapshot.
- [ ] Check the queue, submit one three-element array, and record the exact job ID.
- [ ] After completion, independently compare achieved coordinates and select the weakest-force passing replica.
