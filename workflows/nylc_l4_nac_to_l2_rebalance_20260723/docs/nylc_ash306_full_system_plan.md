# NylC-C18 ASH306 Full-System Preflight Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Produce and audit a neutral full-system ASH306 NylC-C18 434 ps NAC microstate that passes GROMACS and ParmEd preprocessing.

**Architecture:** A pure Python builder verifies immutable inputs, performs chain replacement plus distant-Na removal, emits compact audit metadata, and is invoked by one CPU-only Slurm wrapper. The wrapper runs grompp and ParmEd conversion, records all outcomes, and promotes PASS only after every topology, charge, geometry, and numerical gate succeeds.

**Tech Stack:** Python 3.9, GROMACS 2022.1, ParmEd, pytest, Slurm on SCNet.

## Global Constraints

- Compute only on SCNet; never on Dell.
- Gate residues are 261-266; Thr267 is excluded from the gate group.
- Do not overwrite source trajectories or earlier job directories.
- Keep large topology, coordinate, TPR, prmtop, and rst7 artifacts off GitHub.
- A technical PASS is not a TS, RC, PMF, barrier, or biological conclusion.

---

### Task 1: Pure full-system builder and audit

**Files:**
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/build_nylc_ash306_full_system.py`
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_nylc_ash306_full_system_contract.py`

**Interfaces:**
- Consumes: immutable source GRO/build topology and probe job 61717760.
- Produces: `system.gro`, extracted `topol_Protein_chain_H.itp`, copied lightweight topology inputs, and `build_audit.json`.

- [ ] Write failing tests for immutable SHAs, atom ranges, chain replacement, terminal-oxygen symmetry, farthest-Na selection, NA 144-to-143 update, and compact audit fields.
- [ ] Run the targeted pytest file and confirm failure is caused by the absent builder.
- [ ] Implement only the builder behaviors named in the tests.
- [ ] Run py_compile and targeted pytest; expect all tests PASS.
- [ ] Commit tests and builder to `codex/nylc-l4-nac-to-l2-rebalance`.

### Task 2: CPU-only grompp and ParmEd wrapper

**Files:**
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/slurm/run_nylc_ash306_full_system_preflight.sbatch`
- Modify: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_nylc_ash306_full_system_contract.py`

**Interfaces:**
- Consumes: Task 1 builder.
- Produces: job-scoped `PASS.json` or `FAIL.json`, GROMACS preflight TPR, Amber prmtop/rst7, and appended TSV/JSONL history.

- [ ] Add failing wrapper contract tests for xahcnormal, no mdrun, grompp, ParmEd conversion, absolute job-scoped output, ERR trap, and run-history updates.
- [ ] Run tests and confirm RED.
- [ ] Implement the wrapper with 1 task, 4 CPU cores, 10 GB memory, and a 1-hour limit.
- [ ] Run bash syntax validation and the complete targeted test; expect PASS.
- [ ] Commit wrapper and tests.

### Task 3: Execute and promote evidence

**Files:**
- Modify: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/RUNBOOK.md`
- Create after successful run: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/audit/nylc_C18_trueT267_freeGS.ash306_full_system_preflight.json`

**Interfaces:**
- Consumes: tested Task 1 and Task 2 artifacts.
- Produces: authoritative SCNet job ID and compact GitHub audit.

- [ ] Mirror only tested text files from GitHub to SCNet.
- [ ] Check the live SCNet queue before submission.
- [ ] Submit one preflight job and monitor conditionally until terminal.
- [ ] Inspect PASS/FAIL JSON, GROMACS messages, total charge, atom count, NAC geometry, minimum heavy contact, and run-history entries.
- [ ] If PASS, publish only the compact audit and RUNBOOK update; if FAIL, preserve the exact category and return to root-cause investigation.
