# NylC A1 Nine-Replica NAC Audit Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Produce an auditable scientific comparison of nine fully unrestrained NylC A1 trajectories and select only eligible stable NAC frames.

**Architecture:** Reuse the established NylC primitive generator and per-replica auditor, adding an A1-specific manifest/wrapper that maps the nine completed runs to their authoritative source indices. Run independent array tasks on SCNet, then merge and rank small JSON/TSV outputs.

**Tech Stack:** Python 3.10, GROMACS/DCU 2022.1, Slurm, JSON/JSONL/TSV, GitHub branch `codex/nylc-l4-nac-to-l2-rebalance`.

## Global Constraints

- Compute only on SCNet.
- Gate residues are 261-266; Thr267 is excluded from the gate group.
- NAC is joint distance <=0.35 nm and angle 95-115 degrees.
- Analyze only the fully unrestrained 1 ns `npt300free` windows.
- Preserve all nine denominator slots and every `NOT_EVALUATED_*` reason.
- Do not commit trajectories, topology-size artifacts, credentials, or secrets.

---

### Task 1: Freeze the A1 audit universe

**Files:**
- Create: `manifests/nylc_a1_nac_audit_universe.json`
- Test: `tests/test_nylc_a1_nac_audit_contract.py`

**Interfaces:**
- Consumes: nine `EQUILIBRATION_COMPLETE.json` files and their `npt300free` directories.
- Produces: exactly nine candidate/seed records with run root, source-cycle index, and A1 microstate.

- [ ] Write a failing contract test requiring 3 candidates x 3 seeds, unique slots 0-8, gate 261-266, and fully unrestrained input.
- [ ] Run the test and verify it fails because the manifest is absent.
- [ ] Add the frozen manifest with exact SCNet paths and source candidate indices.
- [ ] Run the test and verify it passes.

### Task 2: Generate and audit primitives independently

**Files:**
- Create: `scripts/run_nylc_a1_nac_audit.sh`
- Create: `slurm/run_nylc_a1_nac_audit.sbatch`
- Modify: `tests/test_nylc_a1_nac_audit_contract.py`

**Interfaces:**
- Consumes: one manifest slot, `run.tpr/run.xtc/run.edr/run.log`, existing `generate_nylc_m1_ensemble_primitives.py`, and `audit_nylc_m1_ensemble_replica.py`.
- Produces: one immutable replica audit directory with primitives, `replica_audit.json`, hashes, and run-history entries.

- [ ] Add failing tests for array size 0-8, immutable code snapshot, numerical scan, joint NAC thresholds, and no-overwrite behavior.
- [ ] Run tests and confirm the new runner/launcher are missing.
- [ ] Implement the minimal wrapper and launcher.
- [ ] Run Python contracts and shell syntax checks until green.
- [ ] Submit the nine-way SCNet audit array and record the exact job ID.

### Task 3: Merge, rank, and document

**Files:**
- Create: `scripts/merge_nylc_a1_nac_audits.py`
- Create: `slurm/run_nylc_a1_nac_merge.sbatch`
- Modify: `RUNBOOK.md`
- Modify: `tests/test_nylc_a1_nac_audit_contract.py`

**Interfaces:**
- Consumes: nine independent audit JSON files.
- Produces: denominator-complete TSV/JSON, candidate summaries, eligible-frame ranking, and final scientific status.

- [ ] Add failing tests that require nine slots and prohibit silently dropping technical failures.
- [ ] Implement merge/ranking with candidate-level and replica-level denominators.
- [ ] Run tests and independent checksum/status verification.
- [ ] Submit merge with `afterany` dependency on the audit array.
- [ ] Update RUNBOOK and run history with technical and scientific outcomes.
