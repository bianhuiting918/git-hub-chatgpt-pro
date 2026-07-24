# NylC Step1 Catalytic-Network GS Screen Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Re-rank the 24 fully unrestrained NylC-C18 PA66-L2 NAC frames by complete Step1 catalytic preorganization before any TS search.

**Architecture:** A Slurm wrapper runs only on SCNet and derives PBC-aware distance/angle series from the immutable 1 ns trajectory. A focused Python analyzer joins those series with the existing strict-NAC and potential-energy data, records every frame, and only nominates a GS if Thr attack geometry and the Tyr146/Asn219 oxyanion-hole geometry coexist.

**Tech Stack:** GROMACS 2022.1 selection tools, Python 3.9 standard library, pytest contract tests, JSON/TSV audit outputs.

## Global Constraints

- Computation runs only on SCNet; no Dell calculation.
- Input trajectory `npt300free/run.xtc`, TPR and EDR are immutable.
- NylC-C18 reactive mapping is Thr267 OG1 8961, L2 C12/O2/N3 10287/10288/10289.
- Gate remains residues 261-266 and excludes Thr267.
- Scientific strict NAC remains distance <= 0.35 nm and attack angle 95-115 degrees.
- GitHub receives scripts, tests, RUNBOOK and compact text/JSON audits only; never trajectories, TPR, topology payloads or secrets.
- A frame lacking oxyanion-hole geometry is not promoted to a Step1 GS solely because its MM potential is low.

---

### Task 1: Complete catalytic-network frame audit

**Files:**
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_nylc_step1_network_screen_contract.py`
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/analyze_nylc_step1_network_nac.py`
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/slurm/run_nylc_step1_network_screen.sbatch`
- Modify: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/RUNBOOK.md`

**Interfaces:**
- Consumes: existing postprocess xvg files for Thr-C distance, attack angle and potential energy; new GROMACS xvg files for Tyr146 N/H-O2, Asn219 ND2/HD21/HD22-O2, Asp306 O-Thr267 N and Asp306-Asp308 O distances.
- Produces: `step1_network_frames.tsv`, `step1_network_audit.json`, and an extracted `source.tmp.gro` promoted to `selected_step1_network_gs.gro` only after time and geometry verification.

- [ ] **Step 1: Write the failing contract test**

The test requires exact active-copy atoms, immutable input paths, no `mdrun`, strict-NAC gates, Tyr146/Asn219 donor geometry, protonation-state caveats for Asp306, mandatory `source.tmp.gro`, and `NOT_EVALUATED_NO_FULL_CATALYTIC_PREORGANIZATION` when no qualified frame exists.

- [ ] **Step 2: Run the test and verify red**

Run:

```bash
python3 -m pytest -q workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_nylc_step1_network_screen_contract.py
```

Expected: failure because the analyzer and Slurm wrapper do not exist.

- [ ] **Step 3: Implement the analyzer and SCNet wrapper**

Use these active-copy atoms:

```text
Tyr146 N/H: 7156/7157
Asn219 ND2/HD21/HD22: 8240/8241/8242
Thr267 N/OG1: 8949/8961
Asp306 OD1/OD2: 9572/9573
Asp308 OD1/OD2: 9591/9592
L2 C12/O2/N3: 10287/10288/10289
```

A full-preorganization frame must first be strict NAC. Oxyanion donor gates are donor-heavy-to-O2 <= 0.35 nm, donor-H-to-O2 <= 0.25 nm and donor-H-O2 angle >= 135 degrees for both Tyr146 and at least one Asn219 hydrogen. Asp distances are reported and ranked but are not hard-gated until the Asp306 protonation microstate is audited.

If qualified frames exist, rank by lower MM potential after all geometry gates. If none exist, write `NOT_EVALUATED_NO_FULL_CATALYTIC_PREORGANIZATION` and do not extract a TS start.

- [ ] **Step 4: Run tests, Python compilation and Bash syntax checks**

Expected: all contract tests pass; `python -m py_compile` and `bash -n` return zero.

- [ ] **Step 5: Submit one CPU-only SCNet analysis job and audit it**

Check `squeue` first, submit on `xahcnormal`, append START and terminal records to TSV/JSONL history, and preserve all failure reasons. A scheduler exit 0 is not enough; inspect `step1_network_audit.json`.

- [ ] **Step 6: Update compact GitHub evidence**

Commit the final compact audit and RUNBOOK outcome. If no full-preorganization frame exists, the next increment is an explicitly bounded NylC recapture/release using oxyanion-hole coordinates, not a TS search and not redocking.
