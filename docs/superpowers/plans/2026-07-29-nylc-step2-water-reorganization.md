# NylC A1 Step2 Water Reorganization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox checkpoints and must preserve technical/scientific gate separation.

**Goal:** Sample natural MM-water exchange around the two frozen Step1 acyl endpoints, identify reproducible preorganized attack-water frames, and only then promote selected complete waters into the Step2 QM region for unrestrained A2 validation.

**Scope:** SCNet computation only. Reuse the frozen Step1 Hamiltonian for sampling: 146 QM atoms, charge 0, 510 electrons, 6 link atoms, no QM water, DFTB3/3OB-3-1, Amber18 sander.MPI. Keep all 40,990 waters as MM during sampling. No PMF, NEB, string, TS, barrier, or mechanism claim.

**Frozen inputs**

- seed26723 restart: `a1_activated_nac_20260726/qmmm/a1_step1_acyl_release_md_continuation/attempt_62216380_0/release_md_endpoint.rst7`
  SHA256 `5d8f76d2c90e3e8c707b640c55a93938f53e18dc30d6f93d54adda467da25f41`
- seed26737 restart: `a1_activated_nac_20260726/qmmm/a1_step1_acyl_release_md_continuation/attempt_62216380_1/release_md_endpoint.rst7`
  SHA256 `4cc60ad4d7be099b3f76040f1ee8c49b91deecebb20ed98570bc192b3d511132`
- prmtop: `a1_activated_nac_20260726/qmmm/a1_unified_core_dftb3_preflight/attempt_62011285/prepared/system.prmtop`
  SHA256 `a61d15bf0bf78675be93275d45f274e808ed6ae450fc1ca21a8e14aee8c12ca0`
- reactive atoms: Nalpha 8949, Thr267 OG1 8960, transferred HG1 8961, L2 C12 10287, O2 10288, N3 10289.
- Gate loop remains residues 261-266; Thr267 is excluded.

---

### Task 1: Add one contract test and observe RED

**Create**

- `workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_nylc_a1_step2_water_reorganization_sampling.py`

- [ ] Assert the production driver exposes a machine-readable `describe` contract with:
  - 2 endpoint seeds × 4 independent velocity replicas = 8 array tasks.
  - 4000 steps at 0.5 fs = 2 ps per task.
  - exactly 8 MPI ranks per task and no array throttle.
  - Step1 contract `146/q0/510e/6link/no_qm_water`.
  - no reactive-coordinate, water-position, or selected-water identity restraint.
  - complete-water universe scan using triclinic minimum image.
  - strict hit thresholds: C12-OW 2.7-3.5 Å; O2-C12-OW 95-125 degrees; donor-H-Nalpha <=2.5 Å; OW-H-Nalpha >=130 degrees; acyl/proton guards pass; at least 3 consecutive saved frames.
  - maximum 2 temporally independent hits retained per endpoint seed.
  - exact restart/prmtop SHA mapping above.
- [ ] Assert the runner uses a unique attempt root, rejects overwrite, keeps full trajectories in scratch, persists only compact hit evidence/final restart/manifest/SHA, and appends `run_history.tsv/jsonl` under `flock`.
- [ ] Assert the sbatch file is `--array=0-7`, `--ntasks=8`, and contains no `%N` throttle.
- [ ] Run on SCNet before production implementation:
  ```bash
  python tests/test_nylc_a1_step2_water_reorganization_sampling.py
  ```
  Expected: RED only because the three production files do not yet exist or the new contract is absent.

### Task 2: Implement the compact sampling/audit driver

**Create**

- `workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/prepare_audit_nylc_a1_step2_water_reorganization_sampling.py`

**Reuse, without changing their scientific thresholds**

- `scripts/prepare_audit_nylc_a1_step2_qmwater_endpoint.py` for triclinic-cell construction, complete-water identity handling, and geometry primitives.
- `scripts/prepare_audit_nylc_a1_step1_acyl_endpoint_stability.py` for NetCDF/Amber trajectory dispatch, engine diagnostics, and Step1 acyl/proton guards.

- [ ] Implement `describe`, `prepare`, `audit`, and `merge-if-ready` subcommands.
- [ ] Map array index as `seed_index = index // 4`, `replica = index % 4`; assign deterministic distinct `ig` values and record them.
- [ ] Generate a Step1-QM sampling mdin with `nstlim=4000`, `dt=0.0005`, periodic restart/trajectory writes sufficient to evaluate consecutive frames, DFTB3, and no DISANG/nmropt/reactive restraint.
- [ ] Parse the actual Amber engine banner and fail technical evaluation unless 146 QM atoms, q0, 510 electrons, 6 links, and no QM water are confirmed.
- [ ] Scan every complete water in every saved frame using the full triclinic cell. For each water use the donor H that gives the better H-Nalpha/OW-H-Nalpha geometry, without changing water identity mid-hit.
- [ ] Apply the strict thresholds and Step1 acyl/proton guards. Collapse consecutive qualifying frames into events; require >=3 saved frames and temporal independence before retaining at most 2 events per seed across its four replicas.
- [ ] Persist only:
  - `SOURCE_MANIFEST.json`
  - `RESULT.json`
  - `WATER_HITS.jsonl` with frame/time/water atom IDs/geometry/source restart and trajectory SHA
  - selected-frame restart(s), if any
  - `SHA256.tsv`
  - `PASS` or `NOT_EVALUATED`
- [ ] Cross-seed gate:
  - `PASS_STEP2_PREORGANIZED_WATER_REPRODUCED` only when both endpoint seeds supply at least one strict event.
  - `FAIL_NO_REPRODUCED_STEP2_PREORGANIZED_WATER` when all 8 tasks are technically complete but the two-seed criterion is unmet.
  - technical failures remain `NOT_EVALUATED_TECHNICAL`.

### Task 3: Implement the runner and Slurm array wrapper

**Create**

- `workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/run_nylc_a1_step2_water_reorganization_sampling.sh`
- `workflows/nylc_l4_nac_to_l2_rebalance_20260723/slurm/run_nylc_a1_step2_water_reorganization_sampling.sbatch`

- [ ] Runner verifies the exact source restart and prmtop SHA before creating work.
- [ ] Persistent output root:
  `a1_activated_nac_20260726/qmmm/a1_step2_water_reorganization_sampling/attempt_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}`; abort if it exists.
- [ ] Scratch root uses `${SLURM_TMPDIR:-/tmp}/nylc_a1_step2_water_reorg_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}`.
- [ ] Run one 8-rank `sander.MPI` trajectory per array task. Full trajectory stays scratch-only; audit before scratch cleanup and copy only compact evidence plus final/selected restarts.
- [ ] Use traps that preserve the first failure and never convert a partial task to PASS.
- [ ] sbatch: one node, 8 tasks, array0-7 without throttle, CPU-only, a bounded wall time suitable for 2 ps DFTB3 sampling.
- [ ] Append one immutable run record under `flock`; do not edit or overwrite prior attempts.

### Task 4: TDD GREEN and immutable SCNet deployment

- [ ] Commit the RED test to branch `codex/nylc-l4-nac-to-l2-rebalance`.
- [ ] Commit only the three production files after RED is observed.
- [ ] Clone the final commit into a new immutable directory below:
  `/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/code_snapshots/`.
- [ ] On SCNet run:
  ```bash
  python tests/test_nylc_a1_step2_water_reorganization_sampling.py
  python -m py_compile scripts/prepare_audit_nylc_a1_step2_water_reorganization_sampling.py
  bash -n scripts/run_nylc_a1_step2_water_reorganization_sampling.sh
  bash -n slurm/run_nylc_a1_step2_water_reorganization_sampling.sbatch
  ```
- [ ] Verify snapshot HEAD, snapshot file hashes, both restart hashes, prmtop hash, array/rank contract, and absence of reaction/water restraints.
- [ ] Any mismatch stops submission as `NOT_EVALUATED_TECHNICAL`; do not patch the immutable snapshot.

### Task 5: Submit once, then gate the next calculation

- [ ] Perform one read-only duplicate check for the new unique job name.
- [ ] Submit exactly once:
  ```bash
  sbatch --job-name=nylc_a2watdyn_<shortcommit> slurm/run_nylc_a1_step2_water_reorganization_sampling.sbatch
  ```
- [ ] Record the exact job ID, commit, snapshot, inputs, resources, and output root in `RUNBOOK.md` and `run_history.tsv/jsonl`.
- [ ] Do not poll in a local loop. A 20-minute heartbeat may perform one discrete queue/status check and immediately audit terminal tasks.
- [ ] After all 8 tasks finish, report the technical denominator and per-seed event denominator separately.
- [ ] Only if `PASS_STEP2_PREORGANIZED_WATER_REPRODUCED`:
  - promote at most 1-2 complete-water events per endpoint seed into the Step2 QM region;
  - verify the new Hamiltonian from the actual engine banner as 149 QM atoms, q0, expected 518 electrons, 6 links;
  - launch two fully unrestrained A2 validation legs per selected event.
- [ ] If the gate fails reproducibly, do not force or pin a water; move to a product-side reverse coupled boundary strategy.
- [ ] No result in this stage is a TS, PMF, activation barrier, or mechanism proof.
