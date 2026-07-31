# NylC A1 Activated-Thr Addition-First Pilot Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build and uniquely submit a two-seed activated-A1 QM/MM pilot that proves source/Hamiltonian integrity before six strictly inherited carbonyl-addition windows.

**Architecture:** A focused Python driver extracts two SHA-frozen activated-A1 MM frames into the frozen Amber topology, performs one reaction-coordinate-free authority minimization, and only then runs up to six serial addition-first windows. A minimal runner and Slurm wrapper execute one seed per task; every stage persists immutable inputs, engine evidence, chemical-integrity checks, and restart hashes.

**Tech Stack:** Python 3, MDAnalysis, ParmEd, Amber18 sander.MPI, Bash, Slurm, unittest.

## Global Constraints

- Chemistry is activated A1: Nalpha-H1/H2/HG1+ and OG1-; HG1 is bonded to Nalpha, not OG1.
- Use seed26723 evt25 frame 240/480 ps and seed26737 evt25 frame 103/206 ps.
- Never reuse output restarts from jobs 62471131, 62477118, or 62500360.
- Preserve 146-QM/q0/510e/6link/DFTB3 and verify the actual engine banner.
- Phase 0 has no reaction-coordinate, proton, angle, or positional restraint.
- Phase 1 biases only OG1-C12 and C12-O2 at 24.0 and 30.0 kcal mol-1 A-2.
- Six serial windows use -0.04 A and +0.03 A target increments relative to the accepted Phase-0 geometry.
- A failed, numerically unhealthy, or chemically invalid restart is never inherited.
- No automatic release, shooting, PMF, NEB, string, Step2 product, TS, barrier, or mechanism claim.

## File map

- Delete obsolete neutral-proton RED contract: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_nylc_a1_step1_mmframe_addition_freeproton.py`.
- Create activated contract: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_nylc_a1_step1_activated_mmframe_addition.py`.
- Create driver: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/prepare_audit_nylc_a1_step1_activated_mmframe_addition.py`.
- Create runner: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/run_nylc_a1_step1_activated_mmframe_addition.sh`.
- Create Slurm wrapper: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/slurm/run_nylc_a1_step1_activated_mmframe_addition.sbatch`.

---

### Task 1: Replace the neutral-proton contract with an activated-A1 RED contract

**Files:**
- Delete: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_nylc_a1_step1_mmframe_addition_freeproton.py`
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_nylc_a1_step1_activated_mmframe_addition.py`
- Test: new file

**Interfaces:**
- Consumes production constants `SOURCES`, `ARRAY_TASKS`, `WINDOWS`, `THR267`, `REACTIVE`, `ATTACK_DELTA_A`, `CARBONYL_DELTA_A`, `FORCE_ATTACK`, and `FORCE_CARBONYL`.
- Produces an executable contract locking sources, hashes, A1 topology, Phase-0 freedom, two Phase-1 restraints, and downstream prohibition.

- [ ] **Step 1: Write the failing contract**

The test must assert:

```python
assert mod.ARRAY_TASKS == 2
assert mod.WINDOWS == 6
assert mod.SOURCES[26723]["frame"] == 240
assert mod.SOURCES[26737]["frame"] == 103
assert mod.SOURCES[26723]["tpr_sha256"] == "c60078a92c2ace51facde4ef64e453f690177fc88b4d6363427f935944fa2e43"
assert mod.SOURCES[26723]["xtc_sha256"] == "1a54f1b5b9f139b746c22d9e0f7e9a4a94eb8154bf2b881888986eedca933d89"
assert mod.SOURCES[26737]["tpr_sha256"] == "dbd19a399547319d10630430ed494d33f6271bab0f30cb0de5a466c6af13ba20"
assert mod.SOURCES[26737]["xtc_sha256"] == "fcba14da98b331368061dcd990f2467628ad77b9b7a9c4ce88090f92e0831b05"
assert mod.describe()["starting_state"] == "NALPHA_H3_PLUS_OG1_MINUS"
assert mod.authority_restraints() == ""
text = mod.reaction_restraints(mod.stage_spec(1, {"attack_A": 2.972, "c12_o2_A": 1.220}))
assert text.count("&rst") == 2
assert str(mod.THR267["hg1"]) not in text
assert str(mod.REACTIVE["n3"]) not in text
```

Also assert no literal `\\n`, exactly two physical `&rst` records, the qPT/angle flags are false, and `automatic_downstream_action == "NONE"`.

- [ ] **Step 2: Commit the RED contract**

Commit message: `test: define activated A1 MM-frame addition pilot`.

- [ ] **Step 3: Observe RED on SCNet**

Deploy the RED commit to a new immutable test snapshot and run:

```bash
/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python tests/test_nylc_a1_step1_activated_mmframe_addition.py
```

Expected: one failure because the activated production driver does not yet exist; remaining driver-dependent tests skip.

---

### Task 2: Implement activated source extraction, authority minimization, serial windows, and audits

**Files:**
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/prepare_audit_nylc_a1_step1_activated_mmframe_addition.py`
- Test: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_nylc_a1_step1_activated_mmframe_addition.py`

**Interfaces:**
- Produces `describe() -> dict`, `task_spec(int) -> dict`, `extract_source(int, Path, str) -> dict`, `authority_restraints() -> str`, `stage_spec(int, Mapping) -> dict`, `reaction_restraints(Mapping) -> str`, `prepare_authority(Path, Path) -> dict`, `audit_authority(Path, Path) -> dict`, `prepare_window(Path, Path, int, Path) -> dict`, `audit_window(Path, Path, int) -> dict`, and `finalize(Path) -> dict`.
- Reuses the proven triclinic PBC, numerical-health, QM-banner, geometry, and A1 covalent-integrity helpers from existing Step1 drivers without changing their old outputs.

- [ ] **Step 1: Implement SHA-bound source constants**

Use the exact paths and hashes from the approved design. `extract_source` verifies all four source-file hashes before MDAnalysis loads the requested frame, transplants coordinates and the six-value triclinic box into the frozen Amber prmtop, persists `source.rst7`, and records its SHA256.

- [ ] **Step 2: Implement source topology and geometry gates**

Validate exact atom identities and require:

```python
required_bonds = {
    frozenset((THR267["nalpha"], THR267["h1"])),
    frozenset((THR267["nalpha"], THR267["h2"])),
    frozenset((THR267["nalpha"], THR267["hg1"])),
    frozenset((THR267["nalpha"], THR267["ca"])),
    frozenset((THR267["ca"], THR267["cb"])),
    frozenset((THR267["cb"], THR267["og1"])),
}
forbidden_bond = frozenset((THR267["og1"], THR267["hg1"]))
```

Reject before calculation unless required bonds exist, the forbidden bond is absent, OG1-C12 <=3.5 A, and O2-C12-OG1 is 95-115 degrees.

- [ ] **Step 3: Implement Phase-0 reaction-coordinate-free input**

Generate the established DFTB3 minimization input with `nmropt=0`, no `DISANG`, no `&wt`, no positional restraint, `maxcyc=2200`, and `ncyc=550`. `authority_restraints()` returns the empty string. Persist authority input/output/mdinfo/restart, engine tail, and SHA256.

- [ ] **Step 4: Implement the Phase-0 gate**

Require normal engine exit, numerical-health pass, actual 146-QM/q0/510e/6link/DFTB3 banner, intact A1/substrate covalent graphs, OG1-C12 <=3.5 A, attack angle 95-115 degrees, and a new output restart SHA. Otherwise classify the exact technical state in the spec and stop without creating window 1.

- [ ] **Step 5: Implement exactly two Phase-1 restraints**

`reaction_restraints(stage)` emits one physical `&rst` line for OG1-C12 and one for C12-O2. It contains no HG1, N3, qPT, or angle atoms and no literal escaped newline.

- [ ] **Step 6: Implement six strict inherited windows**

Window `n` targets `phase0_attack - 0.04*n` and `phase0_c12_o2 + 0.03*n`. The accepted SHA of stage `n-1` equals the input SHA of stage `n`. Stop on the first technical failure, graph failure, attack-angle escape, or missing same-direction response. Never copy a rejected restart forward.

- [ ] **Step 7: Persist compact immutable results**

Write `SOURCE_MANIFEST.json`, `AUTHORITY_RESULT.json`, per-window `WINDOW_MANIFEST.json` and `RESULT.json`, top-level `RESULT.json`, terminal `PASS.json` or `NOT_EVALUATED.json`, and `SHA256.tsv`. Report exact technical and accepted-window denominators.

- [ ] **Step 8: Verify GREEN**

Run on SCNet:

```bash
python tests/test_nylc_a1_step1_activated_mmframe_addition.py
python -m py_compile scripts/prepare_audit_nylc_a1_step1_activated_mmframe_addition.py
python scripts/prepare_audit_nylc_a1_step1_activated_mmframe_addition.py --mode describe
```

Expected: all contract tests pass, compile exits 0, and describe reports two tasks, one authority stage, six windows, two reaction-coordinate restraints, activated A1, and no automatic downstream action.

- [ ] **Step 9: Commit the driver**

Commit message: `feat: add activated A1 MM-frame addition pilot`.

---

### Task 3: Add the minimal runner and Slurm wrapper

**Files:**
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/run_nylc_a1_step1_activated_mmframe_addition.sh`
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/slurm/run_nylc_a1_step1_activated_mmframe_addition.sbatch`

**Interfaces:**
- Runner consumes `A1_ACTIVATED_ADDITION_CODE_ROOT`, `A1_ACTIVATED_ADDITION_GITHUB_COMMIT`, and `SLURM_ARRAY_TASK_ID=0|1`.
- Produces unique `attempt_<job>_<task>` roots under `a1_step1_activated_mmframe_addition`.

- [ ] **Step 1: Implement the runner**

Initialize the source, run authority minimization with `mpirun --bind-to none -np 8 sander.MPI`, audit immediately, and continue only if authority passes. For windows 1..6, prepare, execute, audit, and break on rejection. Use `${SLURM_TMPDIR:-/tmp}`, isolated scratch and pycache paths, a failure trap, and append-only run history.

- [ ] **Step 2: Implement the Slurm wrapper**

Use `#SBATCH --array=0-1`, `-n 8`, `-c 1`, `--mem-per-cpu=2500M`, `--time=24:00:00`, and a new unique job name. Verify `GITHUB_COMMIT` and `SNAPSHOT_SHA256.tsv` before executing the runner.

- [ ] **Step 3: Verify launch contracts**

Run:

```bash
bash -n scripts/run_nylc_a1_step1_activated_mmframe_addition.sh
bash -n slurm/run_nylc_a1_step1_activated_mmframe_addition.sbatch
python scripts/prepare_audit_nylc_a1_step1_activated_mmframe_addition.py --mode describe
```

Expected: both shell checks exit 0 and describe matches the design.

- [ ] **Step 4: Commit launch files**

Commit message: `ops: launch activated A1 MM-frame addition pilot`.

---

### Task 4: Deploy one immutable snapshot and submit once

**SCNet outputs:**
- New immutable snapshot under `code_snapshots/a1_activated_addition_<shortsha>_runtime_20260731`.
- New attempt roots keyed by the new job ID.
- Updated `GITHUB_COMMIT`, `SNAPSHOT_SHA256.tsv`, `RUNBOOK.md`, and append-only run history.

- [ ] **Step 1: Deploy the exact final GitHub commit**

Use old-Git-compatible commands with an explicit working directory. Verify HEAD equals the full expected commit and the snapshot contains no unexpected work files.

- [ ] **Step 2: Run the complete pre-submit gate**

Require targeted tests, Python compile, two `bash -n` checks, describe, snapshot SHA, prmtop SHA `a61d15bf0bf78675be93275d45f274e808ed6ae450fc1ca21a8e14aee8c12ca0`, and all four TPR/XTC hashes to pass.

- [ ] **Step 3: Check for a duplicate unique job name once**

Run one `squeue` query. Reuse an equivalent job bound to the same commit/snapshot; stop without cancellation if the same name has different inputs.

- [ ] **Step 4: Submit exactly once**

Submit array 0-1 with 8 MPI ranks per task, exporting the immutable code root and full commit. Record the exact job ID and append RUNBOOK/run history. Do not query the queue again in the submission turn.

- [ ] **Step 5: Update the heartbeat**

Track only the exact new job ID. Audit each task immediately when it leaves the queue. No downstream calculation is submitted automatically.
