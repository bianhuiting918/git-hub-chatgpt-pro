# Five-candidate PA66-L4 NAC to PA66-L2 rebalance implementation plan

> **For agentic workers:** Execute inline on SCNet. Do not use Dell for computation and do not create a local worktree.

**Goal:** Rebuild five fully unrestrained L4 NAC frames as the same audited PA66-L2 chemistry, then independently equilibrate and compare their unconstrained NAC retention.

**Architecture:** A pure-Python ITP/GRO graph builder discovers candidate-specific source-to-L2 mappings under explicit endpoint-cut constraints. A preparation driver creates isolated candidate directories and machine-readable audits. Slurm jobs run the same staged equilibration protocol for every build that passes preprocessing and EM.

**Tech Stack:** Python 3 standard library, GROMACS-DCU 2022.1, Slurm, GitHub branch `codex/nylc-l4-nac-to-l2-rebalance`.

## Global Constraints

- SCNet root: `/work/home/acshdt1dks/nylon_pa66_scnet_20260708`.
- Task root: `/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723`.
- Audited L2 topology SHA-256: `b0e753c60fd4b71c282d21cc6106a15e73d91d12a20d80e92dd01516162eb301`.
- Preserve all five candidates and continue after per-candidate failures.
- Do not overwrite source trajectories or source branch files.
- Do not commit trajectories, TPR/CPT/EDR/GRO files, large topologies, credentials or secrets.
- NylC gate is residues 261-266; Thr267 is excluded from the gate group.
- Restrained stages are never scientific NAC evidence.

---

### Task 1: Freeze source manifest and extraction audit

**Files:**
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/manifests/candidates.json`
- Create remotely: `inputs/source_manifest.resolved.json`
- Create remotely: `inputs/source_sha256.tsv`

- [ ] Verify every segment, XTC, TPR, source ITP, topology and index exists.
- [ ] Verify each representative cumulative time against geometry.tsv and determine the local 0-100 ps XTC time.
- [ ] Extract one GRO per candidate without modifying source directories.
- [ ] Hash all resolved source inputs and record exact global/local reactive C/O/N indices.

**Expected:** five unique candidate IDs and five extracted coordinate files; no source path changes.

### Task 2: Implement graph mapping with TDD

**Files:**
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_l4_to_l2_graph.py`
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/build_l4_to_l2.py`

- [ ] Write failing tests for ITP parsing, C-O-N bonded-triplet validation, endpoint-only degree mismatch, two genuine cut edges, NylC-C23 regression mapping, NylC-C18 non-C23 mapping, and Nyl12 J1/J2 local numbering.
- [ ] Run the tests on SCNet and capture the expected missing-module failure.
- [ ] Implement minimal graph matching and rerun to green.
- [ ] Add tests for hydrogen parent equivalence, terminal ammonium reconstruction, carboxylate O reconstruction, source C/O/N coordinate identity and GRO atom-count delta.

**Expected:** all pure-Python tests pass and the C23 oracle matches the audited historical mapping.

### Task 3: Build and audit all five candidates

**Files:**
- Create remotely per candidate: `candidates/<id>/build/`
- Create remotely per candidate: `candidates/<id>/audit/build_audit.json`
- Create remotely: `audit/build_summary.tsv`

- [ ] Extract frames into task-owned input paths.
- [ ] Run the builder independently for every candidate, continuing after failures.
- [ ] Copy only required topology components into candidate directories and replace exactly one ligand molecule.
- [ ] Record source-to-L2 mapping, cuts, counts, charge, reactive geometry and minimum nonbonded distance.
- [ ] Run grompp smoke preprocessing for all candidates.

**Expected:** each row is PASS_BUILD or NOT_EVALUATED_BUILD with an exact reason.

### Task 4: Energy minimization gate

**Files:**
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/mdp/em.mdp`
- Create remotely per candidate: `candidates/<id>/em/`
- Update remotely: `audit/em_summary.tsv`

- [ ] Preprocess EM without undocumented warnings.
- [ ] Run EM for every PASS_BUILD candidate.
- [ ] Scan logs for FATAL, NaN, LINCS and SETTLE.
- [ ] Record convergence, maximum force, potential energy, reactive geometry and minimum nonbonded distance.

**Expected:** each candidate becomes PASS_EM or NOT_EVALUATED_EM/FAIL_EM; later candidates still run.

### Task 5: Submit independent staged equilibration

**Files:**
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/run_candidate.sbatch`
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/submit_candidates.sh`
- Create: MDPs for 50 K, 150 K, 300 K, restrained NPT, release NPT and free NPT.
- Update remotely: `run_history.tsv` and `run_history.jsonl`.

- [ ] Query live queue and request one DCU/task with modest CPU and memory per independent candidate.
- [ ] Submit only PASS_EM candidates, each with its own output lineage.
- [ ] Chain stages through afterok dependencies or a single fail-fast per-candidate batch script.
- [ ] Record job IDs, parameters, input/output paths and exit states without secrets.

**Expected:** each eligible candidate is independently queued/running/completed; a failed candidate does not cancel siblings.

### Task 6: Unconstrained comparison and GitHub audit

**Files:**
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/analyze_free_npt.py`
- Create remotely: `audit/free_npt_comparison.tsv`
- Create on GitHub: compact manifests, audit summaries and run history.

- [ ] Analyze only the final >=1 ns fully unrestrained NPT trajectory.
- [ ] Compute PBC-aware distance, angle, joint NAC occupancy and branch-specific gate opening.
- [ ] Extract temperature/pressure statistics and scan numerical warnings.
- [ ] Mark DFTB3/3OB-3-1 preflight eligibility only from unconstrained evidence.
- [ ] Commit only scripts, RUNBOOK, manifests and compact text/JSON/TSV audits.

**Expected:** five-row formal comparison with PASS/FAIL/NOT_EVALUATED reasons and no large simulation artifacts in GitHub.
