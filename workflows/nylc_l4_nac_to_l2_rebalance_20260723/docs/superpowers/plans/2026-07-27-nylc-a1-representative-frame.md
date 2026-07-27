# NylC A1 Representative NAC Frame Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Extract and fail-closed audit the 354 ps full explicit-water A1 NAC coordinate from evt25/seed26723 on SCNet.

**Architecture:** A focused Python auditor validates source-frame uniqueness, extracted-coordinate identity, A1 bonded topology, joint NAC geometry, gate and frame-specific contacts. A shell runner performs GROMACS-native extraction and `grompp -maxwarn 0`, writes immutable outputs/run history, and promotes GRO/PDB only after the auditor passes.

**Tech Stack:** Python 3, MDAnalysis, NumPy, GROMACS 2022.2 CPU, Bash, unittest, Slurm xahcnormal.

## Global Constraints

- Computation and generated coordinates run only on SCNet.
- Source TPR SHA256 is `c60078a92c2ace51facde4ef64e453f690177fc88b4d6363427f935944fa2e43`.
- Source XTC SHA256 is `1a54f1b5b9f139b746c22d9e0f7e9a4a94eb8154bf2b881888986eedca933d89`.
- Selected time is exactly 354.000 ps relative to the 0--1000 ps `npt300free` trajectory.
- Reactive index1 identities are Thr267 OG1 8960 and L2 C12/O2/N3 10287/10288/10289; each must be revalidated against the pinned TPR.
- A1 requires N bonded to H1/H2/HG1 and OG1 not bonded to HG1.
- Joint NAC is C--OG1 <=0.35 nm and O--C--OG1 95--115 degrees.
- Gate is residues 261--266 and excludes Thr267.
- Preserve System, waters, ions and box; no trajectory/topology/checkpoint/secret goes to GitHub.
- Scientific scope is classical fixed-topology NAC preorganization, not proton transfer, TS, PMF or barrier evidence.

---

### Task 1: A1 representative-frame auditor

**Files:**
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/audit_nylc_a1_representative_frame.py`
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_nylc_a1_representative_frame_contract.py`

**Interfaces:**
- Consumes: `audit_frame(tpr, xtc, gro, ndx, time_ps, source_hashes) -> dict`.
- Produces: JSON fields for source hashes/time match, atom mapping, A1 bonds, NAC geometry, gate, exact frame-specific minimum-contact partners, atom/box counts and scientific status.

- [ ] **Step 1: Write failing contract and unit tests**

Test exact constants, one-frame time selection, PBC angle/distance, A1 bond assertions,
contact partner identity, rejection of wrong hashes/atom names/multiple frames, and
`PASS_A1_REPRESENTATIVE_NAC_FRAME` only when every gate passes.

- [ ] **Step 2: Run red tests**

Run:
`python -m unittest -v workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_nylc_a1_representative_frame_contract.py`

Expected: FAIL because the auditor does not yet exist.

- [ ] **Step 3: Implement the minimal auditor**

Use `MDAnalysis.Universe(tpr, xtc)` to require one source timestamp within
0.001 ps and `MDAnalysis.Universe(tpr, gro)` for the extracted frame. Compare
source/extracted coordinates under minimum image with a 0.0015 nm maximum
per-atom displacement tolerance. Validate `atom.index+1`, residue/name,
`bonded_atoms`, 133589 atoms, 79/33 L2 total/heavy atoms, finite box,
NAC geometry, gate opening imported from the existing primitive generator, and
frame-specific argmin ligand-protein and ligand-water heavy-atom pairs.

- [ ] **Step 4: Run green tests and full A1 contract suite**

Run the new test, then:
`python -m unittest discover -v workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests -p 'test_nylc_a1_*contract.py'`

Expected: all PASS.

- [ ] **Step 5: Commit**

Commit message: `feat: audit A1 representative NAC frame`.

### Task 2: Immutable SCNet extraction runner

**Files:**
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/run_nylc_a1_representative_frame.sh`
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/slurm/run_nylc_a1_representative_frame.sbatch`
- Modify: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_nylc_a1_representative_frame_contract.py`

**Interfaces:**
- Consumes the immutable evt25/seed26723 free-run root and A1 build topology.
- Produces an attempt-specific directory containing temporary extraction,
  grompp evidence, audit JSON, hashes, and promoted representative GRO/PDB only
  after PASS.

- [ ] **Step 1: Add failing runner contract tests**

Require `source.tmp.gro`, System selection, `-dump 354`, pinned source hashes,
`grompp -maxwarn 0`, no overwrite, promotion only after PASS, no M1 proton
analysis, and locked run-history STARTED/terminal rows.

- [ ] **Step 2: Run red tests**

Expected: FAIL because runner/Slurm files do not exist.

- [ ] **Step 3: Implement runner and CPU Slurm wrapper**

Use `/public/software/apps/gromacs/2022.2/hpcx-gcc7.3.1/bin/gmx_mpi`.
Write attempt-specific outputs. Verify hashes before extraction. Run
`printf 'System\n' | gmx trjconv -dump 354 -s run.tpr -f run.xtc -o source.tmp.gro`.
Run the Python auditor, `grompp -maxwarn 0` with the A1 topology and
`em_cg_flexible_m1.mdp`, verify the preflight TPR has no restraints, then
promote GRO and create PDB with GROMACS. Preserve NOT_EVALUATED on failure.

- [ ] **Step 4: Run tests and commit**

Expected: full A1 contract suite PASS.
Commit message: `feat: extract A1 representative NAC frame`.

### Task 3: SCNet execution and evidence

**Files:**
- Modify: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/RUNBOOK.md`
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/audit/nylc_a1_representative_frame_result_20260727.json`

- [ ] **Step 1: Deploy immutable GitHub commit to SCNet**

Do not alter the dirty canonical checkout. Create a commit-addressed code
snapshot or isolated checkout in the task root and record its SHA.

- [ ] **Step 2: Check queue and submit one CPU extraction/audit job**

Use `xahcnormal`, modest CPU/memory and no Dell computation. Record exact job ID.

- [ ] **Step 3: Verify technical and scientific output**

Require Slurm 0:0, PASS audit, exact source/output hashes, correct A1 bonds,
joint NAC, acceptable frame-specific contacts, `grompp -maxwarn 0`, and
viewer PDB/GRO hashes. A failure remains NOT_EVALUATED/FAIL with reason.

- [ ] **Step 4: Commit compact evidence and RUNBOOK update**

Do not commit GRO/PDB/TPR/XTC/topologies. Record the remote coordinate path,
job ID, geometry, hashes, contacts and scientific boundary.
