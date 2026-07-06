# PETase Blind QM/MM Mechanism Reproduction Runbook

Updated: 2026-07-06 20:35 CST

This document records the blind PETase QM/MM reproduction work from setup through
TS validation and PMF production. It is written so that another researcher can
reconstruct the calculation logic, locate the exact server artifacts, and repeat
or extend each step. Raw trajectories, restart files, and passwords are not
stored in GitHub.

## Scope and Rule

Goal: derive PETase acylation and deacylation mechanisms from structural and
chemical first principles, then compare against the paper only after our TS
ensemble and preliminary PMFs exist.

Blind rule:

- Allowed inputs: PETase structure, generic serine-hydrolase chemistry,
  docking/pose generation, Amber/Sander QM/MM methodology, ATESA/committor/PMF
  ideas.
- Disallowed inputs before validation: paper TS coordinates, paper umbrella
  windows, paper trajectories, paper PMF values, paper rate-limiting conclusion,
  and paper-specific residue-motion conclusions.
- Methodological learning from the paper is allowed only at the level of
  "use QM/MM, TS/committor, and PMF"; concrete coordinates and results are not
  used as seeds.

## Compute Environment

Primary compute host:

```text
CPU server: 210.73.40.29
Remote project root:
/Dell/Dell14/bianht/petase_blind_qmmm/repo/project01/phase2_blind_petase_qmmm_20260630
```

Do not commit or print server passwords. Use the existing SSH/askpass setup from
the Codex working machine.

Main runtime used for the final formal calculations:

```text
Amber/Sander QM/MM
QM method: DFTB3
Temperature: 310 K
QM/MM timestep used in production: dt=0.0001 ps
Thermostat: Langevin, gamma_ln=5.0
Trajectory write interval: ntwx=20
```

The smaller timestep is deliberate. It was chosen after early QM/MM paths showed
numerical fragility under reactive restraints. This makes the protocol slower
than the paper-like 1 fs setup, but it reduces false failures in this blind
workflow.

## Atom Indices and Key CVs

Amber atom indices are 1-based.

Deacylation:

```text
SerO        1911
His NE2     2992
Acyl C      3843
Acyl O      3844
Water O     3863
Water H2    3865

C_SerO  = distance(3843, 1911)
OW_C    = distance(3863, 3843)
H2_His  = distance(3865, 2992)
eta     = C_SerO - OW_C
angle   = angle(3863, 3843, 3844)
```

Acylation:

```text
SerO        1911
SerH        1912
His NE2     2992
Carbonyl C  3843
Leaving O   3845

SerO_C      = distance(1911, 3843)
C_Oleave    = distance(3843, 3845)
SerH_SerO   = distance(1912, 1911)
SerH_His    = distance(1912, 2992)
eta         = SerO_C - C_Oleave
```

## High-Level Workflow

1. Build a structure-only Michaelis complex hypothesis.
   - Use PETase structure and generic serine hydrolase geometry.
   - Use docking/pose filtering to find plausible ester positioning.
   - Do not import the paper reactant state.

2. Build Amber topology and QM/MM input.
   - Use catalytic residues, substrate atoms, and relevant water in the QM
     region.
   - Keep the first formal protocol conservative: Amber/Sander DFTB3, 310 K,
     dt=0.0001 ps.

3. Search for reactive boundary states.
   - Initial scan and path attempts were used only to find chemically plausible
     seeds.
   - ATESA/aimless shooting and focused committor were used to distinguish TS
     core candidates from product-like or reactant-like structures.

4. Validate TS-like seeds by committor or short commitment tests.
   - Endpoint classification used geometry rather than paper labels.
   - For deacylation, seed 1 / seed 4 focused tests led to a usable TS region.
   - For acylation, cand02 was the accepted TS-core candidate.

5. Build PMF windows from our own paths.
   - Deacylation required multiple rounds of gap and bridge windows before the
     br206 -> r210 region had usable overlap.
   - Acylation first used a 13-window path, then bridge/refine windows to fill
     weak eta regions.

6. Recompute PMFs with exact Amber piecewise restraint energies.
   - A major bug was discovered in early deacylation PMFs: an oversimplified
     harmonic bias and legacy temperature assumptions produced false high
     barriers.
   - Correct bias uses the Amber `&rst` piecewise distance potential with
     linear extension outside r1/r4.

7. Use state free energies as the primary barrier estimate.
   - Projected path-bin PMFs are diagnostic only when bins are sparse or empty.
   - State free energies and neighbor overlap are the main decision basis.

## Deacylation Final State

Final completed production:

```text
/Dell/Dell14/bianht/petase_blind_qmmm/repo/project01/phase2_blind_petase_qmmm_20260630/blind_work/09_umbrella_pmf/deacylation_seed1_attempt040_merged037038_38win_longprod_40000
```

Final analysis:

```text
/Dell/Dell14/bianht/petase_blind_qmmm/repo/project01/phase2_blind_petase_qmmm_20260630/blind_work/09_umbrella_pmf/deacylation_seed1_attempt041_from040_longprod_pymbar
```

Key results:

```text
Windows: 38/38 completed
Frames: 76000 total, 2000 per window
Bad windows: 0
SCC warnings: 0
vlimit exceeded: 0
Hard errors: 0
State barrier from r150: 19.061 +/- 0.073 kcal/mol
Barrier state: br208_mid_from_br206 / br208_mid_from_r210
Block sensitivity:
  first half  = 19.136 +/- 0.103 kcal/mol
  second half = 18.965 +/- 0.103 kcal/mol
Weak MBAR overlap: h119 -> h104 = 0.0198
```

Decision:

- The deacylation barrier is stable around 19.1 kcal/mol.
- The only weak MBAR overlap is on the product-side H-transfer/late region, not
  the barrier region.
- Do not spend compute on h119 -> h104 unless later analysis shows a product-side
  thermodynamic question depends on it.

## Acylation Current State

Validated TS-core candidate:

```text
cand02
Combined committor: forward/product-like 18, reactant-like 13, other 1
pB_committed = 0.580645
pB_all       = 0.5625
```

Short PMF/bridge work completed before the long production:

```text
attempt042: 13-window 4D path production, 20000 steps
attempt043: 7 bridge windows, 20000 steps
attempt045: 5 refine windows, 20000 steps
attempt046: merged attempt042/043/045 preliminary MBAR/overlap
```

Active long production started 2026-07-06 20:26 CST:

```text
/Dell/Dell14/bianht/petase_blind_qmmm/repo/project01/phase2_blind_petase_qmmm_20260630/blind_work/09_umbrella_pmf/acylation_cand02_attempt047_merged042043045_25win_longprod_40000
```

Attempt047 design:

```text
Windows: 25
Steps per window: 40000
MAX_WORKERS: 25
Seeds: final restarts from attempts 042, 043, 045
Restraints: copied per-window umbrella.RST from source windows
Input: input/acyl_4d_path_longprod_40000.in
```

Initial health check at 2026-07-06 20:28 CST:

```text
Running sander processes: 25
Driver: 1
NSTEP: 400 / 40000
Temperature range: 307.89 - 312.82 K
SCC warnings: 0
vlimit exceeded: 0
Hard errors: 0
Estimated remaining time at startup: about 3.8 hours
```

## How to Reproduce the Current Results

### 1. Verify the remote environment

Run on the CPU server:

```bash
cd /Dell/Dell14/bianht/petase_blind_qmmm/repo/project01/phase2_blind_petase_qmmm_20260630
/Dell/Dell14/bianht/micromamba/envs/petase_stage1/bin/python - <<'PY'
import sys
for mod in ["mdtraj", "pymbar", "numpy", "pandas"]:
    m = __import__(mod)
    print(mod, getattr(m, "__version__", "ok"))
print(sys.executable)
PY
```

Expected environment:

```text
Python: /Dell/Dell14/bianht/micromamba/envs/petase_stage1/bin/python
mdtraj: 1.11.0
pandas: 2.2.3
```

### 2. Recheck deacylation final analysis

Read the analysis outputs:

```bash
cd /Dell/Dell14/bianht/petase_blind_qmmm/repo/project01/phase2_blind_petase_qmmm_20260630/blind_work/09_umbrella_pmf/deacylation_seed1_attempt041_from040_longprod_pymbar
cat summary.md
head state_free_energies.tsv
cat block_sensitivity_state_barrier.tsv
sort -k7,7g neighbor_overlap_mbar_matrix.tsv | head
```

Important files:

```text
attempt040_cv_error_summary.tsv
state_free_energies.tsv
block_sensitivity_state_barrier.tsv
neighbor_overlap_mbar_matrix.tsv
neighbor_overlap_s_hist.tsv
pathcv_occupied_bins_pmf.tsv
summary.md
```

### 3. Monitor acylation attempt047

```bash
cd /Dell/Dell14/bianht/petase_blind_qmmm/repo/project01/phase2_blind_petase_qmmm_20260630/blind_work/09_umbrella_pmf/acylation_cand02_attempt047_merged042043045_25win_longprod_40000
ps -u "$USER" -o pid,ppid,stat,etime,cmd | grep -E 'run_acyl_25win_longprod_40000.py|sander' | grep -v grep
wc -l umbrella_run_progress.tsv
find windows -name md.out | wc -l
find windows -name md.rst7 | wc -l
find windows -name md.nc | wc -l
```

Each window should eventually reach `NSTEP = 40000`. Treat these as hard
failure signals:

```text
SCC convergence failure
vlimit exceeded
fatal / segmentation / forrtl errors
missing md.rst7 or md.nc after driver completion
```

Ignore these benign Amber diagnostics when scanning error lines:

```text
Ewald error estimate
TESTING RELATIVE ERROR
```

### 4. Analyze acylation after attempt047 completes

Use the same analysis principles as attempt046:

- compute CV summaries from each `md.nc`;
- parse each `umbrella.RST`;
- evaluate the exact Amber piecewise bias;
- run MBAR at 310 K;
- compute neighbor overlap;
- estimate the barrier from state free energies;
- use eta/path PMFs only as diagnostics if bins are sparse.

The required inputs are already in the attempt047 directory:

```text
complex.prmtop
window_manifest.tsv
input/acyl_4d_path_longprod_40000.in
windows/w*/umbrella.RST
windows/w*/md.nc
windows/w*/md.rst7
```

## Amber Piecewise Bias Formula

For each `&rst` restraint:

```python
def rst_energy(x, b):
    r1, r2, r3, r4 = b["r1"], b["r2"], b["r3"], b["r4"]
    rk2, rk3 = b["rk2"], b["rk3"]
    e = np.zeros_like(x, dtype=float)
    lo = x < r1
    lq = (x >= r1) & (x < r2)
    uq = (x > r3) & (x <= r4)
    hi = x > r4
    e[lq] = rk2 * (x[lq] - r2) ** 2
    e[uq] = rk3 * (x[uq] - r3) ** 2
    e[lo] = rk2 * (r1 - r2) ** 2 + 2 * rk2 * (r1 - r2) * (x[lo] - r1)
    e[hi] = rk3 * (r4 - r3) ** 2 + 2 * rk3 * (r4 - r3) * (x[hi] - r4)
    return e
```

This exact piecewise form is essential. Earlier PMF estimates using a simpler
harmonic approximation produced misleading high barriers.

## Decision Logic Used

The "grill me" questions used before starting or extending a calculation:

1. Does the existing data support the conclusion?
2. What is the biggest way this conclusion could be misleading?
3. Does the next calculation reduce the largest remaining uncertainty?

Examples:

- Deacylation attempt040 was justified because attempt039 fixed the br206 -> r210
  overlap problem but still needed longer production for stable state free
  energies.
- The h119 -> h104 weak overlap was not patched immediately because it is not at
  the barrier, while acylation still lacked a long production comparable to
  deacylation.
- Acylation attempt047 was started because it directly addresses the remaining
  largest uncertainty: the acylation PMF was based on shorter 20000-step windows.

## Current Barrier Summary

As of 2026-07-06 20:35 CST:

```text
Deacylation:
  barrier = 19.061 +/- 0.073 kcal/mol
  status  = stable long-production result

Acylation:
  barrier = pending attempt047 long-production PMF
  status  = running
```

The final paper comparison should wait until acylation attempt047 finishes and
is reanalyzed.

