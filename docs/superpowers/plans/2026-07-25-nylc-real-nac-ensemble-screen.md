# NylC Real-NAC Ensemble Screen Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build, run, and independently audit a 12-conformation × 3-velocity-seed NylC-C18 M1 ensemble screen on SCNet, then advance the top six conformations to a total of 1 ns fully unrestrained validation and emit a defensible QM/MM eligibility manifest.

**Architecture:** Extend the existing `workflows/nylc_l4_nac_to_l2_rebalance_20260723` package with deterministic event selection, parameterized M1 build/equilibration array jobs, replica-level analysis, conformation-level ranking, and an independent final auditor. Legacy single-conformation jobs and outputs remain immutable; every new candidate, replica, and attempt receives a new directory and manifest.

**Tech Stack:** Python 3 standard library + NumPy, GROMACS 2022.1 DCU for MD, GROMACS/CP2K 2023.1 double precision for flexible-water EM where required, Bash/Slurm, `pytest`/unittest, JSON/TSV/JSONL.

## Global Constraints

- Repository: `bianhuiting918/git-hub-chatgpt-pro`.
- Branch: `codex/nylc-l4-nac-to-l2-rebalance`.
- Approved design: `docs/superpowers/specs/2026-07-25-nylc-real-nac-ensemble-screen-design.md`.
- Workflow root in Git: `workflows/nylc_l4_nac_to_l2_rebalance_20260723`.
- SCNet task root: `/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723`.
- Compute, dependencies, caches, trajectories, and large results run or remain on SCNet; never run production MD on Dell.
- GitHub stores only scripts, MDP text, RUNBOOK, small manifests/audits, and checksums; never commit GRO/XTC/TRR/TPR/CPT/EDR, large topology files, credentials, private keys, passwords, or tokens.
- Source M0 trajectory and legacy jobs are read-only. Jobs `61801874` and `61803121` remain cancelled and must not be resumed.
- Scientific NAC gate: distance ≤ 0.35 nm and angle 95–115°.
- NylC gate group is residues 261–266; Thr267 is excluded.
- Microstate is M1: NαH2 / Asp306H / Asp308− / Thr267OH.
- The audited PA66-L2 ITP SHA256 is `b0e753c60fd4b71c282d21cc6106a15e73d91d12a20d80e92dd01516162eb301`.
- Technical restrained stages never count as scientific NAC evidence.
- Stage A is 12 conformations × 3 velocity seeds × 100 ps fully unrestrained NPT.
- Stage B keeps the same replica identities and extends the top six from 100 ps to 1 ns total free time.
- Stage B primary analysis window is 100–1000 ps.
- Conformation PASS requires at least 2/3 replicas with NAC after 100 ps, pooled NAC occupancy ≥ 1%, at least one continuous NAC event ≥ 4 ps, bound substrate, numerical stability, stable thermodynamics, and a NAC cluster reproduced in at least two replicas.
- At most three conformations become QM/MM candidates.
- Fixed-topology MM proton geometry is preorganization evidence only.
- Each remote action appends `run_history.tsv` and `run_history.jsonl` without secrets.
- Per-candidate failures are retained and never block independent candidates.

---

## File Map

### Create

- `workflows/nylc_l4_nac_to_l2_rebalance_20260723/manifests/nylc_C18_m1_real_nac_ensemble.authority.json` — immutable source paths, atom identities, thresholds, seeds, and source hashes.
- `workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/select_nylc_m1_nac_ensemble.py` — event grouping, GRO parsing, Kabsch alignment, local RMSD clustering, deterministic 12-way selection.
- `workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/prepare_nylc_m1_nac_ensemble.py` — candidate extraction/build orchestration and immutable manifests.
- `workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/render_nylc_m1_replica_mdp.py` — renders the three fixed velocity seeds without mutating templates.
- `workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/audit_nylc_m1_ensemble_replica.py` — Stage A/B per-replica numerical, NAC, bound-state, gate, thermo, and proton-network audit.
- `workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/rank_nylc_m1_ensemble.py` — conformation-level 3-replica aggregation, top-six selection, Stage B final gate, top-three QM/MM eligibility.
- `workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/independent_audit_nylc_m1_ensemble.py` — recomputes decisions from primitive XVG/JSON inputs and rejects stale/superseded outputs.
- `workflows/nylc_l4_nac_to_l2_rebalance_20260723/mdp/npt300free_m1_stageA.mdp` — 100 ps fully unrestrained NPT.
- `workflows/nylc_l4_nac_to_l2_rebalance_20260723/mdp/npt300free_m1_stageB_extend.mdp` — 900 ps continuation, total free time 1 ns.
- `workflows/nylc_l4_nac_to_l2_rebalance_20260723/slurm/run_nylc_m1_ensemble_select.sbatch` — CPU selection/extraction job.
- `workflows/nylc_l4_nac_to_l2_rebalance_20260723/slurm/run_nylc_m1_ensemble_build_array.sbatch` — 12-way CPU M1 build/grompp preflight.
- `workflows/nylc_l4_nac_to_l2_rebalance_20260723/slurm/run_nylc_m1_ensemble_em_array.sbatch` — 12-way double-precision flexible-water EM.
- `workflows/nylc_l4_nac_to_l2_rebalance_20260723/slurm/run_nylc_m1_ensemble_stageA_array.sbatch` — 36-way restrained preparation followed by 100 ps free MD.
- `workflows/nylc_l4_nac_to_l2_rebalance_20260723/slurm/run_nylc_m1_ensemble_rank_stageA.sbatch` — independent Stage A aggregation/top-six manifest.
- `workflows/nylc_l4_nac_to_l2_rebalance_20260723/slurm/run_nylc_m1_ensemble_stageB_array.sbatch` — 18-way 900 ps continuation.
- `workflows/nylc_l4_nac_to_l2_rebalance_20260723/slurm/run_nylc_m1_ensemble_final_audit.sbatch` — final independent science audit and QM/MM eligibility.
- `workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_select_nylc_m1_nac_ensemble.py`.
- `workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_prepare_nylc_m1_nac_ensemble.py`.
- `workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_render_nylc_m1_replica_mdp.py`.
- `workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_audit_nylc_m1_ensemble_replica.py`.
- `workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_rank_nylc_m1_ensemble.py`.
- `workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_nylc_m1_ensemble_slurm_contract.py`.
- `workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_independent_audit_nylc_m1_ensemble.py`.

### Modify

- `workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/analyze_nac_series.py` — add complete event list and analysis-window support while retaining schema-v1 compatibility.
- `workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_analyze_nac_series.py` — regression tests for events and time window.
- `workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_mdp_contract.py` — assert Stage A/B lengths and zero restraints.
- `workflows/nylc_l4_nac_to_l2_rebalance_20260723/RUNBOOK.md` — append the approved ensemble procedure, queue-aware submission commands, recovery, and interpretation.
- `.gitignore` — reject remote MD binaries/trajectories and local copies of candidate coordinate artifacts if current patterns do not already cover them.

## Frozen Interfaces

```python
def group_nac_events(rows: list[dict], sample_interval_ps: float = 2.0) -> list[dict]: ...
def read_gro(path: pathlib.Path) -> tuple[list[dict], tuple[float, ...]]: ...
def kabsch_rmsd(reference: "np.ndarray", mobile: "np.ndarray") -> float: ...
def select_twelve(events: list[dict], rmsd_nm: "np.ndarray") -> list[dict]: ...
def build_candidate(source_gro: pathlib.Path, candidate_dir: pathlib.Path, authority: dict) -> dict: ...
def audit_replica(run_root: pathlib.Path, manifest: dict, window: tuple[float, float]) -> dict: ...
def rank_stage_a(candidate_audits: list[dict]) -> dict: ...
def gate_stage_b(candidate_audits: list[dict]) -> dict: ...
```

Candidate IDs are `nac_evtNN_timeTTTTps`; replica IDs are `seed26711`, `seed26723`, and `seed26737`. The same three seeds are used across conformations for paired comparison, but each replica has independent velocities.

---

### Task 1: Extend strict NAC series analysis

**Files:**
- Modify: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/analyze_nac_series.py`
- Modify: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_analyze_nac_series.py`

**Interfaces:**
- Consumes: aligned distance/angle/potential XVG rows.
- Produces: schema version 2 JSON with `nac_events`, `analysis_window_ps`, occupancy, longest event, and optional energy.

- [ ] **Step 1: Write failing event/window tests**

Add tests using 2 ps samples that assert:

```python
audit = audit_series(
    [(0, .30), (2, .31), (4, .40), (6, .32), (8, .33)],
    [(0, 100), (2, 101), (4, 105), (6, 106), (8, 107)],
    .35, 95, 115,
    analysis_start_ps=2.0,
    analysis_end_ps=8.0,
)
assert audit["schema_version"] == 2
assert audit["analysis_window_ps"] == [2.0, 8.0]
assert [(e["start_ps"], e["end_ps"], e["frame_count"]) for e in audit["nac_events"]] == [
    (2.0, 2.0, 1),
    (6.0, 8.0, 2),
]
assert audit["nac_occupancy"] == 3 / 4
```

Also assert a discontinuity larger than 2.000001 ps starts a new event.

- [ ] **Step 2: Run the focused tests and confirm RED**

Run:

```bash
python -m pytest workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_analyze_nac_series.py -q
```

Expected: new tests fail because `audit_series` lacks window arguments and `nac_events`.

- [ ] **Step 3: Implement windowing and explicit event records**

Update the signature exactly to:

```python
def audit_series(
    distance_rows,
    angle_rows,
    distance_max,
    angle_min,
    angle_max,
    potential_rows=None,
    analysis_start_ps=None,
    analysis_end_ps=None,
    sample_interval_ps=2.0,
):
```

Filter only after verifying full-series alignment. Each event record must contain `event_id`, `start_ps`, `end_ps`, `duration_ps`, `frame_count`, `member_times_ps`, mean distance/angle, and minimum potential when present. Add CLI arguments `--analysis-start-ps`, `--analysis-end-ps`, and `--sample-interval-ps`.

- [ ] **Step 4: Run regression and focused tests**

Expected:

```text
test_analyze_nac_series.py: all passed
```

- [ ] **Step 5: Commit**

```bash
git add workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/analyze_nac_series.py workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_analyze_nac_series.py
git commit -m "feat: expose NAC residence events and windows"
```

---

### Task 2: Freeze source authority and select 12 real NAC events

**Files:**
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/manifests/nylc_C18_m1_real_nac_ensemble.authority.json`
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/select_nylc_m1_nac_ensemble.py`
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_select_nylc_m1_nac_ensemble.py`
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/slurm/run_nylc_m1_ensemble_select.sbatch`

**Interfaces:**
- Consumes: job `61780674` TPR/XTC/EDR plus selection XVGs from build job `61801089`.
- Produces: `ensemble_selection.json`, 12 candidate GROs on SCNet, per-frame SHA256, local RMSD matrix, and selection category.

- [ ] **Step 1: Write the immutable authority manifest**

Use these exact scientific fields:

```json
{
  "schema_version": 1,
  "candidate_id": "nylc_C18_trueT267_freeGS",
  "source_job_id": "61780674",
  "selection_audit_job_id": "61801089",
  "source_run_root": "/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/candidates/nylc_C18_trueT267_freeGS/runs/mm_proton_geometry_ext19ns_job_61780674",
  "selection_root": "/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/candidates/nylc_C18_trueT267_freeGS/microstates/M1_nalphaH2_Asp306H/build_job_61801089/selection",
  "source_tpr": "run.tpr",
  "source_xtc": "run.part0002.xtc",
  "source_edr": "run.part0002.edr",
  "analysis_window_ps": [1000.0, 12660.0],
  "sample_interval_ps": 2.0,
  "nac_gate": {"distance_max_nm": 0.35, "angle_min_deg": 95.0, "angle_max_deg": 115.0},
  "central_gate": {"distance_max_nm": 0.32, "angle_min_deg": 100.0, "angle_max_deg": 110.0},
  "expected_frame_count": 5831,
  "expected_nac_frame_count": 85,
  "expected_event_count": 72,
  "old_control_time_ps": 1462.0,
  "reactive_global_index1_m0": {"thr267_og1": 8961, "l2_c": 10287, "l2_o": 10288, "l2_n": 10289},
  "gate_residues": [261, 262, 263, 264, 265, 266],
  "excluded_gate_residues": [267],
  "velocity_seeds": [26711, 26723, 26737]
}
```

The selection job writes actual SHA256 values for source TPR/XTC/EDR and refuses to run if paths or counts disagree.

- [ ] **Step 2: Write failing deterministic-selection tests**

Synthetic tests must cover:

```python
events = group_nac_events(rows, sample_interval_ps=2.0)
assert len(events) == 72
selected = select_twelve(events, rmsd_nm)
assert len(selected) == 12
assert len({x["event_id"] for x in selected}) == 12
assert sum(x["category"] == "three_frame_medoid" for x in selected) == 4
assert sum(x["category"] == "two_frame_medoid" for x in selected) == 4
assert sum(x["category"] == "central_low_energy" for x in selected) == 3
assert [x for x in selected if x["time_ps"] == 1462.0][0]["category"] == "legacy_control"
assert select_twelve(events, rmsd_nm) == select_twelve(events, rmsd_nm)
```

Include a test proving a multi-frame event cannot occupy both the persistence and central-energy categories.

- [ ] **Step 3: Implement GRO parsing, Kabsch RMSD, and selection**

The script must:

1. read the 85 NAC frame times;
2. use `gmx trjconv -dump <time>` to create temporary SCNet GROs;
3. construct one fixed heavy-atom selection: PA66-L2 + Thr267/Asp306/Asp308 + the union of protein residues having any heavy atom within 0.6 nm of PA66-L2 in any NAC frame;
4. align by all protein backbone N/CA/C atoms in the active chain;
5. calculate a symmetric local RMSD matrix using NumPy Kabsch;
6. choose one representative per event;
7. apply the 4 + 4 + 3 + 1 category allocation;
8. break ties by category priority, smaller mean local RMSD, lower within-event potential, smaller distance from (0.30 nm, 105°), then earlier time;
9. emit candidate IDs and SHA256 values.

Temporary GROs remain on SCNet and are never copied into Git.

- [ ] **Step 4: Add a CPU-only selection Slurm job**

Resource contract:

```bash
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem-per-cpu=2500M
#SBATCH -t 02:00:00
```

The job must create a new `ensemble/selection_job_$SLURM_JOB_ID` directory, append START/PASS/FAIL to both run histories with `flock`, and refuse to overwrite.

- [ ] **Step 5: Run tests and a read-only dry run**

Run:

```bash
python -m pytest workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_analyze_nac_series.py workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_select_nylc_m1_nac_ensemble.py -q
python workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/select_nylc_m1_nac_ensemble.py --authority workflows/nylc_l4_nac_to_l2_rebalance_20260723/manifests/nylc_C18_m1_real_nac_ensemble.authority.json --validate-only
```

Expected: tests pass; validate-only prints the three source paths and does not create candidate directories.

- [ ] **Step 6: Commit**

```bash
git add workflows/nylc_l4_nac_to_l2_rebalance_20260723/manifests/nylc_C18_m1_real_nac_ensemble.authority.json workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/select_nylc_m1_nac_ensemble.py workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_select_nylc_m1_nac_ensemble.py workflows/nylc_l4_nac_to_l2_rebalance_20260723/slurm/run_nylc_m1_ensemble_select.sbatch
git commit -m "feat: select twelve real NylC NAC events"
```

---

### Task 3: Parameterize M1 construction and independent preflight

**Files:**
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/prepare_nylc_m1_nac_ensemble.py`
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_prepare_nylc_m1_nac_ensemble.py`
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/slurm/run_nylc_m1_ensemble_build_array.sbatch`

**Interfaces:**
- Consumes: one selected M0 PA66-L2 GRO and the existing audited build assets.
- Produces: one immutable M1 build directory with `system_M1.gro`, topology links/copies on SCNet, `build_audit.json`, and `PASS.json` or exact failure JSON.

- [ ] **Step 1: Write failing build-contract tests**

Test a synthetic chain and manifest:

```python
audit = build_candidate(source_gro, candidate_dir, authority)
assert audit["microstate"] == {
    "Thr267_Nalpha": "NH2_neutral_proxy",
    "Thr267_Ogamma": "OH",
    "Asp306": "ASH_HD2_neutral",
    "Asp308": "ASP_minus",
}
assert audit["charge_delta_e"] == 0.0
assert audit["reaction_geometry"]["joint_nac"] is True
assert audit["gate_definition"] == "NylC residues 261-266; Thr267 excluded"
assert audit["source_gro_sha256"] == sha256(source_gro)
```

Also assert existing target directories raise `FileExistsError`, wrong Thr267/global identities raise `ValueError`, and a source outside the selection manifest is rejected.

- [ ] **Step 2: Extract reusable logic from the legacy build without changing legacy behavior**

Import and reuse:

```python
from prepare_nylc_nalpha_h2_ash306 import (
    drop_gro_atom,
    export_chain_pdb,
    extract_molecule_itp,
    replace_gro_slice,
    rewrite_chain_itp,
)
```

Parameterize candidate ID, source GRO, selected time, selected event, output root, and source SHA. Preserve the audited active chain range 8949–10272 and M1 post-build reactive global indices Thr267 OG1=8960, L2 C/O/N=10287/10288/10289 only after identity validation.

- [ ] **Step 3: Implement all preflight gates**

The build array must verify:

- source candidate is one of exactly 12 selected events;
- source and output hashes;
- M1 active-chain atom count 1324;
- active-chain charge −4.0 e before and after M1 rewrite;
- Thr267 has N/H1/H2/OG1/HG1 and no H3;
- Asp306 is ASH with HD2;
- Asp308 remains ASP without HD2;
- full-system atom count unchanged;
- unchanged atoms outside active chain;
- ordinary active-chain heavy atoms preserved within 0.0011 nm;
- source and built attack geometry remain NAC;
- minimum chain–rest distance ≥ 0.08 nm;
- minimum heavy-atom contact ≥ 0.18 nm;
- `grompp -maxwarn 0` passes.

Write distinct failure labels from the approved design.

- [ ] **Step 4: Add the 12-way CPU build array**

Use `#SBATCH --array=0-11%4`, map array index through `ensemble_selection.json`, and write to:

```text
ensemble/candidates/<candidate_id>/build_job_<jobid>_<arrayid>/
```

No symlink may point into a superseded or cancelled run directory except the read-only old audited build assets explicitly recorded in `build_audit.json`.

- [ ] **Step 5: Run focused tests**

```bash
python -m pytest workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_prepare_nylc_nalpha_h2_ash306.py workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_prepare_nylc_m1_nac_ensemble.py -q
```

Expected: all old and new build tests pass.

- [ ] **Step 6: Commit**

```bash
git add workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/prepare_nylc_m1_nac_ensemble.py workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_prepare_nylc_m1_nac_ensemble.py workflows/nylc_l4_nac_to_l2_rebalance_20260723/slurm/run_nylc_m1_ensemble_build_array.sbatch
git commit -m "feat: parameterize M1 ensemble builds"
```

---

### Task 4: Render replica seeds and implement EM/release arrays

**Files:**
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/render_nylc_m1_replica_mdp.py`
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_render_nylc_m1_replica_mdp.py`
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/mdp/npt300free_m1_stageA.mdp`
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/mdp/npt300free_m1_stageB_extend.mdp`
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/slurm/run_nylc_m1_ensemble_em_array.sbatch`
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/slurm/run_nylc_m1_ensemble_stageA_array.sbatch`
- Modify: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_mdp_contract.py`

**Interfaces:**
- Consumes: 12 PASS builds and seeds 26711/26723/26737.
- Produces: 12 EM PASS/FAIL records and 36 fully released Stage A replicas.

- [ ] **Step 1: Write failing MDP rendering and contract tests**

```python
rendered = render_mdp(template, {"gen-seed": "26711"})
assert parse(rendered)["gen-seed"] == "26711"
assert parse(rendered)["define"] == "-DPOSRES -DPOSRES_L2_1000"
```

Assert:

```python
stage_a = parse_mdp("npt300free_m1_stageA.mdp")
assert float(stage_a["dt"]) * int(stage_a["nsteps"]) == 100.0
assert "define" not in stage_a
assert stage_a["gen-vel"].lower() == "no"

stage_b = parse_mdp("npt300free_m1_stageB_extend.mdp")
assert float(stage_b["dt"]) * int(stage_b["nsteps"]) == 900.0
assert stage_b["continuation"].lower() == "yes"
assert "define" not in stage_b
```

- [ ] **Step 2: Implement strict MDP rendering**

Only keys already present in a template may be replaced. Duplicate keys, missing keys, non-approved seeds, or attempts to inject `define` into a free MDP must fail.

CLI:

```bash
python scripts/render_nylc_m1_replica_mdp.py --template mdp/nvt50_m1.mdp --output run/nvt50.mdp --set gen-seed=26711
```

- [ ] **Step 3: Add the 12-way EM array**

Use the existing double-precision flexible-water EM contract:

```text
GROMACS 2023.1 double precision
integrator=cg
-DPOSRES -DPOSRES_L2_1000 -DFLEXIBLE
Fmax ≤ 500 kJ mol−1 nm−1
no LINCS/SETTLE/NaN/FATAL
```

A failed EM is `FAIL_TECHNICAL_EM`; it does not cancel other array elements.

- [ ] **Step 4: Add the 36-way Stage A array**

Index mapping:

```python
candidate_index = SLURM_ARRAY_TASK_ID // 3
seed_index = SLURM_ARRAY_TASK_ID % 3
seed = [26711, 26723, 26737][seed_index]
```

Resource request before queue-specific adjustment:

```bash
#SBATCH -p xahdnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --gres=dcu:1
#SBATCH --mem-per-cpu=2500M
#SBATCH -t 06:00:00
#SBATCH --array=0-35%8
```

Run `nvt50 → nvt150 → nvt300 → npt300r → npt300rel → npt300free_m1_stageA`. Only `nvt50` receives the rendered seed. Use `-r system_M1.gro` only for restrained stages; the Stage A `.tpr` must contain no position restraints or distance restraints.

Before Stage A mdrun, dump/check the TPR and write `free_tpr_contract.json` containing:

```json
{"position_restraints": 0, "distance_restraints": 0, "scientific_window": "fully_unrestrained_NPT_100ps"}
```

- [ ] **Step 5: Verify tests**

```bash
python -m pytest workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_render_nylc_m1_replica_mdp.py workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_mdp_contract.py -q
```

- [ ] **Step 6: Commit**

```bash
git add workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/render_nylc_m1_replica_mdp.py workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_render_nylc_m1_replica_mdp.py workflows/nylc_l4_nac_to_l2_rebalance_20260723/mdp/npt300free_m1_stageA.mdp workflows/nylc_l4_nac_to_l2_rebalance_20260723/mdp/npt300free_m1_stageB_extend.mdp workflows/nylc_l4_nac_to_l2_rebalance_20260723/slurm/run_nylc_m1_ensemble_em_array.sbatch workflows/nylc_l4_nac_to_l2_rebalance_20260723/slurm/run_nylc_m1_ensemble_stageA_array.sbatch workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_mdp_contract.py
git commit -m "feat: run seeded M1 ensemble release and stage A"
```

---

### Task 5: Audit each free replica and rank Stage A

**Files:**
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/audit_nylc_m1_ensemble_replica.py`
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/rank_nylc_m1_ensemble.py`
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_audit_nylc_m1_ensemble_replica.py`
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_rank_nylc_m1_ensemble.py`
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/slurm/run_nylc_m1_ensemble_rank_stageA.sbatch`

**Interfaces:**
- Consumes: primitive XVGs/logs/GROs for 36 replicas.
- Produces: per-replica audit JSON, per-conformation Stage A summaries, ranked top-six manifest.

- [ ] **Step 1: Write failing replica-audit tests**

Fixtures must cover:

```python
assert audit_replica(pass_root, manifest, (0.0, 100.0))["technical_status"] == "PASS"
assert audit_replica(pass_root, manifest, (0.0, 100.0))["nac"]["nac_frame_count"] > 0
assert audit_replica(lincs_root, manifest, (0.0, 100.0))["scientific_status"] == "NOT_EVALUATED_NUMERICAL_FAIL"
assert audit_replica(unbound_root, manifest, (0.0, 100.0))["scientific_status"] == "FAIL_UNBOUND"
```

The audit must report the actual denominator, sampling interval, longest event, gate opening, thermo mean/stdev/min/max, minimum heavy contact, pocket contacts, and proton preorganization channels.

- [ ] **Step 2: Implement bound-state and proton-network metrics**

Use fixed, pre-registered geometric annotations:

- pocket retained if PA66-L2 has at least 3 heavy-atom contacts ≤ 0.45 nm to the source 0.6 nm pocket-residue union and ligand-pocket COM distance ≤ 1.2 nm;
- severe clash if any ligand–protein heavy-atom distance < 0.18 nm;
- ThrOH donor channel annotated by O-donor/acceptor distance ≤ 0.35 nm and donor–H–acceptor angle ≥ 135°;
- first-shell water oxygen within 0.35 nm of Thr OG1;
- a two-edge water relay requires both hydrogen-bond edges to meet the same distance/angle definition in the same frame.

These are preorganization annotations, not proton-transfer claims.

- [ ] **Step 3: Write failing ranking tests**

Create six synthetic conformations with three replicas each and assert lexicographic ranking by:

```python
key = (
    replica_count_with_nac_after_20ps,
    pooled_nac_occupancy_after_20ps,
    longest_continuous_nac_ps,
    reproduced_cluster_count,
    -unbound_replica_count,
    -technical_failure_count,
)
```

The test must prove selection returns six distinct conformations, never six replicas, and excludes candidates whose three replicas all lack NAC after 20 ps.

- [ ] **Step 4: Implement Stage A ranking**

Write:

- `stageA_replica_audits.json`;
- `stageA_conformation_ranking.tsv`;
- `stageA_top6_manifest.json`;
- `stageA_failures.json`.

If fewer than six candidates remain scientifically evaluable, select all available and mark `NOT_EVALUATED_STAGEB_FEWER_THAN_SIX`; do not fill the set with failed candidates.

- [ ] **Step 5: Add the CPU ranking job**

It must depend on the complete Stage A array with `afterany`, so failed elements are still audited. The job exits nonzero only for audit corruption, missing primitive inputs without explicit failure records, or inconsistent candidate universes—not because some candidates scientifically fail.

- [ ] **Step 6: Run tests and commit**

```bash
python -m pytest workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_audit_nylc_m1_ensemble_replica.py workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_rank_nylc_m1_ensemble.py -q
git add workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/audit_nylc_m1_ensemble_replica.py workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/rank_nylc_m1_ensemble.py workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_audit_nylc_m1_ensemble_replica.py workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_rank_nylc_m1_ensemble.py workflows/nylc_l4_nac_to_l2_rebalance_20260723/slurm/run_nylc_m1_ensemble_rank_stageA.sbatch
git commit -m "feat: audit and rank NylC ensemble stage A"
```

---

### Task 6: Extend top six replicas to 1 ns and enforce the science gate

**Files:**
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/slurm/run_nylc_m1_ensemble_stageB_array.sbatch`
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_nylc_m1_ensemble_slurm_contract.py`
- Modify: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/rank_nylc_m1_ensemble.py`
- Modify: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_rank_nylc_m1_ensemble.py`

**Interfaces:**
- Consumes: exact 18 Stage A checkpoints selected by `stageA_top6_manifest.json`.
- Produces: 18 trajectories totaling 1 ns free time and conformation-level Stage B PASS/FAIL.

- [ ] **Step 1: Write failing Slurm dependency/continuation tests**

Assert the script:

```python
assert "--array=0-17%6" in text
assert "stageA_top6_manifest.json" in text
assert "run.cpt" in text
assert "npt300free_m1_stageB_extend.mdp" in text
assert "gen-vel" not in stage_b_override_text
assert "afterok" not in documented_stage_a_to_audit_dependency
```

Also assert no Stage B output path contains legacy job IDs `61801874` or `61803121`.

- [ ] **Step 2: Implement 18-way continuation**

Use each exact Stage A `.gro/.cpt`; do not regenerate velocities. The Stage B job writes a new sibling `stageB_attempt_<jobid>_<arrayid>` and records its Stage A parent SHA256/checkpoint.

- [ ] **Step 3: Extend conformation gate tests**

Exact Stage B PASS fixture:

```python
decision = gate_stage_b(three_replica_audits)
assert decision["replicas_with_nac_after_100ps"] >= 2
assert decision["pooled_nac_occupancy_100_1000ps"] >= 0.01
assert decision["longest_continuous_nac_ps"] >= 4.0
assert decision["bound_replica_count"] == 3
assert decision["technical_failure_count"] == 0
assert decision["cross_replica_cluster_reproduced"] is True
assert decision["scientific_status"] == "PASS_UNRESTRAINED_M1_ENSEMBLE_NAC"
```

Add one separate failing fixture for each gate and preserve its exact reason list.

- [ ] **Step 4: Implement the 100–1000 ps aggregation**

All occupancies and proton-route comparisons use only frames with `100.0 <= time_ps <= 1000.0`. A continuously maintained NAC qualifies as occurrence; leaving and returning also qualifies. Pooled occupancy denominator is the sum of valid frames across all three replicas, never a mean of percentages.

- [ ] **Step 5: Run tests and commit**

```bash
python -m pytest workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_nylc_m1_ensemble_slurm_contract.py workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_rank_nylc_m1_ensemble.py -q
git add workflows/nylc_l4_nac_to_l2_rebalance_20260723/slurm/run_nylc_m1_ensemble_stageB_array.sbatch workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_nylc_m1_ensemble_slurm_contract.py workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/rank_nylc_m1_ensemble.py workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_rank_nylc_m1_ensemble.py
git commit -m "feat: validate top six M1 conformations to one ns"
```

---

### Task 7: Add independent final audit and QM/MM eligibility manifest

**Files:**
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/independent_audit_nylc_m1_ensemble.py`
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_independent_audit_nylc_m1_ensemble.py`
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/slurm/run_nylc_m1_ensemble_final_audit.sbatch`

**Interfaces:**
- Consumes: frozen universe manifest and primitive Stage A/B metrics.
- Produces: `FINAL_AUDIT.json`, `qmmm_eligible_top3.json`, failure ledger, and compact TSV.

- [ ] **Step 1: Write failing independent-audit tests**

Test rejection of:

- duplicate candidate IDs;
- duplicate seeds;
- missing one of 36 Stage A slots without a failure record;
- a Stage B candidate absent from Stage A top six;
- a PASS whose primitive occupancy recomputes below 1%;
- Thr267 included in the gate group;
- a free TPR declaring restraints;
- legacy/cancelled output paths;
- `COMPLETE.json` without the science gate.

- [ ] **Step 2: Implement independent recomputation**

Do not import `rank_nylc_m1_ensemble.py`. Reimplement the small gate calculation from primitive arrays so agreement is meaningful. Emit:

```json
{
  "technical_status": "PASS",
  "scientific_status": "PASS_UNRESTRAINED_M1_ENSEMBLE_NAC",
  "candidate_universe": 12,
  "stageA_expected_replicas": 36,
  "stageB_expected_replicas": 18,
  "scientific_pass_count": 0,
  "qmmm_eligible_count": 0,
  "qmmm_eligibility": []
}
```

Counts are populated from data; maximum `qmmm_eligible_count` is 3.

- [ ] **Step 3: Define recurrent GS-cluster output**

For each eligible conformation, include the cross-replica NAC cluster medoid, source replica/time, cluster membership by replica, source coordinate checksum on SCNet, reaction geometry, gate opening, proton-route annotation, and the explicit caveat that the medoid is not yet a QM/MM-optimized GS.

- [ ] **Step 4: Add final audit job**

Use CPU-only resources and `afterany` dependency on Stage B. Technical corruption exits nonzero. Scientific failure exits zero after writing a complete failure ledger, because “no eligible NAC” is a scientific result, not a scheduler failure.

- [ ] **Step 5: Run tests and commit**

```bash
python -m pytest workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_independent_audit_nylc_m1_ensemble.py -q
git add workflows/nylc_l4_nac_to_l2_rebalance_20260723/scripts/independent_audit_nylc_m1_ensemble.py workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_independent_audit_nylc_m1_ensemble.py workflows/nylc_l4_nac_to_l2_rebalance_20260723/slurm/run_nylc_m1_ensemble_final_audit.sbatch
git commit -m "feat: independently audit M1 ensemble eligibility"
```

---

### Task 8: Document, secure, and regression-test the workflow

**Files:**
- Modify: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/RUNBOOK.md`
- Modify: `.gitignore`
- Create: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_nylc_m1_ensemble_repository_policy.py`

**Interfaces:**
- Consumes: all new scripts/jobs.
- Produces: exact operator commands, recovery rules, and repository policy checks.

- [ ] **Step 1: Write repository-policy tests**

Scan tracked workflow paths and fail on suffixes:

```python
FORBIDDEN = {".xtc", ".trr", ".tpr", ".cpt", ".edr", ".gro"}
```

Also scan text for PEM private-key headers, `ghp_`, passwords assigned in shell, and raw SSH private key content. Allow literal words such as “password” only in policy prose.

- [ ] **Step 2: Append RUNBOOK sections**

Document exact commands for:

- read-only `squeue` and `sacct`;
- checkout verification at the approved Git commit;
- full test suite;
- selection → build → EM → Stage A → Stage A rank → Stage B → final audit;
- dependency job IDs;
- per-array failure recovery into new attempt directories;
- cancellation without deleting outputs;
- copying only small JSON/TSV summaries back to GitHub;
- scientific interpretation boundaries.

- [ ] **Step 3: Run the complete local-text test suite**

On SCNet checkout:

```bash
python -m pytest workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests -q
```

Expected: all tests pass; record exact test count.

- [ ] **Step 4: Commit**

```bash
git add .gitignore workflows/nylc_l4_nac_to_l2_rebalance_20260723/RUNBOOK.md workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_nylc_m1_ensemble_repository_policy.py
git commit -m "docs: add NylC ensemble runbook and repository guards"
```

---

### Task 9: Deploy and execute on SCNet

**Files:**
- Remote-only generated: `$TASK_ROOT/ensemble/*`, `run_history.tsv`, `run_history.jsonl`.
- GitHub update after each completed checkpoint: scripts/docs/small audit JSON only.

**Interfaces:**
- Consumes: reviewed GitHub commits from Tasks 1–8.
- Produces: exact Slurm job IDs, Stage A/B audit artifacts, and final science status.

- [ ] **Step 1: Read-only queue and repository audit**

Run one narrow SSH session:

```bash
ssh 210.73.40.29 'squeue -u "$USER" -o "%.18i %.12P %.32j %.8T %.10M %.10l %.6D %R"; sacct -X -S 2026-07-25 -u "$USER" --format=JobID,JobName%32,Partition,State,Elapsed,ExitCode -n | tail -80'
```

Then confirm no task job is running and the cancelled jobs remain cancelled. Do not cancel unrelated jobs.

- [ ] **Step 2: Update the SCNet checkout safely**

At:

```text
/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/repo
```

verify branch, status, and remote before fast-forwarding. If SCNet still lacks GitHub authentication, transfer only the reviewed text patch through the secure existing SSH session and record both GitHub commit SHA and remote file SHA256; do not claim the SCNet-only commit is pushed.

- [ ] **Step 3: Run all tests before submission**

```bash
cd /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/repo
python -m pytest workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests -q
```

No production submission if the new tests fail.

- [ ] **Step 4: Submit selection and wait for audited PASS**

```bash
jid_select=$(sbatch --parsable workflows/nylc_l4_nac_to_l2_rebalance_20260723/slurm/run_nylc_m1_ensemble_select.sbatch)
```

Record `jid_select`. Submit build only if `ensemble_selection.json` contains exactly 12 unique events and the source count/hash audit passes.

- [ ] **Step 5: Submit build and EM arrays**

```bash
jid_build=$(sbatch --parsable --dependency=afterok:$jid_select workflows/nylc_l4_nac_to_l2_rebalance_20260723/slurm/run_nylc_m1_ensemble_build_array.sbatch)
jid_em=$(sbatch --parsable --dependency=afterany:$jid_build workflows/nylc_l4_nac_to_l2_rebalance_20260723/slurm/run_nylc_m1_ensemble_em_array.sbatch)
```

The EM array self-skips candidates without build PASS and records `NOT_EVALUATED_BUILD_FAIL`.

- [ ] **Step 6: Submit Stage A and its independent ranker**

After checking current DCU queue, keep or lower the `%8` throttle; do not raise it without available capacity.

```bash
jid_a=$(sbatch --parsable --dependency=afterany:$jid_em workflows/nylc_l4_nac_to_l2_rebalance_20260723/slurm/run_nylc_m1_ensemble_stageA_array.sbatch)
jid_rank_a=$(sbatch --parsable --dependency=afterany:$jid_a workflows/nylc_l4_nac_to_l2_rebalance_20260723/slurm/run_nylc_m1_ensemble_rank_stageA.sbatch)
```

- [ ] **Step 7: Audit Stage A before Stage B submission**

Confirm denominators:

```text
candidate universe = 12
expected Stage A replica slots = 36
completed technical replicas + explicit failed/not-evaluated replicas = 36
top-six conformation IDs are unique
```

Do not advance a candidate merely because Slurm exited 0.

- [ ] **Step 8: Submit Stage B and final audit**

```bash
jid_b=$(sbatch --parsable --dependency=afterok:$jid_rank_a workflows/nylc_l4_nac_to_l2_rebalance_20260723/slurm/run_nylc_m1_ensemble_stageB_array.sbatch)
jid_final=$(sbatch --parsable --dependency=afterany:$jid_b workflows/nylc_l4_nac_to_l2_rebalance_20260723/slurm/run_nylc_m1_ensemble_final_audit.sbatch)
```

Report all exact job IDs and expected array sizes.

- [ ] **Step 9: Monitor without overwriting or premature claims**

Each monitoring cycle:

1. read `squeue`/new `sacct` states;
2. inspect only new task logs and small JSON summaries;
3. append run history;
4. repair only technical failures in new attempt directories;
5. keep scientific failures unchanged;
6. continue independent candidates;
7. distinguish `PASS_TECHNICAL`, `PASS_UNRESTRAINED_M1_ENSEMBLE_NAC`, and `NOT_EVALUATED_*`.

- [ ] **Step 10: Publish compact checkpoint artifacts**

At selection, Stage A, and final audit checkpoints, commit only:

- scripts/MDPs/tests/RUNBOOK changes;
- candidate/time/source-hash manifest;
- compact ranking TSV/JSON;
- independent audit JSON;
- run-history excerpts without secrets.

Never publish coordinate or trajectory files.

---

### Task 10: Hand off only eligible recurrent NAC clusters to QM/MM planning

**Files:**
- Remote/GitHub small artifact: `workflows/nylc_l4_nac_to_l2_rebalance_20260723/audit/nylc_C18_m1_ensemble_qmmm_eligibility.json`
- Next design/plan, created only after results: Step1/Step2/PMF documents.

**Interfaces:**
- Consumes: `qmmm_eligible_top3.json`.
- Produces: a frozen list of zero to three recurrent NAC cluster medoids and explicit next decision.

- [ ] **Step 1: Verify the final gate**

If zero candidates pass, write `FAIL_UNRESTRAINED_M1_ENSEMBLE_NO_QMMM_ENTRY` with all reasons and stop before QM/MM. This does not prove every M1 realization is impossible.

- [ ] **Step 2: Compare proton paths only inside NAC frames**

For eligible candidates, report ThrOH→NαH2, ThrOH→Asp306/Asp308, and water-mediated geometry using the 100–1000 ps NAC denominator. If the denominator is insufficient, write `NOT_EVALUATED_PROTON_ROUTE_INSUFFICIENT_NAC`.

- [ ] **Step 3: Freeze the QM/MM entry universe**

Maximum three entries. Record recurrent-cluster membership, source replica/time, coordinate SHA256 on SCNet, M1 charges, atom mapping, and link-atom boundary candidates.

- [ ] **Step 4: Start a separate incremental Step1 plan**

The next plan must cover QM/MM mask/charge/link-atom preflight, short DFTB3/3OB-3-1 stability smoke, Step1 reaction-coordinate scans, TS/endpoint validation, then only conditionally Step2 and PMF. It must not assume a full mechanism before the ensemble evidence is available.
