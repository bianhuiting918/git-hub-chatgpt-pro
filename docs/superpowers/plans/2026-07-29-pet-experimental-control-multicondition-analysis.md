# PET Experimental-Control Multicondition Analysis Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Reuse the 30 exact-canonical PET `PATCH_PASS` records to test catalytic and external 14 Å stickiness against Nature 2022 product amounts across Table D3 conditions, with Table D6 as low-power sensitivity.

**Architecture:** One dependency-free Python analysis script joins the immutable patch-feature TSV to the authoritative activity long table by `sequence_md5`. It writes condition-level and rank-normalized aggregate statistics plus one hash-bearing audit JSON. The same script contains deterministic self-tests, avoiding additional test files.

**Tech Stack:** Python 3.11 standard library; existing remote surface-screen environment; GitHub Draft PR #2 as code authority; Sugon result storage.

## Global Constraints

- Do not recompute structures, fields, probes, shells, or patches.
- Evaluate 30 exact-canonical PET controls; keep protein 202 as `NONCANONICAL_MAPPING_REQUIRED`.
- Use only `NATURAL2022 / sum aromatic products / mg/L`.
- Keep Table D3, Table D6 substrate forms, and PET/nylon denominators separate.
- Never pool raw mg/L values across conditions.
- Table D4 is excluded because each condition contains only 3–4 records.
- Missing activity is `NOT_EVALUATED_FOR_CONDITION`, not zero or inactive.
- Use deterministic permutation and bootstrap statistics and report BH-adjusted q values.
- Do not submit Slurm work.

---

### Task 1: Add the self-testing multicondition analysis script

**Files:**
- Create: `projects/03-polymer-surface-hotspot-screening/scripts/analyze_pet_control_multicondition_patch.py`

**Interfaces:**
- Consumes:
  - `results/experimental_control31_patch_relative_stickiness_20260729_v1/patch_relative_stickiness_metrics.tsv`
  - `/work/home/acshdt1dks/petase_orbmol_lg1_lg4_layer8343_20260721/inputs/activity/activity_energy_long_authority.tsv`
- Produces:
  - `d3_joined_records.tsv`
  - `d3_condition_associations.tsv`
  - `d3_crosscondition_aggregate.tsv`
  - `d3_crosscondition_associations.tsv`
  - `d3_direction_consistency.tsv`
  - `d6_crosscondition_aggregate.tsv`
  - `d6_form_associations.tsv`
  - `ANALYSIS_PASS.json`

- [ ] **Step 1: Add deterministic self-tests before the analysis entrypoint**

The script must expose and test:

```python
def rankdata(values: list[float]) -> list[float]: ...
def spearman(x: list[float], y: list[float]) -> float: ...
def bh(rows: list[dict], pkey: str, qkey: str) -> None: ...
def parse_condition(value: str) -> tuple[str, str, str, str]: ...
def condition_percentiles(rows: list[dict]) -> dict[str, float]: ...

def self_test() -> None:
    assert rankdata([3, 1, 1, 2]) == [4.0, 1.5, 1.5, 3.0]
    assert abs(spearman([1, 2, 3], [3, 2, 1]) + 1.0) < 1e-12
    assert parse_condition("Table D3|PET_unspecified|50.0C|H7.5") == (
        "Table D3", "PET_unspecified", "50.0C", "H7.5"
    )
    pct = condition_percentiles([
        {"sequence_md5": "a", "activity_value": 0.0},
        {"sequence_md5": "b", "activity_value": 0.0},
        {"sequence_md5": "c", "activity_value": 2.0},
    ])
    assert pct == {"a": 0.25, "b": 0.25, "c": 1.0}
```

- [ ] **Step 2: Verify the self-test initially fails before the script exists**

Run:

```bash
/work/home/acshdt1dks/polymer_surface_hotspot_screen_20260725/envs/surface-screen-py311/bin/python3.11   projects/03-polymer-surface-hotspot-screening/scripts/analyze_pet_control_multicondition_patch.py --self-test
```

Expected before sync: file-not-found failure.

- [ ] **Step 3: Implement the minimal analysis**

Required constants:

```python
PRIMARY_METRICS = (
    "r14_relative_external_sticky_fraction",
    "r14_composite_catalytic_candidate_percentile",
    "r14_mean_catalytic_candidate_percentile_ACOA",
)
ALL_AGGREGATE_METRICS = tuple(
    f"r{radius}_{suffix}"
    for radius in (6, 10, 14)
    for suffix in (
        "relative_external_sticky_fraction",
        "composite_catalytic_candidate_percentile",
        "mean_catalytic_candidate_percentile_ACOA",
    )
)
MIN_CONDITION_N = 10
N_PERM = 20000
N_BOOT = 5000
SEED = 20260729
```

Implementation requirements:

- verify the patch metrics contain exactly 30 unique MD5s;
- verify each activity `(sequence_md5, condition)` key is unique;
- condition-level D3 analysis uses only joined rows and records missing MD5 counts separately;
- per-condition rank percentile uses tied average ranks scaled as `(rank - 1)/(n - 1)`;
- cross-condition aggregation averages percentiles per protein and records `n_conditions`;
- D3 BH correction covers all condition-by-primary-metric tests;
- D3 aggregate BH correction covers the nine radius-by-metric tests;
- D6 BH correction covers all substrate-form-by-primary-metric tests;
- all TSVs use stable sorting and tab delimiters;
- `ANALYSIS_PASS.json` records input hashes, output hashes, exact denominators, excluded Table D4, protein 202 exclusion, random seed, permutation/bootstrap counts, and endpoint limitations.

- [ ] **Step 4: Run self-test and production analysis**

Run:

```bash
PY=/work/home/acshdt1dks/polymer_surface_hotspot_screen_20260725/envs/surface-screen-py311/bin/python3.11
$PY runs/experimental_control31_20260729_v1/analyze_pet_control_multicondition_patch.py --self-test
$PY runs/experimental_control31_20260729_v1/analyze_pet_control_multicondition_patch.py
```

Expected:

```text
SELF_TEST_PASS
MULTICONDITION_ANALYSIS_PASS
```

- [ ] **Step 5: Commit the script on Draft PR #2**

Commit message:

```text
analysis: compare PET patch metrics across assay conditions
```

### Task 2: Document, audit, and report

**Files:**
- Modify: `projects/03-polymer-surface-hotspot-screening/RUNBOOK.md`
- Remote modify: `runs/experimental_control31_20260729_v1/RUNBOOK.md`
- Remote append: `logs/run_history.tsv`

**Interfaces:**
- Consumes: Task 1 output files.
- Produces: independently checkable hashes, denominator summaries, and the scientific interpretation.

- [ ] **Step 1: Add the exact rerun command and claim boundary to both RUNBOOKs**

Document that the activity values are endpoint sums of aromatic products, not `kcat`, barriers, or single-product measurements.

- [ ] **Step 2: Independently audit results without importing the analysis script**

Use a separate one-shot standard-library Python audit to verify:

```python
assert audit["status"] == "MULTICONDITION_PATCH_ANALYSIS_PASS"
assert audit["patch_exact_canonical_n"] == 30
assert audit["protein_202_excluded"] is True
assert audit["table_d4_excluded"] is True
assert all(Path(v["path"]).is_file() for v in audit["outputs"].values())
assert all(sha256(Path(v["path"])) == v["sha256"] for v in audit["outputs"].values())
```

- [ ] **Step 3: Append one run-history record**

Record UTC time, script path, input/output paths, exit status, exact D3 and D6 denominators, and audit SHA256. Do not record credentials.

- [ ] **Step 4: Read the audited output and report**

Report separately:

- D3 cross-condition aggregate associations;
- condition-wise direction counts and adjusted significance;
- H7.5 temperature series;
- D6 substrate-form sensitivity with low-power label;
- what is supported, what is only nominal, and what cannot be inferred.

## Self-review

- Spec coverage: D3, cross-condition ranks, radius sensitivity, D6, D4 exclusion, hashes, and denominator labels are all assigned.
- Placeholder scan: no TBD/TODO or unspecified implementation steps.
- Type consistency: metric names exactly match the existing patch metrics TSV.
- Scope: one analysis script and two documentation updates; no structure or patch recomputation.
