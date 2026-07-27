# Nylonase Pose-Retention, ProteinMPNN, SaProt, and Activity Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build and run an auditable exact-MD5 analysis relating PA66 AD-dimer pose retention, ProteinMPNN, and SaProt compatibility scores to four separate Chem Catalysis 2025 nylonase raster color-span endpoints.

**Architecture:** A pure-Python analysis package first assembles one enzyme-level authority table and explicit `NOT_EVALUATED` ledger, then computes endpoint-specific statistics and sequence-cluster sensitivities, and finally renders static/interactive figures. GitHub stores source, tests, plans, compact input snapshots and compact outputs; the immutable production run, logs, hashes and all generated artifacts live under a versioned Sugon output root.

**Tech Stack:** Python 3.9+, pandas, numpy, scipy, statsmodels, matplotlib, seaborn, plotly, pytest, mmseqs2 when available.

## Global Constraints

- Join all sources only by exact `sequence_md5`; enzyme is the statistical unit.
- Primary pose predictor is strict-NAC pose-retention; loose gate, weight=5 and weight=10 are sensitivities.
- Primary ProteinMPNN score is vanilla `v_48_020`; soluble `v_48_020` is sensitivity.
- Primary SaProt score is exact-canonical `mean_log_likelihood_FULL_PROTEIN` with `status_FULL_PROTEIN=PASS`.
- Preserve missing/failed inputs as reason-specific `NOT_EVALUATED_*`; never fill missing values with zero.
- Keep PA66-L1, PA66-L2, PA6 maximum and PA66 L1 fraction separate.
- Call experimental values raster color-span proxies, not absolute activity.
- Call model outputs pose retention or structure-conditional compatibility, not barriers, Tm or activity predictions.
- Local Project2 remains read-only; computation and generated artifacts are written only to Sugon.
- Production root: `/work/home/acshdt1dks/petase_orbmol_lg1_lg4_layer8343_20260721/outputs/nylon_pose_mpnn_saprot_activity_20260727_v1`.

---

### Task 1: Freeze compact inputs and assemble the exact-MD5 denominator

**Files:**
- Create: `analysis/nylon_pose_mpnn_saprot_activity/assemble.py`
- Create: `analysis/nylon_pose_mpnn_saprot_activity/schema.py`
- Create: `tests/nylon_pose_mpnn_saprot_activity/test_assemble.py`
- Create: `analysis/nylon_pose_mpnn_saprot_activity/inputs/pose_enzyme_retention_summary.tsv`
- Create: `analysis/nylon_pose_mpnn_saprot_activity/inputs/pose_enzyme_retention_by_stratum.tsv`
- Create: `analysis/nylon_pose_mpnn_saprot_activity/inputs/chemcatal2025_nyl95_activity.tsv`

**Interfaces:**
- Consumes: pose summary fields `sequence_md5, strict_nac_rate, retention_rate_loose_gate`; stratum fields `sequence_md5, stratum, retention_rate`; activity fields `sequence_md5, pa66_l1_color_span_px, pa66_l2_color_span_px, pa6_max_color_span_px, sequence`; ProteinMPNN and SaProt schemas verified live in the design.
- Produces: `assemble_authority(...) -> (joined_df, not_evaluated_df, denominator_summary)` and endpoint columns including `l1_fraction`.

- [ ] **Step 1: Write failing tests for uniqueness, exact-MD5 joins and missing values**

```python
def test_assemble_uses_exact_md5_and_never_zero_fills():
    joined, excluded, summary = assemble_authority(
        activity=activity_fixture(),
        pose=pose_fixture(),
        pose_strata=strata_fixture(),
        mpnn=mpnn_fixture(),
        saprot=saprot_fixture(),
    )
    assert summary["activity_authority_unique_md5"] == 3
    assert summary["four_way_intersection"] == 1
    assert joined.loc[joined.sequence_md5 == "a"*32, "l1_fraction"].iat[0] == 0.25
    assert excluded.query("sequence_md5 == @missing_md5").reason.str.startswith("NOT_EVALUATED_").all()
    assert not excluded.select_dtypes("number").fillna(0).eq(0).all(axis=1).any()
```

- [ ] **Step 2: Run the test and verify RED**

Run: `pytest -q tests/nylon_pose_mpnn_saprot_activity/test_assemble.py`  
Expected: import failure because `assemble_authority` does not exist.

- [ ] **Step 3: Implement minimal schema validation and assembler**

```python
def compute_l1_fraction(l1, l2):
    total = l1 + l2
    return np.where(total > 0, l1 / total, np.nan)

def valid_saprot(rows):
    return rows[
        rows["status_FULL_PROTEIN"].eq("PASS")
        & rows["scoring_sequence_scope"].eq("CANONICAL_EXACT")
        & rows["canonical_sequence_md5"].eq(rows["scoring_sequence_md5"])
    ].rename(columns={"canonical_sequence_md5": "sequence_md5",
                      "mean_log_likelihood_FULL_PROTEIN": "saprot_full_mean"})

def valid_mpnn(rows):
    return rows[
        rows["family"].eq("Nylonase")
        & rows["proteinmpnn_status"].eq("PASS")
        & rows["scoring_sequence_scope"].eq("CANONICAL_EXACT")
        & rows["canonical_sequence_md5"].eq(rows["scoring_sequence_md5"])
    ].rename(columns={"canonical_sequence_md5": "sequence_md5",
                      "vanilla_mean_log_likelihood": "mpnn_vanilla_mean",
                      "soluble_mean_log_likelihood": "mpnn_soluble_mean"})
```

- [ ] **Step 4: Run the test and verify GREEN**

Run: `pytest -q tests/nylon_pose_mpnn_saprot_activity/test_assemble.py`  
Expected: all tests pass.

- [ ] **Step 5: Commit the input and assembler slice**

Commit message: `feat: assemble exact-md5 nylon activity denominator`.

### Task 2: Implement endpoint statistics and multiplicity control

**Files:**
- Create: `analysis/nylon_pose_mpnn_saprot_activity/statistics.py`
- Create: `tests/nylon_pose_mpnn_saprot_activity/test_statistics.py`

**Interfaces:**
- Consumes: complete-case endpoint frames with `pose_strict_nac, mpnn_vanilla_mean, saprot_full_mean`.
- Produces: `correlation_table`, `bootstrap_table`, `multivariable_table`, `predictor_correlation_table`, and `model_audit`.

- [ ] **Step 1: Write failing tests for Pearson/Spearman, Holm, bootstrap and VIF gates**

```python
def test_endpoint_statistics_match_known_monotone_fixture():
    out = analyze_endpoint(monotone_fixture(), "pa66_l1_color_span_px", seed=20260727, n_boot=500)
    assert out.correlations.query("predictor == 'pose_strict_nac'").spearman_rho.iat[0] == 1.0
    assert out.correlations.p_holm.between(0, 1).all()
    assert out.bootstrap.n_success.min() >= 450
    assert out.model_audit["n_complete"] == len(monotone_fixture())
```

- [ ] **Step 2: Run and verify RED**

Run: `pytest -q tests/nylon_pose_mpnn_saprot_activity/test_statistics.py`  
Expected: import failure for `analyze_endpoint`.

- [ ] **Step 3: Implement statistics with declared numerical gates**

```python
PRIMARY = ["pose_strict_nac", "mpnn_vanilla_mean", "saprot_full_mean"]
ENDPOINTS = ["pa66_l1_color_span_px", "pa66_l2_color_span_px",
             "pa6_max_color_span_px", "l1_fraction"]
VIF_WITHHOLD = 10.0
CONDITION_WITHHOLD = 1e8

def holm(p):
    order = np.argsort(p)
    adjusted = np.empty(len(p), float)
    running = 0.0
    for rank, idx in enumerate(order):
        running = max(running, min(1.0, (len(p)-rank)*p[idx]))
        adjusted[idx] = running
    return adjusted
```

Use 10,000 deterministic bootstrap resamples in production, standardized OLS coefficients, two-sided 95% percentile intervals, VIF and design condition number. Withhold multivariable coefficients when `max_vif > 10`, condition number exceeds `1e8`, predictor variance is zero, or `n <= p+2`.

- [ ] **Step 4: Run and verify GREEN**

Run: `pytest -q tests/nylon_pose_mpnn_saprot_activity/test_statistics.py`  
Expected: all tests pass.

- [ ] **Step 5: Commit**

Commit message: `feat: add audited nylon endpoint statistics`.

### Task 3: Add sequence-homology sensitivity

**Files:**
- Create: `analysis/nylon_pose_mpnn_saprot_activity/clusters.py`
- Create: `tests/nylon_pose_mpnn_saprot_activity/test_clusters.py`

**Interfaces:**
- Consumes: activity authority exact sequences.
- Produces: deterministic cluster assignments and cluster-aware bootstrap correlation summaries.

- [ ] **Step 1: Write failing tests**

```python
def test_cluster_bootstrap_samples_clusters_not_duplicate_enzymes():
    assignments = {"a": "c1", "b": "c1", "c": "c2"}
    sampled = sample_one_per_cluster(assignments, np.random.default_rng(7))
    assert len(sampled) == 2
    assert set(sampled) <= {"a", "b", "c"}
```

- [ ] **Step 2: Run and verify RED**

Run: `pytest -q tests/nylon_pose_mpnn_saprot_activity/test_clusters.py`  
Expected: import failure.

- [ ] **Step 3: Implement formal and fallback clustering**

Run MMseqs2 `easy-cluster --min-seq-id 0.30 -c 0.80 --cov-mode 0 --cluster-mode 2` when available. If unavailable, mark `NOT_EVALUATED_CLUSTER_TOOL`; do not silently replace it with single-link clustering. Implement one-representative-per-cluster resampling with 10,000 deterministic replicates.

- [ ] **Step 4: Run and verify GREEN**

Run: `pytest -q tests/nylon_pose_mpnn_saprot_activity/test_clusters.py`  
Expected: all tests pass.

- [ ] **Step 5: Commit**

Commit message: `feat: add sequence-cluster sensitivity analysis`.

### Task 4: Render four endpoint-specific figure packages

**Files:**
- Create: `analysis/nylon_pose_mpnn_saprot_activity/plotting.py`
- Create: `tests/nylon_pose_mpnn_saprot_activity/test_plotting.py`

**Interfaces:**
- Consumes: the exact endpoint complete-case frame and companion statistics.
- Produces per endpoint: `*_3d.html`, `*_3d.png`, `*_2d_panels.png`, `*_figure_audit.json`.

- [ ] **Step 1: Write failing tests for row identity and scale separation**

```python
def test_figure_audit_matches_complete_case_md5s(tmp_path):
    audit = render_endpoint_figures(fixture(), "pa66_l2_color_span_px", tmp_path)
    assert audit["endpoint"] == "pa66_l2_color_span_px"
    assert audit["n_points"] == 4
    assert audit["sequence_md5_sha256"] == expected_md5_set_sha256()
    assert "pa6" not in audit["y_axis_label"].lower()
```

- [ ] **Step 2: Run and verify RED**

Run: `pytest -q tests/nylon_pose_mpnn_saprot_activity/test_plotting.py`  
Expected: import failure.

- [ ] **Step 3: Implement figures**

Use pose strict-NAC, vanilla ProteinMPNN and SaProt as 3D axes; color by exactly one endpoint. Produce three 2D predictor panels with regression line and bootstrap band. Hover labels show enzyme and exact MD5; static labels are limited to predeclared endpoint extrema and Cook's-distance influential points.

- [ ] **Step 4: Run and verify GREEN**

Run: `pytest -q tests/nylon_pose_mpnn_saprot_activity/test_plotting.py`  
Expected: all tests pass and PNG/HTML artifacts exist.

- [ ] **Step 5: Commit**

Commit message: `feat: render separate PA66 PA6 nylon figures`.

### Task 5: Build and execute the Sugon production run

**Files:**
- Create: `analysis/nylon_pose_mpnn_saprot_activity/run_analysis.py`
- Create: `analysis/nylon_pose_mpnn_saprot_activity/audit_run.py`
- Create: `analysis/nylon_pose_mpnn_saprot_activity/RUNBOOK.md`
- Create: `tests/nylon_pose_mpnn_saprot_activity/test_audit_run.py`

**Interfaces:**
- Consumes: immutable path/SHA manifest and all modules above.
- Produces: versioned production outputs, `run_history.tsv`, `input_sha256.tsv`, `audit.json`, and root `PASS.json`.

- [ ] **Step 1: Write failing independent-reconstruction tests**

```python
def test_audit_rebuilds_denominators_and_correlations(tmp_path):
    report = audit_run(tmp_path)
    assert report["status"] == "PASS"
    assert report["tests"]["no_missing_filled_zero"]
    assert report["tests"]["four_way_denominator_reconciles"]
    assert report["tests"]["stored_correlations_recompute"]
    assert report["tests"]["figure_md5_sets_match"]
```

- [ ] **Step 2: Run and verify RED**

Run: `pytest -q tests/nylon_pose_mpnn_saprot_activity/test_audit_run.py`  
Expected: import failure.

- [ ] **Step 3: Implement launcher and audit**

Production command:

```bash
python run_analysis.py \
  --pose-summary inputs/pose_enzyme_retention_summary.tsv \
  --pose-strata inputs/pose_enzyme_retention_by_stratum.tsv \
  --activity inputs/chemcatal2025_nyl95_activity.tsv \
  --mpnn /work/home/acshdt1dks/petase_orbmol_lg1_lg4_layer8343_20260721/outputs/proteinmpnn_compatibility_pet8329_nylon4167_20260727_v2_final/proteinmpnn_candidate_scores.tsv \
  --saprot /work/home/acshdt1dks/petase_orbmol_lg1_lg4_layer8343_20260721/outputs/saprot_nylonase_4556_denominator_20260725_v1/nylonase_4556_saprot_denominator.tsv \
  --output . --seed 20260727 --bootstrap 10000
python audit_run.py --run-root .
```

- [ ] **Step 4: Execute smoke, production and independent audit on Sugon**

Expected smoke: all exact-MD5 and endpoint formula tests pass on a balanced subset.  
Expected production: four endpoint packages plus primary and sensitivity statistics.  
Expected audit: root `PASS.json` only after all declared gates pass.

- [ ] **Step 5: Commit**

Commit message: `feat: complete audited nylon pose score activity pipeline`.

### Task 6: Publish compact verified artifacts and report scientific limits

**Files:**
- Create: `analysis/nylon_pose_mpnn_saprot_activity/results/denominator_summary.tsv`
- Create: `analysis/nylon_pose_mpnn_saprot_activity/results/correlation_summary.tsv`
- Create: `analysis/nylon_pose_mpnn_saprot_activity/results/model_summary.tsv`
- Create: `analysis/nylon_pose_mpnn_saprot_activity/results/NOT_EVALUATED_summary.tsv`
- Create: `analysis/nylon_pose_mpnn_saprot_activity/results/audit.json`
- Create: `analysis/nylon_pose_mpnn_saprot_activity/README.md`

**Interfaces:**
- Consumes: only a remote run with root `PASS.json`.
- Produces: compact, checksum-traceable GitHub handoff.

- [ ] **Step 1: Verify remote PASS and hashes**

Run: `python audit_run.py --run-root <versioned-run-root>`  
Expected: exit 0 and every gate true.

- [ ] **Step 2: Publish only compact artifacts**

Exclude structures, docking poses, model weights, caches and credentials. Include the exact remote run root, source commit, input SHA-256 values, denominators and caveats in `README.md`.

- [ ] **Step 3: Final verification**

Run: `pytest -q tests/nylon_pose_mpnn_saprot_activity` and recompute compact-file SHA-256 values.  
Expected: zero failures and hashes matching the remote handoff.

- [ ] **Step 4: Commit**

Commit message: `docs: publish audited nylon correlation results`.
