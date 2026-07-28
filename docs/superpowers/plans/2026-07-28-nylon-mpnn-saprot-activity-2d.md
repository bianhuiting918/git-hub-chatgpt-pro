# Nylonase ProteinMPNN-SaProt 2D Activity Maps Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build audited, versioned Sugon outputs containing three ProteinMPNN-versus-SaProt maps with separate PA66-L1, PA66-L2, and PA6 experimental rank overlays.

**Architecture:** A pure Python join/rank module reads three authoritative TSV files, freezes exact-MD5 background and experimental ledgers, and exposes deterministic tables. A plotting command consumes only those frozen tables, writes interactive Plotly HTML plus static Matplotlib PNG/PDF, then a verifier independently rebuilds counts and hashes before writing `PASS.json`.

**Tech Stack:** Python 3, standard library `csv/hashlib/json`, pandas, numpy, scipy, matplotlib, Plotly, pytest.

## Global Constraints

- All computation, code, dependencies, caches, and outputs stay under the Sugon project.
- Local Project2 is read-only.
- Join only by exact `sequence_md5`; no name or fuzzy sequence joins.
- Background requires exact-canonical, finite vanilla ProteinMPNN and finite full-protein SaProt mean log-likelihood.
- Use 3,985 unique exact-MD5 rows only if live inputs reproduce that count; otherwise stop before plotting and report drift.
- Experimental authority is Nyl01-Nyl95, excluding NylC; missing rows are `NOT_EVALUATED_COMPATIBILITY_MISSING_EXACT_SCORE`, never zero.
- Endpoint values are raster color-span proxies, not absolute activity.
- Maintain `RUNBOOK.md`, append-only `logs/run_history.tsv`, input/output SHA256, and a fully constrained `PASS.json`.

---

### Task 1: Freeze exact-MD5 tables

**Files:**
- Create: `scripts/nylon_mpnn_saprot_2d.py`
- Create: `tests/test_nylon_mpnn_saprot_2d.py`
- Create: `RUNBOOK.md`
- Create: `logs/run_history.tsv`

**Interfaces:**
- Consumes: scope-fixed manifest, final ProteinMPNN score TSV, ChemCatal Nyl95 activity TSV.
- Produces: `build_tables(manifest_path, score_path, activity_path) -> tuple[pandas.DataFrame, pandas.DataFrame, pandas.DataFrame, dict]`.

- [ ] **Step 1: Write the failing join tests**

```python
def test_exact_scope_and_unique_md5_only(tmp_path):
    background, overlay, missing, summary = build_tables(
        fixture_manifest(tmp_path), fixture_scores(tmp_path), fixture_activity(tmp_path)
    )
    assert background["scoring_sequence_md5"].is_unique
    assert set(background["scoring_sequence_scope"]) == {"CANONICAL_EXACT"}
    assert summary["background_n"] == 2

def test_missing_experimental_score_is_not_zero(tmp_path):
    background, overlay, missing, summary = build_tables(
        fixture_manifest(tmp_path), fixture_scores(tmp_path), fixture_activity(tmp_path)
    )
    assert summary["experimental_authority_n"] == 3
    assert summary["experimental_linked_n"] == 2
    assert summary["experimental_not_evaluated_n"] == 1
    assert missing.iloc[0]["technical_status"] == "NOT_EVALUATED_COMPATIBILITY_MISSING_EXACT_SCORE"
    assert not (overlay[["vanilla_mean_log_likelihood", "saprot_mean_log_likelihood"]] == 0).any().any()
```

- [ ] **Step 2: Run tests and confirm RED**

Run: `pytest -q tests/test_nylon_mpnn_saprot_2d.py`  
Expected: collection/import failure because `nylon_mpnn_saprot_2d.py` does not exist.

- [ ] **Step 3: Implement the minimal deterministic join**

```python
ENDPOINTS = {
    "pa66_l1": "pa66_l1_color_span_px",
    "pa66_l2": "pa66_l2_color_span_px",
    "pa6": "pa6_max_color_span_px",
}

def is_finite_number(value):
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False

def build_tables(manifest_path, score_path, activity_path):
    scores = pd.read_csv(score_path, sep="\t", dtype=str)
    nylon = scores[
        (scores.family == "Nylonase")
        & (scores.scoring_sequence_scope == "CANONICAL_EXACT")
        & (scores.proteinmpnn_status == "PASS")
    ].copy()
    keep = nylon.vanilla_mean_log_likelihood.map(is_finite_number)
    keep &= nylon.saprot_mean_log_likelihood.map(is_finite_number)
    background = nylon.loc[keep].drop_duplicates("scoring_sequence_md5", keep=False)
    activity = pd.read_csv(activity_path, sep="\t", dtype=str)
    activity = activity[
        activity.table_s1_id.str.fullmatch(r"Nyl(?:0[1-9]|[1-8][0-9]|9[0-5])")
    ].copy()
    overlay = activity.merge(
        background, left_on="sequence_md5", right_on="scoring_sequence_md5",
        how="inner", validate="one_to_one"
    )
    missing = activity.loc[~activity.sequence_md5.isin(overlay.sequence_md5)].copy()
    missing["technical_status"] = "NOT_EVALUATED_COMPATIBILITY_MISSING_EXACT_SCORE"
    summary = {
        "background_n": len(background),
        "experimental_authority_n": len(activity),
        "experimental_linked_n": len(overlay),
        "experimental_not_evaluated_n": len(missing),
    }
    return background, overlay, missing, summary
```

- [ ] **Step 4: Run tests and confirm GREEN**

Run: `pytest -q tests/test_nylon_mpnn_saprot_2d.py`  
Expected: all tests pass.

### Task 2: Rank each endpoint independently

**Files:**
- Modify: `scripts/nylon_mpnn_saprot_2d.py`
- Modify: `tests/test_nylon_mpnn_saprot_2d.py`

**Interfaces:**
- Consumes: experimental overlay table.
- Produces: `add_endpoint_ranks(overlay: DataFrame) -> DataFrame`.

- [ ] **Step 1: Write failing ranking tests**

```python
def test_endpoint_rank_is_descending_and_ties_are_average():
    frame = pd.DataFrame({"pa66_l1_color_span_px": ["10", "20", "20"]})
    ranked = add_endpoint_ranks(frame)
    assert ranked["pa66_l1_rank_desc"].tolist() == [3.0, 1.5, 1.5]
    assert ranked["pa66_l1_percentile"].between(0, 1).all()

def test_all_three_endpoint_rank_columns_exist(linked_fixture):
    ranked = add_endpoint_ranks(linked_fixture)
    for endpoint in ("pa66_l1", "pa66_l2", "pa6"):
        assert f"{endpoint}_rank_desc" in ranked
        assert f"{endpoint}_percentile" in ranked
```

- [ ] **Step 2: Run the targeted tests and confirm RED**

Run: `pytest -q tests/test_nylon_mpnn_saprot_2d.py -k rank`  
Expected: failure because `add_endpoint_ranks` is absent.

- [ ] **Step 3: Implement ranking**

```python
def add_endpoint_ranks(overlay):
    out = overlay.copy()
    for label, column in ENDPOINTS.items():
        values = pd.to_numeric(out[column], errors="coerce")
        if values.isna().any():
            raise ValueError(f"non-finite endpoint: {column}")
        out[f"{label}_rank_desc"] = values.rank(method="average", ascending=False)
        out[f"{label}_percentile"] = values.rank(method="average", pct=True)
    return out
```

- [ ] **Step 4: Run all tests and confirm GREEN**

Run: `pytest -q tests/test_nylon_mpnn_saprot_2d.py`  
Expected: all tests pass.

### Task 3: Generate identical-axis interactive and static plots

**Files:**
- Modify: `scripts/nylon_mpnn_saprot_2d.py`
- Modify: `tests/test_nylon_mpnn_saprot_2d.py`

**Interfaces:**
- Consumes: frozen background and ranked overlay tables.
- Produces: `make_plots(background, overlay, output_dir) -> list[pathlib.Path]`.

- [ ] **Step 1: Write failing plot-contract test**

```python
def test_plot_contract_uses_gray_background_and_fixed_axes(tmp_path, tables):
    manifest = make_plots(*tables, tmp_path)
    assert {p.suffix for p in manifest} >= {".html", ".png", ".pdf"}
    contract = json.loads((tmp_path / "plot_contract.json").read_text())
    assert contract["background_color"] == "#B8B8B8"
    assert len({tuple(v) for v in contract["axis_ranges"].values()}) == 1
    assert set(contract["endpoints"]) == {"pa66_l1", "pa66_l2", "pa6"}
```

- [ ] **Step 2: Run plot test and confirm RED**

Run: `pytest -q tests/test_nylon_mpnn_saprot_2d.py -k plot_contract`  
Expected: failure because `make_plots` is absent.

- [ ] **Step 3: Implement plots**

Use one globally calculated x/y range for all three endpoints. Draw the complete background first with `color="#B8B8B8"`, then overlay the 80 linked rows with `viridis` percentile color. Plotly hover must include raw proxy, descending rank, percentile, both scores, Nyl ID, and MD5. Matplotlib writes each endpoint PNG/PDF and one 1x3 combined panel.

- [ ] **Step 4: Run all tests and confirm GREEN**

Run: `pytest -q tests/test_nylon_mpnn_saprot_2d.py`  
Expected: all tests pass with non-empty plot files.

### Task 4: Production run and independent audit

**Files:**
- Create: `scripts/verify_nylon_mpnn_saprot_2d.py`
- Create: `tests/test_verify_nylon_mpnn_saprot_2d.py`
- Create at runtime: `PASS.json`, `input_sha256.tsv`, `output_sha256.tsv`, result tables and figures.

**Interfaces:**
- Consumes: immutable inputs and produced artifacts.
- Produces: `verify_run(output_dir, expected_background=3985, expected_authority=95, expected_linked=80) -> dict`.

- [ ] **Step 1: Write a failing audit test**

```python
def test_verifier_rejects_zero_filled_missing_rows(tmp_path, valid_run):
    missing = pd.read_csv(valid_run / "nyl95_not_evaluated.tsv", sep="\t")
    missing["vanilla_mean_log_likelihood"] = 0
    missing.to_csv(valid_run / "nyl95_not_evaluated.tsv", sep="\t", index=False)
    with pytest.raises(AssertionError, match="NOT_EVALUATED"):
        verify_run(valid_run, 3985, 95, 80)
```

- [ ] **Step 2: Run verifier test and confirm RED**

Run: `pytest -q tests/test_verify_nylon_mpnn_saprot_2d.py`  
Expected: import failure because verifier does not exist.

- [ ] **Step 3: Implement independent verification**

The verifier re-reads all TSVs, recomputes unique counts, verifies finite scores, checks linked plus missing equals 95, checks all expected artifacts are non-empty, computes SHA256, and writes `PASS.json` only when every declared gate is true.

- [ ] **Step 4: Run full production and audit**

Run:
```bash
python scripts/nylon_mpnn_saprot_2d.py --manifest ... --scores ... --activity ... --output-dir .
python scripts/verify_nylon_mpnn_saprot_2d.py --output-dir . --expected-background 3985 --expected-authority 95 --expected-linked 80
pytest -q
```

Expected: tests pass; production exits 0; verifier exits 0; `PASS.json` reports all gates true.

- [ ] **Step 5: Append run history**

Append start and completion records with UTC/local timestamps, exact commands, input/output paths, exit codes, counts, and PASS path. Do not record credentials.
