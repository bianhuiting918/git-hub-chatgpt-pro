# Tests

Codex should add tests for the lightweight pipeline.

Required coverage:

```text
test_scoring.py
test_classification.py
test_io.py
test_merge.py
test_folddisco_modes.py
test_layered_motif_schema.py
```

Tests must verify:

- L0/L1/L2/L3 score aggregation;
- Folddisco-after-recall does not create duplicate candidates;
- Folddisco-de-novo unique hits enter C_total;
- T227-like L3 mismatch does not cause L0 failure;
- Pythia/GRASE routing includes strict pass, de-novo hit, partial + strong context, boundary, and negative controls;
- explanations mention sequence, fold, Folddisco mode, layer evidence, pocket, stability, novelty, and penalties.
