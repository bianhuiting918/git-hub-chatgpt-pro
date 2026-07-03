# Source package

Codex should create the Python package here:

```text
src/scaffold_search/
  __init__.py
  schema.py
  io.py
  scoring.py
  reporting.py
  sequence_search.py
  structure_search.py
  motif_search.py
  pocket_rank.py
```

Key implementation requirement: `motif_search.py` must support layered Folddisco evidence:

```text
after_recall
de_novo
L0_catalytic_core
L1_catalytic_environment
L2_substrate_pocket
L3_processing_interface
```
