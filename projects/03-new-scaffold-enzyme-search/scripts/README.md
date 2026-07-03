# Scripts

Codex should create CLI wrappers here. Every script must support `--help`.

Expected scripts:

```text
run_sequence_search.py
run_foldseek_search.py
run_folddisco_after_recall.py
run_folddisco_denovo.py
run_pocket_ranking.py
merge_and_rank_candidates.py
```

`run_folddisco_after_recall.py` validates sequence/Foldseek recalled candidates.

`run_folddisco_denovo.py` scans a broad structure database for new-scaffold motif hits.

The merge script must preserve L0/L1/L2/L3 evidence and explain classification decisions.
