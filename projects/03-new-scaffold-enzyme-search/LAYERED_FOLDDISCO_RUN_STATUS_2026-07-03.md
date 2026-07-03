# Layered Folddisco Run Status

Date: 2026-07-03
Project: enzyme scaffold search v2

## Scope

This status records the first CPU-capped after-recall L0-L3 Folddisco score-only run for three enzyme families:

- PETase / cutinase-like polyester hydrolases
- Nylonase / nylon hydrolases
- Polyurethane / carbamate hydrolases

Laccase was not included in this L0-L3 pass.

## Server Artifacts

Remote project path:

```text
/Dell/Dell14/bianht/enzyme_scaffold_search_v2
```

Created inputs and scripts:

- `configs/layered_folddisco_definitions.tsv`
- `scripts/33_build_layered_folddisco_queries.py`
- `scripts/34_run_layered_folddisco_after_recall.py`
- `scripts/35_score_layered_folddisco_results.py`

Run outputs:

- `results/layered_folddisco/inputs/layered_motif_queries.tsv`
- `results/layered_folddisco/runs/after_recall/query_result_manifest.tsv`
- `results/layered_folddisco/runs/after_recall/scores/layered_folddisco_hit_scores.tsv`
- `results/layered_folddisco/runs/after_recall/scores/layered_folddisco_score_summary.tsv`
- `results/layered_folddisco/runs/after_recall/logs/nohup_layered_after_recall.log`

## CPU Cap

Run command used CPU-capped settings:

```text
python3 scripts/34_run_layered_folddisco_after_recall.py --threads 16 --max-parallel-jobs 1 --max-total-cpu 64 --force
```

This enforces `threads * max_parallel_jobs <= max_total_cpu`, so this run used one 16-thread Folddisco process at a time and stayed below the requested 64-core cap.

## Query Generation

Generated query table:

```text
results/layered_folddisco/inputs/layered_motif_queries.tsv
```

Counts:

| Enzyme class | Total generated queries | After-recall runnable queries |
|---|---:|---:|
| PETase | 35 | 35 |
| nylonase | 19 | 17 |
| polyurethane_carbamate_hydrolase | 16 | 16 |
| Total | 70 | 68 |

The two non-runnable nylonase queries came from the transferred Nyl10 motif because the local seed structure path was unavailable. These are `NOT_EVALUATED`-type setup gaps, not biological failures.

## Folddisco Native Result Summary

Evaluated universe: after-recall local candidate structures from previous sequence/Foldseek/Folddisco recall runs.

Candidate structure counts:

| Enzyme class | Candidate structures |
|---|---:|
| PETase | 4059 |
| nylonase | 4014 |
| polyurethane_carbamate_hydrolase | 4856 |

Manifest result:

| Enzyme class | Runnable queries | Folddisco query status |
|---|---:|---|
| PETase | 35 | 35 OK |
| nylonase | 17 | 17 OK |
| polyurethane_carbamate_hydrolase | 16 | 16 OK |
| Total | 68 | 68 OK |

No query in the completed after-recall run was treated as a biological fail. Folddisco returned native result rows under each query setting and those rows were parsed into project scores.

## Score Outputs

Score parser output:

```text
wrote scores results/layered_folddisco/runs/after_recall/scores scored_rows=150501 summary_rows=68
status {'SCORED': 68}
classes {'PETase': 35, 'nylonase': 17, 'polyurethane_carbamate_hydrolase': 16}
```

Output rows:

| Output file | Rows |
|---|---:|
| `layered_folddisco_score_summary.tsv` | 68 |
| `layered_folddisco_hit_scores.tsv` | 150501 |

## Filtering Policy

This run is score-only.

- `filter_decision = NOT_FILTERED_FOR_PYTHIA`
- `STRICT_PASS` is not used as the main output label in this stage.
- Low Folddisco score, partial motif coverage, or loose geometry is retained for downstream Pythia/GRASE-style scoring.
- Missing seed structure or unavailable input is not called biological failure.

## Important Debug Note

The first attempted full run exposed a runner bug: Folddisco query output subdirectories such as `results/.../results/PETase/` were not being created before calling Folddisco. Folddisco therefore failed with file-create errors. The runner was fixed to create `result_path.parent` before each query, regression-checked with a dry run, and then restarted with `--force`.

The completed run described above is from the fixed runner.

## Next Stage

Use `layered_folddisco_hit_scores.tsv` and `layered_folddisco_score_summary.tsv` as score-only support tables for Pythia-pocket / GRASE-style candidate prioritization. The next integration should preserve:

- Folddisco native hit evidence
- query variant and included layers
- matched count and query count
- min RMSD and loose cutoff
- layer score
- NOT_FILTERED_FOR_PYTHIA status
