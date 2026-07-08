# Reports

This directory should contain lightweight generated summaries only, for example:

```text
candidate_ranking.tsv
candidate_report.md
```

Reports should include:

- sequence and Foldseek recall counts;
- Folddisco-after-recall pass / partial / fail counts;
- Folddisco-de-novo new candidate counts;
- L0/L1/L2/L3 motif score summaries;
- Pythia-Pocket / GRASE routing decisions;
- top candidates with explanations;
- proposed experimental panel design.

Do not commit large raw search outputs, structure databases, docking ensembles, trajectories, or model checkpoints.

Operational runbooks kept here:

- [`FOLDDISCO_STRUCTURE_RECOVERY_RUNBOOK.md`](FOLDDISCO_STRUCTURE_RECOVERY_RUNBOOK.md): how to recover structures from Folddisco online result API when direct AFDB/RCSB/ESMAtlas downloads fail.
