# Config files

Codex should generate lightweight example configs here.

Required examples:

```text
targets.example.yaml
scoring.default.yaml
databases.example.yaml
```

The scoring config must support:

```text
sequence_evidence
fold_evidence
folddisco_score
pocket_embedding_similarity
predicted_stability
substrate_accessibility
novelty
```

Folddisco layer weights must be configurable:

```text
l0_core_score
l1_environment_score
l2_pocket_score
l3_context_score
```

External databases should be referenced by manifest paths only. Do not commit database dumps, model checkpoints, trajectories, or large raw search outputs.
