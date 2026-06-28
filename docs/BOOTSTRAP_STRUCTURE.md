# BOOTSTRAP STRUCTURE — EXECUTABLE CORE (Codex Init Layer)

This document defines the **minimum runnable architecture** required to convert this repository from design-state → executable-state.

It implements a strict rule:

> Project must be runnable even with mocked TS / PLACER / QM/MM.

---

# 🧱 1. FINAL DIRECTORY STRUCTURE (MUST EXIST)

```
git-hub-chatgpt-pro/
│
├── projects/
│   ├── 01-specialized-ts-aware-scorer/
│   │   ├── src/
│   │   │   └── ts_aware_scorer/
│   │   │       ├── __init__.py
│   │   │       ├── geometry.py
│   │   │       ├── electrostatics.py
│   │   │       ├── ensemble.py
│   │   │       ├── models.py
│   │   │       └── featurize.py
│   │   ├── scripts/
│   │   │   ├── train_baseline.py
│   │   │   ├── predict_rank.py
│   │   │   └── prepare_dataset.py
│   │   └── tests/
│   │
│   ├── 02-general-enzyme-prediction/
│   │   ├── src/
│   │   │   └── general_enzyme_predictor/
│   │   │       ├── __init__.py
│   │   │       ├── embeddings.py
│   │   │       ├── predicted_ts.py
│   │   │       ├── placer_screening.py
│   │   │       └── models.py
│   │   ├── scripts/
│   │   │   ├── train_student.py
│   │   │   └── build_features.py
│   │   └── tests/
│
├── models/
│   ├── esm/
│   ├── ts_encoder/        # MOCK
│   └── placer/             # MOCK
│
├── examples/
│   ├── synthetic_project01.csv
│   └── synthetic_project02.csv
│
├── outputs/
├── logs/
├── run_all.sh
└── docs/
    └── BOOTSTRAP_STRUCTURE.md
```

---

# ⚙️ 2. EXECUTION LAYER (RUN RULE)

System MUST support:

```bash
bash run_all.sh
```

This script is the ONLY entrypoint Codex must use.

---

# 🚀 3. run_all.sh (MANDATORY LOGIC)

```
STEP 1: run Project01 baseline
STEP 2: run Project02 student
STEP 3: validate TS leakage
STEP 4: save outputs
```

Expected behavior:

- No manual intervention
- No missing imports
- No external QM/MM dependency

---

# 🧪 4. MINIMAL PIPELINE CONTRACT

## Project 01 (true-TS teacher)

Input:

```
GS + TS + mock QM/MM label
```

Output:

```
ΔG‡ prediction
ranking
```

Model:

```
Ridge Regression (baseline)
```

---

## Project 02 (predicted-TS student)

Input:

```
GS + reaction prior + predicted TS embedding
```

Output:

```
barrier_proxy
ranking
```

Model:

```
MLP or Ridge baseline
```

---

# 🧬 5. REQUIRED MOCK COMPONENTS

## TS Encoder (MUST BE MOCK)

```python
import numpy as np

def encode_ts(gs, reaction_prior):
    return np.random.randn(256)
```

---

## PLACER (MUST BE MOCK)

```python
def generate_conformers(structure):
    return [{"id": 0, "score": 0.5}]
```

---

# 📦 6. SYNTHETIC DATA RULE

If missing:

Project 01:
```
(sample_id, enzyme_id, variant_id, features, delta_G)
```

Project 02:
```
(sample_id, protein_emb, ts_emb, barrier_proxy)
```

---

# 🧪 7. VALIDATION RULES (CRITICAL)

Codex MUST ensure:

## ❌ No leakage
```
projects/02-general-enzyme-prediction MUST NOT access true_TS
```

Check:
```bash
grep -R "true_TS" projects/02-general-enzyme-prediction
```

---

# 📊 8. SUCCESS CRITERIA

System is valid only if:

- Project 01 runs
- Project 02 runs
- run_all.sh executes end-to-end
- outputs/ contains CSV results
- no TS leakage

---

# 🧠 9. DESIGN PRINCIPLE

This system is intentionally:

> “mock-complete but execution-valid”

Meaning:

- Physics is abstracted
- Pipeline is fully runnable
- Interfaces are stable for future QM/MM upgrade

---

# END OF BOOTSTRAP SPEC
