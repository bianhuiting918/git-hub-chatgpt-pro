# Project 01 progress log

Updated: 2026-06-28

## Current execution state

Codex started the Project 01 implementation from the repository task outline.

Important upload scope note:

- This file is a progress summary for ChatGPT Pro and the user.
- The tested scaffold/code/schema/test files were created and validated in a local/CPU temporary workspace.
- They have not all been uploaded to GitHub, per user instruction.
- GitHub should receive selected artifacts only after the user decides what should be kept in the shared repository.

Completed locally / on CPU server in this update:

- Prepared CPU server execution environment for lightweight Python tests:
  - host: `user-PowerEdge-R940`
  - user work root: `/Dell/Dell14/bianht`
  - virtual environment: `~/codex_venvs/enzyme_tasks`
  - Python: 3.11.4
  - installed packages: `pytest`, `jsonschema`, `numpy`, `pandas`, `pyyaml`, `scikit-learn`
- Built a Project 01 scaffold candidate:
  - `configs/`
  - `data_schema/`
  - `scripts/`
  - `src/ts_aware_scorer/`
  - `tests/`
  - `notebooks/README.md`
- Built first-pass JSON schema candidates for Project 01 manifests, labels, feature tables, and teacher export.
- Built deterministic lightweight Python utility candidates for:
  - reaction geometry calculations
  - electrostatic baseline calculations
  - conformer ensemble summaries
  - QM/MM label parsing and validation
  - feature-row assembly from a sample manifest
- Built a first runnable Project 01 baseline candidate:
  - `src/ts_aware_scorer/models.py`
  - `scripts/prepare_dataset.py`
  - `scripts/train_baseline.py`
  - `scripts/predict_rank.py`
  - synthetic Project 01 CSV inputs
- Built unit tests covering the first feature utilities and schema examples.

## Validation completed

CPU server validation command:

```bash
source ~/codex_venvs/enzyme_tasks/bin/activate
cd ~/codex_project01_test
PYTHONPATH=src python -m pytest -q
python scripts/train_baseline.py --features examples/synthetic_project01/features.csv --target delta_delta_G_dagger_vs_reference
python scripts/predict_rank.py --model outputs/project01/models/ridge_baseline.pkl --features examples/synthetic_project01/features.csv --output outputs/project01/predicted_ranking.csv
```

Observed result:

- tests: `19 passed in 0.10s`
- training output: `outputs/project01/models/ridge_baseline.pkl`
- metrics output: `outputs/project01/metrics.json`
- ranking output: `outputs/project01/ranked_variants.csv`
- teacher export output: `outputs/project01/project01_teacher_export.csv`

Synthetic-data training metrics:

```json
{
  "n_samples": 4,
  "rmse": 0.25316307857803194,
  "mae": 0.252256897425643,
  "r2": 0.983657932913439,
  "target": "delta_delta_G_dagger_vs_reference"
}
```

## Scientific boundary preserved

This update does not implement PLACER, QM/MM, DFT, TS search, or protein embedding inference. It only creates the orchestration and validation layer that can consume outputs from those external workflows.

## Next recommended decision

Decide which artifacts should be uploaded to GitHub:

1. progress summary only,
2. minimal runnable scripts and synthetic CSV only,
3. full Project 01 scaffold, schemas, source modules, and tests.
