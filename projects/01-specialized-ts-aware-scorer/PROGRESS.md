# Project 01 progress log

Updated: 2026-06-28

## Current execution state

Codex started the Project 01 implementation from the repository task outline.

Completed in this update:

- Prepared CPU server execution environment for lightweight Python tests:
  - host: `user-PowerEdge-R940`
  - user work root: `/Dell/Dell14/bianht`
  - virtual environment: `~/codex_venvs/enzyme_tasks`
  - Python: 3.11.4
  - installed packages: `pytest`, `jsonschema`, `numpy`, `pandas`, `pyyaml`, `scikit-learn`
- Added Project 01 scaffold directories:
  - `configs/`
  - `data_schema/`
  - `scripts/`
  - `src/ts_aware_scorer/`
  - `tests/`
  - `notebooks/README.md`
- Added first-pass JSON schemas for Project 01 manifests, labels, feature tables, and teacher export.
- Added deterministic lightweight Python utilities for:
  - reaction geometry calculations
  - electrostatic baseline calculations
  - conformer ensemble summaries
  - QM/MM label parsing and validation
  - feature-row assembly from a sample manifest
- Added the first runnable Project 01 baseline layer:
  - `src/ts_aware_scorer/models.py`
  - `scripts/prepare_dataset.py`
  - `scripts/train_baseline.py`
  - `scripts/predict_rank.py`
  - synthetic Project 01 CSV inputs
- Added unit tests covering the first feature utilities and schema examples.

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

## Next recommended step

Proceed to the next Project 01 increment:

- extend teacher export fields and add `scripts/export_teacher_labels.py`
- add baseline comparison documentation
- add a root `run_all.sh` once Project 02 has a matching minimal runnable layer
