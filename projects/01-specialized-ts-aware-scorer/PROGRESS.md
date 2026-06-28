# Project 01 progress log

Updated: 2026-06-29 00:05 CST

## Current execution state

Codex started Project 01 from the repository task outline and has now moved from scaffold validation into a real external GMX-CP2K smoke/batch calculation on the CPU server.

Important upload scope note:

- This file is a progress summary for ChatGPT Pro and the user.
- The tested scaffold/code/schema/test files were created and validated in a local/CPU temporary workspace.
- They have not all been uploaded to GitHub, per user instruction.
- GitHub should receive selected artifacts only after the user decides what should be kept in the shared repository.
- Raw GMX-CP2K outputs, trajectories, `.edr`, `.trr`, wavefunction restart files, and large batch result directories are not uploaded.

Completed locally / on CPU server before the real compute step:

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

## Real GMX-CP2K execution status

User clarified that no PLACER run exists yet; the available real external calculation path is GMX-CP2K.

CPU server: `user-PowerEdge-R940`

WT reference work directory:

```text
/Dell/Dell14/bianht/gromacs_cp2k/PETase_WT_QMMM_vacuum_total
```

Batch work directory:

```text
/Dell/Dell14/bianht/gromacs_cp2k/PETase_QMMM_vacuum_total_batch
```

Active WT GMX-CP2K jobs observed at 2026-06-29 00:03 CST:

- GS: `wt_gs_lg3_qmmm_vacuum_total`, PID `1983512`, status running, wall time about 2h35m.
- TS: `wt_ts_lg4_qmmm_vacuum_total`, PID `1983511`, status running, wall time about 2h35m.
- Both are using the patched `gmx_mpi_d` / GMX-CP2K path and `-ntomp 4 -nsteps 1`.
- Batch worker PID `3529069` is alive and waiting for WT PIDs `1983511,1983512` before launching the remaining mutant calculations.

Current log interpretation:

- Both WT jobs have produced initial GROMACS energy tables and CP2K output.
- `wt_gs_lg3_qmmm_vac.out` is actively updating; observed SCF progression through OT steps around step 20.
- `wt_ts_lg4_qmmm_vac.out` is actively updating; observed one inner SCF loop reaching 30 steps, then entering the next outer SCF iteration.
- The jobs are therefore running real CP2K work, not just dry-run tests.

Important caveat:

- This remains a low-grid, ligand-only QM/MM smoke/batch descriptor workflow. It is useful for first-pass GMX-CP2K-vs-OrbMol correlation screening, not for final production free-energy claims.

## Scientific boundary preserved

The GitHub repository itself still stores orchestration code, schemas, summaries, and lightweight status. External GMX-CP2K calculations run on the CPU server and raw computational outputs remain outside version control.

## Next execution step

Monitor the two WT reference jobs until they exit. If both exit successfully, parse GS/TS energies, record the WT descriptor, then let the batch worker start the remaining PETase mutant jobs and summarize only lightweight results back to GitHub.
