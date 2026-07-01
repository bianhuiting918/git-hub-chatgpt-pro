# Baker reference distal deficit status

Date: 2026-07-02

## Active-pocket counts before deficit expansion

After the Baker-reference shell8 batch:

| Bin | Active-pocket pass | Active target | Status |
|---|---:|---:|---|
| 90 | 18 | 10 | active target met |
| 80 | 11 | 10 | active target met |
| 70 | 16 | 10 | active target met |
| 60 | 0 | 10 | deficit |
| 50 | 1 | 10 | deficit |

These are `FINAL_QUALIFIED_ACTIVE` counts only. Final QMMM-ready qualification still requires the stricter sidechain/ligand-compatibility gate.

## Deficit expansion strategy

Only the deficient 60 and 50 bins were expanded.

Generation script:

- `scripts/generate_baker_reference_distal_deficit_bins.py`

Selection rule:

- Freeze the active-site CA shell.
- Mutate only the farthest positions from the 14-residue Baker motif.

Generated candidates:

| Bin | Shell cutoff | Fixed residues | Mutated positions | Selected candidates |
|---|---:|---:|---:|---:|
| 60 | 10 A | 89 | 64 | 100 |
| 50 | 8 A | 60 | 80 | 100 |

Remote manifest:

- `/data/bht/project01_baker_serhyd_routeB_20260701/manifests/baker_reference_distal_deficit_bins_selected.tsv`

## Running ESMFold job

Run id: `baker_reference_distal_deficit_esmfold_200`

- PID at launch: 193403
- rows: 200
- log: `/data/bht/project01_baker_serhyd_routeB_20260701/logs/baker_reference_distal_deficit_esmfold_200.log`
- status TSV: `/data/bht/project01_baker_serhyd_routeB_20260701/manifests/baker_reference_distal_deficit_esmfold_200_status.tsv`
- summary JSON: `/data/bht/project01_baker_serhyd_routeB_20260701/manifests/baker_reference_distal_deficit_esmfold_200_summary.json`

Next step: run `gate_baker_reference_bins_active.py` on partial or complete outputs, then merge 60/50 active-pocket pass counts with the shell8 pass table.
