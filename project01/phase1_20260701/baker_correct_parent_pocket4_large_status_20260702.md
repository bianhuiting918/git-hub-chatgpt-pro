# Project 01 Baker Correct-Parent 4A Large Candidate Status

Date: 2026-07-02

## Formal Pocket Definition

This status uses only the formal 4.0 A fixed-pocket definition:

- fixed residues: `15,16,19,38,39,42,55,59,73,77,122,148,149,150,151`
- parent: `04_af2/outputs/refined_out_1_bn1_1_20.pdb`
- ligand reference: `03_design/outputs/refined_out_1_bn1_1.pdb`
- ligand: `bn1`

No 8 A or 10 A expanded-pocket candidates are counted.

## Round01 Completed

Remote manifest:

`/data/bht/project01_baker_serhyd_routeB_20260701/manifests/baker_correct_parent_pocket4_ligandmpnn_bins_selected.tsv`

Generation size:

| bin | candidates entering ESMFold |
|---:|---:|
| 90 | 40 |
| 80 | 40 |
| 70 | 40 |
| 60 | 40 |
| 50 | 40 |

ESMFold completed: 200/200 OK.

Round01 active-gate result:

| bin | FINAL_QUALIFIED_ACTIVE | FAIL |
|---:|---:|---:|
| 90 | 3 | 37 |
| 80 | 0 | 40 |
| 70 | 0 | 40 |
| 60 | 0 | 40 |
| 50 | 0 | 40 |

## Large Candidate Pool Generated

Remote manifest:

`/data/bht/project01_baker_serhyd_routeB_20260701/manifests/baker_correct_parent_pocket4_large_ligandmpnn_bins_selected.tsv`

Generation size requested by user and achieved:

| bin | candidates selected |
|---:|---:|
| 90 | 200 |
| 80 | 1000 |
| 70 | 1000 |
| 60 | 1000 |
| 50 | 1000 |

Total selected: 4200.

Remote generation summary:

`/data/bht/project01_baker_serhyd_routeB_20260701/manifests/baker_correct_parent_pocket4_large_ligandmpnn_bins_summary.json`

## Large Batch01 ESMFold Started

Batch01 was deduplicated against round01 sequences.

Remote manifest:

`/data/bht/project01_baker_serhyd_routeB_20260701/manifests/baker_correct_parent_pocket4_large_batch01_selected.tsv`

Batch01 size:

| bin | candidates entering ESMFold |
|---:|---:|
| 90 | 195 |
| 80 | 200 |
| 70 | 200 |
| 60 | 200 |
| 50 | 200 |

Total: 995.

ESMFold run:

- PID at launch: `405163`
- log: `/data/bht/project01_baker_serhyd_routeB_20260701/logs/baker_correct_parent_pocket4_large_batch01_esmfold.log`
- status TSV: `/data/bht/project01_baker_serhyd_routeB_20260701/manifests/baker_correct_parent_pocket4_large_batch01_esmfold_status.tsv`
- summary JSON: `/data/bht/project01_baker_serhyd_routeB_20260701/manifests/baker_correct_parent_pocket4_large_batch01_esmfold_summary.json`

Next action: monitor batch01 ESMFold, run the same 4A active gate on batch01 outputs, merge with round01 passes, and continue with additional large batches for bins below 10.
