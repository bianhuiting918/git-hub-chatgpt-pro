# Baker reference shell8 status

Date: 2026-07-02

## Current target

Final target remains 10 final-qualified sequences per similarity bin:

| Bin | Target | Current active-qualified count before shell8 batch |
|---|---:|---:|
| 90 | 10 | 1 |
| 80 | 10 | 0 |
| 70 | 10 | 0 |
| 60 | 10 | 0 |
| 50 | 10 | 1 |

These are active-pocket gate counts only. Sidechain/ligand compatibility remains a stricter downstream gate before QMMM.

## Why the route was reset

Sample21 de novo parent attempts did not produce a foldable parent under the active-pocket post-ESMFold gate:

- round01 parent: 100/100 ESMFold OK, active-pocket PASS=0/100.
- round02 parent: 100/100 ESMFold OK, active-pocket PASS=0/100.
- Best round02 fail: motif CA RMSD 6.0132 A and motif pair max delta 14.3543 A.

The Baker reference parent was then tested directly:

- parent: `external/serine-hydrolase-design/design_pipeline/04_af2/outputs/super_design_model_4_ptm_seed_0_unrelaxed.pdb`
- ESMFold smoke: 1/1 OK.
- active-pocket gate: PASS, motif CA RMSD 0.2375 A, motif pair max delta 0.7515 A.

## First Baker-reference bin batch

The first Baker-reference 50-sequence batch completed ESMFold:

- ESMFold: 50/50 OK.
- Active gate counts:
  - 90: 1 PASS / 9 FAIL.
  - 80: 0 PASS / 10 FAIL.
  - 70: 0 PASS / 10 FAIL.
  - 60: 0 PASS / 10 FAIL.
  - 50: 1 PASS / 9 FAIL.

Conclusion: fixing only the 14 motif residues is not enough. The active pocket moves after ESMFold.

## Current shell8 batch

The next batch freezes a CA-distance shell around the motif:

- shell rule: fixed residues are all residues within 8 A CA distance of the 14 motif residues.
- fixed residues: 60.
- mutable distal background positions: 100.
- generated candidates: 20 per bin, 100 total.
- ESMFold run id: `baker_reference_shell8_bins_esmfold_100`.
- PID at launch: 167532.
- manifest: `/data/bht/project01_baker_serhyd_routeB_20260701/manifests/baker_reference_shell8_ligandmpnn_bins_selected.tsv`
- status: `/data/bht/project01_baker_serhyd_routeB_20260701/manifests/baker_reference_shell8_bins_esmfold_100_status.tsv`

Decision rule: after ESMFold completes, run `gate_baker_reference_bins_active.py`; only `FINAL_QUALIFIED_ACTIVE` counts toward the active-pocket target, and any bin below 10 must be expanded again.
