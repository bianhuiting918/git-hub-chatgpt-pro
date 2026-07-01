# 2026-07-01 rank021 2D scan sanity path and manual mid-proton diagnostic

## Purpose

After repeated single-endpoint His-mediated TS-seed attempts failed to transfer Ser HG to HID NE2, I tested whether a small relaxed path/scan could reveal a continuous acylation path. This remains a blind first-principles diagnostic, not paper-coordinate replication and not umbrella sampling.

## 2D relaxed QMMM scan scaffold

Directory:

`blind_work/04_qmmm_rank021_reactant_prep/hid237_qmmm_2d_relaxed_scan/`

Coordinates:

- x: Ser132 OG to ligand C005 distance.
- y: approximate proton-transfer coordinate via paired SerH-HisNE2 and SerH-SerOG restraints.

Three sanity-path nodes were run sequentially:

- `n00_reactant`
- `n01_mid`
- `n02_late`

Geometry summary is saved in:

`blind_work/04_qmmm_rank021_reactant_prep/hid237_qmmm_2d_relaxed_scan/sanity_path_geometry.tsv`

| node | SCC fail | DFTBESCF final | restraint | SerO-C | SerH-His | SerH-SerO | PTdiff His-O | O-H-N | attack O006 | attack O004 | O006-MetN | His-O004 | C-O004 |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| start | NA | NA | NA | 3.221 | 1.872 | 1.012 | 0.860 | 166.855 | 67.461 | 106.729 | 3.332 | 4.946 | 1.372 |
| n00 | no | -5107.4087 | 0.7489 | 3.302 | 1.929 | 0.936 | 0.994 | 158.062 | 63.467 | 104.037 | 3.338 | 4.959 | 1.412 |
| n01 | no | 5500.0533 | 3.5628 | 2.970 | 1.829 | 0.992 | 0.836 | 160.655 | 61.969 | 109.570 | 3.345 | 4.878 | 1.399 |
| n02 | no | 5458.8224 | 12.2373 | 2.975 | 1.829 | 0.994 | 0.834 | 160.491 | 61.789 | 109.467 | 3.343 | 4.877 | 1.397 |

Interpretation:

- n01/n02 are not reliable: DFTBESCF became very large positive despite no explicit SCC-failure string.
- Proton transfer coordinate did not progress; SerH-SerO stayed near 0.99 A.
- The 3-node path should not be expanded to 3x3 in this form.

## Manual mid-proton string image

To test whether the minimizer simply could not pull the proton out of the SerO-H well, I manually constructed a mid-proton image from our own rank021 step4 coordinate, placing SerH along the SerO-HisNE2 line. This is not from the paper.

Directory:

`blind_work/04_qmmm_rank021_reactant_prep/hid237_qmmm_string_seed/manual_mid_from_step4/`

Geometry summary is saved in:

`blind_work/04_qmmm_rank021_reactant_prep/hid237_qmmm_string_seed/manual_mid_from_step4/manual_mid_geometry.tsv`

| state | SCC fail | DFTBESCF final | restraint | SerO-C | SerH-His | SerH-SerO | PTdiff His-O | O-H-N | attack O006 | attack O004 | O006-MetN | His-O004 | C-O004 |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| step4_base | NA | NA | NA | 2.768 | 1.548 | 1.055 | 0.493 | 164.818 | 63.644 | 97.583 | 3.147 | 4.628 | 1.353 |
| manual_mid_initial | NA | NA | NA | 2.768 | 1.232 | 1.350 | -0.118 | 180.000 | 63.644 | 97.583 | 3.147 | 4.628 | 1.353 |
| manual_mid_10step | no | -5129.7114 | 1.8193 | 2.939 | 1.444 | 1.030 | 0.414 | 168.074 | 64.094 | 94.028 | 3.143 | 4.676 | 1.370 |

Interpretation:

- A deliberately constructed mid-proton image collapses back toward SerO-H in only 10 QMMM steps.
- This indicates rank021, as currently prepared, does not stabilize a SerH-to-His transfer image.

## Ligand atom assignment check

Checked `ligand.mol2` connectivity in rank021 HID topology prep:

- C005-O006 bond order 2, O006 type `o`: O006 is the carbonyl oxygen.
- C005-O004 bond order 1, O004 type `os`, O004-C003 bond order 1: O004 is the ester/leaving-side oxygen.

Therefore the persistent poor `attack O006` angle is real, not an O004/O006 naming swap.

## Revised root-cause hypothesis

The likely issue is upstream preparation: the HID rebuild used a topology seed that may not preserve the best top5 MM-relaxed stage300 coordinates. The top5 MM-relaxed rank022 had a much better O006 attack angle before HID rebuild; after HID/QMMM preparation, the O006 attack geometry deteriorated.

Next action:

- Build a HID topology while preserving the top5 MM-relaxed coordinates, starting with rank021/rank022 stage300 structures.
- Do not continue expanding the current rank021 scan.
- Do not start umbrella sampling.
