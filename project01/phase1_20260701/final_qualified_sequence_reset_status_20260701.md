# Project 01 final-qualified sequence reset status

Date: 2026-07-01

## Final target

The target is reset to final-qualified sequence counts, not intermediate sequence-bin counts.

Required final panel:

| Similarity bin | Required FINAL_QUALIFIED count | Current FINAL_QUALIFIED count |
|---|---:|---:|
| 90% | 10 | 0 |
| 80% | 10 | 0 |
| 70% | 10 | 0 |
| 60% | 10 | 0 |
| 50% | 10 | 0 |

`FINAL_QUALIFIED` means:

1. sequence is in the requested similarity bin;
2. post-design structure prediction finished successfully;
3. active-pocket motif CA gate passes;
4. sidechain/ligand-compatibility gate passes before any downstream QMMM use.

Intermediate LigandMPNN sequence-bin membership is only a candidate-pool label and is not counted as success.

## Hard evidence so far

Round01 parent search:

- ESMFold completed: 100/100 parent candidates.
- Parent active-pocket gate: 0 PASS / 100 FAIL.
- Best failing parent: motif CA RMSD 5.8851 A, motif pair max delta 13.7388 A.
- Decision: round01 parent candidates are not valid parents for 90/80/70/60/50 sequence-bin generation.

Round02 parent search:

- Generated 100 new parent candidates using diverse LigandMPNN strategies:
  - soluble_mpnn with anti-low-complexity AA bias: 30 candidates;
  - protein_mpnn with anti-low-complexity AA bias: 30 candidates;
  - ligand_mpnn with anti-low-complexity AA bias: 40 candidates.
- ESMFold was relaunched with a compatible manifest after fixing column aliases required by the existing driver.
- Gate remains unchanged: motif CA RMSD <= 1.0 A and motif pair max delta <= 1.0 A.

## Next decision rule

If round02 produces at least one parent PASS, choose the best parent by motif CA RMSD and only then generate 90/80/70/60/50 similarity-bin candidates around that parent.

If round02 produces zero parent PASS, do not enter similarity-bin generation. Generate another parent batch with a stronger fold-preservation strategy.
