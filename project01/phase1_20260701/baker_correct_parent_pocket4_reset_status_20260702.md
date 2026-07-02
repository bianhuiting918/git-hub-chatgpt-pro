# Project 01 Baker Serine Hydrolase Correct-Parent Reset

Date: 2026-07-02

## Decision

The previous Baker-reference sequence batches are not counted toward the final target.

Reason: they used `05_ligand_copy/outputs/super_design_model_4_ptm_seed_0_unrelaxed.pdb` as the parent. That structure does not contain the intended Baker catalytic motif at the target positions:

| file | residue 15 | residue 16 | residue 55 | residue 73 | residue 149 | ligand |
|---|---:|---:|---:|---:|---:|---|
| `03_design/outputs/refined_out_1_bn1_1.pdb` | SER | LEU | HIS | ASP | LEU | `bn1`, 36 atoms |
| `04_af2/outputs/refined_out_1_bn1_1_20.pdb` | SER | LEU | HIS | ASP | LEU | none |
| `05_ligand_copy/outputs/super_design_model_4_ptm_seed_0_unrelaxed.pdb` | ARG | ARG | ARG | ARG | GLU | `bu2`, 32 atoms |

Therefore, any active-pocket pass counts from the wrong-parent `super_design` route are diagnostic only and must not be reported as Project 01 final-qualified sequences.

## Correct Parent And Pocket Definition

Use:

- Sequence-design / structure-prediction parent: `/data/bht/project01_baker_serhyd_routeB_20260701/external/serine-hydrolase-design/design_pipeline/04_af2/outputs/refined_out_1_bn1_1_20.pdb`
- Ligand pocket reference: `/data/bht/project01_baker_serhyd_routeB_20260701/external/serine-hydrolase-design/design_pipeline/03_design/outputs/refined_out_1_bn1_1.pdb`
- Catalytic residues: SER15, HIS55, ASP73
- Oxyanion residues: LEU16, LEU149
- Ligand: `bn1`

Fixed pocket policy for this reset:

- all protein residues with any heavy atom within 4.0 A of `bn1`
- plus SER15, LEU16, HIS55, ASP73, LEU149

The resulting fixed set has 15 residues:

`15,16,19,38,39,42,55,59,73,77,122,148,149,150,151`

This leaves 145 mutable positions in a 160 aa protein, so 90/80/70/60/50 percent similarity bins are all feasible.

## Final Target

The target is not intermediate generation. The target is:

| similarity bin | final-qualified sequences required |
|---:|---:|
| 90% | 10 |
| 80% | 10 |
| 70% | 10 |
| 60% | 10 |
| 50% | 10 |

Intermediate LigandMPNN sequence-bin membership is only a candidate-pool label. It is not counted as success.

## Sequence-Phase Final Gate

For the current sequence-only phase, `FINAL_QUALIFIED_ACTIVE` means:

1. Candidate is in the intended exact mutation bin.
2. Fixed 4 A ligand-pocket residues are unchanged.
3. ESMFold prediction exists.
4. Fixed-pocket CA RMSD to correct AF2 parent is <= 1.0 A.
5. Fixed-pocket pair max distance delta is <= 1.0 A.
6. Catalytic sidechain distance max delta for SER15/HIS55/ASP73 atom pairs is <= 1.0 A.

This is not a PLACER or QMMM success label. PLACER and QMMM are currently deferred by user decision.

## New Lightweight Scripts

- `scripts/generate_baker_correct_parent_pocket4_bins.py`
- `scripts/gate_baker_correct_parent_pocket4_active.py`

Remote output paths:

- candidate manifest: `/data/bht/project01_baker_serhyd_routeB_20260701/manifests/baker_correct_parent_pocket4_ligandmpnn_bins_selected.tsv`
- generation summary: `/data/bht/project01_baker_serhyd_routeB_20260701/manifests/baker_correct_parent_pocket4_ligandmpnn_bins_summary.json`
- ESMFold status: `/data/bht/project01_baker_serhyd_routeB_20260701/manifests/baker_correct_parent_pocket4_esmfold_status.tsv`
- gate output: `/data/bht/project01_baker_serhyd_routeB_20260701/manifests/baker_correct_parent_pocket4_active_gate.tsv`
- gate summary: `/data/bht/project01_baker_serhyd_routeB_20260701/manifests/baker_correct_parent_pocket4_active_gate_summary.json`

## Next Action

Run candidate generation on the GPU host, then run ESMFold with a low memory fraction. After gating, expand any bin with fewer than 10 `FINAL_QUALIFIED_ACTIVE` sequences.
