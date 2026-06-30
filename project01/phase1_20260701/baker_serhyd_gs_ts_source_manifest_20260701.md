# Baker Serine Hydrolase GS/TS Source Manifest

Date: 2026-07-01

Purpose: reset Project 01 serine hydrolase generation around Baker-style motif/theozyme constraints rather than a full ligand-4A pocket lock.

## Source

- Paper: Lauko et al., Science 388, eadu2454 (2025), "Computational design of serine hydrolases".
- Paper data/code reference: Zenodo DOI `10.5281/zenodo.14642126`.
- Code/example repository: `https://github.com/laukoag/serine-hydrolase-design`.
- Local external reference clone, not part of this repo:
  `work/external_refs/serine-hydrolase-design`.

## Available Baker Example Inputs

These files were located locally from the Baker example repository. They are not copied here because some are large; this manifest records the provenance and intended use.

| State / role | Local source file | Residue / ligand code | Use in our reset |
| --- | --- | --- | --- |
| Motif-generation TS-like substrate geometry | `motif_gen/01_sampling_his_stub/inputs/mu1.params` | `mu1` | Rosetta params for near-tetrahedral substrate geometry used for Ser-His theozyme sampling. |
| Motif-generation constraint geometry | `motif_gen/01_sampling_his_stub/inputs/1LNS_mu1.cst` | `mu1` + SER/HIS | Defines Ser attack on substrate carbonyl and Ser-His interaction. |
| Simple theozyme PDB | `motif_gen/02_diffusion/inputs/simple_theozyme.pdb` | `mu1`, SER, HIS | Minimal input motif for RFdiffusion motif expansion. |
| Design-stage substrate params | `design_pipeline/03_design/inputs/bn1.params` | `bn1` | Params for design-stage substrate/theozyme example. |
| Design-stage full theozyme | `design_pipeline/01_diffusion/inputs/theozyme.pdb` | `bn1`, SER, HIS, ASP, oxyanion donors | Input active-site motif for Baker-style scaffold generation. |
| Design-stage enzyme constraints | `design_pipeline/03_design/inputs/theozyme.cst` | `bn1` + SER/HIS/ASP/OXH | Defines Ser attack, oxyanion holes, Ser-His, leaving-group His, and His-Asp geometry. |
| PLACER ground/substrate-bound input | `design_pipeline/06_PLACER/inputs/super_af2_bu2.pdb` | `bu2` | Substrate-bound GS-like input for PLACER state sampling. |
| PLACER substrate-bound ensemble | `design_pipeline/06_PLACER/outputs/substrate_super_af2_bu2_model.pdb` | `bu2` | Example PLACER ensemble for substrate-bound state. |
| PLACER T1 ensemble | `design_pipeline/06_PLACER/outputs/tet1_super_af2_bu2_model.pdb` | `899` on catalytic Ser | Example tetrahedral intermediate 1 ensemble. |
| PLACER acyl-enzyme ensemble | `design_pipeline/06_PLACER/outputs/aei_super_af2_bu2_model.pdb` | `432` on catalytic Ser | Example acyl-enzyme intermediate ensemble. |
| PLACER T2 ensemble | `design_pipeline/06_PLACER/outputs/tet2_super_af2_bu2_model.pdb` | `75I` on catalytic Ser | Example tetrahedral intermediate 2 ensemble. |

## Important Interpretation

For new-backbone generation, Baker does not freeze every residue within ligand 4A. The fixed object is the active-site motif/theozyme:

- substrate or transition-state-like ligand geometry;
- catalytic Ser attack geometry;
- catalytic His geometry;
- catalytic Asp/Glu-His geometry when present;
- oxyanion-hole donor geometry.

The surrounding scaffold is generated around this motif. Ligand 4A pocket residues should be derived after scaffold generation, then used for filtering or later sequence diversification.

## Project Reset Rule

1. Fixed-backbone route: define ligand 4A residues from an existing holo structure and keep those residues fixed while diversifying the non-pocket sequence.
2. New-backbone route: start from Baker-style motif/theozyme inputs (`mu1` or `bn1`) and generate scaffolds around the motif; do not lock the full pre-existing ligand 4A shell.
3. After new-backbone generation: compute ligand 4A pocket on each generated scaffold and filter for catalytic geometry, oxyanion geometry, clashes, and pocket shape.
4. PLACER/QMMM is not the current target. Current target remains sequence/scaffold generation and geometry screening.

## Next Actions

- Use `design_pipeline/01_diffusion/inputs/theozyme.pdb` and `design_pipeline/03_design/inputs/theozyme.cst` as the first Baker-style new-backbone reference.
- Prepare a small run manifest for 90/80/70/60/50 sequence-similarity bins after generated scaffolds pass motif/pocket filters.
- Keep large PDB ensembles, model files, and trajectories off GitHub; sync only manifests, scripts, and concise status reports.
