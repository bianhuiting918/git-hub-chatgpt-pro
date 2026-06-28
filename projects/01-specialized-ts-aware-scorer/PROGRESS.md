# Project 01 progress log

Updated: 2026-06-29 01:20 CST

## Corrected project direction

The working target has been corrected from the earlier GMX-CP2K-first smoke workflow to the intended Baker-style enzyme design pipeline:

1. Generate or search backbones around a serine-hydrolase active-site motif.
2. Fill/design sequences for the generated backbones.
3. Use PLACER to generate local conformer/ligand/side-chain ensembles.
4. Send selected conformers to QM/DFT or QM/MM calculations for barrier labels.
5. Train and evaluate a barrier-learning model for ranking designs.

The earlier GMX-CP2K setup remains useful later as a label-generation backend, but it is not the current front-end pipeline to run through first.

## Upload scope

Only lightweight progress summaries, scripts, manifests, and selected command logs should be committed to GitHub. Large cloned tool repositories, model weights, conda environments, raw PLACER ensembles, trajectories, wavefunction files, and DFT/QM outputs stay on the compute servers.

## GPU server setup completed

Tool work root:

```text
/data/bht/design_tools
```

Downloaded/synchronized tools:

| Tool | Server path | Version checked | Status |
| --- | --- | --- | --- |
| RFdiffusion-AA | `/data/bht/design_tools/src/rf_diffusion_all_atom` | `f913a19` | Repository and paper weights present; runtime environment still needs dependency setup/container decision. |
| ProteinMPNN | `/data/bht/design_tools/src/ProteinMPNN` | `8907e66` | Runs with existing `protrek_gpu` Python environment. |
| LigandMPNN | `/data/bht/design_tools/src/LigandMPNN` | `26ec57a` | Model parameters downloaded; local venv built; Numpy compatibility patch applied to vendored OpenFold files. |
| PLACER | `/data/bht/design_tools/src/PLACER` | `7dc5563` | Official environment created and smoke-tested on CUDA. |

PLACER environment:

```text
/data/bht/design_tools/envs/placer_env
```

PLACER environment validation:

- `torch 2.3.1`
- `dgl 2.4.0+cu121`
- `e3nn 0.5.4`
- `openbabel 3.1.0`
- `networkx 3.4.2`
- `scipy 1.14.1`
- `numpy 2.2.1`
- CUDA device visible: `NVIDIA A100 80GB PCIe`
- `run_PLACER.py --help` starts successfully and resolves default weights at `weights/PLACER_model_1.pt`.

## Smoke tests completed

### ProteinMPNN PETase seed smoke

Input seed:

```text
/data/bht/enzyme_scaffold_search_v2_gpu/seeds/structures/pdb/PET_01__6ILW.pdb
```

Command type: one-sequence CPU smoke, chain A, model `v_48_020`.

Observed output:

```text
/data/bht/design_tools/smoke/proteinmpnn_petase/seqs/PET_01__6ILW.fa
```

Result: completed successfully; generated one designed sequence for the PETase/cutinase-family seed.

### LigandMPNN PETase seed smoke

Input seed:

```text
/data/bht/enzyme_scaffold_search_v2_gpu/seeds/structures/pdb/PET_01__6ILW.pdb
```

Command type: one-batch CPU smoke using downloaded LigandMPNN `protein_mpnn` checkpoint.

Observed outputs:

```text
/data/bht/design_tools/smoke/ligandmpnn_petase/seqs/PET_01__6ILW.fa
/data/bht/design_tools/smoke/ligandmpnn_petase/backbones/PET_01__6ILW_1.pdb
```

Result: completed successfully.

### PLACER serine-hydrolase smoke

Input example:

```text
/data/bht/design_tools/src/PLACER/examples/inputs/denovo_SER_hydrolase.pdb
```

Command type: official PLACER serine-hydrolase example with `-n 1`, mutation `128A:75I`, custom residue JSON, CUDA enabled.

Observed result:

- device: `cuda:0`
- checkpoint: `/data/bht/design_tools/src/PLACER/weights/PLACER_model_1.pt`
- model variables: `5,071,802`
- output PDB: `/data/bht/design_tools/smoke/placer_denovo_ser_hydrolase/denovo_SER_hydrolase_75I_n1_model.pdb`
- output CSV: `/data/bht/design_tools/smoke/placer_denovo_ser_hydrolase/denovo_SER_hydrolase_75I_n1.csv`
- runtime: `4.84 seconds`
- exit status: `0`

Single-sample score line:

```text
fape=1.39487, lddt=0.94432, prmsd=0.90885, plddt=0.96652, plddt_pde=0.93063
```

## Current GPU state

At 2026-06-29 01:18 CST, the A100 is still heavily occupied by another user process:

```text
NVIDIA A100 80GB PCIe, about 4.6 GiB used, about 99-100% GPU utilization
```

Because of that, only one minimal PLACER smoke was run. Full PLACER conformer ensembles and RFdiffusion-AA generation should wait until the GPU is free or an agreed job window is available.

## Immediate next steps

1. Finish RFdiffusion-AA runtime setup: either build a compatible Python environment or decide on Apptainer/Singularity/container route, since Docker is not available to the user without daemon permission.
2. Prepare an active-site motif input from the existing PETase/cutinase seed set, starting from `PET_01__6ILW.pdb` where chain A catalytic residues were verified as Ser160, Asp206, His237.
3. Run a small RFdiffusion-AA backbone-generation smoke around the active-site motif when GPU is available.
4. Feed generated/scaffolded backbones to ProteinMPNN or LigandMPNN for sequence filling.
5. Use PLACER to generate a small conformer ensemble for the first designed scaffold.
6. Select a small number of conformers for DFT/QM/MM label generation; only lightweight manifests and parsed scalar labels should return to GitHub.

## Source repositories

- RFdiffusion-AA: https://github.com/baker-laboratory/rf_diffusion_all_atom
- ProteinMPNN: https://github.com/dauparas/ProteinMPNN
- LigandMPNN: https://github.com/dauparas/LigandMPNN
- PLACER: https://github.com/baker-laboratory/PLACER
