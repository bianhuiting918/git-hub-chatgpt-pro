# Project 01 progress log

Updated: 2026-06-29 07:35 CST

## Corrected project direction

The working target is the intended Baker-style enzyme design pipeline:

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
| RFdiffusion-AA | `/data/bht/design_tools/src/rf_diffusion_all_atom` | `f913a19` | Repository, paper weights, official SIF container, Apptainer env, and manual Python venv are present. Help/config loading works in `rfaa_venv`; full generation has not been launched because GPU is occupied. |
| ProteinMPNN | `/data/bht/design_tools/src/ProteinMPNN` | `8907e66` | Runs with existing `protrek_gpu` Python environment. |
| LigandMPNN | `/data/bht/design_tools/src/LigandMPNN` | `26ec57a` | Model parameters downloaded; local venv built; Numpy compatibility patch applied to vendored OpenFold files. |
| PLACER | `/data/bht/design_tools/src/PLACER` | `7dc5563` | Official environment created and smoke-tested on CUDA. |

Important RFdiffusion-AA runtime details:

- Official container copied to GPU: `/data/bht/design_tools/models/rfdiffusion_aa/rf_se3_diffusion.sif`
- Container size: about `12G`
- Container sha256 on GPU: `aba2b443f9a35a1f409eda30dbef3e810cf6b1fdb163251a0ea377d082c1dc41`
- User-space Apptainer installed at: `/data/bht/design_tools/envs/apptainer_env/bin/apptainer`
- Apptainer version: `1.5.2`
- Apptainer SIF execution currently fails without system/user-namespace permission: `Operation not permitted`
- Manual RFdiffusion-AA venv created at: `/data/bht/design_tools/envs/rfaa_venv`
- `run_inference.py --help` reaches the Hydra config page in `rfaa_venv`, after installing missing packages including `hydra-core`, `omegaconf`, `icecream`, `biopython`, `assertpy`, `fire`, `ipdb`, `rdkit`, and `deepdiff`.

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

### RFdiffusion-AA lightweight runtime smoke

Command type: import/config/help test only, no generation.

Observed result:

- `rfaa_venv` imports `torch`, `dgl`, `e3nn`, `hydra`, `omegaconf`, `icecream`, `rdkit`, `Bio`, `assertpy`, `fire`, `ipdb`.
- After adding `deepdiff`, `run_inference.py --help` reaches the Hydra configuration page.
- No RFdiffusion-AA design has been generated yet.

## 2026-06-29 morning validation

A fresh validation was run at about 07:28 CST through the CPU jump host to the GPU server.

Confirmed:

- GPU host: `dell-PowerEdge-R760`
- RFdiffusion-AA venv exists: `/data/bht/design_tools/envs/rfaa_venv`
- PLACER env exists: `/data/bht/design_tools/envs/placer_env`
- RFdiffusion-AA `run_inference.py --help` still reaches the Hydra configuration page.
- PLACER `run_PLACER.py --help` still starts successfully.

Current blocker:

```text
NVIDIA A100 80GB PCIe, about 4.6 GiB used, 100% GPU utilization
active process: /Dell/Dell9/chenchao23/miniconda3/envs/pxdesign/bin/python
```

No RFdiffusion-AA generation was launched during this validation because the GPU remains fully occupied. After the successful check, subsequent SSH attempts from the local machine to the CPU jump host briefly failed before handshake, so no further remote commands were forced.

## Immediate next steps

1. Prepare an active-site motif input from the existing PETase/cutinase seed set, starting from `PET_01__6ILW.pdb` where chain A catalytic residues were verified as Ser160, Asp206, His237.
2. Run a tiny RFdiffusion-AA backbone-generation smoke around that motif using `/data/bht/design_tools/envs/rfaa_venv` when GPU is free.
3. Feed generated/scaffolded backbones to ProteinMPNN or LigandMPNN for sequence filling.
4. Use PLACER to generate a small conformer ensemble for the first designed scaffold.
5. Select a small number of conformers for DFT/QM/MM label generation; only lightweight manifests and parsed scalar labels should return to GitHub.
6. If Apptainer container execution is preferred, ask the server administrator to enable user namespaces or install a system/setuid Apptainer/Singularity runtime.

## Source repositories

- RFdiffusion-AA: https://github.com/baker-laboratory/rf_diffusion_all_atom
- ProteinMPNN: https://github.com/dauparas/ProteinMPNN
- LigandMPNN: https://github.com/dauparas/LigandMPNN
- PLACER: https://github.com/baker-laboratory/PLACER
