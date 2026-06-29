# Project 01 progress log

Updated: 2026-06-29 13:05 CST

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

No RFdiffusion-AA generation was launched during this validation because the GPU remains fully occupied.

## PLACER n=5 conformer transfer to CPU

At about 2026-06-29 09:13 CST, a very small PLACER conformer ensemble was generated on the GPU and transferred to the CPU server to validate the GPU-to-CPU handoff path.

Input:

```text
/data/bht/design_tools/src/PLACER/examples/inputs/denovo_SER_hydrolase.pdb
```

GPU output directory:

```text
/data/bht/design_tools/runs/placer_denovo_ser_hydrolase_n5_20260629T011253
```

CPU queue directory:

```text
/Dell/Dell14/bianht/project01_conformer_queue/placer_denovo_ser_hydrolase_n5_20260629T011253
```

Transfer archive checksum:

```text
12418688982dafb4f4e69881604ef27b85f92314a62d4be1b68029299286a8a4
```

Generation summary:

- PLACER samples: `5`
- GPU runtime reported by PLACER: `16.08 seconds`
- Output PDB size: about `230K`
- Output CSV size: `721 bytes`
- CPU split models: `5` files under `split_models/model_001.pdb` to `model_005.pdb`
- CPU manifest: `manifest.tsv`

Scores sorted by ascending `prmsd`:

| model_idx | prmsd | plddt | recommended_for_qm |
| --- | --- | --- | --- |
| 4 | 0.8627644777 | 0.9733046293 | yes |
| 1 | 0.8691723943 | 0.9810498953 | yes |
| 5 | 0.8789390326 | 0.9692426920 | yes |
| 2 | 0.8951162100 | 0.9736832380 | no |
| 3 | 0.9116663933 | 0.9686573148 | no |

Interpretation:

This validates the chain `GPU PLACER conformer generation -> archive/checksum -> CPU receive -> model splitting -> manifest/ranking`. It is still not a production DFT/barrier run. Before launching DFT/QM/MM, the next required CPU preparation is to define QM region, charge, multiplicity, protonation state, substrate/TS state, and CP2K/DFT input templates for the selected conformers.

## 2026-06-29 PLACER crop back-mapping and no-min QMMM trial

User clarified that the PLACER outputs should be interpreted as active-site conformer crops, not complete proteins. The full PLACER input was recovered from the GPU examples directory and copied to CPU:

```text
GPU: /data/bht/design_tools/src/PLACER/examples/inputs/denovo_SER_hydrolase.pdb
CPU: /Dell/Dell14/bianht/project01_conformer_queue/placer_denovo_ser_hydrolase_n5_20260629T011253/inputs/denovo_SER_hydrolase_full_input.pdb
```

Checksums on CPU:

```text
8615932565dfd91deed3f7aa8c59f9126e69874785f65ccbe7ce34861c6b7fa2  denovo_SER_hydrolase_full_input.pdb
7d9e2bd92efbd4725383902c424b6d4e403596e83830f2649804835295696d3e  75I.json
```

Full-length no-min mapping directory:

```text
/Dell/Dell14/bianht/project01_conformer_queue/placer_denovo_ser_hydrolase_n5_20260629T011253/full_length_mapped_no_min_20260629_1255
```

Generated for recommended PLACER conformers `model_004`, `model_001`, and `model_005`:

- `full_mapped_no_min.pdb`: full protein coordinates, with PLACER crop coordinates replacing matching active-site atoms.
- `protein_only_mapped_no_min.pdb`: protein-only input for `pdb2gmx`.
- `qex_extra_75I.xyz`: six extra heavy atoms from PLACER 75I reaction fragment.
- `mapping_report.tsv`: atom counts and key distances.

Mapping summary:

| model | protein atoms written | PLACER-replaced protein atoms | QEX atoms | Ser128 OG-C1A A | Ser128 OG-O17 A | Ser126 OG-OAC A |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 004 | 1246 | 575 | 6 | 1.4279 | 2.2505 | 3.2390 |
| 001 | 1246 | 575 | 6 | 1.3311 | 2.3061 | 7.0873 |
| 005 | 1246 | 575 | 6 | 1.4260 | 2.3058 | 3.4526 |

`pdb2gmx` validation succeeded for all three full-length protein-only structures with `amber99sb-ildn` and no water:

```text
160 residues, 2494 atoms after hydrogen addition, total charge -10 e
```

For `model_004`, a no-min GMX preparation was built by adding the six PLACER QEX atoms as an extra GROMACS molecule and centering the system in an 8 nm vacuum box. MM-only one-step steepest descent completed successfully:

```text
Potential Energy = 3.75035417350657e+03 kJ/mol
Maximum force    = 5.15069605455416e+05 on atom 2495
```

GMX-CP2K QMMM smoke status:

- `grompp` succeeds and reports `Number of MM atoms=2494; Number of QM atoms=6`.
- QEX-only QMMM `mdrun` aborts during force evaluation with MPI_ABORT before a CP2K output file is written.
- The same topology without QMMM runs successfully, so the current blocker is the temporary incomplete QEX/QMMM interface, not the full-length back-mapped protein or GROMACS topology.
- PETase successful examples use a complete minimal ligand molecule (`LG3`/`LG4`, 23 atoms including H) as the QM group and have `Bonds removed=0; F_CONNBONDS added=0`. The current QEX model is only six heavy atoms from the 75I reaction fragment and is not yet a complete 75I ligand topology.

Interpretation:

The safe no-min mapping workflow is now validated up to full-length PDB generation and GROMACS MM topology. For production GMX-CP2K labels, the next required step is to build a complete 75I/TS-state ligand representation with hydrogens and a minimal GROMACS topology analogous to PETase `LG3`/`LG4`, then use that full ligand as the QM group. Until that is done, direct QMMM energies from the six-heavy-atom QEX smoke should not be interpreted scientifically.

## Immediate next steps

1. Build a complete 75I ligand/topology from `75I.json`, including hydrogens and the correct atom ordering, instead of the temporary six-heavy-atom QEX molecule.
2. Re-run the `model_004` no-min GMX-CP2K smoke with complete 75I as the QM group.
3. If the smoke succeeds, repeat for `model_001` and `model_005` and parse GROMACS `Potential Energy` only as a technical label candidate.
4. Keep all active-site heavy-atom coordinates from PLACER fixed as the reference geometry; avoid using minimized output coordinates unless a separate geometry-drift report shows the catalytic distances are unchanged.
5. Prepare a tiny RFdiffusion-AA backbone-generation smoke around the catalytic motif when GPU is free, then use ProteinMPNN/LigandMPNN and PLACER on the generated scaffold.
6. If Apptainer container execution is preferred, ask the server administrator to enable user namespaces or install a system/setuid Apptainer/Singularity runtime.

## Source repositories

- RFdiffusion-AA: https://github.com/baker-laboratory/rf_diffusion_all_atom
- ProteinMPNN: https://github.com/dauparas/ProteinMPNN
- LigandMPNN: https://github.com/dauparas/LigandMPNN
- PLACER: https://github.com/baker-laboratory/PLACER
## 2026-06-29 complete 75I QMMM/cluster batch status

The current real-compute scale is one design/mutation setting, `128A:75I`, with 5 PLACER conformers. The first-pass QM/MM comparison batch uses 3 conformers: `model_004`, `model_001`, and `model_005`.

`model_004` now has a completed with-MM QMMM step-0 energy and a no-MM cluster control energy. `model_001` and `model_005` have complete 75I reconstructions and validated QMMM `.tpr` files. Their no-MM cluster jobs are running, and `model_001` with-MM QMMM is running. See `runs/2026-06-29-placer-128a-75i-qmmm-cluster.md`.