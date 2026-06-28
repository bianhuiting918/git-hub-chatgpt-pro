# Project 01 runbook

Updated: 2026-06-29 02:18 CST

This runbook records lightweight, reproducible commands for the current enzyme-design pipeline. It intentionally excludes passwords, raw model weights, cloned repositories, conda environments, and generated structure ensembles.

## Server paths

```text
Tool root: /data/bht/design_tools
Seed structures: /data/bht/enzyme_scaffold_search_v2_gpu/seeds/structures
RFdiffusion-AA venv: /data/bht/design_tools/envs/rfaa_venv
RFdiffusion-AA SIF: /data/bht/design_tools/models/rfdiffusion_aa/rf_se3_diffusion.sif
Apptainer env: /data/bht/design_tools/envs/apptainer_env
PLACER env: /data/bht/design_tools/envs/placer_env
LigandMPNN env: /data/bht/design_tools/envs/ligandmpnn_venv
Existing GPU Python env: /data/bht/miniforge3/envs/protrek_gpu
```

## Verify tool versions

```bash
for d in rf_diffusion_all_atom ProteinMPNN LigandMPNN PLACER; do
  cd /data/bht/design_tools/src/$d
  printf '%s\t' "$d"
  git rev-parse --short HEAD
done
```

Observed on 2026-06-29:

```text
rf_diffusion_all_atom   f913a19
ProteinMPNN             8907e66
LigandMPNN              26ec57a
PLACER                  7dc5563
```

## RFdiffusion-AA container and manual venv

Official SIF container is present and checksum-matched:

```bash
cd /data/bht/design_tools/models/rfdiffusion_aa
ls -lh rf_se3_diffusion.sif
sha256sum rf_se3_diffusion.sif
```

Observed sha256:

```text
aba2b443f9a35a1f409eda30dbef3e810cf6b1fdb163251a0ea377d082c1dc41  rf_se3_diffusion.sif
```

User-space Apptainer is installed:

```bash
/data/bht/design_tools/envs/apptainer_env/bin/apptainer --version
/data/bht/design_tools/envs/apptainer_env/bin/apptainer buildcfg | grep APPTAINER_SUID_INSTALL
```

Observed:

```text
apptainer version 1.5.2
APPTAINER_SUID_INSTALL=0
```

Current Apptainer limitation:

```bash
/data/bht/design_tools/envs/apptainer_env/bin/apptainer exec rf_se3_diffusion.sif python -V
```

Observed failure:

```text
FATAL: while extracting rf_se3_diffusion.sif: root filesystem extraction failed: Operation not permitted
```

So the practical route for now is the manual venv:

```bash
cd /data/bht/design_tools/src/rf_diffusion_all_atom
PYTHONPATH=/data/bht/design_tools/src/rf_diffusion_all_atom:/data/bht/design_tools/src/rf_diffusion_all_atom/lib/rf2aa \
  /data/bht/design_tools/envs/rfaa_venv/bin/python run_inference.py --help | sed -n '1,80p'
```

Expected result: Hydra configuration help is printed. This is an import/config smoke only; it does not generate a backbone.

## ProteinMPNN PETase seed smoke

```bash
cd /data/bht/design_tools/src
CUDA_VISIBLE_DEVICES='' /data/bht/miniforge3/envs/protrek_gpu/bin/python ProteinMPNN/protein_mpnn_run.py \
  --pdb_path /data/bht/enzyme_scaffold_search_v2_gpu/seeds/structures/pdb/PET_01__6ILW.pdb \
  --pdb_path_chains A \
  --path_to_model_weights ProteinMPNN/vanilla_model_weights \
  --model_name v_48_020 \
  --num_seq_per_target 1 \
  --sampling_temp 0.1 \
  --batch_size 1 \
  --seed 7 \
  --out_folder /data/bht/design_tools/smoke/proteinmpnn_petase
```

Expected output:

```text
/data/bht/design_tools/smoke/proteinmpnn_petase/seqs/PET_01__6ILW.fa
```

## LigandMPNN PETase seed smoke

```bash
cd /data/bht/design_tools/src/LigandMPNN
CUDA_VISIBLE_DEVICES='' /data/bht/design_tools/envs/ligandmpnn_venv/bin/python run.py \
  --model_type protein_mpnn \
  --checkpoint_protein_mpnn model_params/proteinmpnn_v_48_020.pt \
  --pdb_path /data/bht/enzyme_scaffold_search_v2_gpu/seeds/structures/pdb/PET_01__6ILW.pdb \
  --out_folder /data/bht/design_tools/smoke/ligandmpnn_petase \
  --seed 7 \
  --batch_size 1 \
  --number_of_batches 1 \
  --temperature 0.1 \
  --verbose 1
```

Expected outputs:

```text
/data/bht/design_tools/smoke/ligandmpnn_petase/seqs/PET_01__6ILW.fa
/data/bht/design_tools/smoke/ligandmpnn_petase/backbones/PET_01__6ILW_1.pdb
```

## PLACER import and CLI validation

```bash
/data/bht/design_tools/envs/placer_env/bin/python - <<'PY'
import torch, dgl, e3nn, openbabel, networkx, scipy, numpy
print('torch', torch.__version__)
print('dgl', dgl.__version__)
print('e3nn', e3nn.__version__)
print('openbabel', openbabel.__version__)
print('cuda_available', torch.cuda.is_available())
if torch.cuda.is_available():
    print('cuda_device', torch.cuda.get_device_name(0))
PY

cd /data/bht/design_tools/src/PLACER
/data/bht/design_tools/envs/placer_env/bin/python run_PLACER.py --help | sed -n '1,80p'
```

## PLACER serine-hydrolase one-sample smoke

Run only when it is acceptable to briefly use the GPU.

```bash
mkdir -p /data/bht/design_tools/smoke/placer_denovo_ser_hydrolase
cd /data/bht/design_tools/src/PLACER
CUDA_VISIBLE_DEVICES=0 timeout 900 /data/bht/design_tools/envs/placer_env/bin/python run_PLACER.py \
  --ifile examples/inputs/denovo_SER_hydrolase.pdb \
  --odir /data/bht/design_tools/smoke/placer_denovo_ser_hydrolase \
  --suffix 75I_n1 \
  -n 1 \
  --mutate 128A:75I \
  --residue_json examples/ligands/75I.json \
  --no-use_sm
```

Observed output on 2026-06-29:

```text
/data/bht/design_tools/smoke/placer_denovo_ser_hydrolase/denovo_SER_hydrolase_75I_n1_model.pdb
/data/bht/design_tools/smoke/placer_denovo_ser_hydrolase/denovo_SER_hydrolase_75I_n1.csv
```

Observed one-sample score:

```text
fape=1.39487, lddt=0.94432, prmsd=0.90885, plddt=0.96652, plddt_pde=0.93063
```

## GPU occupancy check

```bash
nvidia-smi --query-gpu=index,name,memory.used,utilization.gpu --format=csv,noheader
nvidia-smi --query-compute-apps=pid,process_name,used_gpu_memory --format=csv,noheader 2>/dev/null || true
```

Do not launch full RFdiffusion-AA generation or PLACER ensembles while the GPU is at sustained high utilization from another user's process.

## Next runnable target

When the A100 is free, the next practical target is:

1. Extract a compact catalytic motif from PETase/cutinase seed `PET_01__6ILW.pdb` using Ser160, Asp206, His237 on chain A.
2. Run a tiny RFdiffusion-AA backbone generation smoke with `/data/bht/design_tools/envs/rfaa_venv`.
3. Sequence-fill generated backbones with ProteinMPNN or LigandMPNN.
4. Run PLACER on a small designed-scaffold subset before selecting conformers for DFT/QM/MM labels.
