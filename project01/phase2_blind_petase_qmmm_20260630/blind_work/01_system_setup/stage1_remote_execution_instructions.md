# Blind PETase Stage 1 Remote Execution Instructions

Date: 2026-06-30

## Boundary

These instructions are for enabling Stage 1 protein protonation and ligand 3D construction. They do not use paper TS coordinates, paper reaction-coordinate terms, selected CVs, paper trajectories, barriers, or paper mechanism conclusions.

## Current Local State

The GitHub repository already contains:

- structure selection and RCSB download manifest;
- cleaned-structure manifest with SHA256 hashes;
- altloc resolution decisions;
- protonation site scan and hypothesis manifest;
- ligand model definitions and atom-label requirements;
- scripts to reproduce the structure audit and initial cleaned PDB manifests.

The actual generated cleaned PDB coordinate files are local/intermediate artifacts and are reproducible from the RCSB downloads and scripts.

## Remote Login Constraint

Non-interactive key-based SSH was tested with:

```text
ssh -o BatchMode=yes -o ConnectTimeout=10 bianht@210.73.40.29 ...
```

Result: `Permission denied (publickey,password)`.

Do not put passwords in shell commands, scripts, GitHub files, logs, or terminal history. Use an interactive SSH session or configure key-based login.

## Step 1 - Clone Or Update Repo

On the compute server:

```bash
mkdir -p /Dell/Dell14/bianht/petase_blind_qmmm
cd /Dell/Dell14/bianht/petase_blind_qmmm
git clone https://github.com/bianhuiting918/git-hub-chatgpt-pro.git repo
cd repo
```

If the repo already exists:

```bash
cd /Dell/Dell14/bianht/petase_blind_qmmm/repo
git pull --ff-only
```

## Step 2 - Probe Tooling

Run:

```bash
bash project01/phase2_blind_petase_qmmm_20260630/scripts/probe_stage1_compute_environment.sh \
  project01/phase2_blind_petase_qmmm_20260630/blind_work/01_system_setup/environment_probe
```

Expected output files:

```text
project01/phase2_blind_petase_qmmm_20260630/blind_work/01_system_setup/environment_probe/stage1_compute_environment_probe.tsv
project01/phase2_blind_petase_qmmm_20260630/blind_work/01_system_setup/environment_probe/stage1_compute_environment_probe.md
```

## Step 3 - Minimum Tool Gate

Proceed to ligand/protonation execution only if the probe finds at least:

- one ligand builder: RDKit Python module or Open Babel `obabel`;
- one protonation/pKa route: PROPKA, pdb2pqr, H++, Amber reduce, or equivalent;
- one force-field preparation route: AmberTools or an equivalent stack;
- later for QM/MM: Amber/Sander with DFTB3/MM, plus AmberTools topology preparation.

If this gate fails, create an environment such as:

```bash
conda create -n petase_stage1 -c conda-forge rdkit openbabel propka pdb2pqr ambertools biopython -y
conda activate petase_stage1
```

Record the exact command and package versions in the environment probe report.

## Step 4 - Reproduce Initial Structures

Download RCSB structures:

```bash
powershell -ExecutionPolicy Bypass -File project01/phase2_blind_petase_qmmm_20260630/scripts/download_stage1_rcsb_structures.ps1
```

If PowerShell is unavailable on Linux, download the PDB IDs listed in `structure_download_manifest.csv` with `wget` or `curl` into:

```text
work/blind_work/01_system_setup/structures/
```

Then run:

```bash
python3 project01/phase2_blind_petase_qmmm_20260630/scripts/prepare_stage1_initial_structures_v2.py \
  work/blind_work/01_system_setup/structures \
  project01/phase2_blind_petase_qmmm_20260630/blind_work/01_system_setup/retained_water_candidates.tsv \
  work/blind_work/01_system_setup
```

Compare the generated `prepared_structure_manifest.tsv` SHA256 values against the GitHub manifest before continuing.

## Step 5 - Protonation Production Gate

Run a pKa/protonation workflow on `6EQE_chainA_initial_clean_v2.pdb`.

Required record:

- input file SHA256;
- tool name/version;
- command line;
- output file SHA256;
- catalytic `ASP206` state;
- catalytic `HIS237` tautomer/protonation;
- remote His tautomer assignments;
- any changes relative to `protonation_hypothesis_manifest.tsv`.

If catalytic `ASP206` or `HIS237` differs from the primary hypothesis, keep both branches for sensitivity testing.

## Step 6 - Ligand 3D Build Gate

Generate 3D structures for:

- `PET_dimer_capped`;
- `BHET_like`;
- `MHET_like`;
- optional neutral `MHET_like` sensitivity variant only if pKa/protonation logic requires it.

Required record:

- input SMILES;
- tool/version;
- conformer generation command;
- output SDF/MOL2 path;
- atom-label table mapping `Csc_acyl`, `Oox_acyl`, `Olg_ethylene_glycol`, etc. to concrete atom names;
- validation that atom names survive topology conversion.

Do not begin docking until this atom-label gate is satisfied.

