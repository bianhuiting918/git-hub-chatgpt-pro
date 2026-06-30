# Sugon SCNet GROMACS+CP2K QM/MM Entry

This note records the login entry and installed GROMACS+CP2K environment for later QM/MM work on SCNet/Sugon.

Do not commit the SCNet password to this repository. Use the user's secure channel, a local password manager, or an environment variable when automation needs credentials.

## SCNet Platform

- Platform URL: <https://www.scnet.cn/>
- SCNet username: `shuixingzhou`
- Password: not stored in git

Available clusters checked previously:

| Region | Cluster ID | Cluster user | Login node | Notes |
| --- | ---: | --- | --- | --- |
| Northwest 1, Xi'an | `11285` | `acshdt1dks` | `login01` | Use this cluster for GROMACS+CP2K work. Partition `xahdnormal`; accelerator resources are Hygon DCU, not NVIDIA GPU. |
| East China 1, Kunshan | `11250` | `shuixingzhou` | unknown | Previously no normal usable Slurm partition was visible. |

The account observed earlier was bound to team account `acshdt1dks` with account ID `606544861843161088`. Balance and resource counts are time-sensitive; refresh them before submitting expensive jobs.

## Installed Environment

Purpose: GROMACS 2023.1 linked with CP2K 2024.1 for QM/MM.

Build characteristics:

- GROMACS: `2023.1`
- CP2K: `2024.1`
- Precision: double
- MPI: enabled
- OpenMP: enabled
- GROMACS GPU support: disabled
- Intended cluster: Northwest 1, Xi'an

Important paths:

```bash
/work/home/acshdt1dks/opt/gmx-cp2k
/work/home/acshdt1dks/opt/gmx-cp2k/env.sh
/work/home/acshdt1dks/opt/gmx-cp2k/gromacs-2023.1-cp2k-2024.1/bin/gmx_mpi_d
/work/home/acshdt1dks/src/gmx-cp2k/cp2k-2024.1/lib/local/psmp/libcp2k.a
/work/home/acshdt1dks/src/gmx-cp2k
```

The installed environment script sets:

```bash
export GMXBIN=gmx_mpi_d
export CP2K_DATA_DIR=/work/home/acshdt1dks/src/gmx-cp2k/cp2k-2024.1/data
```

## Basic Use

After logging into the Northwest 1 cluster:

```bash
source /work/home/acshdt1dks/opt/gmx-cp2k/env.sh
$GMXBIN --version
```

Prepare a QM/MM run with the CP2K input file:

```bash
source /work/home/acshdt1dks/opt/gmx-cp2k/env.sh
$GMXBIN grompp -f mdp.mdp -c conf.gro -p topol.top -qmi topol-qmmm.inp -o topol.tpr
mpirun -np <N> $GMXBIN mdrun -s topol.tpr -deffnm qmmm
```

Direct executable form:

```bash
gmx_mpi_d grompp ...
gmx_mpi_d mdrun ...
```

## Verification Commands

Use these first when picking up this environment:

```bash
source /work/home/acshdt1dks/opt/gmx-cp2k/env.sh
which "$GMXBIN"
$GMXBIN --version
$GMXBIN grompp -h | grep -i -E 'qmmm|cp2k|qm|qmi'
ldd "$(command -v "$GMXBIN")" | grep 'not found' || echo "no missing shared libraries"
```

Previous verification showed:

- `GMX_CP2K:BOOL=ON` in the GROMACS CMake cache.
- `$GMXBIN --version` runs and reports GROMACS `2023.1`, double precision, MPI, OpenMP.
- `grompp -h` exposes `-qmi [<.inp>]`, the QM input option.
- `ldd` reported no missing shared libraries.

## Slurm Checks

```bash
sinfo -p xahdnormal
squeue -u acshdt1dks
```

Start with a small QM/MM test system before running expensive production jobs.

## Build Notes

The source and build root is:

```bash
/work/home/acshdt1dks/src/gmx-cp2k
```

Build script:

```bash
/work/home/acshdt1dks/src/gmx-cp2k/build.sh
```

Known build details:

- The build filters user Anaconda paths to avoid contaminating compiler/linker state.
- Because of that, the script loads `python/3.8.10` explicitly.
- CP2K toolchain `install/setup` is not compatible with `set -u`; source it with `set +u` temporarily.
- CP2K toolchain arch files must be copied from `tools/toolchain/install/arch/` into the CP2K source `arch/` directory.
- CP2K `libcp2k.pc` needed rpath cleanup so OpenMPI paths are passed as `-Wl,-rpath,/path`, not as bare directory inputs.
- Because this is a double-precision GROMACS build, the executable is `gmx_mpi_d`, not `gmx_mpi`.
