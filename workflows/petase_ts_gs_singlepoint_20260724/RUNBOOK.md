# PETase step1/step2 TS-GS DFTB3 fixed-coordinate diagnostic

## Purpose

Estimate a rapid potential-energy gap for the previously committor-validated PETase TS-core structures:

- step1 acylation: accepted g05 TS core versus three verified Michaelis-like anchors;
- step2 deacylation: accepted seed1 TS core versus four seed1 trajectories that returned to the acyl-enzyme reactant basin.

This is a **fixed-coordinate QM/MM potential-energy diagnostic**:

```
Delta E diagnostic = E(TS snapshot) - E(GS snapshot)
```

It is not an umbrella/MBAR free-energy barrier, does not include entropy or zero-point corrections, and does not prove that either snapshot is a stationary point.

## Model

- Amber `sander` from the existing `petase_stage1` environment.
- QM theory: `DFTB3`, using the installed Amber DFTB/3OB parameterization.
- QM charge: -1.
- step1: the existing 65-atom acylation QM mask, `qm_ewald=0`.
- step2: the existing deacylation QM mask including catalytic water, `qm_ewald=1`.
- The two steps use the same MM topology hash but different QM regions. Compare TS and GS only within the same step.

## Run

On the Dell CPU server:

```bash
cd /path/to/this/workflow
MAX_PARALLEL=4 bash run_ts_gs_singlepoint.sh
```

The script creates:

```
blind_work/10_energy_diagnostics/attempt155_petase_step1_step2_ts_gs_dftb3_singlepoint
```

It records inputs, SHA256 values, per-case Amber outputs, coordinate invariance, hard-error scans, `singlepoint_energies.tsv`, `summary.json`, and `run_history.tsv`.

## Gates

A reportable diagnostic requires:

1. finite energies for every selected TS/GS snapshot;
2. no Amber fatal/error, NaN, or SCC convergence failure;
3. no coordinate change between input and output;
4. a separately reported GS-snapshot spread.

Large GS spread means the raw fixed-coordinate energy gap is dominated by environmental snapshot noise. In that case, do not label the mean as an activation free energy; proceed to matched restrained relaxation or converged PMF.
