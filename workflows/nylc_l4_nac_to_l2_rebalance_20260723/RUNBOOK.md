# RUNBOOK: five-candidate L4 NAC to L2 rebalance

## Boundaries

- Computation runs only on SCNet.
- Remote task root: `/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723`.
- GitHub branch: `codex/nylc-l4-nac-to-l2-rebalance`.
- Source trajectories and branch roots in `manifests/candidates.json` are read-only.
- Do not put GRO/XTC/TRR/TPR/CPT/EDR/large ITP files, credentials or secrets in GitHub.
- Continue later candidates after a per-candidate build or MD failure.

## Environment

```bash
source /work/home/acshdt1dks/opt/gromacs-fastest/env.sh
export GMX=/public/software/apps/Gromacs-DCU2/2022.1/mpi/bin/gmx_mpi
export TASK_ROOT=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
```

The audited L2 ITP must hash to:

```text
b0e753c60fd4b71c282d21cc6106a15e73d91d12a20d80e92dd01516162eb301
```

## Checkout and source audit

Clone the GitHub branch into `$TASK_ROOT/repo`. Do not reuse a Dell checkout.

Resolve and hash every path in `manifests/candidates.json`. For each segment, verify:

- action is FREE_MONITOR;
- source and gate restraint flags are zero at the selected cumulative time;
- the selected local XTC time exists;
- source TPR, XTC, ITP, topology and index exist;
- source reactive C/O/N is a bonded carbonyl-amide triplet.

Append START and terminal states to both `run_history.tsv` and `run_history.jsonl`.

## TDD graph builder

Set the three ITP paths and run:

```bash
export TASK_L2_ITP=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/nylc_gyaq_pa66_l2_nac_qmmm_20260723/inputs/parameterized/PA66_L2_GMX.itp
export NYLC_L4_ITP=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/apo_gate_l4_three_carbonyl_20260715/cyclic_gate_nac_20260719/branches/nylc_gyaq/C23/l4_nfree_GMX.itp
export NYL12_L4_ITP=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/apo_gate_l4_three_carbonyl_20260715/nyl12_ad_l4_j123_20260720/branches/Nyl12/J1/pa66_l4_GMX.itp
python3 -m unittest -v workflows/nylc_l4_nac_to_l2_rebalance_20260723/tests/test_l4_to_l2_graph.py
```

The first run must fail because `build_l4_to_l2.py` does not exist. After implementation, the full test command must pass.

## Build and preprocessing

Use `scripts/prepare_candidates.py` from the repository checkout. It extracts each frame into the task root, calls `build_l4_to_l2.py`, copies the required non-ligand topology includes, installs the audited L2 ITP and position-restraint includes, and writes one audit JSON per candidate.

A build is not eligible for EM until its audit confirms:

- 79 L2 atoms and 33 heavy atoms;
- audited L2 charge;
- exactly two endpoint cut edges;
- exact source C/O/N coordinate preservation within GRO precision;
- correct whole-system atom-count delta;
- finite reactive geometry and minimum nonbonded distance;
- successful grompp topology preprocessing.

## EM and staged equilibration

The common schedule is:

| Stage | Ensemble | Temperature | L2-heavy restraint | Duration |
| --- | --- | ---: | ---: | ---: |
| EM | steep | n/a | 1000 | convergence or 50,000 steps |
| nvt50 | NVT | 50 K | 1000 | 100 ps |
| nvt150 | NVT | 150 K | 500 | 100 ps |
| nvt300 | NVT | 300 K | 100 | 200 ps |
| npt300r | NPT | 300 K | 100 | 200 ps |
| npt300rel | NPT | 300 K | 10 | 200 ps |
| npt300free | NPT | 300 K | none | 1000 ps minimum |

Each candidate uses an independent SCNet job and output directory. Query `squeue` and `sinfo` before submission. Use one DCU task and modest CPU/memory per candidate unless a current cluster probe justifies a different request.

Every stage scans the log for FATAL, NaN, LINCS and SETTLE. A failure records its exact stage and leaves remaining candidate jobs independent.

## Unconstrained scientific audit

Only `npt300free` counts. Use PBC-aware analysis at a common sampling interval to report:

- distance <= 0.35 nm;
- angle 95-115 degrees;
- joint NAC occupancy;
- branch-specific validated gate opening;
- temperature and pressure statistics;
- numerical-warning scan.

Restrained results are diagnostic only. A candidate with no complete free window is NOT_EVALUATED, not inactive. DFTB3/3OB-3-1 preflight eligibility requires the unconstrained audit.

## Recovery

- Never overwrite a completed stage.
- Resume a candidate only from its last valid checkpoint and recorded stage.
- Do not resubmit merely because SSH or scheduler display is transiently unavailable.
- Re-run only the failed candidate and append a new history row.
- Preserve every audit JSON and failure reason.
