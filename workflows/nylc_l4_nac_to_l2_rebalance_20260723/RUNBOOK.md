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
- account for the segment boundary: the first geometry row is the 2 ps XTC sample, so segment start is `first_geometry_time - sample_interval`; never subtract the first geometry time directly;
- recomputed Thr OG1--carbonyl C distance and the Thr OG1--C versus C--O angle agree with the selected `geometry.tsv` row within XTC precision;
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

First run `scripts/extract_sources.py`; use `--force` only for an audited refresh of derived GRO files, never source trajectories. Then run `scripts/prepare_candidates.py`. The latter calls `build_l4_to_l2.py`, copies the required non-ligand topology includes, installs the audited L2 ITP and position-restraint includes, and writes one audit JSON per candidate.

The 2026-07-23 offset audit established corrected XTC local times 14, 44, 52, 28 and 98 ps for NylC-C18, NylC-C23, Nyl50-C18, Nyl12-J1 and Nyl12-J2 respectively. The superseded minus-2-ps derived inputs/builds are preserved under `$TASK_ROOT/superseded/incorrect_minus2ps_20260723`.

A build is not eligible for EM until its audit confirms:

- 79 L2 atoms and 33 heavy atoms;
- audited L2 charge;
- exactly two endpoint cut edges;
- exact source C/O/N coordinate preservation within GRO precision;
- correct whole-system atom-count delta;
- finite reactive geometry and minimum nonbonded distance;
- successful grompp topology preprocessing.

## NylC step1 GS selection and staged release

The original `Fmax <= 500` rule is not a scientific gate. Authentic finite-temperature source frames can have larger instantaneous forces, and over-minimization can move a selected NAC out of its geometric basin. Contact relaxation is acceptable when it finishes, produces coordinates and has no FATAL, NaN, LINCS or SETTLE event; always record, but do not threshold, its maximum force.

Run the NylC-only diagnostic pilot with:

```bash
sbatch workflows/nylc_l4_nac_to_l2_rebalance_20260723/slurm/run_nylc_step1_pilot_array.sbatch
```

Job `61686301` used the corrected C18 and C23 inputs. Both completed restrained contact relaxation, 2 ps at 10 K and 10 ps at 50 K without numerical warnings. Full-trajectory strict NAC statistics were:

| candidate | 10 K occupancy | 50 K occupancy | 50 K longest residence |
| --- | ---: | ---: | ---: |
| NylC-C18 | 0.5025 | 0.3433 | 0.45 ps |
| NylC-C23 | 0.0448 | 0.0000 | 0 ps |

These are restrained diagnostics, not scientific PASS. C23 is retained with its failure evidence but is not extended in the current NylC-first route.

For C18, select only from frames satisfying distance <= 0.35 nm and angle 95-115 degrees. Exclude the initial heating transient when ranking by potential energy. The selected late-window frame is:

```text
50 K trajectory time: 9.55 ps
distance: 0.337 nm
angle: 110.327 degrees
SHA256: c67cbb1d275863606be62628df6829e8ef15fbad39f3b5a51188c311f95ce235
remote file: candidates/nylc_c18_11854ps/selected_step1_gs/c18_50k_late_lowest_potential_nac_9p55ps.gro
```

Run the staged continuation with:

```bash
sbatch workflows/nylc_l4_nac_to_l2_rebalance_20260723/slurm/run_nylc_c18_step1_continuation.sbatch
```

Job `61687591` performs 100 K with weak protein/L2 restraints, 150 K with weaker L2 restraints, 300 K with only 10 kJ mol-1 nm-2 L2 restraints, then a 100 ps fully unrestrained NPT pilot. The 100 ps pilot is not the required final 1 ns window. Audit its complete trajectory before extension.

## Unconstrained scientific audit

Only fully unrestrained NPT counts as scientific NAC evidence. Use `scripts/analyze_nac_series.py` with PBC-aware GROMACS distance/angle series at a common sampling interval and report:

- distance <= 0.35 nm;
- angle 95-115 degrees;
- joint NAC occupancy and longest continuous residence;
- a representative lower-potential NAC frame selected only after thermal equilibration;
- branch-specific gate opening using residues 261-266, excluding Thr267;
- temperature and pressure statistics;
- FATAL, NaN, LINCS and SETTLE counts.

The 100 ps free pilot determines whether to extend C18 to the required >=1 ns free window. Restrained results cannot approve QM/MM. A candidate without a complete free-window audit is `NOT_EVALUATED`.

Step1 DFTB3/3OB-3-1 preparation begins only after a stable C18 reactant/NAC ensemble is identified. Step1 TS, committor and PMF inputs are added incrementally after the reactant basin is validated. Step2 is not started before a defensible step1 acyl-enzyme/product basin exists.

## Recovery

- Never overwrite a completed stage.
- Resume a candidate only from its last valid checkpoint and recorded stage.
- Do not resubmit merely because SSH or scheduler display is transiently unavailable.
- Re-run only the failed candidate and append a new history row.
- Preserve every audit JSON and failure reason.

### 100 ps free-pilot audit and 1 ns extension

Job `61687591` completed with exit code `0:0` and no FATAL, NaN, LINCS or SETTLE event. The restrained 100 K, 150 K and 300 K release stages had NAC occupancies 0.3532, 0.2139 and 0.0579, respectively. These restrained values remain diagnostic only.

The 100 ps fully unrestrained NPT pilot contained 1 NAC frame among 501 frames (occupancy 0.001996). It occurred at 3.6 ps with distance 0.348 nm, angle 97.426 degrees and gate opening 2.631 nm. Mean gate opening was 2.576 nm. Temperature averaged 300.01 K; no numerical instability was detected. This is evidence that the rebuilt L2 can transiently form NAC, but it is not evidence of a stable high-occupancy NAC basin.

The required 1 ns continuation uses the exact 100 ps endpoint and checkpoint:

```bash
sbatch workflows/nylc_l4_nac_to_l2_rebalance_20260723/slurm/run_nylc_c18_step1_free_1ns.sbatch
```

Job `61688773` passed an independent `grompp -maxwarn 0` preflight and is the fully unrestrained 1 ns extension. Its scheduler completion alone is not a scientific PASS; audit the complete trajectory for NAC occupancy/residence, energy-conditioned NAC representatives, gate opening, thermodynamics and numerical events.

