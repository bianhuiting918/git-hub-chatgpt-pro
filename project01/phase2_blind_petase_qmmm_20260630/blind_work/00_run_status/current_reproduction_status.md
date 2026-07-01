# Current Blind Reproduction Status

Updated: 2026-07-02 Asia/Shanghai

Scope: PETase blind first-principles QM/MM mechanism reproduction. The paper is used only for methodology and final calibration. Paper coordinates, trajectories, barriers, and mechanistic conclusions are not used as inputs.

## Current Main Route

The active route is rank022 acylation with Amber/sander DFTB3 QM/MM and ATESA-style aimless shooting. The earlier rank021 route is rejected because short QM/MM equilibration lost the desired catalytic preorganization and old ATESA attempts spawned unstable residual jobs.

Current compute path:

```text
/Dell/Dell14/bianht/petase_blind_qmmm/repo/project01/phase2_blind_petase_qmmm_20260630
```

Active status details are recorded in:

```text
blind_work/00_run_status/20260702_rank022_acylation_ts_seed_and_atesa_smoke.md
```

## Geometry and TS-Search Status

A Michaelis-like rank022 starting structure was generated from structure/docking/reactive-geometry criteria, then refined by restrained and short unrestrained QM/MM. The current independent TS-near seed is not a final TS ensemble.

Key current bracket:

| seed | SerO-C A | C-Oleave A | SerH-SerO A | SerH-His A | SerH-Oleave A | direct shooting result |
|---|---:|---:|---:|---:|---:|---|
| lam26.55 | 1.790 | 2.491 | 1.180 | 2.691 | 1.218 | 3/3 no-basin |
| lam26.55625 | 1.697 | 2.526 | 1.364 | 2.346 | 1.175 | 2/3 product, 1/3 no-basin |
| lam26.5625 | 1.689 | 2.543 | 1.379 | 2.396 | 1.154 | 3/3 product |
| lam26.6 | 1.789 | 2.456 | 1.311 | 2.426 | 1.210 | 7/10 product, 3/10 no-basin |

Interpretation: the proton-transfer coordinates are TS-like, but the acyl-transfer coordinate is still product-biased relative to the paper's final TS reference. This is a TS-near bracket, not a validated transition-state ensemble.

## ATESA Status

ATESA was installed and can run Amber/sander QM/MM jobs on the CPU server. The first failure mode was numerical, not mechanistic: ATESA production jobs read unstable restart velocities and produced `vlimit`/temperature blow-up. That was fixed by using fresh velocities for production (`ntx=1`, `irest=0`), NVT/constant volume, and `dt=0.0001 ps`.

Validated ATESA smoke:

```text
blind_work/05_atesa_acylation/rank022_as_attempt010_lam2655625_freshvel_5moves
```

Result: ran normally for 5 moves with no `vlimit`, `NaN`, or SCC-DFTB errors, but acceptance was `0/5`; `as_raw.out` showed all five attempts classified to the same basin. This indicates the current ATESA center/commitment workflow is still not producing an accepted TS ensemble.

Second ATESA smoke:

```text
blind_work/05_atesa_acylation/rank022_as_attempt011_lam266_freshvel_5moves
```

Result: ran normally for 5 moves with no `vlimit`, `NaN`, or SCC-DFTB errors, but acceptance was again `0/5`; `as_raw.out` showed five `B <-` classifications. Independent endpoint geometry checks show 9/10 fwd/bwd trajectories are product/fwd by the four-distance basin and 1/10 is no-basin. No trajectory committed to the reactant/bwd basin.

Interpretation: fresh-velocity ATESA is numerically usable, but `lam26.6` is still product-biased and does not provide a committor-balanced TS seed. The next target is not umbrella sampling; it is a reactant-side TS-near seed or longer/diagnostic commitment trajectories that actually reach both reactant and product basins.

## GitHub Upload Policy

GitHub should contain the reproducible workflow, scripts, tests, protocols, manifests, and compact status files. Large binary trajectory/restart/topology outputs remain on the CPU server and are excluded by `.gitignore`.

Do not commit passwords, private SSH material, full Amber trajectories, NetCDF files, restart files, or `.prmtop` topology files.

## Next Required Scientific Step

1. Generate a cleaner reactant-side TS-near seed closer to `SerO-C ~1.9 A`, `C-Oleave ~1.9-2.1 A`, and `SerH-Oleave ~1.2 A`; current `lam26.55625/lam26.6` seeds are product-biased.
2. Run paired direct-shooting diagnostics from that seed and require both product/fwd and reactant/bwd commitments, not only product/no-basin.
3. Run ATESA only after the manual classifier shows a plausible two-sided commitment region.
4. Only after committor-balanced TS candidates are obtained should umbrella sampling or PMF windows be started.
