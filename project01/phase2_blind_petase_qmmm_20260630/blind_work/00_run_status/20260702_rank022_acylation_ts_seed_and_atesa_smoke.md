# 2026-07-02 rank022 acylation TS-seed and ATESA smoke status

Purpose: continue blind PETase acylation mechanism reproduction from self-generated rank022 structures, without using paper coordinates as inputs. Paper public geometries are used only as calibration/reference.

## rank021 route

The rank021 ATESA route remains rejected. Residual rank021 sander processes from earlier fake-slurm ATESA attempts were killed again on 2026-07-02 because they kept reappearing from old attempt directories.

## rank022 stage3 short QMMM stability

Directory:

`04_qmmm_rank022_reactant_prep/hid237_qmmm_equilibration_paperlike/`

Start:

`rank022_stage3_qmmm10_start.rst7`

Paper-like QM region for these tests:

- Ser160 sidechain + Asp206 sidechain + His237 sidechain + ligand
- DFTB3, `qmcharge=-1`
- no reaction-coordinate restraint

### 0.05 ps smoke

Input/output:

- `rank022_stage3_equil_0p05ps.in`
- `rank022_stage3_equil_0p05ps.out`
- `rank022_stage3_equil_0p05ps.rst7`

Result: completed normally, SCC warning count 0.

Geometry start -> 0.05 ps:

- SerO-C: 3.142 -> 3.090 A
- SerH-HisNE2: 2.045 -> 2.217 A
- HisNE2-Oleaving: 4.529 -> 4.456 A
- C-Oleaving: 1.371 -> 1.360 A
- attack angle SerO-C-Ocarbonyl: 73.4 -> 79.2 deg

Interpretation: unlike rank021, rank022 stage3 does not immediately lose Ser-His preorganization.

### 0.50 ps continuation

Input/output:

- `rank022_stage3_equil_next0p45ps.in`
- `rank022_stage3_equil_next0p45ps.out`
- `rank022_stage3_equil_next0p45ps.rst7`

Result: completed normally, SCC warning count 0.

Geometry at cumulative 0.50 ps:

- SerO-C 2.160 A
- C-Oleaving 2.708 A
- SerH-HisNE2 1.540 A
- SerH-Oleaving 4.569 A
- HisNE2-Oleaving 5.234 A
- attack angle SerO-C-Ocarbonyl 105.2 deg

Interpretation: stage3 is not a stable equilibrated Michaelis complex. It spontaneously progresses along Ser attack / C-O cleavage but without productive leaving-group proton relay. It is a TS-seed precursor, not a validated MC ensemble.

## Clean restrained TS-seed refinement

Existing `late_net` had better geometry but its original run had 12 SCC warnings. It was rerun from the clean `net_30step` seed with paper-like QM region, `qmcharge=-1`, and DFTB SCF stabilization settings.

Best clean seed selected:

`04_qmmm_rank022_reactant_prep/hid237_qmmm_equilibration_paperlike/rank022_late_net_paperlike_60step.rst7`

SCC warning count across the 20+20+20 step paper-like continuation: 0.

Geometry of selected 60-step seed:

- SerO-C 2.163 A
- C-Oleaving 1.830 A
- SerH-SerO 0.981 A
- SerH-HisNE2 2.396 A
- SerH-Oleaving 1.670 A
- HisNE2-Oleaving 3.421 A
- C=O 1.229 A
- attack angle SerO-C-Ocarbonyl 97.6 deg

The 80-step continuation remained SCC-clean but moved away from this combined geometry, so 60-step is currently preferred.

## ATESA smoke

Directory:

`05_atesa_acylation/rank022_as_attempt006_late60_smoke/`

Issue 1: ATESA/mdtraj could not read the sander-written binary/non-UTF8 restart. Fixed by converting to ASCII restart with ParmEd:

`rank022_late60_ts_seed_ascii.rst7`

Issue 2: with the original ATESA production template (NPT, dt=0.5 fs), forward/backward shooting from `late60` blew up within a few fs with vlimit, huge temperature, and NaN. This is not a valid committor result.

Stabilization test:

`05_atesa_acylation/rank022_as_attempt006_late60_smoke/rank022_late40_md_smoke_nvt_dt0p1fs.in`

Starting from the milder `late40` seed, NVT/constant-volume and dt=0.1 fs completed 100 steps without vlimit/NaN/SCC warning. Temperature relaxed toward 310 K. Geometry changed from:

- late40 start: SerO-C 2.111, C-Oleaving 1.729, SerH-HisNE2 2.262, SerH-Oleaving 1.824, attack 103.3 deg
- after 100 steps: SerO-C 2.017, C-Oleaving 1.738, SerH-HisNE2 1.861, SerH-Oleaving 2.290, attack 109.1 deg

Interpretation: stabilized shooting parameters prevent numerical blow-up, but the proton-relay geometry may drift. The next ATESA attempt should use constant-volume/NVT, dt=0.1 fs, and probably test late40/late60 seeds with short trajectories before scaling up.

## Current next step

Do not claim a TS ensemble yet. Next concrete step is to repair ATESA execution around the stabilized production template or manually launch controlled fwd/bwd trajectories from late40/late60 ASCII seeds, then classify endpoints with the four-distance commitment basins.

## 2026-07-02 continuation: direct shooting smoke and product bracket

Old rank021 ATESA drivers were found still alive and spawning rank021 `sander` jobs. They were killed with a pattern targeting `rank021_find_ts_attempt003`, `rank021_find_ts_attempt004`, and `rank021_as_attempt005`. No rank021/rank022 ATESA or `sander -O` jobs were left running after cleanup.

### Direct-seed shooting smoke

Directory:

`05_atesa_acylation/rank022_manual_shooting_attempt008_direct_seed/`

Reason: ATESA init restart files produced unstable velocities; direct `ntx=1` from ASCII seed with NVT/constant-volume and `dt=0.1 fs` is numerically stable.

Runs:

- `late40_direct_200step.out/rst7/nc`
- `late60_direct_200step.out/rst7/nc`
- `late60_direct_1000step.out/rst7/nc`

All completed without NaN and without SCC warning. Each had an early `vlimit` warning at the first step only, then relaxed to stable temperature. Endpoints did not enter the paper four-distance product/reactant basins:

- late40 200-step: no-basin, SerO-C 2.107, C-Oleaving 1.823, SerH-Oleaving 2.716
- late60 200-step: no-basin, SerO-C 2.183, C-Oleaving 1.979, SerH-Oleaving 2.689
- late60 1000-step: no-basin, SerO-C 2.556, C-Oleaving 2.128, SerH-Oleaving 4.079

Interpretation: the clean late40/late60 TS seeds are not committor-near enough. In unbiased dynamics the proton-relay coordinate drifts away from the leaving oxygen.

### Strong proton-relay restrained seed

Input:

- `rank022_ts_handoff_stronger_proton.RST`
- `rank022_ts_handoff_stronger_proton_20step.in`

Outputs:

- `rank022_ts_handoff_stronger_proton_20step.rst7`
- `rank022_ts_handoff_stronger_proton_40step.rst7`

Both runs had SCC warning count 0.

Geometry:

- strong20: no-basin, SerO-C 1.984, C-Oleaving 2.066, SerH-SerO 1.126, SerH-Oleaving 1.437, attack 108.9 deg
- strong40: no-basin, SerO-C 1.960, C-Oleaving 2.172, SerH-SerO 1.076, SerH-Oleaving 1.336, attack 102.7 deg

Interpretation: these are cleaner later TS-like handoff seeds, closer to product but not product-basin endpoints.

### Product-side bracket

Input:

- `rank022_product_pull.RST`
- `rank022_product_pull_20step.in`

Output:

`rank022_product_pull_from_strong40_20step.rst7`

SCC warning count: 0.

This endpoint satisfies the paper-style product/fwd basin definition:

- SerO-C 1.524 A, product criterion < 1.8 A
- SerH-SerO 1.647 A, product criterion > 1.4 A
- SerH-Oleaving 1.014 A, product criterion < 1.2 A
- C-Oleaving 2.572 A, product criterion > 2.4 A

It also has HisNE2-Oleaving 2.945 A, C=O 1.212 A, attack angle 121.6 deg.

Current bracket for rank022 acylation:

- Reactant side: `04_qmmm_rank022_reactant_prep/hid237_qmmm_equilibration_paperlike/rank022_stage3_qmmm10_start.rst7` is reactant-like by the four-distance basin.
- TS-like seed side: `05_atesa_acylation/rank022_manual_shooting_attempt008_direct_seed/rank022_ts_handoff_stronger_proton_20step.rst7` and `rank022_ts_handoff_stronger_proton_40step.rst7`.
- Product side: `05_atesa_acylation/rank022_manual_shooting_attempt008_direct_seed/rank022_product_pull_from_strong40_20step.rst7`.

Next step: run controlled restrained interpolation/short shooting between strong20/strong40 and product_pull, then repair ATESA to avoid unstable init velocities by starting production trajectories directly from ASCII coordinates with `ntx=1`, NVT/constant-volume, and `dt=0.1 fs`.

## 2026-07-02 01:16 CST update: lam26-lam27 committor bracket

Obstacle/status check requested by user.

No residual `sander -O` / `atesa rank022` processes were observed after the lam26 shooting loop.

Latest direct unrestrained QMMM shooting, NVT/constant volume, `ntx=1`, `dt=0.1 fs`, 500 steps:

- lam25 start: no-basin; 3/3 endpoints no-basin from earlier test.
- lam26 start generated from lam25 with midpoint restraints between lam25 and lam27; minimization ended normally, no SCC/NaN seen. Geometry: SerO-C 1.927, C-Oleave 2.237, SerH-SerO 1.114, SerH-His 2.551, SerH-Oleave 1.396, His-Oleave 3.389, C=O 1.187, attack Ocarb 103.1 deg. Direct endpoints: 3/3 no-basin, drifting back toward reactant side (SerO-C 2.306-2.368, C-Oleave 1.823-1.993, SerH remains on SerO/His side).
- lam27 start: no-basin but product-biased. Geometry: SerO-C 1.750, C-Oleave 2.507, SerH-SerO 1.336, SerH-His 2.780, SerH-Oleave 1.013, His-Oleave 3.182, C=O 1.199, attack Ocarb 108.9 deg. Direct endpoints: 3/3 product/fwd.
- lam29 and lam32 direct endpoints from earlier tests: 3/3 product/fwd each.

Interpretation: the current independent TS-search seed is not yet a TS ensemble. It has a clean product-side/reaction-side bracket between lam26 and lam27. Next single-variable test should be lam26.5 or a short restrained interpolation between lam26 and lam27, followed by expanded direct shooting. ATESA itself still needs restart/velocity initialization repair; direct `ntx=1` shooting is stable enough for the current bracket search.

## 2026-07-02 01:51 CST update: first TS-near committor seed found

Continuation from lam26-lam27 bracket.

- lam26.5 generated from lam26 with midpoint restraints did not land between lam26 and lam27; it relaxed back toward the reactant/proton-on-Ser side: SerO-C 1.984, C-Oleave 2.217, SerH-SerO 1.046, SerH-Oleave 1.580. It was not used for shooting.
- lam26.75 generated by product-side backoff from lam27 overshot into static product/fwd: SerO-C 1.638, C-Oleave 2.524, SerH-SerO 1.491, SerH-Oleave 1.073. It was not used for shooting.
- lam26.6 generated by product-side backoff from lam27 produced a TS-near no-basin seed: SerO-C 1.789, C-Oleave 2.456, SerH-SerO 1.311, SerH-His 2.426, SerH-Oleave 1.210, His-Oleave 3.049, C=O 1.194, attack Ocarb 113.7 deg.
- 10 unrestrained QMMM direct shootings from lam26.6 (`ntx=1`, NVT, dt=0.1 fs, 500 steps) all completed without SCC/NaN. Endpoint counts: 7 product/fwd, 3 no-basin. Product seeds: 26602, 26603, 26605, 26606, 26608, 26609, 26610. No-basin/reactant-side seeds: 26601, 26604, 26607.

Interpretation: lam26.6 is the first independently generated TS-near seed from the blind rank022 route, but it is product-biased (~0.7 product in this short direct-shooting test). Next step is to shift slightly toward the reactant side, likely by a product-side/backoff lam26.55 or lam26.5 variant starting from lam26.6 rather than from lam26.

## 2026-07-02 02:12 CST update: refined lam26.55-lam26.575-lam26.6 bracket

Additional refinement after lam26.6 was found product-biased.

- lam26.55 generated from lam26.6 by reactant-side fine backoff. Geometry: SerO-C 1.790, C-Oleave 2.491, SerH-SerO 1.180, SerH-His 2.691, SerH-Oleave 1.218, His-Oleave 2.986, C=O 1.193, attack Ocarb 106.5 deg. Three direct shootings completed without SCC/NaN; endpoint counts: 3/3 no-basin.
- lam26.575 generated as midpoint between lam26.55 and lam26.6. Geometry: SerO-C 1.585, C-Oleave 2.614, SerH-SerO 1.381, SerH-His 2.676, SerH-Oleave 1.093, His-Oleave 2.915, C=O 1.203, attack Ocarb 112.9 deg. It is statically no-basin only because SerH-SerO remains just below the product cutoff, but direct shooting gave 3/3 product/fwd.
- lam26.6 remains the best current TS-near seed statistically: 7/10 product/fwd, 3/10 no-basin. It is product-biased but provides the first usable transition-region ensemble candidate.

Interpretation: the committor boundary is very sharp in the proton-transfer coordinate. Current local bracket is lam26.55 (reactant/no-basin side) versus lam26.575 (product side), with lam26.6 product-biased but sampled more extensively. Next useful test is a midpoint between lam26.55 and lam26.575 (lam26.5625-like) or, preferably, repairing ATESA initialization so aimless shooting can generate an ensemble around this bracket without manual overfitting.

## 2026-07-02 02:24 CST update: lam26.5625 test

- lam26.5625 generated from lam26.55 toward lam26.575. Geometry: SerO-C 1.689, C-Oleave 2.543, SerH-SerO 1.379, SerH-His 2.396, SerH-Oleave 1.154, His-Oleave 2.945, C=O 1.147, attack Ocarb 117.8 deg. It is statically no-basin but extremely close to the product cutoff.
- Three unrestrained QMMM direct shootings completed without SCC/NaN. Endpoint counts: 3/3 product/fwd.

Updated bracket: lam26.55 is reactant/no-basin side (3/3 no-basin); lam26.5625 is product side (3/3 product/fwd); lam26.6 remains the broader sampled TS-near candidate (7/10 product, 3/10 no-basin). The productive boundary is now narrower than before and dominated by the proton-transfer coordinate.

Next action should be either lam26.55625/lam26.5575 manual bisection or repairing ATESA restart/velocity initialization and launching aimless shooting around lam26.55-lam26.5625.

## 2026-07-02 02:36 CST update: lam26.55625 mixed committor smoke

- lam26.55625 generated from lam26.55 toward lam26.5625. Geometry: SerO-C 1.697, C-Oleave 2.526, SerH-SerO 1.364, SerH-His 2.346, SerH-Oleave 1.175, His-Oleave 2.947, C=O 1.174, attack Ocarb 117.4 deg. It is statically no-basin but close to product threshold.
- Three unrestrained QMMM direct shootings completed without SCC/NaN. Endpoint counts: 2 product/fwd, 1 no-basin. Product seeds: 26556251, 26556253. No-basin seed: 26556252.

Interpretation: lam26.55625 is now the best manually generated TS-near seed in the refined rank022 bracket, closer to a mixed committor than lam26.6 (7/10 product) and more useful than lam26.55 (0/3 product) or lam26.5625 (3/3 product). It should be expanded or used as the center for ATESA once restart/velocity initialization is repaired.

## 2026-07-02 02:38 CST update: ATESA initialization failure evidence

Inspected attempt007 ATESA files:

- `rank022_as_attempt007_late40_nvt_smoke/atesa_input_files/aimless_shooting_init_amber.in` uses `ntx=1`, `irest=0`, `nstlim=1`, `dt=0.00001` to make an init restart.
- `rank022_as_attempt007_late40_nvt_smoke/atesa_input_files/aimless_shooting_prod_amber.in` then uses `ntx=5`, `irest=1`, `dt=0.0001` to continue from that restart.
- Both fwd/bwd prod outputs begin with severe velocity blow-up: `vlimit exceeded for step 0; vmax = 372.1386`, `TEMP(K)=*********`, then >4e5-6e5 K for subsequent steps.

Interpretation: the ATESA smoke failure is not DFTB/SCC failure and not the reaction coordinate itself. It is a restart/velocity initialization problem in the ATESA workflow. The manual direct-shooting runs are stable because they start from coordinates with `ntx=1`, `irest=0`, and assign velocities directly. Next ATESA repair should force production shooting to start from the chosen TS-near coordinates with fresh velocities (`ntx=1`, `irest=0`) or otherwise prevent reading the problematic init restart velocities.

## 2026-07-02 02:45 CST update: ATESA fresh-velocity smoke attempt009

Created `rank022_as_attempt009_lam2655625_freshvel_smoke` centered on lam26.55625. Initial coordinate had to be converted from NetCDF restart to true Amber ASCII restart because MDTraj/ATESA failed on the binary restart with `UnicodeDecodeError`.

Minimal repair tested: `aimless_shooting_prod_amber.in` changed from `ntx=5, irest=1` to `ntx=1, irest=0`, so prod trajectories ignore restart velocities and assign fresh velocities. ATESA then ran normally and exited without BrokenPipe. Prod fwd/bwd no longer show the previous `vmax=372` / >4e5 K blow-up; fwd started around 452 K and cooled to ~305 K by step 10, bwd started around 456 K and cooled to ~309 K by step 10. No SCC/NaN was observed in the smoke logs.

This attempt used `max_moves=1`, so it only validated the initialization repair. Status: acceptance ratio 0/1. Next step: run a slightly larger fresh-velocity ATESA smoke (attempt010) around the same lam26.55625 seed.

## 2026-07-02 continuation: ATESA attempt011 centered on lam26.6

Directory:

`05_atesa_acylation/rank022_as_attempt011_lam266_freshvel_5moves/`

Purpose: single-variable test after attempt010. Keep the repaired fresh-velocity ATESA settings (`ntx=1`, `irest=0`, NVT/constant volume, `dt=0.0001 ps`, 500 production steps) and move the center from `lam26.55625` to the more product-biased direct-shooting seed `lam26.6`.

Run status:

- ATESA exited normally.
- Output count: 10 production outputs (`5 fwd + 5 bwd`).
- Error scan: no `vlimit exceeded`, no `NaN`, no `QMMM SCC-DFTB`, no traceback.
- `work/status.txt`: acceptance ratio `0/5`.
- `work/as_raw.out`: five `B <-` classifications.

Independent endpoint geometry check from the NetCDF last frames:

- 9/10 fwd/bwd endpoints satisfy the product/fwd four-distance basin.
- 1/10 endpoint is no-basin.
- 0/10 endpoints satisfy the reactant/bwd basin.

Interpretation: the repaired ATESA path is numerically stable, but the `lam26.6` center is still product-biased. The one-sided ATESA result is consistent with endpoint geometry and does not yet indicate a final TS ensemble. The next scientific step is to generate a cleaner reactant-side TS-near seed or extend diagnostic commitment trajectories so that both product and reactant basins are actually reached before scaling ATESA or starting umbrella sampling.
