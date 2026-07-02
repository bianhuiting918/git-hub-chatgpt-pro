# PETase pB histogram proxy and next sampling plan (2026-07-02)

Purpose: move from qualitative TS-core geometry toward the paper workflow, where ATESA/LMax identifies a linear RC and validates it with a committor probability histogram centered near pB = 0.5.

## Inputs used

- Existing deacylation committor candidates from attempts 041, 048, 053, 055, and 056.
- Existing acylation accepted TS-core member from attempt028/033 evidence.
- Paper-like deacylation CVs from the Fig. 4 caption: CV1 = d(SerO, water H) - d(water H, water O), CV2 = d(SerO, acyl C) - d(acyl C, water O), CV3 = d(acyl C, water O) - d(water O, water H).

## Generated data files

- `blind_work/06_ts_ensemble/petase_qmmm_stage2_ts_core_milestone_20260702/petase_stage2_pb_histogram_proxy_dataset.tsv`: candidate-level paper-like CVs, proxy RC values, committor counts, and pB excluding undecided outcomes.
- `blind_work/06_ts_ensemble/petase_qmmm_stage2_ts_core_milestone_20260702/petase_stage2_deacylation_pb_bins_proxy.tsv`: pB bins using the paper-coefficient RC proxy.
- `blind_work/06_ts_ensemble/petase_qmmm_stage2_ts_core_milestone_20260702/petase_stage2_deacylation_lmax_proxy.tsv`: small-data LMax/logistic proxy using CV1-CV2-CV3 and current product/reactant counts.
- `blind_work/06_ts_ensemble/petase_qmmm_stage2_ts_core_milestone_20260702/petase_stage2_deacylation_lmax_proxy_coefficients.tsv`: fitted coefficients for the LMax/logistic proxy.
- `blind_work/06_ts_ensemble/petase_qmmm_stage2_ts_core_milestone_20260702/petase_stage2_deacylation_lmax_proxy_bins.tsv`: pB bins in the fitted LMax eta coordinate.

## Key diagnostic result

Directly applying the paper-like deacylation RC form to our minimized candidate structures compresses most candidates into a narrow RC range near the dividing surface. Within that range observed pB still spans reactant-side, TS-like, and product-side outcomes. This means direct coefficient transfer from the paper is not a valid reproduction shortcut for our system.

The fitted LMax/logistic proxy captures the obvious product-side points w241/w243, but it does not fully separate local accepted TS members n01/n05 from neighboring lower-bracket points. This is expected because the current fit uses only minimized start structures and aggregate counts, while the paper LMax uses many aimless-shooting frames and likelihood maximization over trajectory outcomes.

## Current deacylation pB evidence in fitted LMax proxy

| candidate | status | observed pB | predicted pB | eta | product/reactant/undecided |
|---|---|---:|---:|---:|---|
| n04_less_break | undecided_reactant_side | 0.000 | 0.299 | -0.852 | 0/1/3 |
| n02_less_attack | lower_reactant_side | 0.000 | 0.343 | -0.650 | 0/2/2 |
| n05_proton_late | accepted_ts_core | 0.600 | 0.343 | -0.648 | 3/2/3 |
| n00_center | lower_after_extension | 0.200 | 0.345 | -0.641 | 1/4/3 |
| n03_more_break | reactant_side | 0.000 | 0.352 | -0.609 | 0/3/1 |
| n01_more_attack | accepted_ts_core | 0.600 | 0.399 | -0.408 | 3/2/3 |
| bridge_cser244_h150 | reactant_shifted_lower | 0.333 | 0.427 | -0.292 | 2/4/2 |
| bridge_cser240_h145 | accepted_ts_core | 0.600 | 0.516 | 0.064 | 6/4/6 |
| attempt041_w04 | early_reactant_side_bracket | 0.000 | 0.611 | 0.453 | 0/3/5 |
| w239 | current_probe_product_biased | 0.857 | 0.627 | 0.518 | 6/1/5 |
| bridge_cser242_h148 | product_shifted_upper | 0.714 | 0.673 | 0.721 | 5/2/1 |
| w241 | product_side | 1.000 | 0.966 | 3.352 | 4/0/0 |
| w243 | product_side | 1.000 | 1.000 | 8.569 | 4/0/0 |

## Next sampling plan

1. Do not start umbrella sampling yet. The paper does umbrella sampling only after RC discovery/validation; our RC is not yet validated.
2. Build a formal ATESA-style training set from trajectory frames, not only minimized start structures. Use existing committor trajectories around bridge_cser240_h145, n01_more_attack, n05_proton_late, n00_center, bridge_cser244_h150, and w239.
3. For each candidate frame, compute the paper CVs plus extra mechanistic CVs that our blind workflow showed are important: WatH-HisNE2, SerO-C, C-Oleave, and attack angle.
4. Run likelihood-maximization/ATESA RC selection on those frame-level data, then validate bins near eta = 0 with additional unbiased committor shots.
5. Only after pB histogram is centered near 0.5 should umbrella windows be launched along the discovered RC.

## CPU policy

No new QM/MM trajectories were launched for this checkpoint. Future committor补采样 must keep total process threads <=64; recommended batch size is 16-24 Amber jobs with OMP_NUM_THREADS=2.
