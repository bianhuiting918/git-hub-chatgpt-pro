# Project 01 Phase1 strict-pass resampling pilot update

## Current remote run

- run root: `/Dell/Dell14/bianht/project01_similarity_phase1_20260630/placer_resample_strictgate_bin80_rank1_20260630_0312`
- driver PID at launch: `2492441`
- run manifest: `resample_run_manifest.tsv`

## Status

The first strategy has finished and failed before PLACER output generation:

- variant: `fixedlig_n5`
- status: `FAILED_RC_1`
- classification: `FAILED_BEFORE_EVALUATION`, not geometry fail
- root cause from log: when the only ligand is marked fixed with `--fixed_ligand X-bu2-1 --fixed_ligand_noise 0`, PLACER has no predicted ligand index for automatic crop/corruption center selection and raises `IndexError: list index out of range` in `get_crop_around_mol`.

The second strategy is still running at the time of this update:

- variant: `bonded_ser128_c1_n5`
- command concept: `--bonds A-128-75I-OG:X-1-bu2-C1:1.533`
- classification until output exists: `IN_PROGRESS`, not pass/fail

## Next recovery rule

Do not repeat `fixedlig_n5` in the same form. If a fixed-ligand sanity control is needed later, provide an explicit crop/corruption center or use low-noise predicted ligand sampling instead.

If `bonded_ser128_c1_n5` produces PLACER output, continue the normal chain: split conformers, full-bu2 mapping, sidechain completion, fixed-pocket sidechain graft, strict geometry filter. Only strict-pass rows become QMMM-ready.
