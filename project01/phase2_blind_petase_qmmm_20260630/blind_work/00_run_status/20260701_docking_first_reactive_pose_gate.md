# 2026-07-01 docking-first reactive pose gate

Purpose: recover a physically plausible PETase/BHET-like Michaelis starting complex before Amber QM/MM. This does not use paper TS coordinates or paper reaction-coordinate results.

Evidence:
- Existing `gs_pose_manifest.tsv` contains docking outputs and reactive-geometry seeds.
- Existing Vina docking poses: 49 fail by broad serine-hydrolase geometry scoring.
- Existing reactive-geometry seeds: 20 previously pass geometry but all fail the new ligand-protein heavy-atom clash gate at 1.2 A.
- Re-scored pass + no-heavy-clash count from the full manifest: 0.

Focused re-docking tests:
- `POSE_6EQE_BHET_like_E01_exh32_modes40`: Vina exhaustiveness 32, 40 modes, 22 A box centered at Ser160 OG. Result: 40/40 fail geometry. Best SerOG-C around 3.48 A, but attack angle around 54 deg and His relay around 5.6 A.
- `POSE_6EQE_BHET_like_E01_box12_exh64_modes40`: Vina exhaustiveness 64, 40 requested modes, 12 A box centered at Ser160 OG. Result: 32/32 scored modes fail geometry. Best SerOG-C around 3.40 A, but attack angle around 48 deg and His relay around 5.7 A.
- Re-scoring the same poses as BHET-like E02 also fails; E02 can place SerOG-C near 3.15 A in one mode, but attack angle remains around 60 deg and His relay is far.

Interpretation:
- Docking is available and should be the upstream source of no-clash binding poses.
- Unconstrained Vina docking alone is not producing a QMMM-ready reactive Michaelis pose for this focused case.
- Next scientific step is docking-pose-based restrained placement/MM relaxation: start from no-clash docking poses, restrain SerOG-C, carbonyl attack angle, oxyanion contacts, and His relay, then run Amber MM minimization/short MD before Amber QM/MM.
