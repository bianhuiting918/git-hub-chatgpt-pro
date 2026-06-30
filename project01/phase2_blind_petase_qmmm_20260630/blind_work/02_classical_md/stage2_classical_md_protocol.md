# Stage 2 Classical MD Protocol

Purpose: generate a blind ensemble of PETase active-site conformers from accepted Stage 1 ground-state poses.

Boundary: do not use article trajectories, article transition-state frames, article reaction coordinates, or article barrier information.

Gate inputs:

- Stage 1 GS pose row with `pass_fail=pass`.
- Non-pending relaxed structure path with recorded input provenance.
- Ligand atom labels and protonation branch already reviewed.

Required workflow:

1. Minimize, heat, equilibrate, and run independent short MD replicates for each accepted pose.
2. Cluster active-site frames using blind catalytic-geometry descriptors.
3. Score clusters by nucleophile distance, attack angle, His acceptor geometry, and oxyanion stabilization.
4. Write accepted cluster representatives to `productive_conformer_manifest.tsv`.
5. Write rejected clusters or unstable poses to `rejected_pose_manifest.tsv` with explicit reasons.
6. Only productive conformers may feed Stage 4 low-cost QM/MM scans.
