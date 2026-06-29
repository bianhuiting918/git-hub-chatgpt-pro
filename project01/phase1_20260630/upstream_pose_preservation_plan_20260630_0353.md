# Upstream pose-preservation plan if full ligand RMSD remains a hard gate

Current evidence: `strict_all_gate=0/21`, while `reaction_geometry_gate=1/21`. The only reaction-geometry-pass candidate is `bin80_rank1_rec2_model_005`, but it fails full ligand heavy RMSD. If the project keeps full ligand RMSD as a hard launch gate, the next work should preserve the full ligand reference pose before and through PLACER, instead of launching QMMM.

## Option A: reference-pose locked sanity input

Goal: prove the QMMM construction chain can preserve a known good full-ligand pose.

1. Start from the full-length reference/candidate input containing the complete bu2 ligand.
2. Insert fixed active/contact sidechains before PLACER or bypass PLACER for ligand placement.
3. Run the same sidechain/graft/strict filter chain.
4. Passing condition: full ligand RMSD <= 0.5 A, key-distance delta <= 0.35 A, catalytic/fixed-pocket RMSD thresholds pass.

This is a pipeline sanity control, not a conformer ensemble.

## Option B: PLACER with explicit pose anchoring

Goal: let PLACER sample local pocket degrees of freedom while preventing full-ligand drift.

1. Keep `--target_res A-128` or explicit crop center to avoid the fixed-ligand crop error.
2. Add bonded or distance constraints on multiple ligand anchor atoms, not only Ser128OG-C1.
3. Anchor atoms should cover the reactive atom and distal ligand body, for example C1 plus ring/chain atoms selected from bu2.
4. Re-run small n first, then split/map/graft/filter.

Passing condition remains strict all-gate if full ligand RMSD is hard.

## Option C: replace full ligand RMSD with local ligand/contact gate

Goal: keep a meaningful geometric filter without rejecting valid distal ligand ensemble diversity.

1. Define a local ligand atom subset around the reacting carbonyl/ester and oxyanion/contact region.
2. Keep reaction distances as hard constraints.
3. Treat full ligand RMSD as diversity/pose-family metadata rather than a blocker.
4. Use the current `bin80_rank1_rec2_model_005` as the first exploratory +/-MM candidate.

This option changes the project gate and must be labeled exploratory unless approved as the new launch criterion.
