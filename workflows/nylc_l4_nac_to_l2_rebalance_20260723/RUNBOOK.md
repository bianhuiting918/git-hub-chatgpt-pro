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

## Critical atom-identity correction (2026-07-24)

The historical NylC source manifest labeled global atom 8896 as catalytic Thr267 OG1. Direct GRO identity parsing proves that atom 8896 is Thr262 OG1, a gate residue. In the same gate-associated protein chain, catalytic Thr267 OG1 is atom 8961; the other protein copy has Thr267 OG1 at atom 3825.

Therefore all NylC C18/C23 residence labels, restrained pilots, the 100 ps free pilot, and job 61688773 that used atom 8896 are `SUPERSEDED_WRONG_NUCLEOPHILE_IDENTITY`. They are not catalytic NAC evidence. Job 61688773 was cancelled after 9 min 53 s. The corresponding pilot, continuation, 1 ns and 1 ns-audit scripts are hard-blocked with exit code 42 to prevent accidental reuse.

Direct corrected spot checks against atom 8961 gave:

| source representative | true Thr267 distance | true Thr267 angle | result |
| --- | ---: | ---: | --- |
| C18, 11854 ps | 1.661 nm | 124.548 deg | not NAC |
| C23, 29684 ps | 1.426 nm | 103.363 deg | not NAC |

The full corrected C18/C23 residence audit is run by:

```bash
sbatch workflows/nylc_l4_nac_to_l2_rebalance_20260723/slurm/run_true_nylc_thr267_residence_array.sbatch
```

The first array attempt, job 61690618, failed technically because concurrent `gmx angle` processes shared the default `angdist.xvg` output. The repaired version gives every segment a unique angle-distribution path and disables backups; job 61691264 is the retry. No earlier atom-8896 occupancy may be used to select a GS or approve QM/MM.



### Corrected true-Thr267 residence result and recapture pilot

Corrected array job `61691264` completed both branches using catalytic Thr267 OG1 atom 8961:

| branch | all frames | fully unrestrained frames | true NAC frames | closest angle-compatible free frame |
| --- | ---: | ---: | ---: | --- |
| C18 | 16,271 | 1,811 | 0 | segment_000122 at 11956 ps; 1.451 nm, 108.009 deg |
| C23 | 16,221 | 5,021 | 0 | segment_000320 at 31838 ps; 1.310 nm, 97.289 deg |

Thus the legacy L4 C18/C23 trajectories contain no authentic Thr267 NAC. Do not truncate their former representatives directly to L2.

The corrected source extractor is:

```bash
python scripts/extract_true_thr267_recapture_sources.py --all
```

It extracts C18 local 16 ps and C23 local 98 ps from the original segment XTC files through the mandatory temporary name `source.tmp.gro`. It promotes a file to `source.gro` only after time, identity, geometry and SHA256 verification. These are explicitly `VERIFIED_TRUE_THR267_RECAPTURE_START_NOT_NAC`.

Prepare the bounded L4 pilot with:

```bash
python scripts/prepare_true_thr267_recapture_pilot.py --all
```

Both candidates passed `grompp -maxwarn 0`. The pilot is 100 ps NVT at 300 K with freshly generated, recorded velocity seeds, a distance reference moving by -0.004 nm/ps with k=100 kJ mol-1 nm-2, a weak 105 degree angle restraint with k=20, and no gate restraint. Its sole purpose is to test a gentle true-Thr267 approach without severe numerical or steric failure; it cannot establish an unbiased NAC or GS.

Array job `61695008` failed before MD because escaped Slurm variables were interpreted literally. The repaired script passed contract tests and bash syntax validation. Retry job `61695303` runs C18 and C23 independently. After it completes, require `pilot_audit.json` technical PASS before any further approach. A later fully released L4 trajectory must supply the actual NAC residence and lower-potential GS candidate; only then may L4-to-L2 truncation resume.


### True-Thr267 staged response probes

All values below use atom 8961 (Thr267 OG1). Restrained probes are steering diagnostics only.

| job | protocol | C18 end / response / angle | C23 end / response / angle | disposition |
| --- | --- | --- | --- | --- |
| 61695303 | 100 ps, k=100, reference -0.4 nm | 1.5940 / -0.1435 nm / 91.57 deg | 1.4199 / -0.1100 nm / 89.18 deg | numerically stable; insufficient force |
| 61695989 | 100 ps, k=500, reference -0.2 nm | 1.4825 / -0.0320 nm / 89.31 deg | 1.3689 / -0.0589 nm / 95.06 deg | numerically stable; insufficient force; auto-audit filename failed and was rerun manually |
| 61697621 | 25 ps, k=2000, reference -0.1 nm | 1.4312 / +0.0193 nm / 106.99 deg | 1.2951 / +0.0148 nm / 88.98 deg | numerically stable; first correct-direction response, below 0.05 nm gate |

The first two protocols must not be extended from their endpoints. Response3 established a small, positive, numerically stable response. C18 currently has the better attack angle; C23 remains retained independently.

Job `61698403` continues both response3 checkpoints for 100 ps without regenerating velocities. It keeps k=2000 and rate=-0.004 nm/ps, leaves the 261-266 gate unrestrained, and audits each endpoint against its own response3 start. A scheduler exit code or restrained endpoint is not NAC/GS evidence. Only after a controlled approach followed by a fully unrestrained L4 release can a lower-potential stable true-Thr267 NAC be selected for L4-to-L2 truncation.

## True-Thr267 NAC recapture and unbiased GS selection (2026-07-24)

Historical C18/C23 NAC labels that used atom 8896 are superseded because atom
8896 is Thr262 OG1.  The corrected catalytic atom is Thr267 OG1, global atom
8961 for this NylC copy.  The historical fully unrestrained L4 branches contain
no corrected true-Thr267 NAC frames.

The bounded C18 recapture series produced authentic true-Thr267 NAC sampling
without numerical instability.  Release1 (job 61701900; distance restraint
3000 and angle restraint 300) sampled 35/101 NAC frames.  Release2 (job
61702619; distance restraint 1000 and angle restraint 100) sampled 39/101 NAC
frames, 38.61% occupancy and a 4 ps longest continuous run.  Release2 ended at
0.392481 nm and therefore exited the endpoint gate, but had no
LINCS/SETTLE/NaN/FATAL events and a 0.294284 nm minimum heavy-atom
ligand-protein contact.  Both windows remain restrained and cannot establish a
scientific GS.

Run the auditable low-energy NAC selector only on the immutable release2
trajectory:

```bash
sbatch slurm/run_select_true_thr267_low_energy_nac.sbatch
```

Job 61705002 selected the 67.0 ps frame: true-Thr267 distance 0.336447 nm,
attack angle 100.87 degrees, potential energy -1825802.875 kJ/mol, GRO SHA256
`d0f157575193c4662f5463c14f5d9075e817c18142aa74f0f7456b8bc732b4f7`.
This record is a restrained-source starting point, not a scientific GS.

The next gate is a three-seed fully unrestrained pilot:

```bash
sbatch slurm/run_true_thr267_unrestrained_pilot_array.sbatch
```

Array job 61705307 uses velocity seeds 26701-26703.  Each task performs 20 ps
fully unrestrained NVT initialization followed by 100 ps fully unrestrained NPT.
The NPT audit reports corrected NAC occupancy, longest residence, potential
energy, gate opening using residues 261-266 only, thermodynamics, close
contacts, and numerical events.  Compare all three seeds before extending any
seed to at least 1 ns.  A restrained window never counts as the scientific
pass.

### Unrestrained pilot outcome and 1 ns continuation

Array 61705307 completed all three tasks without LINCS/SETTLE/NaN/FATAL.
Seeds 26701 and 26702 sampled no true NAC in 101 NPT frames.  Seed 26703
sampled 6/101 true NAC frames (5.94% occupancy; 1 ps longest continuous run).
Its lowest-potential free NAC was at 8 ps: distance 0.329 nm, angle 114.257
degrees, potential energy -1829950.75 kJ/mol.  Therefore NAC is rare in this
pilot and the 100 ps result is not by itself a stable GS ensemble.

Job 61705692 continues seed 26703 for 1 ns fully unrestrained NPT from its
checkpoint.  The continuation does not restart from the early NAC frame and
does not add restraints.  After completion, select the lowest-potential frame
only from frames satisfying the corrected true-Thr267 NAC definition.  If no
NAC recurs, retain a scientific FAIL/NOT_EVALUATED reason and do not relabel the
restrained frame as GS.



### Rare free NAC L4-to-L2 build and post-EM continuation

The selected L4 source is the fully unrestrained seed-26703 frame at 8 ps
(job 61705307): Thr267 OG1--carbonyl C 0.328757 nm, attack angle 114.257
degrees.  It is an authentic rare free NAC candidate but is explicitly
`NOT_ENSEMBLE_STABLE`; the subsequent 1 ns L4 continuation (job 61705692)
had zero NAC frames.  This source is therefore a working GS hypothesis, not a
claim of a dominant equilibrium basin.

Job 61707444 built the audited PA66-L2 system with exact reactive C/O/N
coordinate preservation, 133589 total atoms, net charge zero, two endpoint cut
edges and `grompp -maxwarn 0` PASS.  The reactive global atoms are Thr267
OG1 8961, L2 carbonyl C 10287 and O 10288.  The gate remains residues 261-266,
excluding Thr267.

Rigid-water steepest descent job 61707615 produced a water-pair collapse and
is rejected.  Single-precision flexible-water CG job 61708159 stopped at
machine precision with Fmax 5844.9 and did not pass the EM gate.  Double-
precision flexible-water CG job 61708647 physically converged in 153 steps to
Fmax 419.978 kJ mol-1 nm-1 and potential energy -2307725.863 kJ mol-1 with no
FATAL, NaN, LINCS or SETTLE event.  Its Slurm state was FAILED only because the
postprocessing Python dictionary contained literal newline escapes.  The
immutable run was independently salvaged by:

```bash
python scripts/audit_em_double_result.py \
  --run-root "$TASK_ROOT/candidates/nylc_C18_trueT267_freeGS/runs/em_cg_flexible_double_retry3" \
  --job-id 61708647
```

The resulting `PASS.json` is `PASS_TECHNICAL_EM`; it is not scientific NAC
evidence.  A standard GROMACS 2022.1 NVT50 preflight from that exact GRO passed
`grompp -maxwarn 0` with TPR SHA256
`4d49d3519939823a15ee3602e1ae311979a661a9326ba0d4092564f6a179f204`.

Continuation job 61708900 was submitted with:

```bash
sbatch slurm/run_nylc_true_thr267_freegs_rebalance_after_em_double.sbatch
```

It starts from the immutable double-precision EM GRO and does not rerun EM.
The stages are restrained NVT 50/150/300 K, restrained NPT, weak-restraint NPT,
then at least 1 ns fully unrestrained NPT.  Only the final unrestrained audit
can decide whether the L2 candidate remains eligible for Step1 DFTB3/3OB-3-1
QM/MM preflight.


#### Independent postprocessing for job 61708900

The first running copy of the continuation script wrote literal escape sequences
to its per-stage JSON and TSV records.  This does not affect the GROMACS
coordinates, checkpoints, energies or trajectories, but those particular
stage JSON files are not trusted.  The branch copy has been corrected and
contract-tested.

A non-destructive dependent postprocessor is queued as job 61709209 with
`afterany:61708900` on the CPU-only `xahcnormal` partition.  It never runs
dynamics.  It independently scans every stage log, preserves a backup before
repairing only the malformed TSV text block, and writes all scientific analysis
into a new `postprocess_job_61709209` directory.  It uses true Thr267 atom
8961, L2 reactive C/O atoms 10287/10288, and the 261-266 gate group excluding
Thr267.  If the final 1 ns artifacts are incomplete or any stage has a numerical
failure, its status remains `NOT_EVALUATED_INCOMPLETE_OR_NUMERICAL_FAIL`.


#### Gated Step1 DFTB3 numerical preflight

Job 61709853 is queued with `afterany:61709209` on the SCNet CPU
`xahcnormal` partition.  It first requires the independent postprocessor
status `PASS_POSTPROCESS_TECHNICAL`.  It then requires at least one fully
unrestrained L2 frame satisfying distance <=0.35 nm and angle 95-115 degrees
and selects the lowest-potential frame only within that qualified NAC set.

If no such frame exists, the job writes
`NOT_EVALUATED_NO_FREE_L2_NAC` and exits without starting Sander.  If the gate
passes, it audits a QM region containing the complete processed Thr267 residue
plus the complete 79-atom PA66-L2 ligand, one protein QM/MM boundary bond,
C/H/N/O 3OB-3-1 coverage, `qmcharge=0`, singlet spin and even electron parity.
It then runs one DFTB3 minimization step plus one 20-step segment.  A numerical
PASS is not a TS, reaction coordinate, barrier, PMF or Step2 endpoint.

The PETase TS workflow is reused only methodologically: choose GS
representatives within an audited reactant/NAC basin, validate any candidate RC
with independent aimless shooting/committor outcomes and endpoint separation,
and require PMF window overlap, block stability, hard-error scans and MBAR
sensitivity.  PETase atom numbers and fitted RC weights are not transferable to
NylC.


The minimal smoke QM region is deliberately not the production Step1 region.
Negoro et al. (FEBS Journal 2023, DOI 10.1111/febs.16755) assign substrate
hydrolysis to the Asp308-Asp306-Thr267 catalytic triad, with Tyr146 and Lys189
as additional nylon-specific catalytic/substrate-binding residues and Asn219
in the Thr267 hydrogen-bond network.  Therefore any Step1 TS/committor/PMF work
must expand at least through Asp306/Asp308 and explicitly test the sensitivity
to including Tyr146/Lys189/Asn219.  The minimal Thr267+L2 smoke result can only
validate file conversion, DFTB3 parameter availability and short-run numerical
entry.


### Fully unrestrained L2 result and clean DFTB3 entry (2026-07-24)

Job 61708900 completed every GROMACS dynamics stage, including 1.000 ns fully
unrestrained NPT. Its Slurm state is FAILED only because the original inline
scientific-audit text had stale escaping. The independent repaired
postprocessor job 61710861 completed 0:0 and is authoritative for the free
window.

For NylC-C18, the branch-specific bonded-graph mapping places the source
reaction amide at L2 N3-C12(O2), not at the C23-only N4-C19(O3) mapping. The
correct global reactive atoms are Thr267 OG1 8961, L2 C12 10287 and O2 10288.
The L2 molecule occupies global atoms 10273-10351.

The 501 fully unrestrained frames contain 24 strict NAC frames (4.7904%
occupancy) in 20 visits. The longest continuous visit is 2 ps. The
lowest-potential strict NAC is 434 ps: 0.339 nm, 114.037 degrees and
-1832730.5 kJ/mol. The minimum all-atom and heavy-atom ligand-protein contacts
are 0.170839 and 0.271442 nm. Gate opening (residues 261-266, excluding
Thr267) is 2.7786 +/- 0.0845 nm. Temperature is 300.035 +/- 0.774 K. All six
stages have zero FATAL, LINCS warning, SETTLE problem and NaN counts. This is
scientific status PASS_UNRESTRAINED_L2_NAC_PRESENT, but it does not imply that
NAC is the dominant equilibrium basin.

The DFTB3 entry attempts are immutable and independently recorded:

- 61710886: missing ParmEd GROMACS include root;
- 61710981: protected prior output directory;
- 61711287: stale L2 global range;
- 61711453: neutral QM charge omitted the +1 processed N-terminal Thr267 and
  link-H parity;
- 61712026: zero-electronic-temperature SCC failure at the one-step gate;
- 61712561: 100 K electronic-temperature rescue passed one step but retained
  one SCC warning in 20 steps, so FAIL_DFTB3_PREFLIGHT;
- 61712692: 200 K sensitivity run completed 0:0. Both the one-step and 20-step
  stages have FINAL RESULTS and Run done with zero SCC, vlimit, SANDER BOMB,
  NaN and FATAL hits.

Job 61712692 uses the minimal numerical smoke region only: complete processed
N-terminal Thr267 plus complete neutral 79-atom L2, one link H, one boundary
bond, qmcharge=+1, singlet and 314 electrons including the link H. Amber18
states that dftb_telec may aid difficult SCC convergence and should be used
only when necessary with careful checking; the 200 K result is therefore a
recorded numerical sensitivity, not a scientific TS/PMF setting. Its PASS.json
SHA256 is
`74a598b392b17998359996d1dd8f6d4468086ac05854778b900a536714e407c9`.

Before Step1 TS/committor/PMF, construct and compare a production core QM region
containing Thr267, Asp306, Asp308 and complete L2 against a sensitivity region
that additionally contains Tyr146, Lys189 and Asn219. Audit formal charge,
link atoms, boundary bonds and electron parity for each region before any
reaction-coordinate scan.


## Step1 GS mechanism correction and retained 434 ps hypothesis (2026-07-24)

Jobs 61713551 and 61713841 tested an assumed direct Tyr146/Asn219-to-substrate-O2
oxyanion-hole gate. Restrained response jobs 61714772 (distance k=500) and
61715211 (distance k=2000) then tested whether those direct contacts could be
recaptured. Both response runs were numerically stable, retained strict
Thr267 attack NAC and had no heavy-atom clash, but the direct donor distances
did not shorten. These four jobs are now
`SUPERSEDED_INVALID_DIRECT_DONOR_ASSUMPTION` and must not gate GS or TS.

The primary NylC study (DOI 10.1111/febs.16755) supports a different
interpretation: Tyr146 contributes to the water/substrate environment, while
Lys189 and Asn219 participate in a Thr267-centered hydrogen-bond network that
includes water. It does not establish Tyr146 and Asn219 as simultaneous direct
hydrogen-bond donors to the substrate carbonyl oxygen. Do not resume the direct
O2 recapture scripts or increase those forces.

Retain the fully unrestrained 1 ns evidence from authoritative postprocessor
job 61710861: 24/501 strict NAC frames (4.7904%), with the lowest-potential
strict NAC at 434 ps (Thr267 OG1-C12 0.339 nm, attack angle 114.037 degrees,
potential -1832730.5 kJ/mol). The minimum heavy ligand-protein distance for the
free trajectory is 0.271442 nm. This 434 ps frame is the current low-potential
strict-NAC GS hypothesis; it is not an equilibrium-dominant state or a TS.

The 434 ps structure contains paper-consistent network evidence at a 0.40 nm
heavy-atom graph cutoff: waters with oxygen atom indices 50165 and 51302 connect
Thr267 OG1 to Asn219, water 81218 connects Asp306 to Asp308, and Lys189 NZ is
0.285 nm from Asn219 OD1. These are structural diagnostics, not proof of a
proton-transfer path. Before Step1 TS, reproduce this water-network audit,
compare Asp306/Asp308 protonation microstates, and construct production QM
regions. The numerical smoke PASS from job 61712692 remains valid but is not a
production TS region.

Compact correction audit:
`audit/nylc_C18_trueT267_freeGS.direct_donor_gate_superseded.json`.


## NylC-C18 neutral ASH306 full-system preflight (2026-07-24)

The isolated active beta-chain protonation probe is authoritative at SCNet job
61717760 (`PROBE_PASS_ASH306_CHAIN`). It regenerated chain H with amber99sb-ildn,
renamed only Asp306 to neutral ASH, added HD2, retained the processed N-terminal
Thr267 atoms, and moved no heavy atom after treating terminal Lys355 OC1/OC2 as
chemically interchangeable. Jobs 61717442 and 61717623 were technical failures
during probe plumbing; job 61717738 reached a scientifically valid chain product
but its first audit compared equivalent terminal oxygen names literally. Keep
all failure records; use job 61717760 as the authoritative probe.

The neutral full-system preflight replaces complete chain H, removes the Na+
farthest from protein plus L2 under the triclinic periodic box, and changes the
topology molecule count from NA 144 to NA 143. The removed ion is source global
atom 133463, residue 41488, 3.950970 nm from the nearest solute heavy atom. This
keeps 133589 atoms and total topology charge effectively zero.

Job 61718597 built the correct neutral system but failed grompp because the
position-restraint reference coordinate was omitted. The wrapper now supplies
`-r system.gro`. Authoritative job 61718715 completed 0:0 and produced
`PASS_ASH306_FULL_SYSTEM_PREFLIGHT`:

- total charge 1.46e-10 e; ASH306 charge approximately 0;
- one 13-atom ASH306 including HD2; one 79-atom PA66-L2;
- Thr267 OG1 to L2 C12 distance 0.339423 nm and attack angle 114.036894 degrees;
- minimum ligand-protein heavy distance 0.280318 nm;
- GROMACS grompp and ParmEd/Amber conversion passed;
- new indices: Thr267 OG1 8961, L2 C12/O2 10288/10289, L2 range 10274-10352.

This is a neutral topology/numerical-entry preflight only. It is not a GS
thermodynamic ranking, reaction coordinate, Step1 TS, Step2 TS, PMF, or barrier.
Next compare the ASH306/Asp308- microstate against the retained all-deprotonated
microstate with the same DFTB3/3ob-3-1 numerical smoke gates before initiating
any transition-path or umbrella workflow.

Compact audit:
`audit/nylc_C18_trueT267_freeGS.ash306_full_system_preflight.json`.


## ASH306/Asp308- production-QM numerical smoke (2026-07-24)

SCNet job 61719422 completed 0:0 and passed both one-step DFTB3/3ob-3-1
numerical smokes at `dftb_telec=200 K`. The core region contains 111 explicit
QM atoms, three link atoms, QM charge 0 and 388 electrons including links. The
network sensitivity region contains 162 explicit QM atoms, seven link atoms,
QM charge +1 and 570 electrons including links. Both regions have FINAL RESULTS
and Run done, with zero SCC convergence, vlimit, SANDER BOMB, segmentation,
forrtl, NaN and FATAL hits.

This PASS establishes numerical feasibility only. Raw absolute energies must
not be compared against the all-deprotonated Asp306-/Asp308- system because the
proton and counterion compositions differ. The mechanistically ready Step1
reactant microstate must be chosen using proton-acceptor geometry and subsequent
reaction-path/committor evidence, not the one-step smoke energy.

Compact audit:
`audit/nylc_C18_trueT267_freeGS.step1_ash306_qm_smoke.json`.
