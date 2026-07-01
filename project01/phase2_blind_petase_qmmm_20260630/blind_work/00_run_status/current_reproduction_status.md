# Current Blind Reproduction Status

Updated: 2026-07-01 Asia/Shanghai

Scope: PETase blind first-principles QM/MM mechanism reproduction. The paper is used only for methodology; no paper coordinates, trajectories, reaction coordinates, barriers, or mechanism conclusions are used as inputs before final validation.

## Server Execution

- CPU host: `210.73.40.29`
- Repository path: `/Dell/Dell14/bianht/petase_blind_qmmm/repo`
- Conda/micromamba environment: `/Dell/Dell14/bianht/micromamba/envs/petase_stage1`
- Active GitHub branch: `main`

## Verified Gates

| Gate | Status | Evidence |
|---|---:|---|
| Stage 1/2 launcher | completed | `launch_blind_stage1_stage2_compute.py` exited `0` on CPU server |
| Reactive pose seed unit test | completed | `test_generate_stage1_reactive_pose_seeds.py`: `Ran 1 test ... OK` |
| Stage 2 manifest unit test | completed | `test_generate_stage2_classical_md_manifests.py`: `Ran 2 tests ... OK` |
| GS pose acceptance | completed | `gs_pose_manifest.tsv`: 76 total rows, 20 accepted rows |
| Accepted pose source | completed | accepted rows all have `generation_method=reactive_geometry_seed` |
| Stage 2 MD queue | completed | `md_replicate_queue.tsv`: 60 jobs, all `not_started` |
| Stage 2 stage labels | completed | 60 jobs have `stage=preacylation_michaelis_complex` |
| Productive conformer selection | not ready | `productive_conformer_manifest.tsv` remains `awaiting_md_replicates` |

## Current Accepted Ground-State Seeds

The accepted structures are blind reactive-geometry Michaelis-complex seeds, not transition states and not production-relaxed conformers. They are suitable for minimization/equilibration and short classical MD screening.

First accepted row checked on server:

| pose_id | Ser-Og-C distance A | attack angle deg | oxyanion count | His acceptor distance A | leaving relay distance A |
|---|---:|---:|---:|---:|---:|
| `REACTIVE_6EQE_BHET_like_E01_001` | 3.300 | 102.95 | 2 | 2.944 | 1.035 |

## Next Required Compute Gate

Run Stage 2 minimization, heating, equilibration, and short replicate MD for the 60 queued jobs. Then cluster active-site frames by blind catalytic geometry and write productive cluster representatives to:

```text
blind_work/02_classical_md/productive_conformer_manifest.tsv
```

Only productive conformer representatives may feed Stage 4 low-cost QM/MM scans. No TS ensemble exists yet.

The current route is Amber/Sander QM/MM. Earlier GMX-CP2K potential-energy work was for a different purpose and is not part of this PETase mechanism-reproduction workflow.

## QM/MM Engine Route

Use the Amber/Sander route for the current PETase mechanism reproduction. The CPU-server environment has AmberTools with `sander`, `tleap`, `antechamber`, and `parmchk2` available in the PETase micromamba environment.

Active Stage 4 input and topology-preparation generators:

```text
scripts/generate_stage4_amber_qmmm_inputs.py
scripts/prepare_stage4_amber_topology_inputs.py
scripts/map_stage4_amber_qmmm_masks_to_topology.py
scripts/prepare_stage4_amber_cleanup_inputs.py
```

The generator prepares blind Amber/Sander DFTB3 QM/MM inputs from accepted seed structures:

- QM/MM minimization with `ifqnt=1` and at least 2000 minimization steps;
- 200 ps DFTB3/MM equilibration at 310 K before TS-search tooling consumes coordinates;
- QM atom selection from catalytic side chains, nearby side chains, and the bound substrate atoms in our own seed structures;
- status explicitly remains `inputs_ready_needs_amber_prmtop_inpcrd` until Amber topology and coordinate files are built or mapped;
- topology preparation now writes ligand extraction, GAFF2/AM1-BCC `antechamber`, `parmchk2`, `tleap`, 15 A TIP3P solvation, and `complex.prmtop`/`complex.inpcrd` output checks;
- QM/MM masks must be remapped after `tleap`, because Amber atom indices differ from the original seed PDB serials after hydrogen/water/ion addition;
- cleanup preparation now stages solvent/ion MM minimization, all-atom MM minimization, short DFTB3/MM minimization, and optional longer DFTB3/MM minimization before production equilibration.

Methodological inspiration from the article is allowed here, but concrete article coordinates, trajectories, reaction-coordinate formulas, selected CVs, windows, barriers, rates, and mechanism conclusions remain blocked until final validation.

## Latest Amber QM/MM Smoke Evidence

Updated: 2026-07-01 Asia/Shanghai

GitHub main records the Amber/Sander Stage 4 route. On the CPU server, clean detached worktrees verified:

- `python -m unittest discover -s project01/phase2_blind_petase_qmmm_20260630/tests`: `Ran 12 tests ... OK`.
- Accepted seed `REACTIVE_6EQE_BHET_like_E01_001` generated Amber QM/MM inputs and AmberTools topology-prep inputs.
- `run_amber_topology_prep.sh` completed with exit `0`, producing `complex.prmtop`, `complex.inpcrd`, `ligand.mol2`, and `ligand.frcmod`.
- `map_stage4_amber_qmmm_masks_to_topology.py` mapped 76 QM atoms from the original seed selection onto the `tleap` atom indices.
- `run_sander_qmmm_smoke.sh` completed with exit `0` after auto-setting `AMBERHOME`; Sander entered DFTB3 QM/MM and completed one minimization cycle.
- After adding staged cleanup inputs at `cd3da25`, the full local test suite reported `Ran 13 tests ... OK`.
- In `/Dell/Dell14/bianht/petase_blind_qmmm/repo_verify_cd3da25_20260701_061947`, 50 solvent/ion MM minimization cycles plus 25 all-atom MM minimization cycles completed and removed the initial `VDWAALS=*************` overflow to finite MM values (`11994.4488`, then `11648.1807`).
- From the cleaned MM restart, a one-step Sander DFTB3/MM minimization completed with `FINAL RESULTS`, `Run done at 06:27:44.185 on 07/01/2026`, `VDWAALS = 11575.5319`, and `DFTBESCF = -7994.4898`.
- A deliberately too-short cleanup test using 10 solvent/ion MM cycles, 5 all-atom MM cycles, and 1 QM/MM cycle crashed in Sander with exit `139` while the starting structure still showed `VDWAALS=*************`. Do not use that ultra-short cleanup as the production route.

Additional science-critical evidence from 2026-07-01:

- The default 200/100/5 cleanup branch was stopped after three DFTB3/MM steps because every QMMM step reported `Convergence could not be achieved in this step`; final observed `VDWAALS = 23525.2070`, `DFTBESCF = 2374.9308`, no `Run done`.
- The earlier 50 solvent/ion MM + 25 all-atom MM branch is the usable Amber start: one QMMM step completed with no SCC warnings, and a four-step continuation also completed with exit `0`, no SCC warnings, `Run done at 11:45:58.264 on 07/01/2026`, final `VDWAALS = 11575.2133`, `DFTBESCF = -8012.6242`.
- However, the QMMM-relaxed geometry is not a productive acylation TS-search start: Ser132 OG to LIG C005 is 4.754 A and Ser132 HG to His209 NE2 is 4.979 A. Direct TS scanning from this frame is blocked until a productive conformer gate supplies a closer Ser-His-substrate arrangement.

No TS/PMF/barrier result exists yet. The next compute gate is productive conformer generation/selection before any Stage 4 reaction-coordinate scan is treated as meaningful.



