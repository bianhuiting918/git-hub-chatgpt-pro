# Current Blind Reproduction Status

Updated: 2026-07-01 Asia/Shanghai

Scope: PETase blind first-principles QM/MM mechanism reproduction. The paper is used only for methodology; no paper coordinates, trajectories, reaction coordinates, barriers, or mechanism conclusions are used as inputs before final validation.

## Server Execution

- CPU host: `210.73.40.29`
- Repository path: `/Dell/Dell14/bianht/petase_blind_qmmm/repo`
- Conda/micromamba environment: `$HOME/micromamba/envs/petase_stage1`
- Active GitHub commit after Stage 2 label fix: `7f30793`

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
- QM/MM masks must be remapped after `tleap`, because Amber atom indices differ from the original seed PDB serials after hydrogen/water/ion addition;`r`n- cleanup preparation now stages solvent/ion MM minimization, all-atom MM minimization, short DFTB3/MM minimization, and optional longer DFTB3/MM minimization before production equilibration.

Methodological inspiration from the article is allowed here, but concrete article coordinates, trajectories, reaction-coordinate formulas, selected CVs, windows, barriers, rates, and mechanism conclusions remain blocked until final validation.
## Latest Amber QM/MM Smoke Evidence

Updated: 2026-07-01 Asia/Shanghai

GitHub main has reached `6526f76` for the Amber/Sander Stage 4 route. On the CPU server, a clean detached worktree at `/Dell/Dell14/bianht/petase_blind_qmmm/repo_verify_6526f76_20260701_060758` verified:

- `python -m unittest discover -s project01/phase2_blind_petase_qmmm_20260630/tests`: `Ran 12 tests ... OK`.
- Accepted seed `REACTIVE_6EQE_BHET_like_E01_001` generated Amber QM/MM inputs and AmberTools topology-prep inputs.
- `run_amber_topology_prep.sh` completed with exit `0`, producing `complex.prmtop`, `complex.inpcrd`, `ligand.mol2`, and `ligand.frcmod`.
- `map_stage4_amber_qmmm_masks_to_topology.py` mapped 76 QM atoms from the original seed selection onto the `tleap` atom indices.
- `run_sander_qmmm_smoke.sh` completed with exit `0` after auto-setting `AMBERHOME`; Sander entered DFTB3 QM/MM and completed one minimization cycle.

This is a technical smoke pass only. The one-step output reports very large energy and `VDWAALS=*************`, so the structure is not yet a chemically relaxed Michaelis complex and no TS/PMF/barrier result exists. The next compute gate is restrained/nonreactive cleanup minimization and staged QM/MM minimization before any TS search or 200 ps equilibration is treated as production evidence.



