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

