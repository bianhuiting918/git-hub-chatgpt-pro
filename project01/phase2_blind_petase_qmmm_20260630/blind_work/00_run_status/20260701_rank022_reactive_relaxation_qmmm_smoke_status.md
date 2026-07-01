# 2026-07-01 rank022 reactive-relaxation and Amber QM/MM smoke status

Main objective: advance blind PETase mechanism reproduction toward a self-generated Michaelis/reactant complex and Amber/Sander QM/MM entry point. GitHub update is deferred until this stage is scientifically cleaner.

## Atom-index correction

The pre-topology Stage1 reactive-relaxation RST files used PDB serials, which do not match Amber atom indices after tleap. A post-topology mapper was deployed on the CPU server:

- `project01/phase2_blind_petase_qmmm_20260630/scripts/map_stage1_amber_reactive_relaxation_restraints_to_topology.py`

For rank022, mapped restraints use Amber topology indices:

- SerOG-C005: `1911,3843`
- attack angle SerOG-C005-O006: `1911,3843,3844`
- leaving O-HisNE2: `3845,2991`
- carbonyl O-MetN: `3844,1915`

This fixes the earlier wrong-Sander-restraint problem.

## Top5 docking-pose MM screening

Top5 no-clash docking poses were built into Amber topologies and screened with 100-step restrained MM minimization:

- topology prep dir: `blind_work/01_system_setup/amber_reactive_relaxation_top5_topology_prep/`
- mapped relaxation dir: `blind_work/01_system_setup/amber_reactive_relaxation_top5_topology_mapped/`
- 100-step geometry summary: `top5_reactive_relaxation_100step_geometry_from_rst7.tsv`

100-step results:

| pose | SerOG-C005 A | attack angle deg | O004-HisNE2 A | O006-MetN A | pass |
|---|---:|---:|---:|---:|---|
| rank022 | 2.784 | 61.412 | 4.401 | 2.814 | no |
| rank027 | 2.969 | 62.753 | 4.473 | 3.048 | no |
| rank002 | 3.169 | 47.732 | 5.743 | 3.005 | no |
| rank021 | 2.861 | 56.462 | 4.200 | 3.076 | no |
| rank031 | 3.538 | 53.045 | 5.923 | 3.375 | no |

## Staged stronger relaxation

Only rank022 and rank021 were carried forward. Stage300 used stronger generic serine-hydrolase restraints, not literature coordinates.

- rank022 stage300: Ser-C 2.779 A, attack 90.100 deg, His 3.113 A, oxyanion 2.936 A, pass yes
- rank021 stage300: Ser-C 2.840 A, attack 68.318 deg, His 3.094 A, oxyanion 2.899 A, pass no

Rank022 was then weak-backed to the original weaker restraints:

- rank022 weakback200: Ser-C 2.940 A, attack 80.339 deg, His 3.009 A, oxyanion 2.964 A, pass yes
- rank022 release100 no NMR restraint: Ser-C 3.462 A, attack 67.078 deg, His 4.084 A, oxyanion 3.196 A, pass no

Interpretation: rank022 is a restrained Michaelis-like candidate suitable for QMMM reactant-prep testing, but not yet an unrestrained stable Michaelis basin.

## Amber/Sander DFTB3 QM/MM smoke

Amber QM/MM input generation and topology qmmask mapping were deployed:

- `generate_stage4_amber_qmmm_inputs.py`
- `map_stage4_amber_qmmm_masks_to_topology.py`

Using rank022 weakback200 coordinates:

1. Default QM region, 64 atoms: LIG + TRP159/SER160/ASP206/HIS237/SER238 side chains.
   - Sander exit: 0
   - DFTB3 entered and read 3ob-3-1 SK files.
   - Warning: SCC-DFTB convergence could not be achieved.

2. Core triad QM region, 41 atoms: LIG + SER160/HIS237/ASP206 side chains.
   - Sander exit: 0
   - Faster than 64-atom test.
   - Warning remains: SCC-DFTB convergence could not be achieved.

Interpretation: Amber/Sander QMMM plumbing is now working from the self-generated rank022 candidate, but DFTB3 SCC settings/QM region/protonation/charge need to be fixed before production QM/MM minimization or TS search.

## Next action

Do not proceed to TS search yet. Next scientific gate is to fix QMMM SCC convergence for rank022 weakback200, then run a short restrained QM/MM minimization/equilibration and test whether the Michaelis geometry remains stable with reduced restraints.
