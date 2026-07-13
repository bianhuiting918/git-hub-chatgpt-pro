# NylC-GYAQ / PA66 L8 slab: current implementation and status

Last verified: 2026-07-13 16:35 (Asia/Shanghai)

## Scope

This record documents the implemented first-pass PA66 slab simulations. It is an execution record for the L8/slab branch of [the tetramer L4/L8 plan](nylc_tetramer_l4_l8_slab_simulation_plan.md), not a claim that the current trajectories establish productive catalysis.

The original scientific plan uses a NylC tetramer and multiple starting orientations. The current production systems are deliberately narrower: one manually placed AB-interface object or one manually placed AD-interface object per trajectory. They are never mixed in the same simulation box.

## Model actually simulated

- Enzyme: NylC-GYAQ model prepared as the processed alpha/beta-chain assembly. The AB and AD systems are separate interface-focused objects selected from the same tetrameric context.
- Polymer: all-atom 100-chain amorphous PA66 slab. The core is PA8-like; selected outer chains were extended to PA16-like lengths to retain longer exposed chain ends without replacing the compact core.
- Polymer topology: provisional GAFF2/AM1-BCC-based fast topology. It is suitable for workflow development and qualitative configuration screening, but not for final binding free energies or catalytic conclusions.
- Starting placement: the protein was manually positioned in ChimeraX at the side of the slab that contains exposed chains. During full-atom remapping, the protein was shifted to remove clashes. The initial nearest heavy-atom distances after remapping were 7.80 A for AB and 5.80 A for AD.

## Implemented workflow

1. Construct and equilibrate the full-atom PA66 slab.
2. Insert either AB or AD next to the exposed-chain side of the slab. Do not include both objects together.
3. Solvate and run 20 ps NVT at 300 K. The NVT final structures are the starts for production.
4. Run 50 ns contact-guard NPT production MD on one DCU per system.

### Contact-guard definition

The production restraint is a GROMACS `flat-bottom-high` pull coordinate between a local protein contact patch and a local slab contact patch.

- Reference distance: 0.600 nm (6 A)
- Force constant: 50 kJ mol-1 nm-2
- At or below 6 A: no pull force; the protein moves freely.
- Above 6 A: a weak restoring force acts only to prevent irreversible departure from the slab.

This is a patch-centre-of-mass approximation, not an exact conditional minimum atom-atom-distance restraint. It supports local contact sampling but must not be interpreted as unbiased adsorption thermodynamics.

## Current production status

| System | Slurm job | State | Last logged step | Last logged time | Contact-guard energy | Interpretation at this point |
| --- | ---: | --- | ---: | ---: | ---: | --- |
| AB | 60765892 | RUNNING | 3,722,600 | 7.4452 ns | 0 kJ mol-1 | Inside the free-contact region at the recorded frame; no catalytic conclusion. |
| AD | 60765894 | RUNNING | 4,561,700 | 9.1234 ns | 0 kJ mol-1 | Inside the free-contact region at the recorded frame; no catalytic conclusion. |

At the same check, both systems were near 300 K and had constraint RMSDs of about 6.2e-6 (AB) and 5.5e-6 (AD). The active Slurm allocation has a 10-day limit, so the 50 ns campaigns require continuation from checkpoint; `run_contact_guard_50ns_autoresume.sbatch` is the continuation entry point and must be checked at allocation end.

## Current visualization export

The trajectories available on 2026-07-13 were exported as solvent-free 50-frame PyMOL-friendly pairs:

| System | Frames | Solvent-free atoms | Files |
| --- | ---: | ---: | --- |
| AB | 50 | 52,513 | `AB_nowater_50frames.pdb` + `AB_nowater_50frames.dcd` |
| AD | 50 | 52,512 | `AD_nowater_50frames.pdb` + `AD_nowater_50frames.dcd` |

The reusable exporter is [`tools/trajectory_export/export_nowater_dcd.py`](../../tools/trajectory_export/export_nowater_dcd.py). It reads a `.tpr` plus `.xtc`, removes water and common monatomic ions by default, samples evenly spaced frames, and writes a reference PDB plus DCD trajectory.

## Remote project layout

Set `SCNET_ROOT` to the project root `nylon_pa66_scnet_20260708` in the relevant SCNet home directory. The following paths are relative to that root:

```text
adsorption_manual_chimerax_ab_ad_v3_shifted_fullatom_slab/
  nylc-gyaq_AB_manual_shifted_fullatom_slab/
    nvt_300K_20ps.gro
    contact_guard_50ns.tpr
    contact_guard_50ns.xtc
    contact_guard_50ns.cpt
    contact_guard_50ns.log
    contact_guard.ndx
    contact_guard_setup_report.tsv
    run_contact_guard_50ns_autoresume.sbatch
  nylc-gyaq_AD_manual_shifted_fullatom_slab/
    nvt_300K_20ps.gro
    contact_guard_50ns.tpr
    contact_guard_50ns.xtc
    contact_guard_50ns.cpt
    contact_guard_50ns.log
    contact_guard.ndx
    contact_guard_setup_report.tsv
    run_contact_guard_50ns_autoresume.sbatch

analysis_visualization_contact_guard_20260713/
  AB_nowater_50frames.pdb
  AB_nowater_50frames.dcd
  AD_nowater_50frames.pdb
  AD_nowater_50frames.dcd
  export_nowater_dcd.py
```

Large coordinates, trajectories, checkpoints, and PyMOL DCD/PDB exports are intentionally not versioned in GitHub. Preserve them on SCNet and the local analysis archive; use the file names above as the manifest.

## Required analysis after production

Do not rank AB versus AD from current elapsed time, distance snapshots, or pull energy alone. After each trajectory reaches the intended sampling length:

1. Verify temperature, pressure, density, constraints, periodic imaging, and checkpoint continuity.
2. Quantify protein-slab contact area, contact residue frequencies, orientation, and residence time.
3. Compare AD-cleft and AB-interface contact fractions separately.
4. Measure the nearest PA66 amide carbonyl C to Thr267 O-gamma distance, attack geometry, and persistence of exposed-chain entry.
5. Classify trajectories as nonproductive adsorption, surface contact, pocket approach, or productive-like candidate. A productive-like candidate still requires independent release replicas and ultimately QM/MM validation.

## Relation to published modelling

The staged contact test is intentionally separated from catalytic-register modelling. The NylC-GYAQ PA66 study manually positioned oligomers near the A/D region and used SMD before short release MD; this can generate catalytic-register candidates but does not demonstrate spontaneous binding. The NylC-HP study used capped/charged PA6 tetramers, incremental docking, and independent 150 ns MD replicas to compare candidate poses, but still did not model spontaneous extraction of a chain from a polymer slab.

References:

- Bocharova et al., *Polym. Chem.* 2025, 16, 1858-1868. https://doi.org/10.1039/D5PY00023H
- Puetz et al., *ChemSusChem* 2025, 18, e202500257. https://doi.org/10.1002/cssc.202500257

