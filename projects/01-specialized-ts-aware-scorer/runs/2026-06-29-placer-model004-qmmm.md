# 2026-06-29 PLACER model_004 QMMM run note

## Scope

This note records the CPU-side no-min GMX-CP2K preparation for the PLACER `model_004` active-site conformer. Large coordinate, trajectory, and CP2K output files remain on the CPU server and are not committed here.

CPU work directory:

```text
/Dell/Dell14/bianht/project01_conformer_queue/placer_denovo_ser_hydrolase_n5_20260629T011253/full_length_mapped_no_min_20260629_1255/model_004
```

## Current data scale

Current PLACER data are from one mutation/design setting, `128A:75I`.

- PLACER conformers generated: `5`
- Recommended first-pass conformers: `model_004`, `model_001`, `model_005`
- Full-length no-min mapped structures prepared so far: `model_004`, `model_001`, `model_005`
- `pdb2gmx` succeeded for all three mapped protein-only structures: `160 residues`, `2494 atoms`, total charge `-10 e`

`model_004` is the active template and first calculation target.

## Structure state

- Full-length protein was reconstructed from the original PLACER full input.
- Matching active-site atoms were replaced by PLACER crop coordinates.
- No free minimization was used before QMMM setup.
- PLACER 75I heavy-atom geometry was preserved as the conformer reference.

PLACER `75I` heavy atoms present in `model_004`:

```text
N CA C O CB OG C1A OAC O17 C18 C19 C20
```

The first technical smoke used only the six extra QEX heavy atoms:

```text
C1A OAC O17 C18 C19 C20
```

That smoke established the GMX-CP2K path but is no longer the main scientific target.

## Technical fixes found

The earlier QMMM aborts were not proof that the PLACER geometry or GROMACS topology was unusable. They were caused by CP2K input/runtime details:

1. `COMMENSURATE` is required under `&DFT / &MGRID` for GAUSS/S-WAVE QMMM coupling.
2. `PREFERRED_DIAG_LIBRARY = SL` is required under `&GLOBAL`; without this, embedded CP2K tried ELPA and aborted with `Setting the ELPA real_kernel failed`.
3. CP2K needs an extended-charge PDB with MM point charges. The current complete-75I file is:

```text
processed_75i_complete_for_cp2k_ext.pdb
```

## Complete 75I reconstruction

Using `75I.json`, the non-leaving 75I/QEXC atoms appended beyond the normal Ser128 protein topology are:

```text
HW C1A OAC O17 C18 C19 C20 H21 H22 H23 H24 H25 H26 H27
```

Excluded as leaving/protein-terminal atoms for the current protein-context model:

```text
OXT HXT H2
```

Kept from the protein/Ser128 topology:

```text
N H CA HA CB HB1 HB2 OG HG C O
```

Ligand-side anchor fit used the six PLACER heavy atoms:

```text
C1A OAC O17 C18 C19 C20
```

Alignment report:

```text
ligand_anchor_rmsd_A=0.542525
qex_complete_count=14
protein_atoms=2494
qex_complete_atoms=14
total_atoms=2508
```

Generated complete-75I files:

```text
processed_75i_complete.gro
processed_75i_complete_for_cp2k_ext.pdb
topol_75i_complete.top
index_ser128_75i_complete.ndx
placer_ser128_75i_complete_qmmm.inp
sp-qmmm-ser128-75i-complete.mdp
sp_qmmm_ser128_75i_complete.tpr
cluster_ser128_75i_complete.xyz
cluster_ser128_75i_complete_noMM.inp
```

## Complete 75I QMMM with MM environment

`grompp` result for the complete-75I QMMM input:

```text
Number of MM atoms=2483; Number of QM atoms=25
Bonds removed=28; F_CONNBONDS added=10
Dihedrals removed=28; Links generated=2
```

QM region:

```text
SER128: N H CA HA CB HB1 HB2 OG HG C O
QEXC:   HW C1A OAC O17 C18 C19 C20 H21 H22 H23 H24 H25 H26 H27
```

Background QMMM job:

```text
pid file: mdrun_ser128_75i_complete.pid
pid at launch: 222526
TPR: sp_qmmm_ser128_75i_complete.tpr
console: mdrun_ser128_75i_complete.console
GROMACS log: sp_qmmm_ser128_75i_complete_run.log
CP2K output: placer_ser128_75i_complete_qmmm.out
exit file: mdrun_ser128_75i_complete.exit
```

Initial CP2K evidence:

```text
GLOBAL| Preferred diagonalization lib. SL
CHARGE_INFO| Total Charge of the Classical System: -10.000000
Number of electrons: 70
OT step 1 energy: -111.4088151124
```

No final QMMM energy had been reached at the time this note was updated.

## No-MM cluster control

A matched no-MM control was started with the same QM atoms only, using standalone CP2K:

```text
cluster_ser128_75i_complete.xyz
cluster_ser128_75i_complete_noMM.inp
cluster_ser128_75i_complete_noMM.out
cluster_ser128_75i_complete_noMM.exit
```

Background no-MM job:

```text
pid file: cluster_ser128_75i_complete_noMM.pid
pid at launch: 222566
```

Initial CP2K evidence:

```text
Number of electrons: 70
outer SCF iter 1 energy: -120.6893765099
outer SCF iter 2 reached at least OT step 9, energy: -120.7811352437
```

No final no-MM energy had been reached at the time this note was updated.

## Superseded technical smoke

The earlier six-heavy-atom `SER128+QEX` QMMM job was stopped after the complete-75I jobs started, to free CPU resources. It had reached the second outer-SCF but was not the final target.

```text
old pid: 150881
old input: sp_qmmm_ser128_qex_prod.tpr
old CP2K output: placer_ser128_qex_qmmm_prod.out
```

## Interpretation

The current pipeline is now at the requested comparison stage for the priority conformer:

1. Complete 75I/Ser128 QMMM with the full MM protein environment is running.
2. A matched no-MM QM cluster control is running.
3. Both calculations use the same PLACER-derived active-site heavy-atom conformer and no free full-system minimization.

Final interpretation must wait for converged CP2K `ENERGY| Total FORCE_EVAL` lines from both jobs.

## Next actions

1. Poll `mdrun_ser128_75i_complete.exit`, `placer_ser128_75i_complete_qmmm.out`, and `cluster_ser128_75i_complete_noMM.exit`.
2. Extract complete-75I QMMM and no-MM cluster final energies when available.
3. Compare the with-MM and no-MM energy behavior as a technical environment-effect control.
4. After `model_004` completes, replicate the validated complete-75I workflow to `model_001` and `model_005`.
