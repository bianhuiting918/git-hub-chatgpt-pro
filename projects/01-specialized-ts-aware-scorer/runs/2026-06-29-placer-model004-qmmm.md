# 2026-06-29 PLACER model_004 QMMM run note

## Scope

This note records the CPU-side no-min GMX-CP2K preparation for the PLACER `model_004` active-site conformer. Large coordinate, trajectory, and CP2K output files remain on the CPU server and are not committed here.

CPU work directory:

```text
/Dell/Dell14/bianht/project01_conformer_queue/placer_denovo_ser_hydrolase_n5_20260629T011253/full_length_mapped_no_min_20260629_1255/model_004
```

## Structure state

- Full-length protein was reconstructed from the original PLACER full input.
- Matching active-site atoms were replaced by PLACER crop coordinates.
- No free minimization was used before QMMM setup.
- Six PLACER 75I/QEX heavy atoms were appended as an extra molecule for the current technical smoke.
- `pdb2gmx` protein topology succeeded with `160 residues`, `2494 atoms`, total charge `-10 e`.

## Technical fixes found

The earlier QMMM aborts were not proof that the PLACER geometry or GROMACS topology was unusable. They were caused by CP2K input/runtime details:

1. `COMMENSURATE` was required under `&DFT / &MGRID` for GAUSS/S-WAVE QMMM coupling.
2. `PREFERRED_DIAG_LIBRARY = SL` was required under `&GLOBAL`; without this, the embedded CP2K tried ELPA and aborted with `Setting the ELPA real_kernel failed`.
3. CP2K needed an extended-charge PDB. Generated:

```text
processed_qex_box8_for_cp2k_ext.pdb
```

Generation summary:

```text
atoms=2500
qm_atoms=17
total_gromacs_charge=-10.000000
mm_charge_written=-10.000000
```

## SER128+QEX QMMM validation

Input files updated:

```text
placer_ser128_qex_qmmm_sp.inp
processed_qex_box8_for_cp2k_ext.pdb
index_ser128_qex.ndx
sp-qmmm-ser128-qex-steep1.mdp
```

`grompp` result:

```text
Number of MM atoms=2483; Number of QM atoms=17
Bonds removed=33; F_CONNBONDS added=15
Dihedrals removed=28; Links generated=2
```

Short 900 s `mdrun` validation reached CP2K SCF output without the previous missing-charge or ELPA aborts:

```text
GLOBAL| Preferred diagonalization lib. SL
Total QM atoms: 17
Number of electrons: 62
OT step 1 energy: -106.1542671308
OT step 2 energy: -108.6588370191
OT step 3 energy: -109.5639576728
```

The 900 s run timed out before full SCF completion, so these are not final energies.

## Background production smoke started

A longer no-min one-step QMMM run was started in the same directory:

```text
pid file: mdrun_ser128_qex_prod.pid
pid at launch: 150881
timeout: 21600 s
threads: 8 OpenMP threads
TPR: sp_qmmm_ser128_qex_prod.tpr
MD log: sp_qmmm_ser128_qex_prod.log
CP2K output: placer_ser128_qex_qmmm_prod.out
console: mdrun_ser128_qex_prod.console
exit file: mdrun_ser128_qex_prod.exit
```

Initial check confirmed the process was running and CP2K output had started:

```text
GLOBAL| Preferred diagonalization lib. SL
CHARGE_INFO| Total Charge of the Classical System: -10.000000
Total QM atoms: 17
```

## Interpretation

The current pipeline is now past the setup blockers and is performing the actual QMMM SCF for the PLACER-derived active site. The remaining issue is wall time/SCF completion, not file-format setup. The current SER128+QEX model is still a technical smoke because QEX contains only six extra heavy atoms; final production labels should use a chemically complete 75I/transition-state representation with hydrogens and an explicit charge/protonation decision.

## Next actions

1. Poll `mdrun_ser128_qex_prod.exit` and `placer_ser128_qex_qmmm_prod.out` for completion or failure.
2. If complete, extract CP2K `ENERGY| Total FORCE_EVAL` and GROMACS potential energy.
3. If not complete within the timeout, keep the input as the validated template and reduce the smoke problem by using a smaller/complete QM test molecule or looser SCF settings only for pipeline testing.
4. Build the chemically complete 75I/TS-state topology before treating energies as scientific labels.
