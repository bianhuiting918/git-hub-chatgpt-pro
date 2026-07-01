# 2026-07-01 rank021 own His-mediated PETase QMMM pose rescreen

## Purpose

Continue the blind, first-principles PETase acylation mechanism workflow without using paper TS coordinates. Rank022 had a stable MC-like QMMM reactant but repeated His-mediated TS-seed attempts did not transfer Ser HG to His NE2. Per systematic debugging, I stopped hammering rank022 and rescreened existing top docking poses.

## Top5 heavy-atom rescreen

Generated: `blind_work/00_run_status/20260701_top5_pose_geometry_rescreen.tsv`.

Key observation from existing top5 MM-relaxed structures:

- rank021 stage300: SerO-C 2.840 A, O006-MetN 2.899 A, SerO-HisNE2 3.065 A, HisNE2-O004 3.094 A, nearest His water 3.276 A.
- rank022 stage300: SerO-C 2.779 A, O006-MetN 2.936 A, SerO-HisNE2 3.297 A, HisNE2-O004 3.113 A, nearest His water 3.219 A.

Interpretation: rank021 is a reasonable independent next candidate because it satisfies nucleophile proximity, oxyanion support, and His access without relying on the paper TS geometry.

## rank021 HID rebuild

Created new working directory:

`blind_work/04_qmmm_rank021_reactant_prep/`

Converted His237 to HID in the rank021 topology seed, changing only His237 residue naming in the PDB before tleap. Topology prep completed successfully:

- `complex.prmtop` and `complex.inpcrd` generated.
- `topology_prep.stderr` is empty.
- `tleap.log`: `Errors = 0`.
- Prmtop check: catalytic residue 209 is `HID` with `HD1` on ND1 and no HE2 on NE2, so NE2 is available as Ser HG acceptor.

## rank021 corrected reactant gate

Corrected gate used only mechanistic criteria:

- Ser OG to BHET carbonyl C.
- Broad carbonyl approach angle.
- Ser HG to HID NE2 preorganization.
- Carbonyl O to Met133 backbone N oxyanion-hole contact.
- Exclude nearest water from HID NE2 site.
- No forced HisNE2-to-leaving-O restraint.

MM gate output:

`blind_work/04_qmmm_rank021_reactant_prep/hid237_reactive_relaxation/mm_rank021_corrected_gate_stage1_300.rst7`

Geometry after MM gate:

| state | SerO-C | attack O006 | attack O004 | O006-MetN | SerH-His | SerO-His | His-O004 | C-O004 | nearWAT-His |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| MM gate | 3.326 | 64.0 | 100.3 | 3.374 | 1.923 | 2.846 | 4.947 | 1.348 | 3.618/res273 |

QMMM reactant gate output:

`blind_work/04_qmmm_rank021_reactant_prep/hid237_qmmm_acylation_reactant_gate/qmmm_rank021_reactant_gate_stage1_10step.rst7`

- No `QMMM SCC-DFTB: SCC failed` warning found.
- `FINAL RESULTS` present.

Geometry after QMMM reactant gate:

| state | SerO-C | attack O006 | attack O004 | O006-MetN | SerH-His | SerH-SerO | O-H-N angle | His-O004 | C-O004 | nearWAT-His |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| reactant_qmmm | 3.221 | 67.5 | 106.7 | 3.332 | 1.872 | 1.012 | 166.9 | 4.946 | 1.372 | 3.589/res273 |

Interpretation: rank021 gives a stable QMMM MC-like candidate with clean His preorganization and no water competition, but attack angle to O006 remains low.

## rank021 own His-mediated TS-seed tests

Started from rank021 QMMM reactant gate. These are own-mechanism tests, not paper-coordinate pulls.

Outputs:

- `qmmm_rank021_own_his_step1_30step.rst7`
- `qmmm_rank021_own_his_step2_30step.rst7`

Both outputs have `FINAL RESULTS`; no `QMMM SCC-DFTB: SCC failed` warning found by grep.

Geometry progression:

| state | SerO-C | attack O006 | attack O004 | O006-MetN | SerH-His | SerH-SerO | O-H-N angle | His-O004 | C-O004 | nearWAT-His |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| reactant_qmmm | 3.221 | 67.5 | 106.7 | 3.332 | 1.872 | 1.012 | 166.9 | 4.946 | 1.372 | 3.589/res273 |
| own_step1 | 2.925 | 58.1 | 106.3 | 3.298 | 1.718 | 1.009 | 160.2 | 4.762 | 1.365 | 3.676/res273 |
| own_step2 | 3.155 | 66.0 | 102.4 | 3.125 | 1.551 | 1.015 | 160.8 | 4.719 | 1.443 | 3.847/res273 |

Interpretation:

- rank021 improves over rank022 as a screened pose because SerO-C can approach and His remains preorganized without water competition.
- However, even stronger SerH-His / SerO-H elongation restraints do not produce Ser HG transfer: SerH-SerO remains about 1.01 A.
- The issue is therefore not simply rank022 pose strain. Current evidence points to a missing mechanistic degree of freedom or electronic environment, likely one of:
  1. Asp206 needs to be included in the QM region for the rank021 path.
  2. His/Asp orientation or protonation network still does not stabilize imidazolium formation.
  3. The chosen SerO-C/SerH-His coordinate is insufficient; an angle/collective coordinate or carbonyl polarization coordinate may be required.
  4. A different top pose may still be better, but rank021 should be tested with Asp-QM before being rejected.

## Next action

Run a minimal rank021 Asp-in-QM smoke/step test from `own_step2`, because rank022 already showed Asp-QM is stable and rank021 now has a better MC-like screened geometry. Do not start umbrella sampling yet; no validated TS seed exists.
