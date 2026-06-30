# Blind PETase Stage 3/4 Mechanism And QM/MM Exploration Protocol

Date: 2026-06-30

## Boundary

This protocol uses accepted Stage 1/2 structures, substrate chemistry, generic serine-hydrolase chemistry, and pH/protonation evidence. It must not use the PETase paper's concrete TS coordinates, reaction-coordinate formulas, selected CVs, trajectories, sampling windows, barrier values, rate constants, or mechanism conclusions.

## Generated Files

Run:

```bash
python3 project01/phase2_blind_petase_qmmm_20260630/scripts/generate_stage3_mechanism_tree.py \
  --out-root project01/phase2_blind_petase_qmmm_20260630/blind_work
```

Expected outputs:

```text
blind_work/03_mechanism_tree/mechanism_hypotheses.yaml
blind_work/03_mechanism_tree/candidate_cv_sets.tsv
blind_work/04_qmmm_exploration/path_screening_table.tsv
blind_work/04_qmmm_exploration/ts_like_guess_manifest.tsv
```

## Stage 3 Grill Gate

Before any QM/MM scan starts:

1. Each path must map to at least one accepted Stage 1 pose or a documented sensitivity branch.
2. Each CV candidate must be traceable to atom labels from `qm_atom_labels/*_atoms.tsv` and protein residue labels from the prepared structure.
3. Protonation-branch paths can only be active if the protonation gate supports the branch.
4. Water-network paths can only be active if classical equilibration shows persistent waters in the relevant geometry.
5. A path with no accepted starting pose remains `not_started` and must not be converted into a TS claim.

## Stage 4 Low-Cost QM/MM Gate

For each active path:

1. Build a QM region from the reacting ester/acyl group, Ser side chain, catalytic His, catalytic Asp if proton-coupled behavior is tested, and any relay water/leaving group atoms required by the path.
2. Run low-cost exploration with DFTB3/MM, xTB/MM, or low-cost DFT/MM relaxed scans, string images, or constrained optimizations.
3. Record every attempted path in `path_screening_table.tsv`, including failures.
4. Extract TS-like guesses only from scan maxima, string images, or saddle-like constrained structures produced by this blind workflow.
5. Populate `ts_like_guess_manifest.tsv` only after a structure exists and has been checked for endpoint chemistry.

## Rejection Rules

Reject a path before TS refinement if:

- the starting pose fails Stage 1 geometry filters;
- required atom labels are ambiguous or missing;
- protonation state is unsupported by the protonation gate;
- low-cost scan follows unphysical chemistry;
- the path requires a paper-specific coordinate, CV, or conclusion to proceed.

## Output Is Not A Mechanism

The mechanism tree and CV table are hypotheses and screening inputs. A mechanism claim requires TS refinement, vibrational/saddle evidence where applicable, committor or dynamical validation, and later PMF/free-energy analysis.
