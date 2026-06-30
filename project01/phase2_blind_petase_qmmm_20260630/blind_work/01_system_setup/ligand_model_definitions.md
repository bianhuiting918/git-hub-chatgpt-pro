# Blind PETase Ligand Model Definitions

Date: 2026-06-30

## Boundary

These substrate models are defined from PET ester chemistry and generic serine-hydrolase requirements only. They do not use the PETase paper's specific TS structures, reaction-coordinate terms, selected CVs, sampling windows, or mechanistic conclusions.

## Purpose

Stage 1 needs substrate models that are chemically close enough to PET to test acylation and deacylation, but small enough for docking, classical pre-equilibration, and later QM/MM sampling.

The initial blind workflow therefore keeps three non-covalent substrate models plus one covalent deacylation precursor definition.

## Model Set

### `PET_dimer_capped`

Purpose: primary acylation substrate model.

Working SMILES:

```text
OCCOC(=O)c1ccc(C(=O)OCCOC(=O)c2ccc(C(=O)OCCO)cc2)cc1
```

Interpretation: hydroxyethyl-capped PET dimer with two terephthalate units and an ethylene glycol linker. This is the first production substrate because it preserves PET-like aromatic/ester repetition while staying tractable.

Formal charge: 0.

Required atom labels:

- `Csc_acyl`: scissile ester carbonyl carbon attacked by Ser O-gamma.
- `Oox_acyl`: carbonyl oxygen that can occupy an oxyanion-hole geometry.
- `Olg_ethylene_glycol`: ester leaving oxygen.
- `Ser160_OG`: enzyme nucleophile; not part of ligand topology but required in pose filters.

### `BHET_like`

Purpose: small neutral control substrate for pose generation and docking sanity checks.

Working SMILES:

```text
OCCOC(=O)c1ccc(C(=O)OCCO)cc1
```

Interpretation: bis(2-hydroxyethyl) terephthalate-like model. Use to test whether the active site can support a chemically plausible ester attack geometry without the extra conformational burden of a PET dimer.

Formal charge: 0.

Required atom labels match the same scissile ester convention as `PET_dimer_capped`.

### `MHET_like`

Purpose: product-side and deacylation reference fragment.

Working pH-7 anionic SMILES:

```text
OCCOC(=O)c1ccc(C(=O)[O-])cc1
```

Interpretation: mono(2-hydroxyethyl) terephthalate-like fragment with an explicit carboxylate. A neutral carboxylic-acid variant may be generated only as a protonation-state sensitivity test.

Formal charge: -1 for the default pH-7 model.

### `MHET_like_acyl_enzyme_precursor`

Purpose: deacylation starting model after enzyme acylation.

This is not a standalone small-molecule SMILES model. It is a protein-covalent model in which Ser160 O-gamma is ester-linked to the acyl carbonyl of an MHET-like fragment. The distal carboxylate state must be documented as either anionic or neutral before QM/MM setup.

Required atom labels:

- `Cacyl_ser_ester`: acyl carbonyl carbon linked to Ser160.
- `Oox_acyl`: acyl carbonyl oxygen.
- `Ser160_OG`: covalent leaving oxygen during deacylation.
- `Wat_nuc_O`: nucleophilic water oxygen selected after classical pre-equilibration.

## Parameterization Route

For MM-only preparation and docking/pre-equilibration:

1. Generate 3D conformers from the listed SMILES.
2. Assign stable atom names before parameterization.
3. Use GAFF2/AM1-BCC or CGenFF consistently across the model set.
4. Record atom-name mapping in `qm_atom_labels/*.tsv`.
5. Validate that the named scissile ester atoms survive topology conversion.

For QM/MM:

- Do not rely on MM bonded terms for the reacting atoms.
- The QM region must include the scissile ester atoms, Ser160 side-chain atoms involved in reaction, catalytic His ring, Asp carboxylate atoms if proton-coupled behavior is being tested, and the relevant water/leaving group atoms.
- The exact QM region is selected later from the mechanism hypothesis tree, not from the paper's selected CVs.

## Grill Gates

Before any docking or substrate placement is accepted:

1. Can the chosen substrate model answer the PETase mechanism question, or is it too small to represent PET-like binding?
2. Are both PET-like repetition and small-substrate controls represented?
3. Are atom labels stable across SMILES-to-3D-to-topology conversion?
4. Is the deacylation precursor clearly a covalent enzyme model rather than an ambiguous free ligand?
5. Are protonation variants explicit rather than silently mixed?
