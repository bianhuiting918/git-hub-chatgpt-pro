# Project 01 always-questioning protocol

Updated: 2026-06-29

## Purpose

Use this protocol as a project-local "skill-like" checklist before and after every Project 01 action. Its job is to force active questioning, find weak assumptions early, and prevent the workflow from quietly drifting away from the scientific goal:

```text
learn how protein-background environment changes the active pocket and stabilizes or destabilizes the ligand / transition-state-like fragment
```

This protocol is not a replacement for calculations. It is a standing audit layer for deciding what to calculate, whether a result is label-quality, and whether an embedding/feature really corresponds to the intended protein environment.

## When to apply

Apply this protocol whenever doing any of the following:

- choosing enzyme systems or scaffolds
- defining fixed/designable residues
- generating sequence panels
- selecting PLACER poses
- reconstructing full protein-ligand structures
- running QMMM or no-MM cluster calculations
- extracting embeddings or geometry features
- building labels or training data
- uploading progress to GitHub

## The five standing questions

Ask these five questions at every step:

1. **What scientific variable is being changed?**
   - sequence background
   - pocket side-chain identity
   - active-site pose
   - ligand geometry
   - scaffold/backbone
   - calculation settings

2. **What must stay fixed for the comparison to mean anything?**
   - catalytic motif
   - ligand/reaction template
   - QM atom definition
   - protonation and charge state
   - embedding layer definition
   - QMMM/no-MM comparison protocol

3. **Does this step produce a label-quality result or only a technical smoke result?**
   - label-quality requires documented convergence and consistent setup
   - smoke results can validate the pipeline but should not train final models

4. **Could this result be explained by a confounder instead of protein-background effect?**
   - ligand pose changed too much
   - active-site geometry changed
   - QM region changed
   - protonation changed
   - SCF failed or used different settings
   - sequence identity bin is mislabeled

5. **What evidence proves the step succeeded?**
   - exact output file
   - command log
   - convergence line
   - energy line
   - residue mask table
   - sequence identity table
   - embedding feature table
   - GitHub progress note

## Stage-specific audit

### 1. Enzyme-system selection

Required questions:

- Does the enzyme have a clearly defined active-site motif?
- Is the ligand / substrate / TS-like state representable with a stable QM region?
- Is the expected rate-limiting chemistry mainly controlled by the active pocket?
- Is there a meaningful difference between no-MM cluster and with-MM QMMM?

Current selected systems:

| System | Role | Reason |
| --- | --- | --- |
| Baker serine hydrolase | Main design-enzyme system | Tests de novo design, covalent-like serine chemistry, oxyanion-hole/pocket effects. |
| KSI | Mechanistic control system | Tests electrostatic / hydrogen-bond stabilization of enolate or TS-like ligand. |

Reject or defer systems when:

- protonation or covalent intermediates are too ambiguous for a first-pass label
- ligand state is hard to define consistently
- active pocket is not the main determinant of the modeled step

### 2. Fixed/designable residue masks

Required outputs:

```text
residue_ligand_distance_layers.tsv
fixed_residues.txt
design_shell_residues.txt
distal_residues.txt
mask_summary.json
```

Required checks:

- Catalytic residues are fixed.
- Direct ligand-contact residues are fixed or explicitly marked as lightly sampled.
- Second/third shell positions are the main sequence-design target.
- Distal positions are not allowed to dominate identity-bin definitions.
- Mask was computed from a ligand-bound structure, not only sequence indices.

Default layer definitions:

```text
active_site = catalytic residues and essential reaction geometry residues
pocket_4A = residues with any heavy atom <= 4 A from ligand/reactive fragment
shell_4_8A = residues with min distance > 4 A and <= 8 A
shell_8_12A = residues with min distance > 8 A and <= 12 A
global = all residues
```

Question to ask before accepting a mask:

```text
If this residue mutates, could it directly change the reaction definition rather than only the protein background?
```

If yes, fix it unless the experiment is explicitly about active-site mutation.

### 3. Sequence-panel generation

Required questions:

- Which positions were fixed?
- Which positions were allowed to vary?
- What is the target identity bin?
- Is identity measured globally and by pocket/shell layer?
- Are the designed sequences chemically plausible around the ligand?

Record identity as:

```text
identity_active_site
identity_pocket_4A
identity_shell_4_8A
identity_shell_8_12A
identity_global
```

Do not accept a sequence panel if:

- catalytic motif changed unexpectedly
- pocket identity bin is unknown
- only global identity was computed
- generated sequences collapse into a single similarity regime

### 4. PLACER pose generation and selection

Required questions:

- How many poses were generated per sequence?
- Were poses selected only by score, or was geometry diversity retained?
- Does each selected pose preserve catalytic geometry?
- Does the ligand/reactive fragment remain comparable across sequences?

Minimum pose metadata:

```text
sequence_id
pose_id
PLACER score fields
ligand RMSD / pose cluster ID
Ser OG - reactive C distance
oxyanion-hole donor distances
key catalytic residue distances/angles
```

Reject or mark as low confidence when:

- ligand pose is incompatible with the intended reaction
- active-site distances are outside the acceptable range
- pose quality score is poor
- reconstruction cannot recover the complete ligand/QM fragment

### 5. QMMM and no-MM cluster calculations

Required questions:

- Is the QM atom set identical between with-MM and no-MM calculations?
- Is charge/multiplicity identical and documented?
- Is the geometry unminimized or minimized by a documented protocol?
- Did SCF converge?
- Is the energy line present?
- Are warnings acceptable for a smoke result only, or label-quality?

Required label fields:

```text
E_cluster_noMM
cluster_scf_status
E_QMMM_withMM
qmmm_scf_status
delta_E_env = E_QMMM_withMM - E_cluster_noMM
qm_atoms
mm_atoms
charge
multiplicity
geometry_source
```

Label-quality rule:

```text
Do not use a pair as final training label unless both QMMM and no-MM calculations are converged or explicitly accepted under a documented convergence threshold.
```

Known current warning points:

| System/model | Issue | Action |
| --- | --- | --- |
| serine hydrolase `model_004` | no-MM cluster SCF warning | Rerun with stricter SCF before final label use. |
| serine hydrolase `model_001` | QMMM SCF warning | Rerun with stricter SCF before final label use. |

### 6. Embedding extraction

Required questions:

- Was the layer mask computed from the exact pose used for the label?
- Are embeddings residue-aligned to the structure sequence?
- Are ligand-contact residues pooled separately from second/third shell residues?
- Is whole-protein embedding included only as a control, not the only feature?
- Are delta embeddings computed relative to the correct reference?

Required embedding groups:

```text
E_active_site
E_pocket_4A
E_shell_4_8A
E_shell_8_12A
E_global
delta_E_active_site
delta_E_pocket_4A
delta_E_shell_4_8A
delta_E_shell_8_12A
```

Do not accept an embedding table if:

- residue numbering is not traceable
- ligand-distance layers are missing
- sequence and structure residue order are inconsistent
- embeddings were pooled over the whole protein only

### 7. Model training and interpretation

Required baselines:

```text
geometry only
embedding only
geometry + embedding
```

Required splits:

```text
pose split
sequence/background split
scaffold split when scaffold diversity exists
```

Question to ask before claiming success:

```text
Does geometry + embedding outperform geometry-only on held-out sequence/background splits?
```

If no, do not scale blindly. Recheck label quality, embedding layer definitions, and sequence/background diversity.

## GitHub upload discipline

Upload:

- Markdown plans and progress notes
- small TSV/JSON manifests
- scripts needed to reproduce masks or feature tables
- selected command summaries

Do not upload:

- large raw PLACER ensembles
- trajectories
- wavefunction/restart files
- model weights
- conda/venv directories
- full QMMM output batches unless explicitly curated and small

Question before every commit:

```text
Does this commit help another agent reproduce the decision or understand the current state without storing large compute artifacts?
```

## Required progress-note template

Every major compute/update note should include:

```text
What changed:
What stayed fixed:
Why these structures/sequences were selected:
Evidence files:
Convergence / warning status:
What is label-quality vs smoke-only:
Next concrete action:
```

## Stop conditions

Stop and ask for review before scaling when any of these occur:

- catalytic motif changes unexpectedly
- ligand/reactive state cannot be reconstructed consistently
- QMMM and no-MM QM regions diverge
- SCF warnings become common
- embedding layer assignment is ambiguous
- sequence panel does not span intended identity bins
- geometry-only model performs as well as geometry + embedding

## Current next actions under this protocol

1. For Baker serine hydrolase, finish the fixed/designable mask and verify catalytic/pocket fixed residues.
2. Generate a small sequence panel with controlled identity bins.
3. For each sequence, compute layer-wise identity and reject sequences that violate fixed motif constraints.
4. Run PLACER only on accepted sequence backgrounds.
5. Select a small number of diverse poses for QMMM/no-MM.
6. Rerun warning calculations from the first pilot before using them as final labels.
7. Start the same protocol for KSI only after the serine hydrolase sequence-panel path is stable.
