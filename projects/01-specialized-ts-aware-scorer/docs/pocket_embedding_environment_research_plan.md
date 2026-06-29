# Project 01 research plan: pocket-layer embeddings for enzyme-background effects

Updated: 2026-06-29

## 1. Core research question

Project 01 is now technically able to run the intended label-generation path:

```text
designed protein / sequence background
-> PLACER active-site ligand conformer ensemble
-> reconstruct full protein + complete ligand/QM fragment
-> QMMM full-enzyme calculation with MM environment
-> no-MM active-site cluster calculation
-> compare protein-background stabilization effects
```

The next scientific question is no longer whether the pipeline can run. The next question is:

```text
Which protein-background structures should be generated and labeled so that a model can learn
how protein environment changes stabilize or destabilize the ligand / transition-state-like fragment?
```

The working hypothesis is:

```text
Protein embeddings contain useful information about the protein environment, but the signal should be
extracted at the active-pocket and shell levels rather than only as a whole-protein embedding.
```

The primary learning target should therefore be the environment effect:

```text
delta_E_env = E_QMMM_with_MM - E_cluster_no_MM
```

This is not a rigorous reaction barrier by itself. It is a practical label for how the protein background changes the energy of the same active-site / ligand state relative to an isolated active-site cluster.

## 2. What must be varied

The catalytic motif and ligand reaction template should remain fixed in the first learning phase. Variation should be introduced in a controlled hierarchy:

| Tier | Variation type | Purpose |
| --- | --- | --- |
| A | Same sequence or near-identical sequence, different PLACER poses | Measures local conformer / pose noise. |
| B | Same scaffold, catalytic motif fixed, first shell mostly fixed, second/third shell varied | Measures pocket-environment control without changing reaction definition. |
| C | Same scaffold, catalytic motif fixed, broader sequence identity changes | Measures stronger background sequence effects while keeping fold comparable. |
| D | Different scaffold/backbone, catalytic motif retained | Measures true background diversity and long-range structural effects. |

The immediate next phase should start with tiers A-C. Tier D is useful later, after the label format and embedding extraction are stable.

## 3. How to generate different sequence-similarity backgrounds

### 3.1 Fixed and designable positions

For each protein design/backbone, define position groups from the ligand-bound structure:

| Group | Definition | Default action |
| --- | --- | --- |
| Catalytic motif | Ser/His/Asp or current motif atoms needed for chemistry | Fixed. |
| Direct ligand-contact shell | Residues with any heavy atom within 4-5 A of ligand / reactive fragment | Fixed or lightly sampled. |
| Second shell | 5-8 A from ligand / reactive fragment | Main design target. |
| Third shell | 8-12 A from ligand / reactive fragment | Main design target. |
| Distal background | Outside 12 A | Used to create global background diversity. |

The key point is that "sequence similarity" should not be measured only by full-length identity. For this project, similarity should be measured at multiple layers:

```text
identity_active_site
identity_pocket_4A
identity_shell_4_8A
identity_shell_8_12A
identity_global
```

This allows us to ask whether a large global sequence change matters if the local pocket remains similar, or whether a small shell change can dominate the environment effect.

### 3.2 Sequence-generation tiers

Generate sequences in controlled bins:

| Bin | Global identity target | Pocket identity target | Shell identity target | Use |
| --- | ---: | ---: | ---: | --- |
| A1 | 95-100% | 95-100% | 95-100% | Pose/conformer noise baseline. |
| B1 | 80-95% | 90-100% | 60-90% | Main second-shell/third-shell effect set. |
| C1 | 60-80% | 80-100% | 40-80% | Strong same-scaffold background set. |
| C2 | 40-60% | 70-100% | 30-70% | Stress test for background diversity. |

The first production-sized pilot should not start with all bins. A practical pilot is:

```text
1 scaffold
20-50 designed sequences
3 bins: A1, B1, C1
10-20 PLACER poses per sequence
2-3 selected poses per sequence for QMMM/no-MM labels
```

### 3.3 Practical generation route

Use the existing design tools in this order:

1. Start from a catalytic-motif scaffold or Baker-style generated scaffold.
2. Define fixed positions from the active-site motif and ligand-contact shell.
3. Use ProteinMPNN or LigandMPNN to sample sequences with different temperatures / fixed-position masks.
4. Filter generated sequences into identity bins.
5. Predict or retain structures for each sequence.
6. Run PLACER to generate ligand / active-site conformer ensembles.
7. Select representative poses for QMMM/no-MM labeling.

The selection should avoid simply taking the top PLACER score every time. Keep diversity in:

```text
Ser OG - ligand reactive C distance
oxyanion-hole donor distances
His/Asp orientation if present
ligand pose RMSD
PLACER plddt/prmsd
pocket heavy-atom RMSD
```

## 4. How to extract embeddings at different pocket layers

### 4.1 Embedding sources

Use at least two embedding families if available:

| Embedding type | What it captures | Limitation |
| --- | --- | --- |
| Sequence model embeddings | Sequence background, evolutionary/biophysical context | Does not directly know the PLACER pose. |
| Structure-conditioned embeddings | Local fold, geometry, residue environment | Requires reliable structure input. |
| Physics/geometric features | Distances, hydrogen bonds, electrostatics proxies | Hand-designed, may miss nonlinear sequence effects. |

The first robust baseline should combine:

```text
sequence/residue embeddings
+ structure-derived pocket/shell masks
+ geometry and electrostatic proxy features
```

Do not rely on whole-protein embedding alone. It should be included as a control, not as the main representation.

### 4.2 Layer definitions

For each PLACER pose or reconstructed full protein structure, define residue sets by distance to ligand / reactive fragment:

```text
active_site = catalytic residues and covalent/near-covalent reaction atoms
pocket_4A = residues with any heavy atom <= 4 A from ligand / reactive fragment
shell_4_8A = residues with min distance > 4 A and <= 8 A
shell_8_12A = residues with min distance > 8 A and <= 12 A
global = all residues
```

The residue masks must be computed from the pose structure, not from sequence indices alone. This is important because PLACER can move ligand and side-chain conformations.

### 4.3 Pooling strategies

For each layer, compute several pooled representations:

```text
mean_pool(layer_embeddings)
max_pool(layer_embeddings)
distance_weighted_pool(layer_embeddings, weight = 1 / (distance_to_ligand + epsilon))
attention_pool(optional, learned later)
```

For matched designs, also compute difference embeddings:

```text
delta_embedding_layer = embedding_variant_layer - embedding_reference_layer
```

Useful derived terms:

```text
E_pocket_4A
E_shell_4_8A
E_shell_8_12A
E_global
delta_E_pocket_4A
delta_E_shell_4_8A
delta_E_shell_8_12A
E_pocket_4A - E_global
E_shell_4_8A - E_active_site
```

These features directly test whether local pocket embedding or broader shell embedding explains the QMMM/no-MM difference.

## 5. Label design

Each selected pose should have a label record with:

```text
E_cluster_noMM
cluster_scf_status
E_QMMM_withMM
qmmm_scf_status
delta_E_env = E_QMMM_withMM - E_cluster_noMM
geometry_features
placer_scores
sequence_identity_bins
embedding_layer_features
```

The first successful real-compute pilot already shows why SCF status must be explicit:

| Model | QMMM with MM | no-MM cluster | Use as label |
| --- | ---: | ---: | --- |
| model_004 | -119.764414147390838 a.u., SCF converged | -120.786870692322040 a.u., SCF warning | Technical label only until cluster rerun. |
| model_001 | -119.301135074270945 a.u., SCF warning | -121.359969245670726 a.u., SCF converged | Technical label only until QMMM rerun. |
| model_005 | -119.224325197741592 a.u., SCF converged | -121.139529659314576 a.u., SCF converged | First clean paired label. |

Clean labels require both QMMM and no-MM calculations to be converged or accepted under a documented threshold.

## 6. Modeling plan

Start with small, interpretable baselines before deep models.

### 6.1 Baseline 1: geometry only

Features:

```text
Ser OG - ligand reactive C distance
oxyanion-hole donor distances
His/Asp catalytic geometry if present
ligand RMSD / pose cluster ID
pocket residue distance statistics
simple electrostatic proxy from charged/polar residues
PLACER plddt/prmsd/fape
```

Purpose:

```text
Can simple geometry explain the label?
```

### 6.2 Baseline 2: embedding only

Features:

```text
active_site embedding
pocket_4A embedding
shell_4_8A embedding
shell_8_12A embedding
global embedding
delta embeddings relative to reference sequence
```

Purpose:

```text
Do learned protein embeddings contain background-environment information?
```

### 6.3 Baseline 3: geometry + embedding

Features:

```text
all geometry features
+ all layer embeddings
+ identity-bin metadata
```

Purpose:

```text
Does embedding add predictive power beyond direct geometric descriptors?
```

Embedding is useful only if Baseline 3 improves over Baseline 1 on held-out sequence/background splits.

## 7. Validation splits

Avoid random pose-level splits as the only evaluation. They can overestimate performance because poses from the same sequence/background are highly correlated.

Use three split types:

| Split | Held-out unit | What it tests |
| --- | --- | --- |
| Pose split | PLACER pose | Local conformer interpolation. |
| Sequence split | Designed sequence/background | Generalization to new sequence backgrounds. |
| Scaffold split | Backbone/scaffold | True out-of-background generalization. |

The first pilot should report pose and sequence splits. Scaffold split becomes meaningful only after Tier D data exists.

## 8. First executable pilot

The next concrete pilot should be small enough to run, but diverse enough to test the idea.

### 8.1 Inputs

```text
one catalytic-motif scaffold
one ligand / reaction template
20-50 designed sequences
3 sequence-similarity bins: A1, B1, C1
10-20 PLACER poses per sequence
2-3 selected poses per sequence for QMMM/no-MM
```

### 8.2 Required outputs

For each selected pose:

```text
full reconstructed protein-ligand structure
complete ligand/QM fragment
QMMM input and tpr
no-MM cluster input
QMMM output energy and SCF status
cluster output energy and SCF status
embedding layer feature record
geometry feature record
teacher-label row
```

### 8.3 Minimum useful dataset size

The first useful dataset should target:

```text
20 sequences x 2 poses = 40 paired labels
```

A stronger first dataset:

```text
50 sequences x 3 poses = 150 paired labels
```

The 40-label dataset is enough to test whether the labels and features are coherent. The 150-label dataset is the first reasonable scale for model comparison.

## 9. Immediate next actions

1. Rerun the two warning calculations from the first QMMM pilot:
   - `model_004` no-MM cluster with stricter SCF.
   - `model_001` QMMM with stricter SCF.
2. Define fixed/designable position masks for the current denovo serine-hydrolase example.
3. Generate a small sequence panel with controlled identity bins.
4. Run PLACER for each sequence/background.
5. Implement layer-mask extraction from ligand-bound structures:
   - active site
   - pocket <= 4 A
   - shell 4-8 A
   - shell 8-12 A
   - global
6. Extract sequence/structure embeddings per layer.
7. Build a teacher-label table linking:
   - sequence/background ID
   - pose ID
   - layer embeddings
   - geometry features
   - QMMM/no-MM labels
8. Train three baseline models:
   - geometry only
   - embedding only
   - geometry + embedding

## 10. Decision rule for continuing

Proceed to larger-scale labeling only if:

```text
geometry + embedding model > geometry-only model
```

on held-out sequence/background splits.

If embedding does not improve prediction, the next move should be to improve label quality, feature definitions, or structure diversity before scaling QMMM.

If embedding does improve prediction, the next move should be:

```text
increase sequence/background diversity
increase paired QMMM/no-MM labels
add different scaffold/backbone backgrounds
```

The scientific goal remains:

```text
learn how protein-background environment changes the active pocket and stabilizes or destabilizes the ligand / transition-state-like fragment.
```
