# Project 01: Specialized TS-aware enzyme scorer

Created: 2026-06-28

## 1. Project goal

This project builds a **TS-aware enzyme scorer（过渡态感知酶打分器）** for a fixed enzyme/reaction system.

The target setting is:

```text
known enzyme family + known substrate + known GS + known TS
→ evaluate whether a mutant or generated enzyme stabilizes TS relative to GS
→ output ΔG‡ / activity ranking
```

Here:

- TS（transition state，过渡态）is allowed as an input because the reaction is fixed and mechanism-known.
- GS（ground state，基态）is the bound substrate/reactant complex.
- ΔG‡（activation free energy，活化自由能/反应能垒）is the main prediction target.
- DFT（density functional theory，密度泛函理论）or QM/MM（quantum mechanics / molecular mechanics，量子力学/分子力学）provides supervision labels.

This project is the first controlled stage before attempting general enzyme prediction.

## 2. Scientific hypothesis

For a fixed reaction, enzyme activity improvement can be decomposed into several physical components:

```text
activity improvement
≈ TS relative stabilization
+ catalytic geometry quality
+ productive conformer population
- strain penalty
- product inhibition / release penalty
```

The core physical quantity is:

```text
ΔG‡ = G_TS - G_GS
```

The model should learn why one enzyme variant or generated design has a lower TS-GS barrier than another.

## 3. Initial enzyme systems

Priority systems:

1. PETase（PET水解酶）
2. PUase / polyurethane hydrolase（聚氨酯水解酶）
3. nylonase（尼龙酶）
4. CALB（Candida antarctica lipase B，水解酶基准体系）

Optional later systems:

5. laccase（漆酶，氧化还原体系）
6. carbonic anhydrase（碳酸酐酶，金属酶/质子转移基准）

## 4. Model inputs

For each enzyme variant or generated design:

### 4.1 Protein identity input

Use a general protein model embedding（蛋白嵌入表示）such as ESM-style embedding.

```text
sequence or structure
→ h_protein
```

This represents:

- mutation context（突变背景）
- sequence constraints（序列约束）
- fold prior（折叠先验）
- long-range residue context（远程残基上下文）

### 4.2 Reaction state / conformation input

Use PLACER-style embedding（蛋白-配体构象集合嵌入）or extracted geometric features from PLACER-generated conformers.

```text
GS complex ensemble
TS complex ensemble
→ h_GS_ensemble, h_TS_ensemble
```

This represents:

- ligand pose（配体姿态）
- side-chain rotamers（侧链旋转异构体）
- catalytic geometry（催化几何）
- conformational accessibility（构象可达性）

### 4.3 Explicit TS/GS physical features

Because this is a fixed-reaction scorer, TS can be supplied explicitly.

Recommended features:

```text
TS coordinates
GS coordinates
forming bond length
breaking bond length
attack angle
proton transfer distance
oxyanion-hole distances
metal coordination distances, if applicable
hydrogen-bond network descriptors
```

### 4.4 DFT / QM/MM supervision features

DFT/QM/MM outputs should be treated as supervision and calibration signals:

```text
ΔG‡
ΔE_TS-GS
q_TS
q_GS
Δq = q_TS - q_GS
reference electric potential / field, if available
```

## 5. Model output

Primary output:

```text
predicted ΔG‡
```

Secondary outputs:

```text
TS stabilization score
GS binding score
reaction-projected electrostatic score
productive conformer probability
variant ranking
```

## 6. Recommended architecture

```text
Protein sequence/structure
        ↓
Protein encoder, e.g. ESM
        ↓
h_protein

GS/TS complexes + PLACER conformer ensembles
        ↓
Conformation encoder
        ↓
h_GS, h_TS

Known reaction TS template / DFT-derived state descriptors
        ↓
TS encoder / reaction-state encoder
        ↓
z_TS, z_GS, Δ-state features

Fusion layer
        ↓
physics-aware scorer
        ↓
ΔG‡ / activity ranking
```

Important: for Project 01, TS can be an input because the task is mechanism-conditioned evaluation.

## 7. Suggested physical decomposition

Use a modular score before training any complex model:

```text
S_total =
  w_field   * S_field
+ w_geom    * S_geometry
+ w_ens     * S_ensemble
- w_strain  * S_strain
- w_prod    * S_product
```

Where:

- `S_field`: reaction-projected electric-field score（反应投影电场分数）
- `S_geometry`: catalytic geometry score（催化几何分数）
- `S_ensemble`: productive conformer population（有效催化构象比例）
- `S_strain`: conformational or ligand strain penalty（构象/配体应变惩罚）
- `S_product`: product inhibition / release penalty（产物抑制/释放惩罚）

The reaction-projected electrostatic term can be written as:

```text
S_field = Σ_i Δq_i * φ(r_i)
```

Where:

- `Δq_i = q_i_TS - q_i_GS`
- `φ(r_i)` is the enzyme electrostatic potential at the reaction atom/probe point.

## 8. Data schema

Each sample should correspond to one enzyme variant/design and one reaction state pair.

```yaml
enzyme_id: PETase_example
variant_id: WT_or_mutant_name
reaction_id: PET_hydrolysis
substrate_id: PET_model_fragment

protein:
  sequence: path_or_string
  structure_pdb: path
  embedding_path: path

states:
  GS:
    complex_pdb: path
    placer_ensemble_dir: path
    charges_path: path
  TS:
    complex_pdb: path
    placer_ensemble_dir: path
    charges_path: path
  Product:
    complex_pdb: optional_path

labels:
  delta_G_dagger: float_or_null
  delta_E_TS_GS: float_or_null
  kcat: float_or_null
  KM: float_or_null
  kcat_over_KM: float_or_null

features:
  geometry_json: path
  electrostatic_json: path
  ensemble_summary_json: path
```

## 9. Training strategy

Start simple.

### Stage 1: feature baseline

Use manually computed physical features and train:

- linear regression（线性回归）
- ridge regression（岭回归）
- random forest（随机森林）
- small MLP（小型多层感知机）

Target:

```text
ΔG‡ or ΔΔG‡ relative to WT
```

### Stage 2: embedding fusion

Add ESM embedding and PLACER-derived embedding to the feature baseline.

### Stage 3: contrastive ranking

Use pairwise mutant ranking:

```text
variant A lower barrier than variant B
```

Loss:

```text
L_rank = max(0, margin + score_A - score_B)
```

where lower score means lower predicted barrier.

## 10. Compute plan

### GPU tasks

- ESM embedding generation
- PLACER inference / ensemble generation, if available on GPU
- neural model training

### CPU tasks

- DFT/QM/MM TS-GS calculations
- restrained optimization
- reaction coordinate scans
- classical electrostatic calculations

With ~1000 CPU cores, target scale for the first paper-level dataset:

```text
3-5 enzyme systems
100-500 variants/designs total
10^5-10^6 PLACER conformers
500-2000 DFT/QM/MM anchor points
```

Do not attempt exhaustive QM/MM potential energy surfaces.

## 11. Evaluation

Compare against:

1. embedding-only activity prediction
2. geometry-only catalytic scoring
3. docking / binding-score style scoring
4. PLACER ensemble score without TS-GS physics
5. full proposed TS-aware score

Metrics:

```text
Pearson / Spearman correlation with ΔG‡ or activity
pairwise ranking accuracy
top-K enrichment
ability to identify known beneficial mutations
```

## 12. Deliverable

A working Project 01 implementation should produce:

```text
input variant structures + known TS/GS
→ feature table
→ trained scorer
→ ranked enzyme variants/designs
→ decomposition report showing why variants are predicted better/worse
```
