# Project 01: Enzyme-specific true-TS barrier scorer

Created: 2026-06-28  
Updated: 2026-06-28

## 1. Project goal

This project builds a **true-TS-aware barrier scorer（真实过渡态感知势垒打分器）** for a specific enzyme, enzyme family, or well-defined catalytic class.

The target setting is:

```text
specific enzyme/family + substrate/reaction step + known GS + true/refined TS
+ QM/MM or DFT barrier labels
→ learn TS-GS barrier changes
→ predict ΔG‡ / ΔΔG‡ / variant ranking for new mutants or generated designs
```

Here:

- TS（transition state，过渡态）is explicitly allowed as an input.
- GS（ground state，基态）is the bound substrate/reactant complex.
- **true TS** means a TS or TS-like reaction state constructed and refined by a consistent protocol: QM/MM optimization, DFT optimization, reaction-coordinate scan, constrained optimization, or a validated TS template.
- `ΔG‡`, `ΔE_TS-GS`, `ΔΔG‡`, and `ΔΔE‡` are the primary prediction labels.
- Experimental `kcat` or `kcat/KM` is not the first-stage training target. It may be used later as weak external validation.
- Project 01 is the **teacher / oracle-TS model** for Project 02.

## 2. Scientific hypothesis

For a fixed enzyme system or catalytic class, computed barrier reduction can be decomposed into physical components:

```text
computed barrier improvement
≈ TS relative stabilization
+ catalytic geometry quality
+ productive conformer population
+ favorable reaction-projected electrostatics
- conformational / ligand strain penalty
- nonproductive binding penalty
```

The core physical quantity is:

```text
ΔG‡ = G_TS - G_GS
```

The most robust learning target is often the relative barrier:

```text
ΔΔG‡ = ΔG‡_variant - ΔG‡_reference
```

The model should learn why one mutant or generated design stabilizes the true/refined TS relative to GS better than another.

## 3. First-round enzyme systems

Priority systems for the first round:

1. serine hydrolases（丝氨酸水解酶）
2. PETase / cutinase-like PET hydrolases（PET 水解酶 / 角质酶类似 PET 水解酶）
3. CALB / lipase / esterase benchmark systems（南极假丝酵母脂肪酶 B / 脂肪酶 / 酯酶）

Possible later systems:

4. PUase / polyurethane hydrolase（聚氨酯水解酶）
5. nylonase（尼龙酶）
6. metal hydrolases（金属水解酶）
7. laccase / redox enzymes（漆酶 / 氧化还原酶，后期跨催化类型探索）
8. carbonic anhydrase（碳酸酐酶，后期金属酶/质子转移基准）

Recommended first reaction step:

```text
serine_hydrolase_acylation_TS1
```

This step includes:

- Ser-Oγ nucleophilic attack on the substrate carbonyl carbon.
- proton transfer through His.
- oxyanion-hole stabilization of the developing negative charge.
- formation of a tetrahedral intermediate-like or acylation TS-like state.

## 4. Role of PLACER in Project 01

Project 01 should use **true TS** as the high-confidence state definition, but PLACER can still be used as a conformer generation and screening interface.

Recommended workflow:

```text
true/refined TS template + GS complex
→ PLACER-generated or PLACER-screened GS/TS-compatible conformer ensemble
→ catalytic geometry filtering
→ clustering / representative conformer selection
→ restrained MM, DFT, or QM/MM relaxation
→ QM/MM or DFT barrier label
→ true-TS barrier model training
```

Important boundary:

```text
PLACER = conformer candidate generation / active-site pose screening
QM/MM or DFT = energy and barrier label source
true/refined TS = Project 01 supervision condition
```

Do not treat raw PLACER poses as final quantum truth. The repository should store PLACER metadata and links to external results, but not large raw PLACER ensemble dumps.

## 5. Model inputs

For each enzyme variant, generated design, or conformer-level state pair:

### 5.1 Protein identity input

Use protein model embedding（蛋白嵌入表示）or structure-derived embedding.

```text
sequence or structure
→ h_protein
```

This represents:

- mutation context（突变背景）
- sequence constraints（序列约束）
- fold prior（折叠先验）
- long-range residue context（远程残基上下文）

### 5.2 True GS/TS state input

Use known or refined state-pair structures:

```text
GS complex
true/refined TS complex
→ h_GS, h_TS_true
```

State structures can come from:

- experimental structures plus modeled substrate.
- hand-built reaction templates.
- TS analogs or intermediate-like templates.
- QM/MM constrained optimization.
- DFT cluster optimization.
- reaction-coordinate scans.
- PLACER-screened conformers refined by external quantum workflows.

### 5.3 Reaction geometry features

Recommended true-TS/GS features:

```text
forming bond length
breaking or weakening bond length
nucleophilic attack angle
proton transfer distance
oxyanion-hole hydrogen-bond distances
catalytic triad geometry
leaving-group alignment
substrate carbonyl orientation
metal coordination distances, if applicable
hydrogen-bond network descriptors
```

### 5.4 Electrostatic and charge features

Recommended features:

```text
q_TS
q_GS
Δq = q_TS - q_GS
electrostatic potential at reaction atoms
reaction-projected electric field / potential
active-site dipole or field proxy
```

The reaction-projected electrostatic term can be written as:

```text
S_field = Σ_i Δq_i * φ(r_i)
```

Where:

- `Δq_i = q_i_TS - q_i_GS`
- `φ(r_i)` is the enzyme electrostatic potential at the reaction atom/probe point.

### 5.5 Conformer ensemble features

Use PLACER or external conformer sampling outputs to summarize:

```text
productive conformer fraction
soft-min ensemble score
active-site side-chain rotamer state
ligand pose cluster
GS/TS-compatible conformer count
geometry-filter pass rate
```

### 5.6 QM/MM or DFT supervision labels

External quantum calculations should provide:

```text
E_GS
E_TS
G_GS, if available
G_TS, if available
ΔE_TS-GS
ΔG‡
ΔΔE‡ vs reference
ΔΔG‡ vs reference
label_tier
protocol_id
```

Every label must record its `protocol_id`, because barrier values can depend on QM region, basis, functional, embedding, sampling, constraints, protonation state, and reaction coordinate.

## 6. Model outputs

Primary outputs:

```text
predicted ΔG‡
predicted ΔΔG‡ vs WT/reference
within-template pairwise ranking
```

Secondary outputs:

```text
TS stabilization score
GS stabilization / binding score
reaction-projected electrostatic score
geometry feasibility score
productive conformer probability
variant or design ranking
prediction uncertainty
recommended next QM/MM samples for active learning
```

## 7. Recommended architecture

```text
Protein sequence/structure
        ↓
Protein encoder, e.g. ESM or structure embedding
        ↓
h_protein

GS complex + PLACER/QM/MM conformer metadata
        ↓
GS state encoder
        ↓
h_GS

true/refined TS complex + PLACER/QM/MM conformer metadata
        ↓
true-TS state encoder
        ↓
h_TS_true

Reaction template + physical features
        ↓
reaction-state feature encoder
        ↓
z_rxn, f_phys

Δ-state fusion
        ↓
[h_TS_true - h_GS, h_TS_true, h_GS, h_protein, z_rxn, f_phys]
        ↓
TrueTSBarrierScorer
        ↓
ΔG‡ / ΔΔG‡ / ranking / mechanism components
```

Important: for Project 01, true/refined TS can be an input because the task is mechanism-conditioned, high-confidence barrier learning.

## 8. Suggested physical decomposition

Use a modular score before training any complex model:

```text
S_total =
  w_field   * S_field
+ w_geom    * S_geometry
+ w_ens     * S_ensemble
- w_strain  * S_strain
- w_nonprod * S_nonproductive_binding
```

Where:

- `S_field`: reaction-projected electric-field or electrostatic score（反应投影电场/静电分数）
- `S_geometry`: catalytic geometry score（催化几何分数）
- `S_ensemble`: productive conformer population（有效催化构象比例）
- `S_strain`: conformational or ligand strain penalty（构象/配体应变惩罚）
- `S_nonproductive_binding`: nonproductive pose penalty（非反应性结合惩罚）

`S_product` or product-release terms can be added later, but they are not required for first-stage TS-GS barrier learning.

## 9. Data schema

Each sample should correspond to one conformer-level state pair:

```text
enzyme × variant/design × reaction step × conformer × GS/TS pair
```

Example manifest:

```yaml
sample_id: PETase_S238F_conf_017_acylation_TS1

enzyme:
  enzyme_id: PETase
  catalytic_class: serine_hydrolase
  family_id: PETase_like_cutinase
  variant_id: S238F
  sequence: path_or_string
  structure_pdb: path
  protein_embedding_path: path_or_null

reaction:
  reaction_template_id: serine_hydrolase_acylation_TS1
  reaction_step: acylation
  substrate_id: PET_model_fragment
  forming_bond: [Ser_Ogamma, substrate_carbonyl_C]
  breaking_bond: [carbonyl_C, leaving_group_O_or_N]
  proton_transfer: [Ser_Ogamma_H, His_N]

conformer:
  conformer_id: conf_017
  conformer_source: PLACER_screened_true_TS
  placer_run_id: path_or_string_or_null
  placer_score: float_or_null
  geometry_filter_passed: true
  cluster_id: cluster_03

states:
  GS:
    raw_complex_pdb: path_or_null
    relaxed_complex_pdb: path_or_null
    qmmm_structure_pdb: path_or_null
    charges_path: path_or_null
  TS:
    true_ts_source: qmmm_refined_or_dft_refined_or_template
    raw_ts_pdb: path_or_null
    relaxed_ts_pdb: path_or_null
    qmmm_ts_pdb: path_or_null
    charges_path: path_or_null

labels:
  delta_G_dagger: float_or_null
  delta_E_TS_GS: float_or_null
  delta_delta_G_dagger_vs_reference: float_or_null
  reference_variant_id: WT_or_other_reference
  label_tier: QM_MM_scan_or_DFT_cluster_or_xTB_or_proxy
  protocol_id: qmmm_protocol_v001

features:
  geometry_json: path
  electrostatic_json: path
  ensemble_summary_json: path
  physical_feature_table_row_id: string_or_null
```

## 10. Training strategy

Start simple and preserve interpretability.

### Stage 1: physical feature baseline

Use manually computed physical features and train:

- linear regression（线性回归）
- ridge regression（岭回归）
- random forest（随机森林）
- gradient boosting, optional
- small MLP（小型多层感知机）, optional

Targets:

```text
ΔG‡
ΔE_TS-GS
ΔΔG‡ relative to WT/reference
pairwise ranking within the same reaction template
```

### Stage 2: true-TS embedding fusion

Add protein embedding, GS state embedding, true-TS state embedding, and PLACER-derived ensemble summaries to the physical baseline.

### Stage 3: contrastive ranking

Use pairwise mutant/design ranking:

```text
variant A lower barrier than variant B
```

Loss:

```text
L_rank = max(0, margin + score_A - score_B)
```

where lower score means lower predicted barrier.

### Stage 4: uncertainty and active learning

Add uncertainty outputs so the model can recommend new QM/MM calculations:

```text
high predicted benefit + high uncertainty
→ priority sample for external QM/MM or DFT calculation
```

## 11. Compute plan

### GPU tasks

- protein embedding generation
- PLACER inference / ensemble generation or screening, if available on GPU
- state embedding extraction
- neural model training

### CPU / external compute tasks

- QM/MM conformational sampling
- QM/MM or DFT restrained optimization
- reaction-coordinate scans
- transition-state refinement or TS-like constrained optimization
- DFT cluster calculations
- classical electrostatic calculations

First-round target scale:

```text
1-3 initial enzyme systems: serine hydrolase / PETase / CALB
20-100 variants/designs per system
100-1000 PLACER or external conformers per variant before filtering
3-10 representative GS/TS state pairs per variant for external quantum labeling
100-1000 high-quality QM/MM or DFT labels for the first model round
```

Do not attempt exhaustive QM/MM free-energy surfaces in the first version.

## 12. Evaluation

Compare against:

1. embedding-only barrier prediction
2. GS-only scoring without true TS
3. geometry-only catalytic scoring
4. docking/binding-score-style scoring
5. PLACER ensemble score without TS-GS physics
6. full true-TS-aware score

Metrics:

```text
Pearson / Spearman correlation with computed ΔG‡ or ΔΔG‡
pairwise ranking accuracy
top-K enrichment for low-barrier variants/designs
leave-variant-out performance
leave-enzyme-out performance within a catalytic class
uncertainty calibration, if available
```

Experimental activity is not the main metric for the first stage. It can be reported later as optional external validation when assay conditions are comparable.

## 13. Deliverable

A working Project 01 implementation should produce:

```text
input variant/design structures + GS + true/refined TS + protocol metadata
→ conformer-level feature table
→ true-TS barrier scorer
→ ΔG‡ / ΔΔG‡ predictions
→ ranked enzyme variants/designs
→ mechanism decomposition report
→ exported teacher labels for Project 02
```

The key scientific claim should be:

> With true/refined TS supervision and controlled QM/MM or DFT barrier labels, enzyme-specific or catalytic-class-specific models can learn how protein environments stabilize TS relative to GS, enabling computed barrier prediction and variant/design ranking.
