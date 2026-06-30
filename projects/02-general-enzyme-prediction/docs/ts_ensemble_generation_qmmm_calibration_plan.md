# TS ensemble generation and QM/MM calibration research plan

Created: 2026-06-30  
Scope: `one enzyme × one substrate × one reaction step` first-stage plan, with a later path to predicted-TS screening and catalytic-class transfer.  
Owner context: this document refines the existing Project 01 / Project 02 plan into a concrete TS-generation work package.

---

## 0. Executive summary

The immediate research goal is not to predict all enzyme activity. The immediate goal is narrower and more testable:

```text
For one fixed enzyme, one fixed substrate, and one fixed reaction step,
learn whether a mutant active site preferentially stabilizes the transition-state ensemble
relative to the ground-state ensemble,
and use that learned relationship to rank new mutants by computed barrier.
```

The proposed first-stage pipeline is:

```text
MPNN-generated mutants
+ fixed enzyme/substrate/reaction template
+ QM/MM-derived reference TS ensemble
+ TS-generator-derived predicted TS ensemble
→ constrained placement / docking into each mutant pocket
→ matched GS–TS state-pair construction
→ geometry filtering and ensemble aggregation
→ QM/MM labels for representative state pairs
→ small interpretable barrier model
→ mutant ranking and active-learning selection
```

The most important scientific claim should be:

```text
A generated TS conformer ensemble can be used for fast screening only if it preserves
the QM/MM-derived TS ensemble's information about mutant ranking and TS-over-GS stabilization.
```

Therefore this work package has two tightly coupled subprojects:

```text
Phase 1A: QM/MM-derived TS teacher ensemble
Phase 1B: TS-generator ensemble calibration and replacement test
```

The output is not a generic docking score. The output is a barrier-relevant score:

```text
predicted ΔG‡ / ΔΔG‡ / pairwise mutant ranking
```

The central physical quantity is:

```text
ΔG‡ = G_TS-bound - G_GS-bound
```

or, in a first computational round:

```text
ΔE‡ = E_TS-bound - E_GS-bound
```

The model must learn **preferential TS stabilization relative to GS**, not absolute TS binding.

---

## 1. First-principles starting point

### 1.1 Rate is exponentially controlled by the activation free energy

For a fixed reaction step under comparable conditions, catalytic ranking is primarily controlled by the activation free energy:

```text
k ∝ exp(-ΔG‡ / RT)
```

A small barrier error is large in rate space. At room temperature, roughly:

```text
1.36 kcal/mol ≈ 10-fold rate change
```

Therefore the model should not be trained first on broad experimental activity labels unless assay conditions are clean and comparable. The first-stage label should be computed and protocol-controlled:

```text
ΔG‡, ΔE‡, ΔΔG‡, ΔΔE‡
```

### 1.2 Enzyme catalysis is differential stabilization

A mutant can bind the TS strongly and still not lower the barrier if it also stabilizes the GS. The relevant effect is:

```text
differential stabilization = stabilization(TS) - stabilization(GS)
```

Equivalently:

```text
barrier lowering = ΔG‡_reference - ΔG‡_mutant
```

The first-stage model must therefore consume matched GS and TS state pairs whenever possible.

### 1.3 A transition state is not a normal ligand pose

A TS-like object has reaction-coordinate content:

```text
forming bonds
breaking or weakening bonds
proton transfers
charge redistribution
orbital alignment / attack geometry
strained substrate conformations
active-site electric-field alignment
```

A ligand-only docking model may find a bound pose that is stable but nonproductive. This is especially dangerous for enzyme design because many mutations improve binding while destroying reaction geometry.

### 1.4 An enzyme-specific TS ensemble is not one structure

For a fixed enzyme and substrate, the transition state should be treated as an ensemble of reaction-compatible structures:

```text
TS ensemble = {TS_i}
```

Different members can differ in:

```text
substrate dihedral
forming-bond distance
leaving-group distance
proton-transfer coordinate
oxyanion orientation
catalytic residue rotamers
local water arrangement
nearby polar residue geometry
```

For mutant screening, the useful question is not:

```text
Does this mutant stabilize one hand-picked TS pose?
```

The useful question is:

```text
Does this mutant support a productive population of TS-like poses,
and does this population lower TS relative to GS?
```

---

## 2. Definitions and boundaries

### 2.1 Ground state, GS

In this document, GS means the bound reactant complex for the chosen reaction step:

```text
protein_mutant + substrate/reactant pose + catalytic residues before bond making/breaking
```

For serine hydrolase acylation TS1, GS usually means:

```text
Ser-Oγ is positioned for nucleophilic attack but has not formed the covalent acyl bond.
The carbonyl carbon is near planar.
His is positioned to accept the Ser proton.
The oxyanion hole may be preorganized but the oxyanion is not fully developed.
```

### 2.2 Transition state, TS

TS means the reaction-state geometry around the activation barrier. In practice the first-stage plan may use several tiers:

```text
true TS stationary point
barrier-top QM/MM scan snapshot
constrained TS-like geometry
tetrahedral intermediate-like proxy
TS analog / reaction-template-derived proxy
```

The exact tier must be recorded for every sample.

### 2.3 TS ensemble, TSE

TSE means a set of TS-like conformers for the same reaction template:

```text
E_QMMM = QM/MM-derived reference TS ensemble
E_gen  = TS-generator-derived predicted TS ensemble
E_wrong = negative-control TS ensemble
```

The first project goal is to determine whether `E_gen` can reproduce the screening signal of `E_QMMM`.

### 2.4 TS-bound complex

A TS-bound complex is:

```text
protein mutant + TS-like reaction-state object placed in active site
```

This object may include only substrate atoms, or it may include a reactive fragment that contains catalytic residues. For serine hydrolases, a ligand-only TS is often insufficient. A better object includes:

```text
substrate carbonyl fragment
Ser-Oγ nucleophile geometry
developing Ser-Oγ–C bond
oxyanion oxygen orientation
His proton-transfer geometry
leaving-group vector
```

### 2.5 Predicted TS generator

A TS generator is any external model or workflow that proposes TS-like structures:

```text
reaction prior + reactant/product + atom mapping
→ predicted TS conformer ensemble
```

The TS generator is not assumed to be correct. It is a proposal mechanism. It becomes useful only after calibration against QM/MM-derived TS ensembles and QM/MM barrier labels.

### 2.6 Project boundary

This document does not implement quantum chemistry, docking engines, TS generator internals, or large trajectory storage inside the repository.

The repository should store:

```text
schemas
manifests
feature extractors
small examples
training scripts
calibration reports
protocol documentation
```

External compute should store:

```text
raw QM/MM outputs
large docking outputs
large TS-generator outputs
trajectories
model checkpoints
```

---

## 3. The first-stage scientific question

### 3.1 Primary question

```text
For a fixed enzyme, fixed substrate, and fixed reaction template,
can TS ensemble docking/placement features predict QM/MM ΔΔG‡ across MPNN-generated mutants?
```

### 3.2 Secondary question

```text
Can a TS-generator-derived ensemble replace the QM/MM-derived reference TS ensemble
for fast mutant screening in this fixed system?
```

### 3.3 What would count as success?

Minimum success:

```text
1. The model ranks mutants better than GS-only, docking-only, and geometry-only baselines.
2. The generated TS ensemble gives similar top-K enrichment to the QM/MM-derived TS ensemble.
3. Wrong-TS and shuffled-reaction controls perform substantially worse.
4. The model identifies high-uncertainty or high-promise mutants for new QM/MM labels.
```

Strong success:

```text
1. E_gen covers the major E_QMMM reaction-center clusters.
2. E_gen-derived features correlate with E_QMMM-derived features across mutants.
3. E_gen-derived model predictions correlate with QM/MM ΔΔG‡ labels.
4. Top-K mutants selected by E_gen are enriched for low computed barriers.
5. Failure modes are interpretable by reaction geometry, electrostatics, or ensemble collapse.
```

---

## 4. “Grill me” adversarial review

This section states the objections a hard reviewer would raise and how the plan must answer them.

### 4.1 Objection: You are learning binding, not catalysis

Risk:

```text
The model may learn that polar/tight pockets bind charged TS-like structures.
That does not prove barrier lowering.
```

Required answer:

```text
Every barrier label must use matched GS–TS pairs.
Feature design must include TS-GS difference features.
Baselines must include TS-only and GS-only models.
```

Kill criterion:

```text
If TS-only features perform as well as TS-GS differential features,
and wrong-TS controls also perform well,
the model is not learning TS-specific catalysis.
```

### 4.2 Objection: The generated TS ensemble is chemically wrong

Risk:

```text
The TS generator may produce conformers that are plausible small-molecule TS geometries
but incompatible with the enzyme reaction coordinate.
```

Required answer:

```text
Compare E_gen against E_QMMM in reaction-coordinate space before docking.
Use geometric coverage, cluster recall, and wrong-cluster filtering.
Do not allow E_gen to enter screening unless it passes coverage gates.
```

Kill criterion:

```text
If E_gen does not cover the main E_QMMM clusters,
the generator cannot be used as a replacement in this system.
```

### 4.3 Objection: Docking creates nonreactive poses

Risk:

```text
Docking may place a TS-like conformer in the pocket with good shape complementarity
but wrong Ser/His/oxyanion/leaving-group geometry.
```

Required answer:

```text
Use catalytic-constraint docking or post-docking geometry filters.
Score productive fraction, not just best pose.
```

Kill criterion:

```text
If most top-scoring poses fail reaction geometry filters,
the docking protocol is not usable for training.
```

### 4.4 Objection: QM/MM labels are inconsistent

Risk:

```text
Different protonation states, QM regions, constraints, or methods can change barriers.
```

Required answer:

```text
Every label must carry protocol_id.
First-stage training should use one controlled protocol.
Cross-protocol labels must be tiered and not mixed blindly.
```

Kill criterion:

```text
If repeated calculations under the same protocol are not reproducible within the target tolerance,
do not train the model on those labels yet.
```

### 4.5 Objection: Mutant train/test leakage inflates performance

Risk:

```text
Conformers from the same mutant appear in train and test.
The model memorizes the mutant rather than learning barrier physics.
```

Required answer:

```text
Split by variant_id, not conformer_id.
Report leave-variant-out performance.
```

Kill criterion:

```text
If conformer-random split is good but variant-level split fails,
the model is not ready for mutant screening.
```

### 4.6 Objection: The model cannot generalize beyond this enzyme/substrate

Answer:

```text
Correct. Phase 1 does not claim broad generalization.
It only claims within-system mutant ranking.
Generalization requires reaction-template conditioning, multiple substrates, multiple enzyme families,
and anchor QM/MM labels for every new regime.
```

---

## 5. “Superpower” framing

The intended “superpower” is not a bigger neural network. It is a controlled feedback loop:

```text
TS generator proposes chemically plausible reaction-state ensembles.
QM/MM tells us which proposed states actually correspond to low barriers in the enzyme pocket.
A small model learns the pocket features that explain the QM/MM barriers.
The model then selects the next mutants and poses where QM/MM is most valuable.
```

This creates a compounding advantage:

```text
better TS ensemble proposals
→ better QM/MM label efficiency
→ better barrier surrogate
→ better active-learning selection
→ fewer wasted quantum calculations
```

The first useful superpower is narrow:

```text
For this one enzyme and one substrate,
screen many MPNN mutants without QM/MM on every mutant.
```

The later superpower is broader:

```text
Learn when a predicted TS ensemble can replace a QM/MM-derived TS ensemble
within a reaction template.
```

---

## 6. Phase 1A — QM/MM-derived TS teacher ensemble

### 6.1 Objective

Generate a reference TS ensemble for the chosen enzyme/substrate/reaction step.

```text
WT/reference enzyme + fixed substrate + reaction template
→ QM/MM-derived TS ensemble E_QMMM
```

This ensemble is the high-confidence reference for the first model and for calibrating generated TS ensembles.

### 6.2 Input requirements

Minimum input:

```yaml
enzyme:
  enzyme_id: string
  reference_variant_id: WT_or_design_reference
  structure_pdb: path
  sequence: path_or_string
  catalytic_residues:
    nucleophile: residue_id
    general_base: residue_id
    acid_or_hbond_partner: residue_id_or_null
    oxyanion_hole_donors: list[residue_id_or_atom_id]

substrate:
  substrate_id: string
  structure_file: path
  atom_names: consistent
  protonation_state: documented
  charge_state: documented

reaction_template:
  reaction_template_id: string
  reaction_step: string
  atom_mapping: path
  forming_bonds: list
  breaking_or_weakening_bonds: list
  proton_transfers: list
  required_geometry_features: list
  expected_delta_q: path_or_null
```

### 6.3 Recommended serine hydrolase TS1 reaction definition

For serine hydrolase acylation TS1:

```yaml
reaction_template_id: serine_hydrolase_acylation_TS1
reaction_step: acylation
forming_bonds:
  - [Ser_Ogamma, substrate_carbonyl_C]
breaking_or_weakening_bonds:
  - [substrate_carbonyl_C, leaving_group_O_or_N]
proton_transfers:
  - [Ser_Ogamma_H, His_Nepsilon_or_Ndelta]
stabilizing_interactions:
  - oxyanion_hole_Hbond_to_carbonyl_O
  - His_Asp_or_His_Glu_Hbond
required_geometry_features:
  - Ser_Ogamma_C_distance
  - Oγ_C_O_attack_angle
  - carbonyl_tetrahedralization
  - carbonyl_O_oxyanion_hole_distances
  - SerH_HisN_distance
  - leaving_group_bond_distance
  - leaving_group_alignment_angle
```

### 6.4 QM/MM protocol tiers

The first-stage plan should allow label tiers, but training should initially use one tier or explicitly model tier effects.

```text
Tier 0: geometry proxy only
  no QM/MM label
  used only for filtering

Tier 1: QM/MM single-point ΔE_TS-GS
  cheap
  useful for rough active learning

Tier 2: restrained QM/MM optimization
  recommended first training tier
  GS and TS-like state locally optimized under consistent constraints

Tier 3: reaction-coordinate scan / NEB / string / umbrella
  higher confidence
  used to validate top candidates and calibrate Tier 1/2

Tier 4: free-energy PMF / ensemble barrier
  expensive
  not required for first version
```

### 6.5 Reaction-coordinate scan design

For serine hydrolase acylation TS1, use collective variables such as:

```text
CV1 = d(Ser_Oγ, C_carbonyl)
CV2 = d(C_carbonyl, leaving_group_O_or_N)
CV3 = d(Ser_H, His_N)
CV4 = carbonyl_tetrahedralization_coordinate
CV5 = oxyanion_O_to_oxyanion_hole_donor_distance(s)
```

Candidate scan coordinates:

```text
RC1 = d(Ser_Oγ, C_carbonyl) - d(C_carbonyl, leaving_group)
RC2 = d(Ser_H, His_N) - d(Ser_Oγ, Ser_H)
```

The final protocol should specify:

```yaml
qmmm_protocol_id: qmmm_serhyd_acyl_TS1_v001
qm_region:
  - substrate_reactive_fragment
  - Ser_side_chain_from_CB_or_Ogamma
  - His_side_chain
  - Asp_or_Glu_side_chain_if_part_of_charge_relay
  - optional_oxyanion_hole_backbone_fragments
mm_region:
  - remaining_protein
  - crystallographic_or_modeled_waters
  - ions
link_atoms: documented
method:
  qm_method: DFT_or_semiempirical_or_ab_initio
  basis_or_param: documented
  mm_forcefield: documented
  electrostatic_embedding: true_or_false
constraints:
  - reaction_coordinate_constraints
  - backbone_restraints_if_any
solvent_boundary: documented
protonation_states: documented
outputs:
  - optimized_GS
  - optimized_TS_like
  - energies
  - charges_GS
  - charges_TS
  - delta_q
```

### 6.6 Extracting E_QMMM

From scans or constrained optimizations:

```text
1. Identify barrier-top or TS-like windows.
2. Extract snapshots near the top.
3. Remove chemically invalid structures.
4. Cluster in reaction-coordinate space, not only Cartesian RMSD.
5. Select representative TS conformers.
```

Recommended clustering features:

```text
Ser_Oγ_C_distance
C_leaving_group_distance
SerH_HisN_distance
Oγ_C_O_attack_angle
carbonyl_tetrahedralization
oxyanion-hole H-bond distances
substrate key dihedrals
His rotamer descriptor
```

Output:

```yaml
ts_ensemble_id: E_QMMM_reference_v001
source: qmmm_scan_or_constrained_optimization
reaction_template_id: serine_hydrolase_acylation_TS1
n_raw_snapshots: int
n_filtered_snapshots: int
n_clusters: int
representatives:
  - conformer_id: TSQ_001
    structure_pdb: path
    reaction_geometry_json: path
    qmmm_energy: float_or_null
    cluster_population: int
    protocol_id: qmmm_serhyd_acyl_TS1_v001
```

---

## 7. Phase 1B — TS-generator ensemble calibration

### 7.1 Objective

Test whether a TS generator can replace or approximate the QM/MM-derived reference TS ensemble for this fixed enzyme/substrate/reaction step.

```text
reaction prior + substrate + product/intermediate + atom mapping
→ TS generator
→ E_gen
→ compare with E_QMMM before and after docking
```

### 7.2 Required TS-generator metadata

Every generated ensemble should store:

```yaml
predicted_ts_ensemble:
  ensemble_id: E_gen_modelX_v001
  ts_prediction_model_id: string
  model_version: string
  input_reaction_template_id: string
  input_reactant_structure: path
  input_product_or_intermediate_structure: path_or_null
  atom_mapping: path
  n_requested: int
  n_generated: int
  n_valid_after_chemistry_filter: int
  confidence_scores_path: path_or_null
  generation_parameters:
    seed: int
    temperature_or_noise_scale: float_or_null
    constraints: object_or_null
  output_structures_dir: external_path
```

### 7.3 Pre-docking chemistry filters

Before protein docking, remove generated structures that violate reaction chemistry.

For serine hydrolase TS1:

```text
filter_1: Ser/nucleophile-compatible attack geometry if Ser fragment is included
filter_2: carbonyl carbon tetrahedralization within allowed range
filter_3: forming bond distance within TS-like range
filter_4: leaving-group distance not reactant-like unless intended
filter_5: proton-transfer coordinate plausible if proton transfer is modeled
filter_6: no impossible valence or atom identity mismatch
filter_7: substrate stereochemistry and atom mapping preserved
```

If the TS generator outputs ligand-only conformers, the pre-docking filter should still check:

```text
carbonyl fragment geometry
leaving-group vector
key substrate dihedrals
oxyanion oxygen orientation proxy
```

### 7.4 Geometry coverage against E_QMMM

Compare `E_gen` and `E_QMMM` in reaction-coordinate space.

Suggested distance:

```text
D(TS_a, TS_b) =
  w1 * zscore(|d_forming_a - d_forming_b|)
+ w2 * zscore(|d_breaking_a - d_breaking_b|)
+ w3 * zscore(|angle_attack_a - angle_attack_b|)
+ w4 * zscore(|tetrahedralization_a - tetrahedralization_b|)
+ w5 * zscore(dihedral_distance)
+ w6 * zscore(proton_transfer_distance_difference)
```

Metrics:

```text
cluster_recall:
  fraction of E_QMMM clusters with at least one E_gen conformer within threshold

cluster_precision:
  fraction of E_gen clusters near any E_QMMM cluster

nearest_neighbor_distance:
  for each E_QMMM representative, nearest E_gen distance

MMD_or_Wasserstein_distance:
  distribution-level distance in reaction-coordinate space

outlier_fraction:
  fraction of E_gen conformers that are chemically invalid or far from all E_QMMM clusters
```

Gate:

```text
E_gen should not proceed as a replacement unless it covers the dominant E_QMMM clusters.
If it only produces a single collapsed TS mode while E_QMMM has multiple productive modes,
it may still be usable for limited screening but not as a general substitute.
```

### 7.5 Negative-control TS ensembles

Construct at least three negative controls:

```text
E_wrong_reaction:
  TS ensemble from the wrong reaction step or wrong atom mapping.

E_wrong_substrate:
  TS ensemble from a related but mismatched substrate.

E_GS_like:
  reactant/ground-state conformers used as if they were TS conformers.

E_random_pose:
  random or weakly constrained placement into the pocket.
```

A valid TS-generator signal must beat these controls.

---

## 8. Mutant placement / docking without PLACER

### 8.1 Objective

For each mutant and each TS ensemble, place TS conformers into the active site and generate reaction-compatible state pairs.

```text
mutant pocket + TS ensemble
→ constrained placement / docking / local relaxation
→ geometry-filtered TS-bound poses
→ matched GS poses
```

This plan does not require PLACER. Any method can be used if it supports the necessary catalytic constraints and returns traceable pose metadata.

### 8.2 Minimum constraints

For serine hydrolase TS1:

```text
constraint 1: Ser_Oγ to carbonyl_C distance
constraint 2: Ser_Oγ-C-O attack angle
constraint 3: oxyanion O to oxyanion-hole donor distances
constraint 4: His_N to Ser_H or proton donor/acceptor geometry
constraint 5: His-Asp/Glu charge-relay distance
constraint 6: leaving-group vector alignment
constraint 7: no severe clashes around reaction center
```

### 8.3 Placement outputs

For every mutant and ensemble:

```yaml
placement_result:
  variant_id: string
  enzyme_id: string
  ensemble_id: string
  reaction_template_id: string
  n_input_conformers: int
  n_raw_poses: int
  n_geometry_pass: int
  n_selected_for_qmmm: int
  pose_summary_path: path
  representative_pose_ids:
    - pose_id
```

Pose-level row:

```yaml
pose_id: string
variant_id: string
ensemble_id: string
ts_conformer_id: string
placement_method: constrained_docking_or_minimization
raw_pose_pdb: external_path
relaxed_pose_pdb: external_path_or_null
placement_score: float_or_null
reaction_geometry:
  forming_bond_distance: float
  breaking_bond_distance: float_or_null
  attack_angle: float
  oxyanion_hbond_distances: list[float]
  proton_transfer_distance: float_or_null
  leaving_group_alignment: float_or_null
filters:
  geometry_pass: bool
  clash_pass: bool
  strain_pass: bool
  protonation_pass: bool
cluster_id: string_or_null
```

### 8.4 Ensemble aggregation

For each mutant, compute:

```text
productive_TS_fraction = n_geometry_pass / n_input_conformers
topk_mean_geometry_score
softmin_TS_score
pose_cluster_entropy
dominant_cluster_population
best_productive_pose_score
mean_productive_pose_score
number_of_distinct_productive_clusters
```

Soft-min:

```text
S_softmin = -τ * log(sum_i exp(-S_i / τ))
```

Use the same aggregation for `E_QMMM`, `E_gen`, and controls.

---

## 9. Matched GS–TS state-pair construction

### 9.1 Why matched pairs are required

A TS pose without a comparable GS pose cannot define a barrier. The label must represent:

```text
E_TS - E_GS
```

not:

```text
E_TS alone
```

### 9.2 Matched GS construction strategies

Allowed strategies:

```text
Strategy A: scan-derived pair
  GS and TS come from the same QM/MM reaction-coordinate path.

Strategy B: reverse-mapped pair
  TS pose is converted back to reactant geometry while preserving the binding frame,
  followed by local relaxation.

Strategy C: common-frame docking pair
  GS and TS are both placed using the same reaction-center frame and then locally relaxed.

Strategy D: ensemble-matched cluster pair
  GS and TS representatives are matched by substrate orientation, active-site rotamers,
  and reaction-center atom correspondence.
```

Disallowed for first-stage primary labels:

```text
Unmatched global GS docking pose paired with unrelated TS pose.
```

### 9.3 State-pair manifest

```yaml
sample_id: enzyme_variant_reaction_conformer_pair_id
enzyme_id: string
variant_id: string
substrate_id: string
reaction_template_id: string
reaction_step: string
ts_ensemble_id: E_QMMM_or_E_gen
ts_conformer_id: string
pose_id: string
state_pair_id: string

GS:
  construction_strategy: scan_derived_or_reverse_mapped_or_common_frame
  raw_complex_pdb: external_path
  relaxed_complex_pdb: external_path_or_null
  qmmm_input_path: external_path_or_null

TS:
  ts_source: qmmm_reference_or_ts_generator_or_control
  raw_complex_pdb: external_path
  relaxed_complex_pdb: external_path_or_null
  qmmm_input_path: external_path_or_null

labels:
  delta_E_TS_GS: float_or_null
  delta_G_dagger: float_or_null
  delta_delta_E_vs_reference: float_or_null
  delta_delta_G_vs_reference: float_or_null
  label_tier: geometry_proxy_or_qmmm_sp_or_qmmm_opt_or_qmmm_scan_or_pmf
  protocol_id: string_or_null

features:
  geometry_features_path: path
  electrostatic_features_path: path_or_null
  ensemble_summary_path: path
  pocket_features_path: path_or_null
```

---

## 10. QM/MM labeling plan for mutants

### 10.1 Sampling strategy

Do not run exhaustive QM/MM on every pose.

Recommended first round:

```text
Mutants:
  50–100 MPNN variants
  include WT/reference
  include sequence-similarity bins
  include active-site and second-shell variants
  include predicted good/bad variants if available

Per mutant:
  generate placements for E_QMMM
  generate placements for E_gen
  generate placements for controls
  select 3–10 representative GS–TS pairs for QM/MM labeling

Total first-round labels:
  100–500 useful labels for initial model
  500–1000 if compute allows
```

### 10.2 Representative selection

Select state pairs using:

```text
1. best productive geometry
2. dominant productive cluster
3. alternative productive cluster
4. borderline pose near geometry threshold
5. high-disagreement pose between E_QMMM and E_gen
6. high-uncertainty pose from current model
```

Avoid selecting only the best-looking poses. The model needs both good and bad examples.

### 10.3 Label fields

Every label row must include:

```text
sample_id
variant_id
enzyme_id
substrate_id
reaction_template_id
ts_ensemble_id
ts_source
state_pair_id
E_GS
E_TS
G_GS if available
G_TS if available
delta_E_TS_GS
delta_G_dagger if available
delta_delta_E_vs_reference
delta_delta_G_vs_reference if available
label_tier
protocol_id
qmmm_status
failure_reason if failed
```

### 10.4 Failure labels are useful

Failed QM/MM or failed geometry cases should not simply disappear. Record:

```text
failed_to_converge
proton_transfer_broke
substrate_left_pocket
Ser_C_distance_collapsed
His_rotamer_flipped
oxyanion_hole_lost
water_invaded_reaction_center
severe_clash_after_relaxation
```

These failure modes can train a feasibility classifier.

---

## 11. Feature design

### 11.1 Geometry features

For each state pair:

```text
forming_bond_distance_TS
forming_bond_distance_GS
breaking_bond_distance_TS
breaking_bond_distance_GS
attack_angle_TS
attack_angle_GS
tetrahedralization_TS
tetrahedralization_GS
proton_transfer_distance_TS
proton_transfer_distance_GS
oxyanion_hbond_min_distance_TS
oxyanion_hbond_min_distance_GS
oxyanion_hbond_count_TS
oxyanion_hbond_count_GS
His_Asp_distance_TS
His_Asp_distance_GS
leaving_group_alignment_TS
leaving_group_alignment_GS
substrate_key_dihedrals_TS
substrate_key_dihedrals_GS
```

### 11.2 Differential geometry features

```text
delta_forming_bond_distance = TS - GS
delta_attack_angle = TS - GS
delta_oxyanion_hbond_count = TS - GS
delta_proton_transfer_distance = TS - GS
delta_substrate_strain = TS - GS
```

### 11.3 Electrostatic features

At minimum:

```text
enzyme electrostatic potential at reaction atoms in GS
enzyme electrostatic potential at reaction atoms in TS
delta_q = q_TS - q_GS
reaction_projected_field_score = sum_i(delta_q_i * phi_i)
local charged residue counts
distance-weighted charge near oxyanion
dipole alignment along reaction vector
```

### 11.4 Pocket and environment features

```text
pocket volume around reaction center
solvent exposure of oxyanion
buried water count near reaction center
hydrophobic contact count
polar contact count
clash score
side-chain rotamer identity for catalytic residues
second-shell residue distances
active-site backbone RMSD to reference
```

### 11.5 Ensemble features

For each mutant and each ensemble source:

```text
productive_fraction
geometry_pass_rate
n_productive_clusters
dominant_cluster_fraction
pose_entropy
softmin_geometry_score
softmin_energy_proxy
topk_mean_score
best_score
mean_score
score_variance
E_gen_vs_E_QMMM_cluster_overlap
```

### 11.6 Protein identity features

For Phase 1, keep these simple:

```text
variant_id
mutation list
active-site mutation indicator
distance of mutated residues to reaction center
MPNN sequence identity bin
optional protein embedding path
```

Do not let protein identity features dominate without reaction-state controls.

---

## 12. Modeling strategy

### 12.1 First model should be small and interpretable

Recommended order:

```text
1. linear / ridge regression on physical features
2. random forest or gradient boosting
3. small MLP only after baseline behavior is understood
4. pairwise ranking loss after regression baseline exists
```

### 12.2 Primary targets

```text
delta_E_TS_GS
delta_G_dagger if available
delta_delta_E_vs_reference
delta_delta_G_vs_reference if available
pairwise ranking within same reaction template
```

### 12.3 Auxiliary targets

```text
geometry_feasibility
QM/MM convergence probability
productive_TS_fraction
field_score
TS_generator_reliability_score
```

### 12.4 Recommended losses

Regression:

```text
L_barrier = MSE(predicted_delta_delta_G, qmmm_delta_delta_G)
```

Ranking:

```text
L_rank = max(0, margin + score_A - score_B)
```

where lower score means lower barrier.

Feasibility classification:

```text
L_feas = binary_cross_entropy(predicted_feasible, qmmm_or_geometry_feasible)
```

Multi-task:

```text
L_total =
  λ_barrier * L_barrier
+ λ_rank    * L_rank
+ λ_feas    * L_feas
+ λ_field   * L_field_proxy
+ λ_geom    * L_geometry_proxy
+ λ_ens     * L_ensemble_proxy
```

### 12.5 Split policy

Never split by conformer for main evaluation.

Required splits:

```text
variant-level split
leave-variant-out
sequence-similarity-bin split
optional active-site-mutation-held-out split
```

Forbidden as primary metric:

```text
random conformer-level split
```

---

## 13. Calibration experiment: E_QMMM vs E_gen

### 13.1 Experiment matrix

For each mutant:

```text
score mutant with E_QMMM
score mutant with E_gen
score mutant with E_wrong_reaction
score mutant with E_GS_like
score mutant with random-pose control
construct matched GS features
compare to QM/MM labels
```

### 13.2 Metrics

Geometry coverage:

```text
E_gen cluster recall against E_QMMM
E_gen cluster precision
nearest-neighbor distance in reaction-coordinate space
outlier fraction
```

Mutant score consistency:

```text
Spearman(S_E_gen, S_E_QMMM)
KendallTau(S_E_gen, S_E_QMMM)
topK_overlap(E_gen, E_QMMM)
pairwise_ranking_agreement(E_gen, E_QMMM)
```

Barrier relevance:

```text
Spearman(pred_E_gen, ΔΔG‡_QMMM)
KendallTau(pred_E_gen, ΔΔG‡_QMMM)
topK_enrichment_low_barrier
pairwise_accuracy_against_QMMM_labels
```

Control separation:

```text
performance(E_gen) > performance(E_wrong)
performance(E_gen) > performance(E_GS_like)
performance(E_gen) > performance(docking_score_only)
performance(E_gen) > performance(GS_only)
```

### 13.3 Suggested go/no-go thresholds

These thresholds can be adjusted after the first data round.

Proceed with E_gen as screening input if:

```text
1. E_gen covers dominant E_QMMM clusters.
2. E_gen-vs-E_QMMM mutant score Spearman is positive and practically meaningful.
3. E_gen model beats GS-only and docking-only baselines on variant-level split.
4. Wrong-TS and shuffled controls are clearly worse.
5. Top-K enrichment is useful for selecting QM/MM candidates.
```

Pause and repair TS generator or docking protocol if:

```text
1. E_gen fails reaction-center coverage.
2. E_gen produces many chemically invalid conformers.
3. E_gen correlates with wrong-TS controls more than E_QMMM.
4. Model performance disappears under variant-level split.
5. Best-pose score works but ensemble productive fraction fails, suggesting false-positive docking.
```

---

## 14. Active-learning loop

### 14.1 First loop

```text
initial QM/MM labels
→ train baseline model
→ score unlabeled mutants
→ select new QM/MM candidates
```

Selection should include:

```text
top predicted low-barrier mutants
high uncertainty mutants
high disagreement between E_QMMM and E_gen
high disagreement between model and physical heuristic
representative mutants from under-sampled sequence-similarity bins
```

### 14.2 Acquisition score

Example:

```text
acquisition =
  w1 * predicted_barrier_improvement
+ w2 * uncertainty
+ w3 * ensemble_disagreement
+ w4 * diversity_score
- w5 * obvious_geometry_failure_penalty
```

### 14.3 Stop condition

Stop adding QM/MM labels when:

```text
top-K enrichment plateaus
model uncertainty is low for the candidate pool
new labels do not change ranking
or compute budget is exhausted
```

---

## 15. Required controls and ablations

### 15.1 Controls

```text
GS-only model
TS-only model
TS-GS differential model
docking-score-only model
geometry-only model
electrostatics-only model
ensemble-only model
wrong-TS ensemble
GS-like ensemble used as TS
random-pose control
shuffled reaction template
shuffled labels
```

### 15.2 Ablation questions

```text
Does adding TS features improve over GS-only?
Does adding TS-GS difference improve over TS-only?
Does adding electrostatics improve over geometry-only?
Does ensemble aggregation improve over best-pose?
Does E_gen approach E_QMMM performance?
Does wrong-TS destroy performance?
Does variant-level split remain valid?
```

### 15.3 Minimum reporting table

```text
model_name
input_features
ensemble_source
split_type
n_train_variants
n_test_variants
n_qmmm_labels
Spearman_delta_delta_G
KendallTau
pairwise_accuracy
topK_enrichment
wrong_TS_gap
notes
```

---

## 16. Data products

### 16.1 Manifests

Recommended files:

```text
data/manifest_variants.csv
data/manifest_reaction_template.yaml
data/manifest_ts_ensembles.csv
data/manifest_placements.csv
data/manifest_state_pairs.csv
data/qmmm_labels.csv
data/features_conformer_level.csv
data/features_variant_level.csv
data/model_predictions.csv
```

Large raw files stay outside GitHub. Manifests store paths and metadata.

### 16.2 TS ensemble manifest schema

```yaml
ensemble_id: string
ensemble_source: qmmm_reference_or_ts_generator_or_control
enzyme_id: string_or_null
substrate_id: string
reaction_template_id: string
generator_model_id: string_or_null
qmmm_protocol_id: string_or_null
n_conformers: int
n_valid_conformers: int
structure_dir: external_path
feature_table_path: path
metadata:
  atom_mapping_path: path
  reactant_path: path
  product_or_intermediate_path: path_or_null
  generation_seed: int_or_null
  notes: string
```

### 16.3 Variant-level feature row

```yaml
variant_id: string
enzyme_id: string
substrate_id: string
reaction_template_id: string
ensemble_source: E_QMMM_or_E_gen_or_control
productive_fraction: float
softmin_score: float
topk_mean_score: float
pose_entropy: float
best_geometry_score: float
mean_field_score: float_or_null
mean_delta_TS_GS_score: float_or_null
n_productive_clusters: int
qmmm_delta_delta_G_label: float_or_null
split_group: string
```

---

## 17. Implementation tasks for the repository

### 17.1 Suggested new package modules

For Project 02 or a shared package:

```text
src/general_enzyme_predictor/
  ts_ensemble.py
  ts_coverage.py
  ts_placement_features.py
  matched_state_pairs.py
  qmmm_calibration.py
  calibration_metrics.py
```

### 17.2 Suggested scripts

```text
scripts/build_ts_ensemble_manifest.py
scripts/compare_ts_ensembles.py
scripts/summarize_ts_placements.py
scripts/build_matched_state_pairs.py
scripts/merge_qmmm_labels.py
scripts/train_ts_generator_calibration_model.py
scripts/compare_qmmm_vs_generated_ts.py
```

### 17.3 Suggested docs

```text
docs/ts_ensemble_generation_qmmm_calibration_plan.md
docs/ts_generator_input_output_contract.md
docs/qmmm_reference_tse_protocol.md
docs/matched_gs_ts_state_pair_protocol.md
docs/ts_generator_go_no_go_gates.md
```

### 17.4 Suggested tests

```text
tests/test_ts_ensemble_manifest.py
tests/test_ts_coverage_metrics.py
tests/test_softmin_ensemble_score.py
tests/test_matched_state_pair_schema.py
tests/test_calibration_metrics.py
tests/test_no_conformer_leakage.py
```

### 17.5 Minimum implementation milestone

Milestone 1 should not require real QM/MM outputs. It should use a synthetic mini dataset:

```text
3 variants
2 TS ensemble sources: E_QMMM_synthetic and E_gen_synthetic
1 wrong-TS control
5 conformers per ensemble
fake geometry features
fake labels
```

The synthetic example should prove:

```text
coverage metrics run
productive fraction runs
variant-level split runs
baseline model runs
wrong-TS control is reported
```

---

## 18. First executable analysis plan

### 18.1 Week-zero data audit

Checklist:

```text
choose enzyme_id
choose substrate_id
choose reaction_template_id
confirm catalytic residues
confirm atom names
confirm substrate protonation state
confirm product/intermediate representation
confirm QM region draft
confirm MPNN variant list
confirm external paths for large files
```

### 18.2 First QM/MM teacher build

```text
1. Prepare WT/reference GS complex.
2. Run initial QM/MM constrained scan.
3. Extract TS-like barrier-top structures.
4. Cluster into E_QMMM.
5. Build 10–50 representative TS conformers.
6. Build matched GS for each representative.
```

### 18.3 First TS-generator calibration

```text
1. Generate 50–200 TS conformers from the TS generator.
2. Apply chemistry filters.
3. Compare E_gen to E_QMMM in reaction-coordinate space.
4. Keep E_gen clusters that pass coverage and chemistry checks.
5. Build wrong-TS and GS-like controls.
```

### 18.4 First mutant screening round

```text
1. Select 20–50 mutants across MPNN similarity bins.
2. Place E_QMMM, E_gen, and controls into every mutant.
3. Compute geometry and ensemble features.
4. Select 3–5 state pairs per mutant for QM/MM.
5. Compute first 100–250 labels.
6. Train baseline models.
7. Evaluate variant-level split and controls.
```

### 18.5 First decision

Possible decisions:

```text
Decision A:
  E_gen passes gates.
  Use E_gen for larger mutant screening.

Decision B:
  E_gen partially passes.
  Use E_gen only with strict geometry filters and active-learning uncertainty.

Decision C:
  E_gen fails.
  Use E_QMMM-derived templates for Phase 1 and improve TS generator/docking before Project 02.
```

---

## 19. Claim discipline

### 19.1 Allowed Phase 1 claim

```text
Within a fixed enzyme–substrate–reaction template,
TS ensemble placement features trained on QM/MM labels can rank mutants by computed barrier.
```

### 19.2 Allowed Phase 1B claim

```text
For this fixed system, the TS-generator ensemble preserves enough of the QM/MM-derived TS ensemble signal
to support fast mutant screening, provided geometry coverage and negative-control gates are passed.
```

### 19.3 Not allowed yet

```text
The model predicts experimental kcat for all variants.
The TS generator predicts true enzyme TS structures without QM/MM calibration.
The model generalizes to unrelated substrates.
The model generalizes to unrelated enzyme mechanisms.
Docking score alone is a barrier.
```

---

## 20. Expansion path after Phase 1

### 20.1 Same enzyme, multiple related substrates

Add:

```text
substrate descriptors
TS conformer descriptors
leaving-group descriptors
substrate-specific Δq
substrate-specific E_QMMM anchor labels
```

Required evaluation:

```text
leave-substrate-out
few-shot new-substrate adaptation
```

### 20.2 Same catalytic class, multiple enzyme families

Add:

```text
reaction-centered coordinate frames
family_id
active-site motif descriptors
structure/sequence embeddings
enzyme-family splits
```

Required evaluation:

```text
leave-enzyme-out
leave-family-out
within-catalytic-class ranking
```

### 20.3 Cross catalytic class

Do not use one undifferentiated model head. Use:

```text
shared pocket encoder
+ reaction-template encoder
+ reaction-specific heads
```

Each catalytic class needs its own TS-generation and QM/MM calibration protocol.

---

## 21. Final practical summary

The TS-generation part should be treated as a calibrated proposal system:

```text
QM/MM-derived E_QMMM = teacher reference
TS-generator-derived E_gen = fast proposal
docking/placement = pocket compatibility test
matched GS–TS QM/MM = barrier label
small model = learned surrogate
negative controls = proof that the model uses reaction-state information
```

The first-stage deliverable is a working, falsifiable pipeline for one enzyme and one substrate. If it succeeds, it becomes the foundation for a predicted-TS student model and later reaction-template-conditioned generalization.
