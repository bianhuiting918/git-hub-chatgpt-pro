# CODEX_TASKS: 同功能新骨架酶搜索 pipeline

本文件给 Codex 执行。目标是把 `README.md` 和 `MOTIF_STRATEGY.md` 中的科学流程变成可运行、可扩展、可测试的轻量 pipeline。

当前实现重点：**Sequence / Foldseek 负责召回；Folddisco 负责分层 active-site / pocket constellation 证据；Folddisco 同时支持 after-recall 验证和 de-novo 扫库；Pythia-Pocket / GRASE 做口袋、活性潜力和稳定性重排序。**

不要下载或提交大型数据库。

---

## Phase 0 — 初始化项目骨架

### Task 0.1 创建目录

在 `projects/03-new-scaffold-enzyme-search/` 下创建：

```text
configs/
data/
motifs/
src/scaffold_search/
scripts/
tests/
reports/
```

### Task 0.2 创建 Python 包

创建：

```text
src/scaffold_search/__init__.py
src/scaffold_search/schema.py
src/scaffold_search/io.py
src/scaffold_search/scoring.py
src/scaffold_search/reporting.py
src/scaffold_search/sequence_search.py
src/scaffold_search/structure_search.py
src/scaffold_search/motif_search.py
src/scaffold_search/pocket_rank.py
```

要求：

- Python 版本目标：3.10+；
- 尽量只依赖标准库、`pandas`、`pyyaml`；
- 如果要加依赖，先在项目内写明用途；
- 不添加 notebook 作为核心实现。

---

## Phase 1 — 配置文件

### Task 1.1 `configs/targets.example.yaml`

创建示例 target 配置：

```yaml
targets:
  - target_id: petase_cutinase
    function_label: PET/cutin/polyester ester bond hydrolysis
    enzyme_family_hint: alpha_beta_hydrolase
    substrate_class: polyester
    folddisco_modes: [after_recall, de_novo]
    motif_layers: [L0_catalytic_core, L1_catalytic_environment, L2_substrate_pocket, L3_processing_interface]
    novelty_goal: remote_or_new_scaffold

  - target_id: nylc_ntn_hydrolase
    function_label: nylon oligomer / polyamide amide bond hydrolysis
    enzyme_family_hint: N-terminal nucleophile hydrolase
    substrate_class: polyamide
    folddisco_modes: [after_recall, de_novo]
    motif_layers: [L0_catalytic_core, L1_catalytic_environment, L2_substrate_pocket, L3_processing_interface]
    novelty_goal: new_scaffold_or_remote_ntn
```

### Task 1.2 `configs/scoring.default.yaml`

创建默认评分权重：

```yaml
weights:
  sequence_evidence: 0.15
  fold_evidence: 0.15
  folddisco_score: 0.30
  pocket_embedding_similarity: 0.20
  predicted_stability: 0.10
  substrate_accessibility: 0.05
  novelty: 0.05

folddisco_layer_weights:
  l0_core_score: 0.45
  l1_environment_score: 0.25
  l2_pocket_score: 0.20
  l3_context_score: 0.10

folddisco_evidence_merge:
  strategy: max
  after_recall_weight: 1.0
  de_novo_weight: 1.0

penalties:
  missing_l0_core_when_confident: 0.50
  low_active_site_confidence: 0.25
  incompatible_oligomeric_state: 0.15
  bad_expression_or_solubility_proxy: 0.10

classification_thresholds:
  conservative_positive:
    min_sequence_evidence: 0.70
    min_total_score: 0.65
  remote_homolog_or_remote_fold:
    max_sequence_evidence: 0.60
    min_fold_evidence: 0.60
    min_l0_core_score: 0.60
    min_total_score: 0.60
  novel_scaffold_explorer:
    max_sequence_evidence: 0.35
    max_fold_evidence: 0.50
    min_denovo_l0_core_score: 0.65
    min_l1_environment_score: 0.45
    min_pocket_embedding_similarity: 0.60
    min_total_score: 0.55
  needs_review:
    min_total_score: 0.45
    allow_l0_partial: true
    allow_low_active_site_confidence: true
  negative_or_boundary_control:
    max_total_score: 0.45
```

### Task 1.3 `configs/databases.example.yaml`

只写 manifest，不写真实大库路径：

```yaml
databases:
  sequence:
    uniprot_reference_proteomes: /path/to/uniprot/ref_proteomes
    metagenome_proteins: /path/to/metagenome/proteins.faa
  structure:
    pdb: /path/to/pdb/mmcif
    alphafold_db_subset: /path/to/afdb/subset
    esm_atlas_subset: /path/to/esm_atlas/subset
  candidate_structures:
    recalled_candidate_manifest: /path/to/recalled_candidate_structures.tsv
  annotations:
    enzyme_annotation_table: /path/to/enzyme_annotations.tsv
```

---

## Phase 2 — 示例 seed 数据

### Task 2.1 `data/seed_enzymes.example.tsv`

列：

```text
target_id	seed_id	protein_name	organism	accession	sequence	structure_id_or_path	known_activity_note
```

至少放 4 行 toy examples。可以使用占位序列，但必须标注为 toy。

### Task 2.2 `data/catalytic_residues.example.tsv`

列：

```text
target_id	seed_id	residue_label	chain	residue_number	residue_name	role	motif_layer	required
```

示例：

- PETase/cutinase L0：`Ser_nucleophile`, `His_general_base`, `Asp_or_Glu_acid`；
- PETase/cutinase L1：`oxyanion_hole_residue`；
- PETase/cutinase L2：`aromatic_clamp`, `hydrophobic_groove_contact`；
- NylC L0：`Thr_nucleophile` 或 seed-specific catalytic-core labels；
- NylC L1：`Asp_or_Glu_network`, `amide_polarizer`；
- NylC L2：`polyamide_groove_contact`；
- NylC L3：`processing_or_interface_context`, including T227-like environmental residues where applicable。

不要把 T227-like 环境残基默认标为 L0 strict required。

---

## Phase 3 — Schema

### Task 3.1 `schema.py`

实现 dataclass 或 TypedDict：

```python
Candidate
SequenceHit
FoldseekHit
LayeredMotifHit
PocketHit
ScoredCandidate
```

`LayeredMotifHit` 必须包含：

```text
candidate_id
target_id
motif_id
motif_version
folddisco_mode                 # after_recall | de_novo
motif_source                   # validated_recall | de_novo | both | none; merge 后使用
core_motif_status              # STRICT_PASS | PARTIAL | NO_HIT | NOT_EVALUATED
l0_core_score
l1_environment_score
l2_pocket_score
l3_context_score
folddisco_score
core_rmsd
required_match_fraction
context_match_fraction
pocket_context_score
best_motif_layer               # L0 | L1 | L2 | L3
residue_mapping
missing_required_labels
partial_match_labels
active_site_confidence
hit_status                     # pass | partial | fail | low_confidence | not_run
notes
```

`ScoredCandidate` 至少包含：

```text
candidate_id
source_id
source_database
sequence_identity_to_best_seed
sequence_score
foldseek_score
foldseek_coverage
folddisco_after_recall_score
folddisco_denovo_score
folddisco_score
folddisco_mode
motif_source
core_motif_status
l0_core_score
l1_environment_score
l2_pocket_score
l3_context_score
best_motif_layer
motif_id
motif_residue_mapping
pocket_score
pocket_id
stability_score
substrate_accessibility_score
active_site_confidence
novelty_score
penalties
total_score
candidate_class
recommended_panel
explanation
```

---

## Phase 4 — Layered motif YAML 模板

### Task 4.1 `motifs/petase_cutinase.layered.example.yaml`

创建 PETase/cutinase layered motif 示例：

```yaml
motif_id: petase_cutinase_layered_active_site_pocket
version: toy_v1
target_id: petase_cutinase
intended_modes: [after_recall, de_novo]

layers:
  L0_catalytic_core:
    hard_filter: true
    weight: 0.45
    residues:
      - label: Ser_nucleophile
        allowed_residue_names: [SER]
        role: nucleophile
        required: true
      - label: His_general_base
        allowed_residue_names: [HIS]
        role: general_base
        required: true
      - label: Asp_or_Glu_acid
        allowed_residue_names: [ASP, GLU]
        role: acid
        required: true
    geometry_constraints:
      - pair: [Ser_nucleophile, His_general_base]
        distance_angstrom: [2.5, 4.2]
      - pair: [His_general_base, Asp_or_Glu_acid]
        distance_angstrom: [2.5, 4.5]

  L1_catalytic_environment:
    hard_filter: false
    weight: 0.25
    residues:
      - label: Oxyanion_hole_1
        allowed_residue_names: [SER, THR, GLY, ASN, GLN, TYR, BACKBONE_NH]
        role: oxyanion_stabilization
        required: false

  L2_substrate_pocket:
    hard_filter: false
    weight: 0.20
    residues:
      - label: Aromatic_clamp_1
        allowed_residue_names: [TRP, TYR, PHE]
        role: polymer_binding
        required: false
      - label: Hydrophobic_groove_1
        allowed_residue_names: [LEU, ILE, VAL, MET, PHE, TYR, TRP]
        role: hydrophobic_groove
        required: false

  L3_processing_interface:
    hard_filter: false
    weight: 0.10
    residues:
      - label: Lid_or_surface_context_1
        allowed_residue_names: [GLY, ALA, SER, THR, LEU, ILE, VAL, TYR, TRP, PHE]
        role: dynamics_or_surface_binding
        required: false
```

### Task 4.2 `motifs/nylc_ntn_hydrolase.layered.example.yaml`

创建 NylC layered motif 示例。不要硬编码真实编号；用 seed annotation 映射。

```yaml
motif_id: nylc_layered_polyamide_pocket_constellation
version: toy_v1
target_id: nylc_ntn_hydrolase
intended_modes: [after_recall, de_novo]

layers:
  L0_catalytic_core:
    hard_filter: true
    weight: 0.45
    residues:
      - label: Thr_nucleophile_or_seed_core_1
        allowed_residue_names: [THR, SER, CYS]
        role: nucleophile_or_seed_core
        required: true
      - label: Core_acid_base_or_seed_core_2
        allowed_residue_names: [ASP, GLU, HIS, SER, THR]
        role: acid_base_network_or_seed_core
        required: true
      - label: Core_acid_base_or_seed_core_3
        allowed_residue_names: [ASP, GLU, HIS, SER, THR]
        role: acid_base_network_or_seed_core
        required: false
    geometry_constraints:
      - pair: [Thr_nucleophile_or_seed_core_1, Core_acid_base_or_seed_core_2]
        distance_angstrom: [2.5, 5.5]

  L1_catalytic_environment:
    hard_filter: false
    weight: 0.25
    residues:
      - label: Amide_polarizer_1
        allowed_residue_names: [ASN, GLN, SER, THR, TYR, HIS, BACKBONE_NH]
        role: amide_polarization
        required: false
      - label: Acidic_network_context_1
        allowed_residue_names: [ASP, GLU]
        role: acidic_network_context
        required: false

  L2_substrate_pocket:
    hard_filter: false
    weight: 0.20
    residues:
      - label: Polyamide_groove_contact_1
        allowed_residue_names: [LEU, ILE, VAL, MET, PHE, TYR, TRP, ASN, GLN, SER, THR]
        role: substrate_groove
        required: false
      - label: Polyamide_groove_contact_2
        allowed_residue_names: [LEU, ILE, VAL, MET, PHE, TYR, TRP, ASN, GLN, SER, THR]
        role: substrate_groove
        required: false

  L3_processing_interface:
    hard_filter: false
    weight: 0.10
    residues:
      - label: Processing_or_interface_context_1
        allowed_residue_names: [GLY, ALA, SER, THR, ASP, GLU, ASN, GLN]
        role: processing_or_interface
        required: false
      - label: T227_like_environmental_context
        allowed_residue_names: [THR, SER, ALA, GLY, ASP, GLU, ASN, GLN]
        role: environmental_context_not_strict_core
        required: false
```

---

## Phase 5 — 外部工具 wrappers

所有 wrapper 必须支持两种模式：

1. 真实运行外部工具；
2. 读取已有 TSV/JSON 结果并继续 pipeline。

如果二进制不可用，必须抛出明确错误：

```text
ExternalToolUnavailable: foldseek executable not found. Provide --precomputed-hits to continue.
```

### Task 5.1 `sequence_search.py`

实现接口：

```python
def run_mmseqs_search(query_fasta: str, db_path: str, out_tsv: str, *, min_identity: float | None = None) -> str: ...
def run_hmmer_search(hmm_path: str, db_fasta: str, out_tsv: str) -> str: ...
def load_sequence_hits(path: str) -> list[SequenceHit]: ...
def score_sequence_hit(hit: SequenceHit) -> float: ...
```

### Task 5.2 `structure_search.py`

实现接口：

```python
def run_foldseek_search(query_structure: str, structure_db: str, out_tsv: str, *, mode: str = "global") -> str: ...
def load_foldseek_hits(path: str) -> list[FoldseekHit]: ...
def score_foldseek_hit(hit: FoldseekHit) -> float: ...
```

必须支持 `mode="global"` 和 `mode="local_active_domain"`。

### Task 5.3 `motif_search.py`: 双模式 + 分层 Folddisco

实现以下接口：

```python
def load_layered_motif_yaml(path: str) -> dict: ...

def run_folddisco_after_recall(
    motif_yaml: str,
    candidate_structure_manifest: str,
    out_tsv: str,
) -> str: ...

def run_folddisco_denovo(
    motif_yaml: str,
    structure_db: str,
    out_tsv: str,
) -> str: ...

def load_layered_motif_hits(path: str, *, folddisco_mode: str | None = None) -> list[LayeredMotifHit]: ...

def score_layered_motif_hit(hit: LayeredMotifHit, layer_weights: dict) -> float: ...

def merge_folddisco_evidence(
    after_recall_hits: list[LayeredMotifHit],
    denovo_hits: list[LayeredMotifHit],
) -> dict[str, LayeredMotifHit]: ...
```

`run_folddisco_after_recall` 的输入是 sequence/Foldseek 已召回候选的结构 manifest，不是全库。

`run_folddisco_denovo` 的输入是 broad structure database，用于发现 sequence/Foldseek 漏掉的新骨架。

评分公式：

```text
folddisco_score =
  0.45 × l0_core_score
+ 0.25 × l1_environment_score
+ 0.20 × l2_pocket_score
+ 0.10 × l3_context_score
```

合并规则：

```text
folddisco_after_recall_score = best score among after_recall hits
folddisco_denovo_score = best score among de_novo hits
folddisco_score = max(folddisco_after_recall_score, folddisco_denovo_score)

if both scores exist:
    motif_source = both
elif after_recall score exists:
    motif_source = validated_recall
elif de_novo score exists:
    motif_source = de_novo
else:
    motif_source = none
```

### Task 5.4 `pocket_rank.py`

实现接口：

```python
def run_pythia_pocket(structure_path: str, out_json: str) -> str: ...
def run_grase_rank(candidate_table: str, out_tsv: str) -> str: ...
def load_pocket_hits(path: str) -> list[PocketHit]: ...
def score_pocket_hit(hit: PocketHit) -> float: ...
```

必须允许没有 Pythia/GRASE 时读取预计算结果。

Pocket hit 字段建议：

```text
candidate_id
pocket_id
pocket_center
pocket_residues
pocket_embedding_similarity
active_site_overlap_fraction
substrate_accessibility_score
predicted_stability
model_name
model_version
```

---

## Phase 6 — 候选合并和去重

### Task 6.1 `io.py`

实现：

```python
def read_tsv(path: str): ...
def write_tsv(rows, path: str): ...
def read_yaml(path: str): ...
def write_json(obj, path: str): ...
```

### Task 6.2 candidate merge

实现：

```python
def merge_evidence(
    sequence_hits,
    foldseek_hits,
    folddisco_after_recall_hits,
    folddisco_denovo_hits,
    pocket_hits,
) -> list[Candidate]: ...
```

去重逻辑：

1. 优先按 accession / UniProt ID / source ID；
2. 其次按 sequence hash；
3. 结构-only 候选按 structure ID + chain + domain range；
4. 同一候选来自多个来源时合并证据，不覆盖证据；
5. Folddisco-de-novo 独有候选必须进入 `C_total`，即使没有 sequence/Foldseek 证据；
6. Folddisco-after-recall 只作为 recalled candidate 的机制证据，不应创建与原 recalled candidate 重复的新 candidate。

### Task 6.3 novelty score

实现：

```python
def compute_novelty_score(sequence_score: float, fold_score: float, known_family_flag: bool) -> float:
    ...
```

建议逻辑：

```text
base = 1 - max(sequence_score, 0.7 * fold_score)
if known_family_flag:
    base -= 0.15
return clamp(base, 0, 1)
```

---

## Phase 7 — 可解释评分和分类

### Task 7.1 `scoring.py`

实现：

```python
def clamp(x: float, low: float = 0.0, high: float = 1.0) -> float: ...
def combine_folddisco_scores(after_recall_score: float | None, denovo_score: float | None, strategy: str = "max") -> tuple[float, str]: ...
def weighted_score(candidate: Candidate, weights: dict, penalties: dict) -> float: ...
def classify_candidate(candidate: Candidate, thresholds: dict) -> str: ...
def build_explanation(candidate: Candidate) -> str: ...
def should_send_to_pocket_ranker(candidate: Candidate) -> bool: ...
```

评分必须满足：

- 所有输入分数先 clamp 到 `[0, 1]`；
- penalty 不允许让总分小于 0；
- 缺失证据不能自动当 0，除非该证据为该类别必需；
- explanation 必须说明进入类别的原因；
- explanation 必须指出 Folddisco 证据来自 `after_recall`、`de_novo`、`both` 或 `none`；
- explanation 必须展示 L0/L1/L2/L3 简要证据。

### Task 7.2 classification 规则

实现五类：

```text
conservative_positive
remote_homolog_or_remote_fold
novel_scaffold_explorer
needs_review
negative_or_boundary_control
```

默认逻辑：

1. `L0 NO_HIT` 且 `active_site_confidence` 高 → 优先 `negative_or_boundary_control`；
2. sequence high + total high → `conservative_positive`，但若 Folddisco-after-recall 明确 L0 冲突，应降为 `needs_review` 或 `negative_or_boundary_control`；
3. sequence low/moderate + fold high + L0/L1 支持 + total high → `remote_homolog_or_remote_fold`；
4. sequence low + fold low + Folddisco-de-novo L0 high + L1/L2 支持 + pocket high → `novel_scaffold_explorer`；
5. L0 partial 但 L1/L2/pocket 强，或 active-site confidence 低 → `needs_review`；
6. total low 或机制明显不合理 → `negative_or_boundary_control`。

不要实现为：

```text
full motif fail → delete
```

### Task 7.3 Pythia/GRASE routing

`should_send_to_pocket_ranker` 应返回 true 的候选包括：

1. L0 strict pass；
2. Folddisco-de-novo 直接命中；
3. L0 partial 但 L1/L2 强；
4. sequence / Foldseek 强召回但 Folddisco 不完整的边界候选；
5. 少量 negative / boundary controls。

---

## Phase 8 — CLI scripts

创建以下脚本。所有脚本都必须有 `--help`。

### Task 8.1 `scripts/run_sequence_search.py`

参数：

```text
--query-fasta
--db
--out
--tool mmseqs|hmmer|precomputed
--precomputed-hits
```

### Task 8.2 `scripts/run_foldseek_search.py`

参数：

```text
--query-structure
--structure-db
--out
--mode global|local_active_domain|precomputed
--precomputed-hits
```

### Task 8.3 `scripts/run_folddisco_after_recall.py`

参数：

```text
--motif-yaml
--candidate-structure-manifest
--out
--precomputed-hits
```

输出必须包含：

```text
folddisco_mode=after_recall
l0_core_score
l1_environment_score
l2_pocket_score
l3_context_score
core_motif_status
```

### Task 8.4 `scripts/run_folddisco_denovo.py`

参数：

```text
--motif-yaml
--structure-db
--out
--precomputed-hits
```

输出必须包含：

```text
folddisco_mode=de_novo
l0_core_score
l1_environment_score
l2_pocket_score
l3_context_score
core_motif_status
```

### Task 8.5 `scripts/run_pocket_ranking.py`

参数：

```text
--candidate-table
--structures-dir
--out
--tool pythia-pocket|grase|precomputed
--precomputed-hits
```

### Task 8.6 `scripts/merge_and_rank_candidates.py`

参数：

```text
--sequence-hits
--foldseek-hits
--folddisco-after-recall-hits
--folddisco-denovo-hits
--pocket-hits
--scoring-config
--out-ranking
--out-report-md
```

输出：

```text
reports/candidate_ranking.tsv
reports/candidate_report.md
```

---

## Phase 9 — 报告生成

### Task 9.1 `reporting.py`

实现：

```python
def render_candidate_table(scored_candidates) -> str: ...
def render_top_candidates(scored_candidates, top_n: int = 50) -> str: ...
def render_layered_motif_summary(scored_candidates) -> str: ...
def render_pocket_ranker_routing(scored_candidates) -> str: ...
def render_panel_design(scored_candidates) -> str: ...
def write_markdown_report(scored_candidates, path: str) -> None: ...
```

报告必须包括：

1. 输入文件摘要；
2. `C_sequence` 数量；
3. `C_foldseek` 数量；
4. Folddisco-after-recall pass / partial / fail 数量；
5. Folddisco-de-novo 新候选数量；
6. L0/L1/L2/L3 分数分布摘要；
7. 进入 Pythia/GRASE 的候选数量和理由；
8. 每类候选数量；
9. Top candidates；
10. novel scaffold explorer 列表；
11. needs_review 列表；
12. negative/boundary controls；
13. 每个 top 候选 explanation；
14. 下一步实验建议。

### Task 9.2 `reports/README.md`

说明报告目录只保存轻量 summary，不保存大型原始结果。

---

## Phase 10 — 测试

创建 `tests/`：

```text
tests/test_scoring.py
tests/test_classification.py
tests/test_io.py
tests/test_merge.py
tests/test_folddisco_modes.py
tests/test_layered_motif_schema.py
```

测试必须覆盖：

- score clamp；
- penalty；
- layered motif YAML parsing；
- L0/L1/L2/L3 score aggregation；
- conservative_positive 分类；
- remote_homolog_or_remote_fold 分类；
- novel_scaffold_explorer 分类；
- needs_review 分类；
- L0 no-hit + high confidence 进入 negative/boundary；
- L0 partial + L1/L2 strong 进入 needs_review 或 exploration；
- 多来源证据合并；
- Folddisco-after-recall 不创建重复 candidate；
- Folddisco-de-novo 独有候选进入 `C_total`；
- `motif_source=both` 的合并逻辑；
- Pythia/GRASE routing；
- explanation 非空且提到 Folddisco mode 和 L0/L1/L2/L3 evidence。

---

## Phase 11 — Toy end-to-end demo

添加一个最小 demo：

```text
data/toy_sequence_hits.tsv
data/toy_foldseek_hits.tsv
data/toy_folddisco_after_recall_hits.tsv
data/toy_folddisco_denovo_hits.tsv
data/toy_pocket_hits.tsv
```

要求包含至少 10 个候选：

1. 高同源 PETase-like，L0/L1 通过；
2. 高同源但缺 catalytic Ser 的 PETase-like 负控；
3. 低序列但 Foldseek 高的 cutinase-like，L0/L1 通过；
4. Foldseek 高但 L0 no-hit 且 active-site confidence 高的负控；
5. Foldseek 高但 L0 partial、L1/L2 强的 needs_review；
6. Folddisco-de-novo L0/L1/L2 高、sequence/fold 低的新骨架候选；
7. Folddisco-de-novo L0 高但 L2 pocket 不合理的边界候选；
8. NylC-like Ntn hydrolase 远缘候选，L0/L1 通过；
9. NylC motif partial hit 且 active-site confidence 低的 needs_review；
10. T227-like L3 mismatch 但 L0/L1/L2 合理的候选，验证不会被误判为 L0 fail。

运行命令示例写入 `README.md` 或 `reports/README.md`：

```bash
python scripts/merge_and_rank_candidates.py \
  --sequence-hits data/toy_sequence_hits.tsv \
  --foldseek-hits data/toy_foldseek_hits.tsv \
  --folddisco-after-recall-hits data/toy_folddisco_after_recall_hits.tsv \
  --folddisco-denovo-hits data/toy_folddisco_denovo_hits.tsv \
  --pocket-hits data/toy_pocket_hits.tsv \
  --scoring-config configs/scoring.default.yaml \
  --out-ranking reports/candidate_ranking.tsv \
  --out-report-md reports/candidate_report.md
```

---

## Phase 12 — 后续真实数据接入

第一轮完成后再接入真实数据。

真实接入顺序建议：

1. 准备 seed enzymes 和 catalytic residue annotation；
2. 跑 MMseqs2/HMMER，生成 sequence hits；
3. 准备 seed structures 和 candidate structures；
4. 跑 Foldseek global/local；
5. 重新定义 layered Folddisco motif ensemble；
6. 将 `C_sequence ∪ C_foldseek` 结构化为 recalled candidate structure manifest；
7. 跑 Folddisco-after-recall，对 recalled candidates 做 L0/L1/L2/L3 机制筛选；
8. 跑 Folddisco-de-novo，直接扫 PDB / AFDB subset / ESMAtlas subset / metagenome structure set；
9. 将已有 sequence / Foldseek / Folddisco 结果映射到新 schema；
10. 选择进入 Pythia-Pocket/GRASE 的候选；
11. 跑 Pythia-Pocket/GRASE 或读取预计算结果；
12. 合并排序；
13. 手动审查 top 50；
14. 输出实验板。

---

## Definition of Done

第一轮完成标准：

- [ ] 项目目录完整；
- [ ] 示例配置和 toy 数据存在；
- [ ] layered motif YAML 示例存在；
- [ ] schema 支持 L0/L1/L2/L3；
- [ ] wrapper 接口存在且有清晰错误；
- [ ] `run_folddisco_after_recall.py` 存在；
- [ ] `run_folddisco_denovo.py` 存在；
- [ ] candidate merge 可运行；
- [ ] Folddisco-after-recall 和 Folddisco-de-novo 证据能合并；
- [ ] Pythia/GRASE routing 可解释；
- [ ] scoring 可运行；
- [ ] classification 可解释；
- [ ] toy demo 可生成 `reports/candidate_ranking.tsv` 和 `reports/candidate_report.md`；
- [ ] tests 通过；
- [ ] 没有提交大型数据库、checkpoint、trajectory、docking ensemble 或临时下载。

---

## 注意事项

- Sequence search 和 Foldseek 的结果必须进入 Folddisco-after-recall 做机制筛选。
- Folddisco 不应只作为后置过滤器；还必须有 de-novo 扫库分支。
- Folddisco motif 必须分 L0/L1/L2/L3，不能只定义一个 strict full motif。
- L0 是 catalytic core；L2/L3 是口袋和上下文，不要混用。
- T227-like 环境残基不能默认作为 L0 strict required。
- Pythia-Pocket/GRASE 是 ranker，不是唯一真值。
- NylC 要按 Ntn hydrolase 和 polyamide groove 处理，不要简单套 Ser-His-Asp esterase 模式。
- Folddisco fail 不等于直接删除；需要结合 L0/L1/L2/L3、active-site confidence、pocket/GRASE 分数和人工审查策略。
- 每个 top candidate 必须保留证据链：sequence、fold、Folddisco mode、L0/L1/L2/L3、pocket、stability、novelty、penalty、explanation。
