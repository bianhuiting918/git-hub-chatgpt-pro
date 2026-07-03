# CODEX_TASKS: 同功能新骨架酶搜索 pipeline

本文件给 Codex 执行。目标是把 `README.md` 中的科学流程变成可运行、可扩展、可测试的轻量 pipeline。

第一轮实现重点：**数据结构、wrapper 接口、候选合并、可解释评分、报告生成**。不要下载或提交大型数据库。

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

## Phase 1 — 定义输入配置和 schema

### Task 1.1 `configs/targets.example.yaml`

创建示例 target 配置：

```yaml
targets:
  - target_id: petase_cutinase
    function_label: PET/cutin/polyester ester bond hydrolysis
    enzyme_family_hint: alpha_beta_hydrolase
    substrate_class: polyester
    required_motif_groups:
      - ser_his_asp_or_glu
      - oxyanion_hole
      - polymer_binding_groove
    novelty_goal: remote_or_new_scaffold

  - target_id: nylc_ntn_hydrolase
    function_label: nylon oligomer / polyamide amide bond hydrolysis
    enzyme_family_hint: N-terminal nucleophile hydrolase
    substrate_class: polyamide
    required_motif_groups:
      - n_terminal_thr_nucleophile
      - acidic_residue_network
      - amide_binding_groove
    novelty_goal: new_scaffold_or_remote_ntn
```

### Task 1.2 `configs/scoring.default.yaml`

创建默认评分权重：

```yaml
weights:
  sequence_evidence: 0.15
  fold_evidence: 0.20
  motif_geometry: 0.25
  pocket_embedding_similarity: 0.20
  predicted_stability: 0.10
  substrate_accessibility: 0.05
  novelty: 0.05

penalties:
  missing_required_catalytic_residue: 0.50
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
    min_total_score: 0.60
  novel_scaffold_explorer:
    max_sequence_evidence: 0.35
    max_fold_evidence: 0.50
    min_motif_geometry: 0.65
    min_pocket_embedding_similarity: 0.60
    min_total_score: 0.55
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
  annotations:
    enzyme_annotation_table: /path/to/enzyme_annotations.tsv
```

### Task 1.4 示例 seed 数据

创建 `data/seed_enzymes.example.tsv`，列：

```text
target_id	seed_id	protein_name	organism	accession	sequence	structure_id_or_path	known_activity_note
```

至少放 4 行 toy examples。可以使用占位序列，但必须标注为 toy。

创建 `data/catalytic_residues.example.tsv`，列：

```text
target_id	seed_id	residue_label	chain	residue_number	residue_name	role	required
```

示例：

- PETase/cutinase：`Ser_nucleophile`, `His_general_base`, `Asp_or_Glu_acid`, `oxyanion_hole_residue`；
- NylC：`Thr_nucleophile`, `Asp_network_1`, `Asp_network_2`, `processing_site_context`。

### Task 1.5 `schema.py`

实现 dataclass 或 TypedDict：

```python
Candidate
SequenceHit
FoldseekHit
MotifHit
PocketHit
ScoredCandidate
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
motif_score
motif_name
motif_residue_mapping
pocket_score
pocket_id
stability_score
substrate_accessibility_score
novelty_score
penalties
total_score
candidate_class
recommended_panel
explanation
```

---

## Phase 2 — 外部工具 wrappers

所有 wrapper 必须支持两种模式：

1. 真实运行外部工具；
2. 读取已有 TSV/JSON 结果并继续 pipeline。

如果二进制不可用，必须抛出明确错误：

```text
ExternalToolUnavailable: foldseek executable not found. Provide --precomputed-hits to continue.
```

### Task 2.1 `sequence_search.py`

实现接口：

```python
def run_mmseqs_search(query_fasta: str, db_path: str, out_tsv: str, *, min_identity: float | None = None) -> str: ...
def run_hmmer_search(hmm_path: str, db_fasta: str, out_tsv: str) -> str: ...
def load_sequence_hits(path: str) -> list[SequenceHit]: ...
def score_sequence_hit(hit: SequenceHit) -> float: ...
```

评分建议：

```text
sequence_score = max(
  normalized_identity,
  normalized_hmm_bitscore,
  conserved_required_residue_fraction
)
```

### Task 2.2 `structure_search.py`

实现接口：

```python
def run_foldseek_search(query_structure: str, structure_db: str, out_tsv: str, *, mode: str = "global") -> str: ...
def load_foldseek_hits(path: str) -> list[FoldseekHit]: ...
def score_foldseek_hit(hit: FoldseekHit) -> float: ...
```

必须支持 `mode="global"` 和 `mode="local_active_domain"`。

评分建议：

```text
fold_evidence = weighted_mean(
  alignment_score_norm,
  coverage_norm,
  active_site_alignment_bonus,
  active_site_confidence
)
```

### Task 2.3 `motif_search.py`

实现接口：

```python
def run_folddisco_search(motif_yaml: str, structure_db: str, out_tsv: str) -> str: ...
def load_motif_hits(path: str) -> list[MotifHit]: ...
def score_motif_hit(hit: MotifHit) -> float: ...
```

Motif hit 必须保存：

```text
motif_name
candidate_id
residue_mapping
rmsd_or_geometry_error
required_residue_match_fraction
optional_context_match_fraction
active_site_confidence
```

评分建议：

```text
motif_geometry = 0.45 * required_residue_match_fraction \
               + 0.25 * geometry_score \
               + 0.20 * optional_context_match_fraction \
               + 0.10 * active_site_confidence
```

### Task 2.4 `pocket_rank.py`

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

## Phase 3 — Motif 模板

### Task 3.1 `motifs/petase_cutinase.example.yaml`

示例内容：

```yaml
motif_id: petase_cutinase_ser_his_asp_oxyanion_aromatic_clamp
version: toy_v1
description: Toy PETase/cutinase active-site constellation. Coordinates are placeholders.
required_residues:
  - label: Ser_nucleophile
    allowed_residue_names: [SER]
    role: nucleophile
  - label: His_general_base
    allowed_residue_names: [HIS]
    role: general_base
  - label: Asp_or_Glu_acid
    allowed_residue_names: [ASP, GLU]
    role: acid
  - label: Oxyanion_hole_1
    allowed_residue_names: [SER, THR, GLY, ASN, GLN, TYR, BACKBONE_NH]
    role: oxyanion_stabilization
optional_context_residues:
  - label: aromatic_clamp_1
    allowed_residue_names: [TRP, TYR, PHE]
  - label: hydrophobic_groove_1
    allowed_residue_names: [LEU, ILE, VAL, MET, PHE, TYR, TRP]
geometry_constraints:
  - pair: [Ser_nucleophile, His_general_base]
    distance_angstrom: [2.5, 4.2]
  - pair: [His_general_base, Asp_or_Glu_acid]
    distance_angstrom: [2.5, 4.5]
```

### Task 3.2 `motifs/nylc_ntn_hydrolase.example.yaml`

示例内容：

```yaml
motif_id: nylc_ntn_thr_asp_asp_amide_groove
version: toy_v1
description: Toy NylC/Ntn hydrolase motif. Coordinates are placeholders.
required_residues:
  - label: Thr_nucleophile
    allowed_residue_names: [THR]
    role: n_terminal_nucleophile
  - label: Acidic_network_1
    allowed_residue_names: [ASP, GLU]
    role: acid_or_base_network
  - label: Acidic_network_2
    allowed_residue_names: [ASP, GLU]
    role: acid_or_base_network
  - label: Amide_binding_polar_1
    allowed_residue_names: [ASN, GLN, SER, THR, TYR, HIS, BACKBONE_NH]
    role: amide_binding
optional_context_residues:
  - label: processing_site_context
    allowed_residue_names: [GLY, ALA, SER, THR, ASP, GLU]
  - label: oligomer_groove_hydrophobic_1
    allowed_residue_names: [LEU, ILE, VAL, MET, PHE, TYR, TRP]
geometry_constraints:
  - pair: [Thr_nucleophile, Acidic_network_1]
    distance_angstrom: [2.5, 5.0]
  - pair: [Thr_nucleophile, Amide_binding_polar_1]
    distance_angstrom: [3.0, 7.0]
```

---

## Phase 4 — 候选合并和去重

### Task 4.1 `io.py`

实现：

```python
def read_tsv(path: str): ...
def write_tsv(rows, path: str): ...
def read_yaml(path: str): ...
def write_json(obj, path: str): ...
```

### Task 4.2 candidate merge

在 `scoring.py` 或单独文件中实现：

```python
def merge_evidence(sequence_hits, foldseek_hits, motif_hits, pocket_hits) -> list[Candidate]: ...
```

去重逻辑：

1. 优先按 accession / UniProt ID / source ID；
2. 其次按 sequence hash；
3. 结构-only 候选按 structure ID + chain + domain range；
4. 同一候选来自多个来源时合并证据，不覆盖证据。

### Task 4.3 novelty score

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

## Phase 5 — 可解释评分和分类

### Task 5.1 `scoring.py`

实现：

```python
def clamp(x: float, low: float = 0.0, high: float = 1.0) -> float: ...
def weighted_score(candidate: Candidate, weights: dict, penalties: dict) -> float: ...
def classify_candidate(candidate: Candidate, thresholds: dict) -> str: ...
def build_explanation(candidate: Candidate) -> str: ...
```

评分必须满足：

- 所有输入分数先 clamp 到 `[0, 1]`；
- penalty 不允许让总分小于 0；
- 缺失证据不能自动当 0，除非该证据为该类别必需；
- explanation 必须说明进入类别的原因。

### Task 5.2 classification 规则

实现四类：

```text
conservative_positive
remote_homolog_or_remote_fold
novel_scaffold_explorer
negative_or_boundary_control
```

默认逻辑：

1. 缺失 required catalytic residue → 优先 negative_or_boundary_control；
2. sequence high + total high → conservative_positive；
3. sequence low/moderate + fold high + total high → remote_homolog_or_remote_fold；
4. sequence low + fold low + motif high + pocket high → novel_scaffold_explorer；
5. 其他 → boundary_control 或 needs_review。

可以额外加入 `needs_review`，但最终报告必须能映射到实验板类别。

---

## Phase 6 — CLI scripts

创建以下脚本。所有脚本都必须有 `--help`。

### Task 6.1 `scripts/run_sequence_search.py`

参数：

```text
--query-fasta
--db
--out
--tool mmseqs|hmmer|precomputed
--precomputed-hits
```

### Task 6.2 `scripts/run_foldseek_search.py`

参数：

```text
--query-structure
--structure-db
--out
--mode global|local_active_domain|precomputed
--precomputed-hits
```

### Task 6.3 `scripts/run_folddisco_search.py`

参数：

```text
--motif-yaml
--structure-db
--out
--precomputed-hits
```

### Task 6.4 `scripts/run_pocket_ranking.py`

参数：

```text
--candidate-table
--structures-dir
--out
--tool pythia-pocket|grase|precomputed
--precomputed-hits
```

### Task 6.5 `scripts/merge_and_rank_candidates.py`

参数：

```text
--sequence-hits
--foldseek-hits
--motif-hits
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

## Phase 7 — 报告生成

### Task 7.1 `reporting.py`

实现：

```python
def render_candidate_table(scored_candidates) -> str: ...
def render_top_candidates(scored_candidates, top_n: int = 50) -> str: ...
def render_panel_design(scored_candidates) -> str: ...
def write_markdown_report(scored_candidates, path: str) -> None: ...
```

报告必须包括：

1. 输入文件摘要；
2. 每类候选数量；
3. Top candidates；
4. novel scaffold explorer 列表；
5. negative/boundary controls；
6. 每个 top 候选 explanation；
7. 下一步实验建议。

### Task 7.2 `reports/README.md`

说明报告目录只保存轻量 summary，不保存大型原始结果。

---

## Phase 8 — 测试

创建 `tests/`：

```text
tests/test_scoring.py
tests/test_classification.py
tests/test_io.py
tests/test_merge.py
```

测试必须覆盖：

- score clamp；
- penalty；
- conservative_positive 分类；
- remote_homolog_or_remote_fold 分类；
- novel_scaffold_explorer 分类；
- missing required catalytic residue 进入 negative/boundary；
- 多来源证据合并；
- explanation 非空。

---

## Phase 9 — Toy end-to-end demo

添加一个最小 demo：

```text
data/toy_sequence_hits.tsv
data/toy_foldseek_hits.tsv
data/toy_motif_hits.tsv
data/toy_pocket_hits.tsv
```

要求包含至少 6 个候选：

1. 高同源 PETase-like；
2. 低序列但 Foldseek 高的 cutinase-like；
3. Folddisco + pocket 高但 sequence/fold 低的新骨架候选；
4. 缺 catalytic Ser 的负控；
5. NylC-like Ntn hydrolase 远缘候选；
6. NylC motif 命中但 pocket 不合理的边界候选。

运行命令示例写入 `README.md` 或 `reports/README.md`：

```bash
python scripts/merge_and_rank_candidates.py \
  --sequence-hits data/toy_sequence_hits.tsv \
  --foldseek-hits data/toy_foldseek_hits.tsv \
  --motif-hits data/toy_motif_hits.tsv \
  --pocket-hits data/toy_pocket_hits.tsv \
  --scoring-config configs/scoring.default.yaml \
  --out-ranking reports/candidate_ranking.tsv \
  --out-report-md reports/candidate_report.md
```

---

## Phase 10 — 后续真实数据接入

第一轮完成后再接入真实数据。

真实接入顺序建议：

1. 准备 seed enzymes 和 catalytic residue annotation；
2. 跑 MMseqs2/HMMER，生成 sequence hits；
3. 准备 seed structures 和 candidate structures；
4. 跑 Foldseek global/local；
5. 从已知结构和复合物构建 Folddisco motif ensemble；
6. 跑 Folddisco；
7. 跑 Pythia-Pocket/GRASE 或读取预计算结果；
8. 合并排序；
9. 手动审查 top 50；
10. 输出实验板。

---

## Definition of Done

第一轮完成标准：

- [ ] 项目目录完整；
- [ ] 示例配置和 toy 数据存在；
- [ ] wrapper 接口存在且有清晰错误；
- [ ] candidate merge 可运行；
- [ ] scoring 可运行；
- [ ] classification 可解释；
- [ ] toy demo 可生成 `reports/candidate_ranking.tsv` 和 `reports/candidate_report.md`；
- [ ] tests 通过；
- [ ] 没有提交大型数据库、checkpoint、trajectory、docking ensemble 或临时下载。

---

## 注意事项

- 新骨架搜索宁可并联召回后再重排序，不要早期用单一阈值裁剪。
- Folddisco motif 必须包含催化残基之外的口袋/底物定位信息，否则假阳性会很多。
- Pythia-Pocket/GRASE 是 ranker，不是唯一真值。
- NylC 要按 Ntn hydrolase 和 polyamide groove 处理，不要简单套 Ser-His-Asp esterase 模式。
- 每个 top candidate 必须保留证据链：sequence、fold、motif、pocket、stability、novelty、penalty、explanation。
