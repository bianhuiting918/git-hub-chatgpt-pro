# Project 03: 同功能新骨架酶搜索计划

本项目面向 **在现有功能酶基础上搜索具有相同功能的新骨架酶**。典型目标包括：

- 塑料降解相关酶：PETase、cutinase、lipase/esterase、polyurethane urethanase 等；
- NylC / nylon hydrolase / 6-aminohexanoate oligomer hydrolase 等 polyamide 降解酶；
- 未来可扩展到其他明确催化机制、底物和活性位点的酶类。

核心思想：**序列搜索、Foldseek、Folddisco 并联召回；GRASE / Pythia-Pocket 统一重排序；实验结果回流到 motif、HMM 和 ranker。**

不要把流程做成严格的单向漏斗。真正的新骨架候选可能在序列搜索或整体 Foldseek 中分数很低，如果早期阈值过严会被提前丢弃。

---

## 1. 目标定义

### 1.1 项目目标

给定一组已知功能酶，构建一个可复用 pipeline：

1. 从已知酶出发扩展序列同源空间；
2. 用整体结构相似性搜索低序列相似候选；
3. 用局部三维催化 motif 搜索可能的收敛进化或新骨架候选；
4. 用口袋级模型和稳定性/活性推荐模型重排序；
5. 输出实验优先级候选板，包括保守阳性、远缘候选、新骨架探索候选和负控。

### 1.2 非目标

本项目第一阶段不直接保证：

- 精确预测 `kcat` 或 `kcat/KM`；
- 替代 docking、MD、QM/MM 或湿实验；
- 自动提交大规模结构库、模型 checkpoint、原始轨迹或大型 docking 输出。

第一阶段的产物是 **可解释、可追踪、可实验验证的候选排序系统**。

---

## 2. 工具分工

| 层级 | 工具 | 主要职责 | 输出证据 |
|---|---|---|---|
| 序列同源 | BLAST / MMseqs2 / HMMER / jackhmmer | 找近缘和远缘同源，建立 MSA、HMM、保守残基图谱 | `sequence_score`、保守催化残基、家族分支、同源正负例 |
| 整体结构相似 | Foldseek | 搜索低序列相似但整体 fold 或 domain 相似的候选 | `fold_score`、coverage、TM-score-like 指标、domain alignment |
| 局部三维 motif | Folddisco | 搜索不连续催化基序、金属配位、Ntn hydrolase 核心几何等 | `motif_geometry_score`、残基映射、几何偏差 |
| 口袋/活性/稳定性 | GRASE / Pythia-Pocket / Pythia | 识别更可能具有目标底物口袋、可表达稳定且活性的候选 | `pocket_embedding_similarity`、`active_site_confidence`、`stability_score` |
| 物理合理性检查 | docking / PLACER / MD / 几何规则 | 检查底物 pose、亲核攻击距离、氧阴离子洞、入口开放性 | `substrate_pose_score`、`mechanistic_sanity_score` |

---

## 3. 候选召回策略

候选来源必须并联：

```text
C_total = C_sequence ∪ C_foldseek ∪ C_folddisco
```

之后统一去重、聚类、打分和排序。

### 3.1 序列搜索线

用途：

- 构建已知功能酶家族边界；
- 标注高度保守的催化残基；
- 产生可靠正例和负例；
- 给 Pythia-Pocket / GRASE 校准提供训练或验证集合。

不要只保留高 identity hit。应保存两类结果：

1. **高置信同源候选**：用于实验阳性和 family anchor；
2. **边界候选**：低 identity 但保留核心 motif，用于结构和口袋重排序。

### 3.2 Foldseek 结构线

用途：

- 找低序列相似但 domain fold 相近的候选；
- 构建结构邻域和 fold cluster；
- 区分“同 fold 新序列”和“真正不同 fold 的新骨架”。

Foldseek 应执行两种搜索：

1. **global/domain search**：关注整体 domain 或全长结构相似；
2. **local active-domain search**：只用携带活性位点的 domain，避免多结构域蛋白的无关区域干扰。

记录至少以下字段：

```text
query_id
target_id
query_coverage
target_coverage
alignment_score
tm_like_score_or_equivalent
evalue
aligned_catalytic_residues
active_site_rmsd_if_available
active_site_plddt_mean
pae_or_domain_confidence
source_database
```

### 3.3 Folddisco 局部 motif 线

用途：

- 找整体 fold 不明显相似但局部催化几何相似的候选；
- 搜索收敛进化或新骨架中的功能位点；
- 对 Foldseek 漏掉的候选补召回。

Folddisco query 不应只包含最小 catalytic triad。必须设计为 **active-site constellation**，包含催化残基、底物定位残基、氧阴离子洞或 amide-binding groove 等。

建议为每个酶类准备多个 motif ensemble：

```text
motifs/
  petase_cutinase/
    ser_his_asp_oxyanion_v1.yaml
    ser_his_asp_oxyanion_aromatic_clamp_v2.yaml
    substrate_bound_pose_ensemble_v3.yaml
  nylon_nylc/
    ntn_thr_asp_asp_v1.yaml
    ntn_thr_processing_site_groove_v2.yaml
    oligomer_binding_groove_v3.yaml
```

每个 motif 命中都要保存残基映射，而不是只保存分数。

---

## 4. PET/PUR 塑料降解酶特化规则

### 4.1 PETase / cutinase / esterase 类

优先关注：

- catalytic Ser-His-Asp/Glu；
- oxyanion hole；
- 疏水浅槽；
- 芳香或疏水夹持残基；
- 可接近的 polymer-binding surface；
- active-site loop/lid 的构象可信度；
- 高温或溶剂条件下稳定性。

主要风险：普通 α/β hydrolase 中 Ser-His-Asp/Glu motif 极多。必须用口袋、底物可达性和 docking 几何过滤假阳性。

### 4.2 Polyurethane urethanase 类

优先关注：

- urethane bond 定位能力；
- glycolysis / solvent 条件兼容性；
- 高稳定性和可表达性；
- 反应口袋对商业 polyurethane 片段的可接近性。

GRASE 类模型在这里应作为主要 ranker 之一，但不能替代几何和实验验证。

---

## 5. NylC / nylon hydrolase 特化规则

NylC 不应被当作普通丝氨酸水解酶处理。它更接近 Ntn hydrolase 语境，需要关注成熟链 N-terminal nucleophile 和自加工/组装相关因素。

优先关注：

- N-terminal Thr 等价亲核残基；
- Asp-Asp-Thr 或等价酸性残基网络；
- 自切割/成熟加工位点；
- polyamide oligomer binding groove；
- 是否能容纳 6-aminohexanoate oligomer；
- scissile amide bond 是否能被放置到正确亲核攻击几何；
- 可能的寡聚状态和界面稳定性。

NylC 的 Folddisco motif 应至少包含：

```text
nucleophile_thr_equivalent
acidic_residue_1_equivalent
acidic_residue_2_equivalent_or_general_base_network
processing_site_context
amide_binding_groove_residues
optional_oligomer_interface_contacts
```

NylC 的 Pythia-Pocket/GRASE 重排序应强调：

- 口袋长度；
- 极性/疏水分布是否适合 polyamide；
- amide chain 是否能以 productive pose 进入；
- 活性位点附近 pLDDT / PAE 是否可信；
- 组装状态是否可能影响活性。

---

## 6. 推荐评分框架

第一版采用可解释加权评分。后续可替换为学习型 ranker。

```text
S_total =
  0.15 × sequence_evidence
+ 0.20 × fold_evidence
+ 0.25 × motif_geometry
+ 0.20 × pocket_embedding_similarity
+ 0.10 × predicted_stability
+ 0.05 × substrate_accessibility
+ 0.05 × novelty
- penalty_missing_required_catalytic_residue
- penalty_low_active_site_confidence
- penalty_incompatible_oligomeric_state
- penalty_bad_expression_or_solubility_proxy
```

按目标调整权重：

- 快速得到可测活性：提高 `sequence_evidence`、`fold_evidence`；
- 发现新骨架：提高 `motif_geometry`、`pocket_embedding_similarity`、`novelty`；
- 工业条件：提高 `predicted_stability`、溶剂/温度相关 proxy。

候选最终分四类输出：

| 类别 | 定义 | 用途 |
|---|---|---|
| conservative_positive | 高序列/结构相似，高稳定性 | 确保实验板有阳性 |
| remote_homolog_or_remote_fold | 低序列相似但 Foldseek 强 | 远缘扩展 |
| novel_scaffold_explorer | Folddisco + pocket 高分，但序列和整体 fold 弱 | 主攻新骨架 |
| negative_or_boundary_control | 缺关键残基或口袋不合理 | 校准模型和实验背景 |

---

## 7. 推荐仓库结构

Codex 实现时优先建立如下轻量结构。不要提交大型数据库、checkpoint 或原始轨迹。

```text
projects/03-new-scaffold-enzyme-search/
  README.md
  CODEX_TASKS.md
  configs/
    targets.example.yaml
    scoring.default.yaml
    databases.example.yaml
  data/
    seed_enzymes.example.tsv
    catalytic_residues.example.tsv
  motifs/
    petase_cutinase.example.yaml
    nylc_ntn_hydrolase.example.yaml
  src/
    scaffold_search/
      __init__.py
      io.py
      schema.py
      sequence_search.py
      structure_search.py
      motif_search.py
      pocket_rank.py
      scoring.py
      reporting.py
  scripts/
    run_sequence_search.py
    run_foldseek_search.py
    run_folddisco_search.py
    run_pocket_ranking.py
    merge_and_rank_candidates.py
  reports/
    README.md
```

External databases should be referenced by manifest files only.

---

## 8. 最小可运行产物

Codex 第一轮不要追求完整跑通真实 AFDB / PDB / metagenome 全库。先完成：

1. schema；
2. config examples；
3. wrapper interfaces；
4. toy input；
5. deterministic scoring；
6. candidate merge and report generation；
7. tests for scoring and schema validation。

真实工具调用通过 wrappers 接入：

```text
MMseqs2 / HMMER wrappers
Foldseek wrapper
Folddisco wrapper
Pythia-Pocket / Pythia / GRASE wrapper
Docking / geometry sanity wrapper
```

如果外部工具不可用，wrapper 必须给出清晰错误，并允许读取预先生成的 TSV/JSON 结果继续后续步骤。

---

## 9. 输出文件标准

最终 ranking 表至少包含：

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

每个候选必须有可解释理由。例如：

```text
Candidate X is classified as novel_scaffold_explorer because it lacks detectable sequence homology to seed PETase enzymes, has weak global Foldseek similarity, but matches the Ser-His-Asp + oxyanion-hole + aromatic-clamp motif and has a high Pythia-Pocket similarity score for the shallow hydrophobic PET-binding groove.
```

---

## 10. 实验板设计建议

每轮实验建议 96-well 或 384-well 板中至少包含：

- 20–40 个 conservative_positive；
- 20–40 个 remote_homolog_or_remote_fold；
- 20–80 个 novel_scaffold_explorer；
- 10–20 个 negative_or_boundary_control；
- 已知阳性酶、失活突变体、空载体、底物空白。

实验结果必须回流：

```text
assay_results.tsv
  candidate_id
  substrate
  condition
  assay_signal
  normalized_activity
  expression_level
  soluble_fraction
  stability_assay
  active_or_inactive_label
  notes
```

回流后更新：

1. HMM profiles；
2. motif ensemble；
3. scoring weights；
4. pocket ranker calibration；
5. negative set。

---

## 11. Codex 工作原则

1. 先读本 `README.md`，再读 `CODEX_TASKS.md`。
2. 第一轮只实现轻量 pipeline 和数据格式，不下载大数据库。
3. 所有外部工具都通过 wrapper 隔离。
4. 每一步必须输出可追踪中间文件。
5. 所有 ranking 都必须能解释，不允许只有黑箱分数。
6. 不提交大型二进制、checkpoint、结构库、trajectory、docking ensemble 或 scratch notebook。
7. 新增代码需有基本测试，尤其是 schema validation、score aggregation、candidate classification。
