# Project 03: 同功能新骨架酶搜索计划

本项目面向 **在现有功能酶基础上搜索具有相同功能的新骨架酶**。典型目标包括：

- 塑料降解相关酶：PETase、cutinase、lipase/esterase、polyurethane urethanase 等；
- NylC / nylon hydrolase / 6-aminohexanoate oligomer hydrolase 等 polyamide 降解酶；
- 未来可扩展到其他明确催化机制、底物和活性位点的酶类。

当前版本的核心思想：

> **Sequence search 和 Foldseek 主要负责召回；Folddisco 不再是单个 motif 硬筛选器，而是多层 active-site / pocket constellation 证据系统；Folddisco 同时承担 after-recall 机制验证和 de-novo 新骨架扫库；Pythia-Pocket / GRASE 对 Folddisco 支持的候选做口袋环境、活性潜力和稳定性重排序。**

详细 motif 规范见：[`MOTIF_STRATEGY.md`](MOTIF_STRATEGY.md)。

---

## 1. 总体流程

不要把流程做成单一路径。推荐数据流是：

```text
Known functional enzymes
  │
  ├─ Sequence search: MMseqs2 / BLAST / HMMER
  │    → C_sequence
  │
  ├─ Foldseek global/local structure search
  │    → C_foldseek
  │
  └─ Layered Folddisco motif ensemble construction
       → M_target = {L0, L1, L2, L3}

C_recalled = C_sequence ∪ C_foldseek

Layer A: Folddisco-after-recall
  Folddisco(M_target, structures_of_C_recalled)
    → E_folddisco_after_recall

Layer B: Folddisco-de-novo
  Folddisco(M_target, broad_structure_database)
    → C_folddisco_de_novo

C_total = C_recalled ∪ C_folddisco_de_novo

C_total + all evidence
  → merge / cluster / deduplicate
  → Pythia-Pocket / GRASE ranking
  → docking / geometry sanity check
  → experimental panel
```

当前阶段 **不需要重新跑全部 sequence search 和 Foldseek**。优先工作是：

1. 重新定义 layered Folddisco motif ensemble；
2. 将已有 sequence / Foldseek / Folddisco 结果映射到新 evidence schema；
3. 对已有 `C_sequence ∪ C_foldseek` 补跑缺失的 `after_recall` 层级；
4. 对大结构库或代表性 subset 跑 `de_novo`；
5. 重新合并、分类、排序。

---

## 2. 工具分工

| 层级 | 工具 | 主要职责 | 输出证据 |
|---|---|---|---|
| 召回 | BLAST / MMseqs2 / HMMER / jackhmmer | 找近缘和远缘同源，建立 MSA、HMM、保守残基图谱 | `sequence_score`、家族分支、同源正负例 |
| 召回 | Foldseek | 搜索低序列相似但整体 fold 或 domain 相似的候选 | `fold_score`、coverage、TM-score-like 指标、domain alignment |
| 机制验证 | Folddisco-after-recall | 对 `C_sequence ∪ C_foldseek` 做分层 active-site / pocket motif 验证 | `folddisco_mode=after_recall`、L0–L3 分数、残基映射 |
| 新骨架召回 | Folddisco-de-novo | 用 layered motif 直接扫大结构库，寻找 sequence/Foldseek 漏掉的新骨架 | `folddisco_mode=de_novo`、`motif_source=de_novo` |
| 口袋/活性/稳定性重排序 | GRASE / Pythia-Pocket / Pythia | 对 Folddisco 支持或边界保留候选做口袋级和稳定性排序 | `pocket_score`、`active_site_confidence`、`stability_score` |
| 物理合理性检查 | docking / PLACER / MD / 几何规则 | 检查底物 pose、亲核攻击距离、过渡态稳定、入口开放性 | `substrate_pose_score`、`mechanistic_sanity_score` |

---

## 3. Layered Folddisco motif / pocket 定义

Folddisco motif 应从“单个 rigid motif”升级为 **多层 active-site / pocket constellation**。

| 层级 | 名称 | 内容 | 作用 | 是否硬筛 |
|---|---|---|---|---|
| L0 | catalytic core | 必需催化残基、亲核残基、general acid/base、关键金属配位残基等 | 判断是否具备最小反应机制 | 通常硬筛，但允许 seed-specific exception |
| L1 | catalytic environment | oxyanion hole、acidic network、amide-binding polar residues、过渡态稳定残基 | 判断反应几何是否完整 | 半硬筛 / 强加权 |
| L2 | substrate pocket | 底物 4–6 Å 接触残基、polymer/polyamide groove、hydrophobic/aromatic clamp、口袋入口 | 判断底物是否可进入和被正确定位 | 加权，不直接淘汰 |
| L3 | processing / interface / dynamics | NylC processing site、可能影响成熟加工的残基、寡聚界面、lid loop、远端 allosteric residues | 辅助解释表达、成熟、组装和活性状态 | 辅助证据，不硬筛 |

关键规则：

```text
L0/L1 = 机制核心与反应环境
L2/L3 = 口袋、底物定位、processing、interface、dynamics 上下文
```

不要把 L2/L3 残基当成 L0 strict required residue。比如 T227 这类可能影响口袋、processing 或界面的残基，应进入 L2/L3，而不应与 catalytic core 一样作为 strict 必需点。

---

## 4. Folddisco-after-recall

输入：

```text
C_recalled = C_sequence ∪ C_foldseek
structures_of_C_recalled
M_target = layered motif ensemble
```

目的：

- 验证 sequence / Foldseek hit 是否保留目标催化几何；
- 区分“同源但可能失活”、“同 fold 但底物谱不对”和“机制上仍然可信”的候选；
- 给 Pythia-Pocket / GRASE 提供更干净的候选池；
- 形成可解释证据，而不是一刀切过滤。

候选处理规则：

```text
L0 STRICT_PASS:
  保留并强化机制证据。

L0 PARTIAL + L1/L2 strong:
  保留为 needs_review 或 exploration candidate。

L0 NO_HIT + active_site_confidence high:
  降级为 negative_or_boundary_control。

L0 NO_HIT + active_site_confidence low:
  标记为 needs_review，不直接删除。
```

禁止实现为：

```text
full motif fail → delete
```

必须实现为：

```text
layered motif evidence → classify / downgrade / review
```

---

## 5. Folddisco-de-novo

输入：

```text
M_target = layered motif ensemble
broad_structure_database = PDB / AFDB subset / ESMAtlas subset / metagenome structure set
```

目的：

- 直接搜索大结构库；
- 找 sequence / Foldseek 低分但局部机制几何强的新骨架；
- 产生 `novel_scaffold_explorer` 候选。

典型目标模式：

```text
sequence_score low
fold_score low
l0_core_score high
l1_environment_score medium/high
l2_pocket_score medium/high
pocket_score high after Pythia/GRASE
```

如果 Folddisco 只作为 `C_sequence ∪ C_foldseek` 的后置筛选器，这类候选会被漏掉。因此 de-novo 分支必须保留。


### 5.1 Folddisco 命中结构找回

当 de-novo Folddisco 命中存在，但 AFDB/RCSB/ESMAtlas 直接 URL 下载失败时，不要直接判定为缺结构。先按 [`reports/FOLDDISCO_STRUCTURE_RECOVERY_RUNBOOK.md`](reports/FOLDDISCO_STRUCTURE_RECOVERY_RUNBOOK.md) 使用 Folddisco result API 通过 `ticket + database + dbkey/id` 找回结构，并对 API 返回的 PDB-like 文件做 fixed-column 标准化后再进入 Pythia-Pocket、RMSD 或 docking 分析。

---

## 6. 推荐 Folddisco 输出字段

最终 ranking 表和中间 evidence 表都应支持以下字段：

```text
candidate_id
target_id
motif_id
motif_version
folddisco_mode                 # after_recall | de_novo | both
motif_source                   # validated_recall | de_novo | both | none
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
best_motif_layer
residue_mapping
missing_required_labels
partial_match_labels
active_site_confidence
hit_status                     # pass | partial | fail | low_confidence | not_run
notes
```

第一版 Folddisco 分层评分：

```text
folddisco_score =
  0.45 × l0_core_score
+ 0.25 × l1_environment_score
+ 0.20 × l2_pocket_score
+ 0.10 × l3_context_score
```

---

## 7. 总评分框架

第一版采用可解释加权评分。后续可替换为学习型 ranker。

```text
S_total =
  0.15 × sequence_evidence
+ 0.15 × fold_evidence
+ 0.30 × folddisco_score
+ 0.20 × pocket_embedding_similarity
+ 0.10 × predicted_stability
+ 0.05 × substrate_accessibility
+ 0.05 × novelty
- penalty_missing_l0_core_when_confident
- penalty_low_active_site_confidence
- penalty_incompatible_oligomeric_state
- penalty_bad_expression_or_solubility_proxy
```

`folddisco_score` 由 `after_recall` 和 `de_novo` 两类证据合并：

```text
folddisco_score = max(
  folddisco_after_recall_score,
  folddisco_denovo_score
)
```

并保留：

```text
motif_source = validated_recall | de_novo | both | none
```

---

## 8. 进入 Pythia-Pocket / GRASE 的候选规则

进入 Pythia-Pocket / GRASE 的候选不应只限于 full strict pass。建议包括：

1. L0 strict pass 的候选；
2. Folddisco-de-novo 直接扫库命中的候选；
3. L0 partial 但 L1/L2 很强的候选；
4. sequence / Foldseek 强召回但 Folddisco 不完整的边界候选；
5. 少量 negative / boundary controls，用于校准模型和实验背景。

这保证 GRASE/Pythia-Pocket 不会只看到“完美 motif”，也能学习或校准边界情况。

---

## 9. PET/PUR 塑料降解酶特化规则

### PETase / cutinase / esterase

L0 catalytic core：

- Ser nucleophile；
- His general base；
- Asp/Glu acid；
- Ser-His-Asp/Glu 距离和方向约束。

L1 catalytic environment：

- oxyanion hole；
- backbone NH 或 side-chain donor；
- tetrahedral intermediate stabilization environment；
- catalytic loop 局部可信度。

L2 substrate pocket：

- 疏水浅槽；
- aromatic clamp；
- polymer-binding surface；
- 底物 scissile bond 定位残基；
- 入口开放性。

L3 processing / dynamics / interface：

- lid loop 或 flexible loop；
- 热稳定相关远端区域；
- 可能影响界面吸附的 surface patch；
- 多聚体或晶体接触导致的口袋遮挡风险。

### Polyurethane urethanase

重点增加：

- urethane bond 定位；
- glycolysis / solvent 条件兼容性；
- 高稳定性和可表达性；
- 商业 polyurethane 片段可接近性。

GRASE 类模型可作为强 ranker，但不能替代 L0/L1 机制几何检查。

---

## 10. NylC / nylon hydrolase 特化规则

NylC 不应被当作普通 Ser-His-Asp esterase 处理。NylC residue labels 必须来自 seed-specific annotation，不要在代码中硬编码某个编号。不同 NylC-like 蛋白的成熟链编号、processing site 和结构域边界可能不同。

L0 catalytic core：

- N-terminal Thr 或等价亲核残基；
- 直接参与亲核攻击或 general acid/base 网络的核心残基；
- 用户或 seed annotation 指定的最小 catalytic core，例如某些项目中可用 `A189/A219/A267` 这类 seed-specific label 表示。

L1 catalytic environment：

- Asp/Glu acidic network；
- amide bond polarization residues；
- oxyanion-like stabilization environment；
- processing-site-proximal polar residues。

L2 substrate pocket：

- polyamide / 6-aminohexanoate oligomer binding groove；
- 4–6 Å 底物接触残基；
- 适合 amide chain 的极性/疏水分布；
- oligomer 入口和通道宽度。

L3 processing / interface：

- processing site；
- T227 这类可能影响局部口袋、processing 或界面的残基；
- 寡聚界面；
- 成熟链稳定性和组装状态相关区域。

---

## 11. 候选分类

最终候选至少分为：

| 类别 | 定义 | 用途 |
|---|---|---|
| conservative_positive | sequence/Foldseek 强，Folddisco-after-recall 支持或不冲突 | 确保实验板有阳性 |
| remote_homolog_or_remote_fold | sequence 较弱，Foldseek 强，L0/L1 支持 | 远缘扩展 |
| novel_scaffold_explorer | sequence/Foldseek 弱，但 Folddisco-de-novo + pocket 高分 | 主攻新骨架 |
| needs_review | L0 partial/fail 但 L1/L2/pocket 支持，或 active-site confidence 低 | 重新建模、人工审查、小比例探索 |
| negative_or_boundary_control | 缺关键 L0 且局部结构可信，或口袋/机制明显不合理 | 校准模型和实验背景 |

---

## 12. 推荐仓库结构

Codex 实现时优先建立如下轻量结构。不要提交大型数据库、checkpoint 或原始轨迹。

```text
projects/03-new-scaffold-enzyme-search/
  README.md
  CODEX_TASKS.md
  MOTIF_STRATEGY.md
  configs/
    targets.example.yaml
    scoring.default.yaml
    databases.example.yaml
  data/
    seed_enzymes.example.tsv
    catalytic_residues.example.tsv
    toy_sequence_hits.tsv
    toy_foldseek_hits.tsv
    toy_folddisco_after_recall_hits.tsv
    toy_folddisco_denovo_hits.tsv
    toy_pocket_hits.tsv
  motifs/
    petase_cutinase.layered.example.yaml
    nylc_ntn_hydrolase.layered.example.yaml
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
    run_folddisco_after_recall.py
    run_folddisco_denovo.py
    run_pocket_ranking.py
    merge_and_rank_candidates.py
  reports/
    README.md
    FOLDDISCO_STRUCTURE_RECOVERY_RUNBOOK.md
```

External databases should be referenced by manifest files only.

---

## 13. 实验板设计建议

每轮实验建议 96-well 或 384-well 板中至少包含：

- 20–40 个 conservative_positive；
- 20–40 个 remote_homolog_or_remote_fold；
- 20–80 个 novel_scaffold_explorer；
- 5–20 个 needs_review；
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
2. layered motif ensemble；
3. L0/L1/L2/L3 权重；
4. Folddisco-after-recall 阈值；
5. Folddisco-de-novo motif ensemble；
6. pocket ranker calibration；
7. negative set。

---

## 14. Codex 工作原则

1. 先读本 `README.md`，再读 `MOTIF_STRATEGY.md`，最后读 `CODEX_TASKS.md`。
2. 第一轮只实现轻量 pipeline 和数据格式，不下载大数据库。
3. 所有外部工具都通过 wrapper 隔离。
4. Folddisco 必须支持 `after_recall` 和 `de_novo` 两种模式。
5. Folddisco evidence 必须支持 L0/L1/L2/L3 分层。
6. 不把 T227 或类似环境残基默认当作 L0 strict required。
7. 不因 full motif fail 直接删除候选。
8. 每一步必须输出可追踪中间文件。
9. 所有 ranking 都必须能解释，不允许只有黑箱分数。
10. 不提交大型二进制、checkpoint、结构库、trajectory、docking ensemble 或 scratch notebook。
11. 新增代码需有基本测试，尤其是 schema validation、score aggregation、candidate classification、Folddisco evidence merge 和 layered motif parsing。
