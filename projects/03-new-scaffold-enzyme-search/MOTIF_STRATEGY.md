# Layered Folddisco motif / pocket strategy

本文件定义 Project 03 中 Folddisco 的新定位：**Folddisco 不再是单个 motif 硬筛选器，而是多层 active-site / pocket constellation 证据系统**。

Sequence search 和 Foldseek 的主要职责是高召回。Folddisco 有两个职责：

1. `after_recall`：对 sequence / Foldseek 已召回候选做机制几何验证；
2. `de_novo`：直接扫大结构库，寻找 sequence / Foldseek 漏掉的新骨架。

Pythia-Pocket / GRASE 不替代 Folddisco，而是在 Folddisco 给出机制几何证据后，对口袋环境、活性潜力、稳定性和实验优先级做重排序。

---

## 1. 多层 motif 原则

每个 target 不再只定义一个严格 motif，而是定义一个 motif ensemble。每个 motif 按层级拆分为：

| 层级 | 名称 | 内容 | 作用 | 是否硬筛 |
|---|---|---|---|---|
| L0 | catalytic core | 必需催化残基、亲核残基、general acid/base、关键金属配位残基等 | 判断是否具备最小反应机制 | 通常硬筛，但允许 seed-specific exception |
| L1 | catalytic environment | oxyanion hole、acidic network、amide-binding polar residues、过渡态稳定残基 | 判断反应几何是否完整 | 半硬筛 / 强加权 |
| L2 | substrate pocket | 底物 4–6 Å 接触残基、polymer/polyamide groove、hydrophobic/aromatic clamp、口袋入口 | 判断底物是否可进入和被正确定位 | 加权，不直接淘汰 |
| L3 | processing / interface / dynamics | NylC processing site、可能影响成熟加工的残基、寡聚界面、lid loop、远端 allosteric residues | 辅助解释表达、成熟、组装和活性状态 | 辅助证据，不硬筛 |

关键原则：**L0 与 L1 更接近机制必需项，L2 与 L3 是口袋和上下文证据，不应与 catalytic core 混为同等严格条件。**

例如用户提到的 T227 类残基，如果它主要影响口袋环境、底物定位、processing 或界面，而不是直接 catalytic core，就应放入 L2/L3，不能作为 L0 strict required residue。

---

## 2. Folddisco 的两个运行模式

### 2.1 `after_recall`

输入：

```text
C_recalled = C_sequence ∪ C_foldseek
structures_of_C_recalled
motif ensemble M_target = {L0, L1, L2, L3}
```

目的：

- 验证 sequence / Foldseek hit 是否保留目标机制几何；
- 从同源或同 fold 候选中剔除明显失活、底物谱不兼容或口袋错误的条目；
- 产生可解释证据链，而不是简单 pass/fail。

输出：

```text
folddisco_mode = after_recall
core_motif_status
core_rmsd
required_match_fraction
context_match_fraction
pocket_context_score
best_motif_layer
layer_scores
residue_mapping
hit_status
```

### 2.2 `de_novo`

输入：

```text
broad_structure_database = PDB / AFDB subset / ESMAtlas subset / metagenome structures
motif ensemble M_target = {L0, L1, L2, L3}
```

目的：

- 直接搜索大结构库；
- 找 sequence / Foldseek 低分但局部机制几何强的新骨架；
- 产生 `novel_scaffold_explorer` 候选。

输出必须标记：

```text
folddisco_mode = de_novo
motif_source = de_novo
```

---

## 3. 推荐字段

每个 Folddisco hit 至少输出：

```text
candidate_id
target_id
motif_id
motif_version
folddisco_mode                 # after_recall | de_novo
motif_source                   # validated_recall | de_novo | both | none, merge 后确定
core_motif_status              # STRICT_PASS | PARTIAL | NO_HIT | NOT_EVALUATED
l0_core_score
l1_environment_score
l2_pocket_score
l3_context_score
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

`layer_scores` 可以用 JSON 字段保存，例如：

```json
{
  "L0": 0.95,
  "L1": 0.80,
  "L2": 0.55,
  "L3": 0.20
}
```

---

## 4. 推荐评分

### 4.1 单个 Folddisco hit 层级评分

第一版采用可解释加权：

```text
folddisco_score =
  0.45 × l0_core_score
+ 0.25 × l1_environment_score
+ 0.20 × l2_pocket_score
+ 0.10 × l3_context_score
```

其中：

```text
l0_core_score:
  catalytic core 残基匹配和几何约束。

l1_environment_score:
  oxyanion hole / acidic network / amide-binding polar residues 等局部反应环境。

l2_pocket_score:
  底物接触残基、groove、clamp、入口开放性等。

l3_context_score:
  processing site、oligomer interface、lid loop、远端上下文。
```

### 4.2 硬筛逻辑

只允许 L0 在明确情况下触发强降级：

```text
L0 STRICT_PASS:
  保留并增强证据。

L0 PARTIAL + L1/L2 strong:
  保留为 needs_review 或 exploration candidate。

L0 NO_HIT + active_site_confidence high:
  降级为 negative_or_boundary_control。

L0 NO_HIT + active_site_confidence low:
  标记为 needs_review，不直接删除。
```

不要执行：

```text
full motif fail → delete
```

必须执行：

```text
layered motif evidence → classify / downgrade / review
```

---

## 5. 进入 Pythia-Pocket / GRASE 的候选

进入 Pythia-Pocket / GRASE 的候选不应只限于 full strict pass。建议包括：

1. L0 strict pass 的候选；
2. Folddisco-de-novo 直接扫库命中的候选；
3. L0 partial 但 L1/L2 很强的候选；
4. sequence / Foldseek 强召回但 Folddisco 不完整的边界候选；
5. 少量 negative / boundary controls，用于校准模型和实验背景。

---

## 6. PETase / cutinase / esterase motif 层级

### L0 catalytic core

- Ser nucleophile；
- His general base；
- Asp/Glu acid；
- 核心 Ser-His-Asp/Glu 距离和方向约束。

### L1 catalytic environment

- oxyanion hole；
- backbone NH 或 side-chain donor；
- 支持 tetrahedral intermediate 的极性环境；
- catalytic loop 局部可信度。

### L2 substrate pocket

- 疏水浅槽；
- aromatic clamp；
- polymer-binding surface；
- 底物 scissile bond 定位残基；
- 入口开放性。

### L3 processing / dynamics / interface

- lid loop 或 flexible loop；
- 热稳定相关远端区域；
- 可能影响界面吸附的表面 patch；
- 多聚体或晶体接触导致的口袋遮挡风险。

---

## 7. NylC / nylon hydrolase motif 层级

NylC 的 residue labels 必须来自 seed-specific annotation，不要在代码中硬编码某个编号。不同 NylC-like 蛋白的成熟链编号、processing site 和结构域边界可能不同。

### L0 catalytic core

- N-terminal Thr 或等价亲核残基；
- 直接参与亲核攻击或 general acid/base 网络的核心残基；
- 用户或 seed annotation 指定的最小 catalytic core，例如某些项目中可用 `A189/A219/A267` 这类 seed-specific label 表示。

### L1 catalytic environment

- Asp/Glu acidic network；
- amide bond polarization residues；
- oxyanion-like stabilization environment；
- processing-site-proximal polar residues。

### L2 substrate pocket

- polyamide / 6-aminohexanoate oligomer binding groove；
- 4–6 Å 底物接触残基；
- 适合 amide chain 的极性/疏水分布；
- oligomer 入口和通道宽度。

### L3 processing / interface

- processing site；
- 例如 T227 这类可能影响局部口袋、processing 或界面的残基；
- 寡聚界面；
- 成熟链稳定性和组装状态相关区域。

---

## 8. Motif YAML 推荐格式

每个 motif YAML 应支持层级：

```yaml
motif_id: nylc_layered_pocket_constellation
version: toy_v1
target_id: nylc_ntn_hydrolase
intended_modes: [after_recall, de_novo]

layers:
  L0_catalytic_core:
    hard_filter: true
    weight: 0.45
    residues:
      - label: Thr_nucleophile
        allowed_residue_names: [THR]
        role: n_terminal_nucleophile
        required: true
      - label: Core_acid_base_1
        allowed_residue_names: [ASP, GLU, HIS, SER, THR]
        role: acid_base_network
        required: true
    geometry_constraints:
      - pair: [Thr_nucleophile, Core_acid_base_1]
        distance_angstrom: [2.5, 5.5]

  L1_catalytic_environment:
    hard_filter: false
    weight: 0.25
    residues:
      - label: Amide_polarizer_1
        allowed_residue_names: [ASN, GLN, SER, THR, TYR, HIS, BACKBONE_NH]
        role: amide_polarization
        required: false

  L2_substrate_pocket:
    hard_filter: false
    weight: 0.20
    residues:
      - label: Polyamide_groove_contact_1
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
```

---

## 9. 推荐执行顺序

当前阶段不需要重跑全部 sequence / Foldseek。推荐顺序：

1. 重新定义 layered Folddisco motif ensemble；
2. 将已有 sequence / Foldseek / Folddisco 结果映射到新 evidence schema；
3. 对已有 `C_sequence ∪ C_foldseek` 补跑缺失的 `after_recall` 层级；
4. 对大结构库或代表性 subset 跑 `de_novo`；
5. 重新合并候选；
6. 将 L0 strict pass、de_novo hit、L0 partial + L1/L2 strong、强召回边界候选和负控送入 Pythia-Pocket / GRASE；
7. 输出新的 ranking 和实验板。

---

## 10. Codex 实现要求

Codex 实现时必须：

- 支持 layered motif YAML；
- 在 schema 中区分 L0/L1/L2/L3 分数；
- 在报告中展示每个 top candidate 的分层证据；
- 不把 T227 或类似环境残基作为默认 L0 strict required；
- 不因 full motif fail 直接删除候选；
- 允许不同 target 用不同层级权重；
- 允许已有结果重新映射，不强制重跑 sequence / Foldseek。
