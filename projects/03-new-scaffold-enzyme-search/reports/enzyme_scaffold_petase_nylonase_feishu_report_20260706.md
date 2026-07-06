# PETase 与尼龙酶骨架搜索筛选方法及初步统计

生成日期：2026-07-06  
项目服务器路径：`/Dell/Dell14/bianht/enzyme_scaffold_search_v2`  
适用范围：PETase 与尼龙酶（NylC / Nyl50 / Nyl10）骨架搜索、口袋分层统计、Pythia-pocket 评分前的交接报告。

## 1. 研究目标

本项目的目标不是只找到与已知酶序列相似的同源蛋白，而是尽可能系统地召回可能具有相似催化口袋或可迁移功能环境的新骨架。当前策略以“高召回、分层描述、后续再排序”为原则：

1. 用序列搜索扩大同源和远缘同源候选。
2. 用 Foldseek 从整体结构相似性角度召回已知或预测结构。
3. 用 Folddisco 从活性残基和周围口袋 motif 的结构几何角度召回候选。
4. 将 sequence / Foldseek / Folddisco 的候选合并去重，避免过早用单一方法排除潜在新骨架。
5. 对合并后的候选计算全局序列相似性、整体结构 RMSD、Pythia-pocket embedding cosine similarity，再进行多指标排序和人工审查。

关键原则：Folddisco、Foldseek、BLAST/MMseqs2 都是召回工具。工具返回 hit 不能直接等同于功能确认；下载失败、缺结构、缺序列、没有 predicted pocket 也不能写成生物学失败，应标记为 `NOT_EVALUATED` 或对应技术状态。

## 2. 总体筛选流程

| 阶段 | 工具/方法 | 目的 | 当前解释口径 |
|---|---|---|---|
| 序列大库召回 | BLASTP / MMseqs2 / HMMER | 从 UniRef90、UniProtKB/TrEMBL 等大库召回同源和远缘同源序列 | tool-native hit，只作为候选来源 |
| 结构召回 | Foldseek | 通过整体或局部结构相似性召回结构候选 | raw alignment 与 unique target 分开统计 |
| 口袋/基序召回 | Folddisco | 用催化残基和周围 residue layer 定义 motif，检索结构库 | Folddisco hit 是检索命中，不是最终功能 pass |
| Layer 统计 | L0、L0+L1、L0+L2、L0+L1+L2 | 描述候选对核心催化残基及周围口袋的支持程度 | 先统计，不作为进入 Pythia 前的硬筛 |
| Pocket embedding | Pythia-pocket | 由模型预测 pocket residues，并计算 pocket embedding cosine | 当前使用 `prob_threshold = 0.6` |
| 全局序列相似性 | Needleman-Wunsch global alignment | 计算全长序列差异 | 使用 `PID_global = identical_residues / alignment_columns_including_gaps` |
| 整体骨架差异 | Foldseek/TM-align 类 RMSD 统计 | 评估候选整体骨架是否与 seed 接近或偏离 | 用于散点图和排序解释 |

## 3. Layer 与口袋定义

当前 layer 是项目定义的分层口袋描述，不是 Folddisco 自带的唯一标准。Folddisco 自身负责按给定 motif 和数据库参数返回结构命中；我们另外统计不同 layer 下候选对 motif 的支持情况。

| Layer | 定义 | 用途 |
|---|---|---|
| L0 | 核心催化残基或最关键的反应相关残基 | 作为 Pythia-pocket 对比和结构 motif 的核心 anchor |
| L1 | L0 周围更近的口袋支撑残基，通常直接参与底物定位、催化几何或局部环境 | 用来描述核心催化环境是否完整 |
| L2 | 更外层的口袋环境或二级支撑残基，可能影响底物通道、形状和稳定性 | 用来描述口袋上下文是否接近 |
| L0+L1 | 核心催化残基加近邻口袋 | 比 L0 更严格，但仍避免过度收缩 |
| L0+L2 | 核心催化残基加外层环境 | 用来捕获可能保留核心和外层形状但近邻发生替换的候选 |
| L0+L1+L2 | 当前最完整的项目口袋层 | 用于后续优先级解释，不作为唯一进入 Pythia 的门槛 |

这里的“通过 layer”应理解为“该候选在对应 layer 下有结构映射/几何支持”，而不是功能已确认。后续 Pythia-pocket cosine、整体 RMSD 和全局序列差异会一起用于排序。

## 4. PETase 种子与搜索

### 4.1 PETase seed 来源

当前 PETase 有 5 个结构 seed。早期序列大库搜索使用了 4 个有明确 UniProt accession 的 seed；后续已补跑 PET_05_3VIS 的 PDB-derived sequence。五个 seed 不是互相替代关系，而是用于覆盖不同 PET hydrolase / cutinase scaffold。

| Seed | Accession / 来源 | 代表结构 | 说明 | 状态 |
|---|---|---|---|---|
| PET_01 | A0A0K8P6T7 | 6ILW | IsPETase，典型 mesophilic PETase | accepted |
| PET_02 | G9BY57 | 4EB0 | LCC，高温 PET hydrolase/cutinase | accepted |
| PET_03 | Q6A0I4 | 4CG1 | TfCut2 / Thermobifida cutinase scaffold | accepted |
| PET_04 | W0TJ64 | 4WFI / 4WFJ | Cut190，Ca2+ regulated PET-active scaffold | accepted |
| PET_05 | PDB-derived sequence | 3VIS | Est119/TaCut-like cutinase，用于补充 thermophilic actinobacterial cutinase 多样性 | accession 待最终确认 |

### 4.2 PETase 序列搜索统计

| 搜索集合 | 方法/数据库 | 初步统计 |
|---|---|---:|
| 早期 PETase broad recall，4 个 accession seed | BLASTP / MMseqs2 / HMMER against UniRef90 + TrEMBL | 7,306 unique candidates |
| PET_05_3VIS 补充搜索 | BLASTP + MMseqs2 against UniRef90 + TrEMBL | 6,358 unique targets |
| 已有候选 FASTA 中可用于 seed 相似性计算的 PETase 序列 | candidate FASTA vs seed FASTA | 25,367 candidate sequences；25,258 有 seed hit |

PET_05_3VIS 补充搜索的细分结果：

| 方法 | 数据库 | Raw hits | Unique targets |
|---|---|---:|---:|
| BLASTP | UniRef90 | 3,433 | 1,730 |
| BLASTP | TrEMBL | 4,043 | 2,033 |
| MMseqs2 | UniRef90 | 1,936 | 1,936 |
| MMseqs2 | TrEMBL | 2,218 | 2,218 |
| 合并去重 | UniRef90 + TrEMBL | - | 6,358 |

### 4.3 PETase 结构候选统计

当前 PETase 结构候选来自 Foldseek 与 direct Folddisco 的合并结果。RMSD 报告中每个 PETase seed 对应 10,089 个候选行。

| 来源 | 行数 |
|---|---:|
| Foldseek-only | 3,959 |
| Folddisco-only | 6,030 |
| Foldseek and Folddisco overlap | 100 |
| 每个 PETase seed 可评估结构候选行 | 10,089 |

解释：这里的行数是结构/RMSD 可评估集合，不等同于序列搜索全集。PETase 后续需要把 PET_05_3VIS 补充搜索结果合并到 PETase broad recall 后再去重，并将新增候选进入结构预测和 embedding 队列。

## 5. 尼龙酶种子与搜索

### 5.1 尼龙酶 seed 来源

尼龙酶当前以 NylC、Nyl50、Nyl10 三类为核心。NylC 的 WT / GYAQ / HP 多数属于同一大骨架背景；Nyl50 和 Nyl10 是更需要关注的新 scaffold / 功能方向。

| Seed | Scaffold group | 用途 | 说明 |
|---|---|---|---|
| NylC_WT_Q79F77 / NylC_WT_5XYO | NylC-like | sequence + structure | PA6-biased baseline NylC |
| NylC_GYAQ_5Y0M | NylC-like | sequence + structure | thermostable GYAQ background |
| NylC_HP_from_GYAQ | NylC-like | sequence + backbone model | HP mutation background，主要用于 backbone retrieval |
| Nyl50_9CXR_9DYS | Nyl50-like | sequence + structure | PA66-selective new scaffold candidate |
| Nyl10_A0A1M5P6R3 | Nyl10-like | sequence + structure/model | PA66-selective purified enzyme；AFDB v6 model 用于结构搜索 |

### 5.2 尼龙酶序列搜索统计

| 搜索集合 | 方法/数据库 | 初步统计 |
|---|---|---:|
| 早期 nylonase broad recall | BLASTP / MMseqs2 / HMMER against UniRef90 + TrEMBL | 33,675 unique candidates |
| PA6/PA66 扩展 seed：NylC WT/GYAQ/HP、Nyl50 | BLASTP + MMseqs2 against UniRef90 + TrEMBL | 已完成 |
| Nyl10 追加搜索 | BLASTP + MMseqs2 against UniRef90 + TrEMBL | 已完成 |
| 已有候选 FASTA 中可用于 seed 相似性计算的尼龙酶序列 | candidate FASTA vs seed FASTA | 15,083 candidate sequences；14,803 有 seed hit |

Nyl10 不是完全搜不到。它在大库序列搜索中有 self-hit，Nyl50 与 NylC 搜索也能命中 Nyl10，但相似度层级不同：Nyl50 到 Nyl10 的 BLAST pident 约 40.6%，query coverage 约 93%；NylC 到 Nyl10 的 BLAST pident 约 30-31%，query coverage 约 76%。这说明 Nyl10 与 NylC-like seed 的序列距离较远，不能只用序列相似性判断是否同类。

### 5.3 尼龙酶结构合并集合与 layer 统计

当前重新定义 layer 后的 nylonase merged candidate universe：

| 指标 | 数量 |
|---|---:|
| structure_merged_candidates_nylonase | 21,211 |
| with existing available structure | 4,035 |
| redefined Folddisco evidence rows | 43,306 |
| redefined Folddisco unique accessions | 24,010 |
| predicted-structure manifest 当前行数 | 3,192 rows / 3,166 unique candidates |
| predicted all-seed embedding 当前行数 | 6,200 rows / 3,106 unique candidates |

Layer 支持统计：

| Layer | Any support | Full support |
|---|---:|---:|
| L0 | 15,284 | 5,383 |
| L0+L1 | 6,701 | 159 |
| L0+L2 | 7,320 | 265 |
| L0+L1+L2 | 5,301 | 51 |

候选来源 route 统计：

| Route | Candidates |
|---|---:|
| direct_folddisco | 17,134 |
| foldseek | 3,555 |
| foldseek;direct_folddisco | 163 |
| sequence;direct_folddisco | 63 |
| sequence;foldseek | 275 |
| sequence;foldseek;direct_folddisco | 21 |

解释：上表中 direct_folddisco 是 Folddisco 数据库直接检索返回的候选；foldseek 是结构相似性检索返回的候选；sequence 代表序列搜索召回后进入结构相关流程的候选。各 route 的候选后续需要按 accession / sequence / structure ID 合并去重。

## 6. Pythia-pocket embedding 与相似性

当前 Pythia-pocket 的使用方式是：不直接用我们手工定义的 L0-L3 residue 作为最终 pocket，而是让 Pythia 按模型预测 pocket residues，然后提取 predicted-pocket embedding。我们使用 seed 的核心口袋作为比较对象，计算候选与 seed pocket embedding 的 cosine similarity。

当前参数与解释：

| 项目 | 当前设置 |
|---|---|
| Pocket residue 来源 | Pythia-pocket predicted pocket |
| Probability threshold | 0.6 |
| 相似性指标 | cosine similarity |
| 没有 predicted pocket | 标记为 `NO_PREDICTED_POCKET`，不等同功能失败 |
| 缺结构或解析失败 | 标记为 `NOT_EVALUATED` |

尼龙酶 all-seed predicted-pocket embedding 当前状态：

| 状态 | 数量 |
|---|---:|
| OK | 5,570 |
| NO_PREDICTED_POCKET | 630 |

按 seed 的 all-seed predicted-pocket rows：

| Seed | Rows |
|---|---:|
| Nyl10_A0A1M5P6R3 | 3,106 |
| Nyl50_9DYS | 2,790 |
| NylC_WT_5XYO | 304 |

已读取到的 NylC 与 Nyl10/Nyl50 的 Pythia-pocket cosine：

| 比较 | Cosine | 口径 |
|---|---:|---|
| NylC_WT_5XYO vs Nyl10_A0A1M5P6R3 | 0.860980 | Pythia predicted pocket, prob_threshold 0.6 |
| NylC_WT_5XYO vs Nyl50_9DYS | 0.925020 | Pythia predicted pocket, prob_threshold 0.6 |
| NylC_WT_5XYO vs NylC_WT_5XYO | 1.000000 | self |

## 7. 当前结果应如何解读

1. PETase 的 5 个结构 seed 已经在结构/口袋层面使用；早期序列大库层面原本只有 4 个明确 UniProt accession seed，PET_05_3VIS 已补跑。
2. 尼龙酶中 NylC、Nyl50、Nyl10 均已经进行大库序列搜索。Nyl10 不是大库中不存在，而是与 NylC-like 的序列距离较远，且结构/embedding 可评估集合与原始序列搜索全集不是同一个 denominator。
3. Folddisco 的结果要区分：
   - `FOLDDISCO_HIT`：Folddisco 在给定 motif 和数据库下返回了 hit。
   - `STRICT_PASS`：我们定义的完整 motif 覆盖和 RMSD 阈值满足。
   - `PARTIAL_OR_LOOSE_HIT`：有命中但不满足完整严格标准。
   - `NO_HIT`：结构被评估但没有命中。
   - `NOT_EVALUATED`：缺结构、下载失败、解析失败或尚未进入队列。
4. 目前不建议在进入 Pythia-pocket 前用 L0+L1+L2 做硬筛，因为这样可能丢掉有保守核心口袋但周围环境发生变化的新 scaffold。
5. 后续应统一输出 candidate-level 主表：candidate_id、enzyme_class、route、best_seed、sequence_available、structure_available、global_PID、global_RMSD、Pythia-pocket cosine、layer support、status label。

## 8. 下一步建议

| 优先级 | 任务 | 目的 |
|---|---|---|
| P0 | 合并 PET_05_3VIS 补充序列搜索结果到 PETase broad recall，并按 accession / sequence 去重 | 保证 PETase 五个 seed 在序列层面完整覆盖 |
| P0 | 对 PETase 与尼龙酶合并候选统一计算 `PID_global` | 避免 BLAST 局部相似性导致“高相似但 coverage 低”的误读 |
| P0 | 对已有结构和预测结构统一计算 Pythia-pocket cosine | 建立最终排序核心指标 |
| P1 | 对可评估结构计算整体 backbone RMSD，并与 pocket cosine 联合绘制散点图 | 识别“口袋相似但骨架不同”的候选 |
| P1 | 将 sequence / Foldseek / Folddisco 来源统一为 candidate-level 表 | 方便后续人工筛选、聚类和去重 |
| P2 | 对 top candidates 进行 90% 或 95% 序列去重后再看 scaffold 多样性 | 避免同一家族冗余占据优先级 |

## 9. 推荐状态标签

| 标签 | 含义 |
|---|---|
| TOOL_HIT | 工具原生检索命中，例如 Foldseek/Folddisco/BLAST/MMseqs 返回 hit |
| STRUCTURE_AVAILABLE | 已有可用于 RMSD/Pythia 的结构 |
| SEQUENCE_AVAILABLE | 已有可用于全局序列比对的序列 |
| PREDICTED_POCKET_OK | Pythia-pocket 在阈值下得到 pocket embedding |
| NO_PREDICTED_POCKET | Pythia-pocket 阈值下无 pocket residue，不等于生物学失败 |
| NOT_EVALUATED | 缺结构、缺序列、下载/解析失败或尚未进入队列 |
| PROJECT_PRIORITY | 后续根据多指标综合排序得到的候选优先级，不是工具原生标签 |
