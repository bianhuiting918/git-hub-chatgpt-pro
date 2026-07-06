# PETase 与尼龙酶骨架搜索筛选方法及初步统计

生成日期：2026-07-06  
更新日期：2026-07-06  
项目服务器路径：`/Dell/Dell14/bianht/enzyme_scaffold_search_v2`  
适用范围：PETase 与尼龙酶（NylC / Nyl50 / Nyl10）骨架搜索、口袋分层统计、Pythia-pocket 评分、cluster 筛选、二聚体/四聚体界面口袋筛选前的交接报告。

## 1. 研究目标

本项目的目标不是只找到与已知酶序列相似的同源蛋白，而是尽可能系统地召回可能具有相似催化口袋或可迁移功能环境的新骨架。当前策略以“高召回、分层描述、后续再排序”为原则：

1. 用序列搜索扩大同源和远缘同源候选。
2. 用 Foldseek 从整体结构相似性角度召回已知或预测结构。
3. 用 Folddisco 从活性残基和周围口袋 motif 的结构几何角度召回候选。
4. 将 sequence / Foldseek / Folddisco 的候选合并去重，避免过早用单一方法排除潜在新骨架。
5. 对合并后的候选计算全局序列相似性、整体结构 RMSD、Pythia-pocket embedding cosine similarity，再进行多指标排序和人工审查。
6. 对尼龙酶候选增加 oligomer-aware 判断：二聚体/四聚体倾向、界面口袋、跨链 substrate tunnel 与 PA6/PA66 productive docking。

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
| 寡聚体界面 | ColabFold-Multimer / AlphaFold-Multimer / template assembly | 判断尼龙酶候选是否形成 NylC/Nyl50-like functional interface | 不问“任意二聚体”，而问“正确界面是否形成” |
| 界面口袋/通道 | P2Rank / fpocket / CAVER / POVME / docking | 判断二聚体/四聚体界面是否构成 PA6/PA66 可进入的 catalytic pocket/tunnel | 放在 monomer 初筛之后、实验候选之前 |

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

这里的“通过 layer”应理解为“该候选在对应 layer 下有结构映射/几何支持”，而不是功能已确认。后续 Pythia-pocket cosine、整体 RMSD、全局序列差异、cluster 与 oligomer-aware pocket geometry 会一起用于排序。

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

### 6.1 Nylonase-specific pocket cosine 分层

GRASE 在 Aes72-centered urethanase discovery 中使用的是 key catalytic residue embedding cosine 1.00-0.97 四个 tier。但 nylonase 已知有效 seed 之间本身更分散：NylC vs Nyl50 为 0.925020，NylC vs Nyl10 为 0.860980。因此不能把 0.97-1.00 直接作为 nylonase 硬阈值，否则会把 Nyl50-like 和 Nyl10-like 已知有效方向排除。

尼龙酶建议使用 multi-seed max cosine，而不是只用 NylC seed：

```text
max_pocket_cosine = max(
  cosine_to_NylC,
  cosine_to_Nyl50,
  cosine_to_Nyl10
)
```

建议分层：

| Tier | 阈值 | 解释 | 处理建议 |
|---|---:|---|---|
| Nyl-pCOS-high | `max_pocket_cosine >= 0.93` | NylC-like / Nyl50-like 近口袋候选 | 高优先级，但仍需 cluster 去冗余和二聚体界面检查 |
| Nyl-pCOS-mid | `0.86 <= max_pocket_cosine < 0.93` | Nyl10-like / Nyl50-like 远缘有效候选区间 | 必须保留；不能按 GRASE 阈值丢弃 |
| Nyl-pCOS-low-but-salvageable | `0.80 <= max_pocket_cosine < 0.86` | 探索区间 | 需要 strong motif、Folddisco、二聚体界面或 CAVER tunnel 支持 |
| Nyl-pCOS-exploratory | `< 0.80` | 低相似探索区间 | 一般不进首轮实验，除非 cluster novelty 或几何证据很强 |

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
5. 后续应统一输出 candidate-level 主表：candidate_id、enzyme_class、route、best_seed、sequence_available、structure_available、global_PID、global_RMSD、Pythia-pocket cosine、layer support、cluster ID、oligomer/interface status、status label。

## 7A. 文献筛选策略对本项目的约束

### 7A.1 GRASE 查找聚氨酯酶：高 embedding 阈值适合 Aes72-centered urethanase，不宜直接套到 nylonase

GRASE 的 template 是 Aes72。作者先实验测试 14 个文献或专利来源的 PURase/urethane hydrolase，只有 Aes72、GatA、SP1 对 N-aryl carbamates 有明显活性，其中 Aes72 对 2,4-TDA-DEG 的活性最高，因此被选为后续搜索模板。

GRASE 前先做了两个 baseline：

1. **传统序列/结构搜索。** 用 Aes72 序列通过 BLASTp 搜 NCBI nonredundant protein database；用 Aes72 结构通过 Foldseek 搜 AlphaFold Protein Structure Database。按 sequence identity 或 structural identity 选 top-ranked candidates，最终 B1-B5 和 F1-F5 上实验。只有 4 个能产生 TDA monomer，F4/F5 不能产生 TAC 或 TDA，说明全局序列/结构相似性不足以判断 urethanase 活性。
2. **cluster + complex prediction。** 从 126 个 clusters 中采样 216 条 representative sequences，预测每条与 2,4-TDA-DEG 的 complex，按 Rosetta energy 选 12 个 CP candidates 上实验。只有 2 个产生 TDA monomer，说明 docking/complex energy 可以辅助，但不能单独作为主筛选指标。

最终有效路线是 GRASE：Pythia-Pocket 对 Aes72 和候选结构的 key catalytic residues 生成 residue-level embeddings，计算 cosine similarity；Pythia 同时用 NLL 估计稳定性。作者把候选与 Aes72 的 Pythia-Pocket cosine similarity 分成 **1.00-0.97 的四个 tier**，每个 tier 内选 Pythia-predicted stability 最高的 3 条，最终得到 GB1-GB12 和 GF1-GF12 共 24 个 GRASE candidates。24 个候选中 21 个表达且有 urethanase activity，17 个能产生 TDA monomers，8 个优于 Aes72，GB6/AbPURase 最强。

对本项目的关键判断：**GRASE 的 0.97-1.00 threshold 不能直接套到 nylonase。** 当前 nylonase 统计中 NylC-WT vs Nyl10 的 predicted-pocket cosine 为 0.860980，NylC-WT vs Nyl50 为 0.925020。如果使用 `cosine >= 0.97` 作为硬阈值，Nyl10 和 Nyl50 这两个已知有效方向都会被排除。nylonase 应使用 multi-seed tier，而不是单 seed 高阈值过滤。

### 7A.2 Metal-coordination mining：机制几何 motif 先定义功能类型，cluster 和底物假说决定实验序列

该 Nature 工作研究的是 **Fe(II)/αKG radical halogenase**，更准确说是“卤化酶发现”，不是一般意义的卤代烃脱卤酶。其核心不是 embedding，而是机制驱动的金属配位几何筛选。

作者先从 cupin / Fe(II)-αKG structural space 出发，整理得到约 530,814 个 AF2 structural models。搜索 motif 是：

```text
2His metal-binding site
+
缺失 nearby Asp/Glu 酸性第三配体
+
Gly/Ala 小残基形成开放配位位点
```

原因是 radical halogenase 需要开放金属配位位点让 Cl- / Br- 等卤离子结合；普通 Fe/αKG hydroxylase 通常是 2His-1Asp/Glu facial triad。几何筛选后得到 putative radical halogenase candidates，再通过 UniProt annotation filtering、30% identity SSN、人工去除明显假阳性 cluster，得到 curated cluster atlas。

关键实验序列不是“每个 cluster 都挑一个”，而是：

```text
几何 motif 大规模筛选
↓
30% SSN 识别 family / cluster
↓
选择 genomic context 最有解释力的 Cluster X
↓
对 Cluster X 做 BLAST 扩展与 >50% identity SSN
↓
结合 genome neighborhood / PaperBLAST 推断底物
↓
从有明确底物假说的 subcluster 中选 AspX 和 BtnX
```

AspX 选择逻辑：Cluster X 的某个 subcluster 邻近 amino acid transporter、ATP-grasp / amino-acid-modifying genes，因此推测底物可能是 free amino acids；作者筛 20 种 canonical L-amino acids，发现 AspX 氯化 L-aspartic acid。

BtnX 选择逻辑：另一个 subfamily 与 CurA-like ACP-dependent halogenase 有远缘关系，但缺少典型 ACP-dependent gene context；PaperBLAST / genome context 指向 ECF biotin transporter 和 killer plasmid context，因此推测底物是 biotin，最终验证为 2R-chlorobiotin 生成。

对本项目的可迁移逻辑：

```text
motif / pocket / geometry
负责判断“可能有目标反应能力”

30% cluster
负责判断“属于哪个 family / scaffold group”

50% cluster
负责判断“哪个 subfamily 值得选代表”

annotation / closest characterized homolog / genome context
负责判断“实验应该测 PA6、PA66 还是更具体的 oligomer/product panel”
```

### 7A.3 NylC / Nyl10 / Nyl50 文献：不是 GRASE 式高 embedding 筛选，而是 diversity panel + 实验读出

Nyl10 / Nyl50 文献的路线和 GRASE 完全不同。作者没有先用 pocket embedding 或 docking 筛选，而是以 NylC-GYAQ 氨基酸序列作为 query，在 ColabFold 的 MMseqs2-based sequence search 中搜索 UniRef100，得到 2,839 条 homologous sequences；随后用 HHfilter 要求与 NylC-GYAQ 的 coverage >= 50%，留下 2,643 条。

为了构建 96-well plate，作者从 2,643 条中选择 95 条 maximally diverse homologs，核心标准是 lowest maximum pairwise sequence identity，而不是最像 NylC。6 条有 fragment、不以 Met 起始或含 unknown amino acids 的序列用 UniProtKB 中最近 BLASTp match 替换，最终 panel 为 89 条 UniRef100 + 6 条 UniProtKB，两两 identity 在 30%-50%。随后全部合成表达，用 PA6 / PA66 粗筛直接读出活性和产物谱，最终识别 Nyl10、Nyl12、Nyl50 等代表。

这说明 NylC/Nyl10/Nyl50 的策略不是：

```text
选 pocket embedding 最像 NylC 的序列
```

而是：

```text
在 NylC homolog space 中最大化序列多样性
↓
用实验直接读出 PA6/PA66 活性与产物谱
```

对本项目的含义：当前 Nyl10 与 NylC 的 BLAST pident 只有约 30-31%，Nyl50 与 Nyl10 约 40.6%，而 NylC vs Nyl10 pocket cosine 也只有 0.861。因此 Nyl10-like scaffold 很容易被高 identity 或高 pocket cosine 筛法排除。必须保留 30%-50% identity 的远缘 cluster 代表。

### 7A.4 Nylonase 需要引入二聚体/四聚体-aware 筛选

NylC/Nyl50-like nylonase 与 PETase 不同。PETase 通常可按 monomer active-site pocket 处理；NylC/Nyl50-like 系统可能依赖 oligomeric interface 形成完整底物通道或活性口袋。

Nyl50 结构分析显示，Nyl50 的 putative substrate access tunnel 是跨两个 monomers 的 U-shaped tunnel；入口半径约 2.5 Å，active-site 附近 bottleneck radius 约 1.5 Å，tunnel volume 约 631 Å³，并且一个 large subpocket 位于 monomer interface，含有来自两个亚基的残基。文献还指出 NylC-GYAQ 和 Nyl50 的 putative active sites 跨越 A/D dimer interface，提示 oligomerization 对活性很可能重要。

因此，本项目后续不应只做：

```text
monomer Pythia-pocket embedding
monomer CAVER
monomer docking
```

这些可以作为高召回初筛，但不能作为最终实验排序依据。nylonase 最终应新增：

```text
oligomer-aware pocket / tunnel / interface screening
```

核心判断不是“是否形成任意二聚体”，而是：

```text
是否形成 NylC/Nyl50-like 正确二聚体/四聚体界面；
该界面是否参与 catalytic pocket / substrate access tunnel。
```

## 8. 下一步筛选决策

### 8.1 必须进行 cluster，建议使用 90%、50%、30% 三层 identity

当前 nylonase 候选规模已经较大：broad recall 33,675 unique candidates，candidate FASTA 15,083 sequences，structure merged candidates 21,211，predicted all-seed embedding 3,106 unique candidates。这个规模下如果不 cluster，近重复序列会占据 top rank，实验候选会缺乏 family/scaffold 多样性。

建议三层 cluster：

| Identity | 用途 | 解释 |
|---:|---|---|
| 90% | 工程去冗余 | 删除近重复、strain-level redundancy、几乎相同序列 |
| 50% | subfamily 代表选择 | 对应 NylC 文献 diversity panel 的 30%-50% identity 区间，适合挑实验候选 |
| 30% | family / scaffold 级别 SSN | 对应 metal-coordination mining 的 family-level cluster 逻辑，用于识别远缘 scaffold 和 novel cluster |

可以加 70% identity 作为过渡层，但第一版建议先做：

```text
cluster90 + cluster50 + cluster30
```

### 8.2 不建议直接使用 GRASE 的 0.97-1.00 embedding 阈值

GRASE 的 0.97-1.00 是 Aes72-centered urethanase discovery 的经验阈值。nylonase 已知有效 seed 更分散，因此应使用：

```text
max cosine to NylC / Nyl50 / Nyl10
```

而不是只看 NylC seed。

建议 candidate-level 主表新增：

```text
cosine_to_NylC
cosine_to_Nyl50
cosine_to_Nyl10
max_pocket_cosine
best_pocket_seed
pocket_cosine_tier
```

### 8.3 引入二聚体/四聚体预测，但不要全量跑 AlphaFold-Multimer

建议漏斗：

```text
Step 1
sequence QC:
去 fragment、去 X、去长度异常、保留 Ntn-hydrolase / Nyl-like motif

Step 2
MMseqs2 cluster:
90%, 50%, 30%

Step 3
每个 50% cluster 选 1-3 条代表：
优先 Pythia-pocket OK
优先 layer support 强
优先 structure pLDDT 高
优先多 route 支持 sequence + Foldseek + Folddisco

Step 4
对 top 500-1000 representative candidates 做 homodimer prediction
必要时做 homotetramer prediction

Step 5
计算 oligomer confidence:
ipTM
interface PAE
pDockQ / pDockQ2
buried SASA
interface contact number
multi-seed consistency
A/D or A/B contact-map similarity to Nyl50/NylC

Step 6
只对 oligomer confidence 高，或 confidence 中等但 motif/tunnel 很强的候选，做 interface-pocket analysis
```

### 8.4 口袋几何匹配应从 monomer pocket 升级到 interface-pocket geometry

建议新增一个 `interface_pocket_geometry_score`，并与 Pythia-pocket cosine 并列，而不是替代它。建议指标：

```text
1. catalytic Thr / Ntn-hydrolase autocleavage motif 是否完整
2. catalytic residue RMSD to NylC/Nyl50/Nyl10
3. interface pocket 是否包含来自两个 chain 的 residues
4. candidate 是否有跨链 substrate tunnel
5. CAVER bottleneck radius
6. CAVER tunnel volume
7. tunnel-lining residues composition
8. PA6 / PA66 oligomer fragment docking 是否形成 productive pose
9. PA66-like / PA6-like substrate channel 是否有差异
```

Nyl50 可作为 tunnel 参考：

```text
Nyl50 reference tunnel:
entrance radius ≈ 2.5 Å
bottleneck radius ≈ 1.5 Å
tunnel volume ≈ 631 Å³
tunnel spans both monomers
one large subpocket is at monomer interface
```

### 8.5 首轮实验候选应按 group 分层抽样，而不是按单一总分 top N

如果首轮实验规模为 24-48 条，建议：

| Group | 选择逻辑 | 建议数量 |
|---|---|---:|
| A. NylC-like 高把握候选 | max cosine high，motif 完整，二聚体/四聚体界面合理 | 6-10 |
| B. Nyl50-like / PA66-selective 候选 | cosine_to_Nyl50 高，interface tunnel 接近 Nyl50，PA66 docking pose 合理 | 6-10 |
| C. Nyl10-like 远缘候选 | cosine 不一定高，但 sequence / structure / motif 指向 Nyl10-like | 6-10 |
| D. low-cosine but strong-geometry candidates | cosine 0.80-0.86，但 Folddisco motif、CAVER tunnel、dimer interface 强 | 4-8 |
| E. cluster novelty representatives | 30% cluster 中远离 NylC/Nyl50/Nyl10，但 motif 和结构可信 | 4-8 |

## 9. Candidate-level 主表建议字段

后续统一主表不应只包括 route、RMSD 和 pocket cosine。建议新增 cluster、multi-seed pocket、oligomer 与 interface pocket 字段：

```text
candidate_id
enzyme_class
route
sequence_available
structure_available
cluster90_id
cluster50_id
cluster30_id
best_seed_by_sequence
best_seed_by_structure
best_seed_by_pocket_cosine
global_PID
global_RMSD
cosine_to_NylC
cosine_to_Nyl50
cosine_to_Nyl10
max_pocket_cosine
best_pocket_seed
pocket_cosine_tier
L0_support
L0L1_support
L0L2_support
L0L1L2_support
oligomer_model_status
predicted_oligomer_state
dimer_ipTM
dimer_interface_PAE
dimer_pDockQ2
interface_BSA
interface_contact_count
interface_contactmap_similarity_to_Nyl50
CAVER_tunnel_found
CAVER_bottleneck_radius
CAVER_tunnel_volume
CAVER_cross_chain
interface_pocket_residue_count
dock_PA6_productive_pose
dock_PA66_productive_pose
dock_best_substrate
dock_reactive_distance
dock_attack_angle
final_priority_group
status_label
```

## 10. 推荐状态标签

| 标签 | 含义 |
|---|---|
| TOOL_HIT | 工具原生检索命中，例如 Foldseek/Folddisco/BLAST/MMseqs 返回 hit |
| STRUCTURE_AVAILABLE | 已有可用于 RMSD/Pythia 的结构 |
| SEQUENCE_AVAILABLE | 已有可用于全局序列比对的序列 |
| PREDICTED_POCKET_OK | Pythia-pocket 在阈值下得到 pocket embedding |
| NO_PREDICTED_POCKET | Pythia-pocket 阈值下无 pocket residue，不等于生物学失败 |
| NOT_EVALUATED | 缺结构、缺序列、下载/解析失败或尚未进入队列 |
| CLUSTER_REPRESENTATIVE | 该候选是 90%/50%/30% cluster 的代表序列之一 |
| OLIGOMER_MODEL_OK | 二聚体/四聚体模型置信度可接受 |
| INTERFACE_POCKET_OK | oligomer interface 上存在与 catalytic site/tunnel 相关的 pocket |
| CROSS_CHAIN_TUNNEL_OK | CAVER/MOLE 显示存在跨链 substrate access tunnel |
| PRODUCTIVE_DOCKING_POSE | PA6/PA66 fragment docking 形成可反应构象 |
| PROJECT_PRIORITY | 后续根据多指标综合排序得到的候选优先级，不是工具原生标签 |

## 11. 当前结论

1. **GRASE 的 0.97-1.00 pocket cosine 阈值不能直接用于 nylonase。** NylC vs Nyl50 为 0.925020，NylC vs Nyl10 为 0.860980；直接套用会漏掉已知有效 scaffold。
2. **Nylonase 应使用 multi-seed pocket cosine。** 以 NylC、Nyl50、Nyl10 三类 seed 的 max cosine 分层，保留 0.86-0.93 的远缘有效候选区间。
3. **必须 cluster。** 建议同时做 90%、50%、30% identity：90% 去冗余，50% 选 subfamily 代表，30% 做 family/scaffold SSN。
4. **monomer-only pocket 只能作为初筛。** NylC/Nyl50-like enzyme 很可能依赖 A/D dimer interface 或四聚体中的跨链 tunnel，最终排序必须引入 oligomer-aware pocket / tunnel / docking。
5. **实验候选不应按单一总分 top N。** 应按 NylC-like、Nyl50-like、Nyl10-like、low-cosine strong-geometry、cluster novelty 分组抽样。
