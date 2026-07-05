# 塑料降解酶与尼龙酶研究课题框架

_Last updated: 2026-07-05_

本文件整理两个互相关联但科学问题不同的课题方向。

- **课题一：通用塑料降解酶新骨架发现与 QM/MM 能量代理模型**
- **课题二：NylC/Nyl10/Nyl50-like 尼龙水解酶对 PA6/PA66 的界面识别、链段进入、催化机制和产物选择性**

两个课题共享部分方法，例如 docking、MD、enhanced sampling、QM/MM、reaction-state ensemble、MPNN/LigandMPNN 和轻量级预测模型。但它们的定位不同：

```text
课题一：方法学平台
从已知催化机制出发，发现“相似催化口袋、不同整体骨架”的新塑料降解酶，并建立 QM/MM-supervised reaction-state energy predictor。

课题二：具体体系机制
以 NylC/Nyl10/Nyl50-like 尼龙水解酶为对象，解析 PA6/PA66 表面吸附、链段进入、酰胺键水解和产物选择性机制。
```

---

## 课题一：反应状态引导的塑料降解酶新骨架发现与序列优化

### 1. 核心科学问题

塑料降解酶研究已经较多，单纯在已有 PETase、cutinase、PUase 或其他已知骨架上做点突变，创新性有限。本课题的核心创新点是：

```text
不是继续优化一个已知骨架，
而是以已知塑料降解酶的催化口袋和 QM/MM 反应状态为约束，
寻找或生成整体骨架不同、但能够支持相同催化机制的新 scaffold。
```

也就是：

```text
known catalytic mechanism
        ↓
minimal catalytic motif / theozyme
        ↓
reaction-state ensemble: ES / NAC / TS1 / TI / TS2 / PS
        ↓
scaffold search or scaffold generation
        ↓
new fold / new backbone supporting similar catalytic geometry
        ↓
experimental validation of weak activity
        ↓
fixed-backbone sequence design and QM/MM-supervised energy model
```

### 2. 关键概念边界

#### 2.1 固定的不是整个 pocket，而是最小催化 motif

不建议把整个口袋完全刚性固定。更合理的是固定：

```text
1. 催化残基的关键原子几何；
2. 亲核攻击距离和角度；
3. oxyanion hole 或等效过渡态稳定几何；
4. general acid/base 或关键水分子位置；
5. 底物 / 过渡态类似物的关键反应坐标；
6. 必要的 access channel 方向。
```

这个最小口袋可以称为：

```text
minimal catalytic motif
reaction-state motif
theozyme-like motif
```

#### 2.2 MPNN 不是新骨架发现工具

需要区分两个阶段：

```text
新骨架发现：
structure search / motif scaffolding / RFdiffusion-like backbone generation

骨架内序列优化：
ProteinMPNN / LigandMPNN / fixed-backbone sequence design
```

ProteinMPNN 或 LigandMPNN 主要解决“给定 backbone，生成更合适的序列”的问题；它们本身不负责发现全新的 global scaffold。

### 3. 研究路线

#### Step 1：复现并扩展已有 QM/MM 反应路径

以已有文献报道的塑料降解酶为起点，复现其 QM/MM 反应路径，得到多个反应状态：

```text
ES：enzyme-substrate / Michaelis complex
NAC：near-attack conformation
TS1：第一过渡态
TI：tetrahedral intermediate 或相应中间体
TS2：第二过渡态
PS：product state / product complex
```

重点不是得到单一 TS 结构，而是得到每个状态的 **reaction-state ensemble**：

```text
ES ensemble = {ES_1, ES_2, ..., ES_n}
TS1-like ensemble = {TS1_1, TS1_2, ..., TS1_m}
TS2-like ensemble = {TS2_1, TS2_2, ..., TS2_p}
```

这些 ensemble 可以来自：

```text
QM/MM umbrella sampling
QM/MM metadynamics
reaction-coordinate constrained sampling
multiple minimum-energy paths
transition-state-like windows
NAC-filtered MD snapshots
```

#### Step 2：基于 reaction-state motif 搜索或生成新骨架

输入：

```text
minimal catalytic motif
TS-like ligand / intermediate-like ligand
critical catalytic geometry constraints
```

候选来源：

```text
1. PDB / AlphaFold DB structure search；
2. Foldseek / DALI / TM-align 几何搜索；
3. motif scaffolding / RFdiffusion-like backbone generation；
4. 结构数据库中未注释或远缘功能蛋白。
```

筛选目标：

```text
1. 能承载相同催化 motif；
2. 整体 fold / topology 与已知塑料降解酶不同；
3. 口袋具有进入底物或聚合物链段的通道；
4. 不严重破坏 protein foldability；
5. catalytic geometry 可在局部优化后保持。
```

#### Step 3：多状态初步打分，排除明显不合理骨架

对候选新骨架进行初筛，不直接追求精确 kcat，而是进行 mechanism-aware triage。

建议评分项：

```text
A. Foldability
- AlphaFold / ESMFold pLDDT
- PAE
- designed backbone recovery
- global collapse risk

B. Catalytic motif retention
- catalytic atom RMSD
- nucleophile-carbonyl distance
- attack angle
- general acid/base geometry

C. ES compatibility
- substrate binding geometry
- pocket clash
- access channel openness

D. NAC population proxy
- near-attack distance and angle
- substrate strain

E. TS-like stabilization
- oxyanion hole H-bonds
- electrostatic stabilization
- key proton-transfer geometry

F. Intermediate compatibility
- tetrahedral geometry compatibility
- acyl-enzyme or covalent intermediate accommodation

G. Product release
- product overbinding penalty
- blocked exit channel penalty
```

核心原则：

```text
不要只选择 substrate ground-state binding 最强的骨架。
真正有催化意义的是 selective transition-state stabilization：
TS-like state 相对于 ES 得到更好的稳定。
```

#### Step 4：实验验证新骨架是否具有弱催化功能

建议第一轮实验验证 48–96 个候选骨架。

验证层次：

```text
一级：可溶模型底物 / 低聚物底物
- 快速判断是否有基本催化能力。

二级：真实塑料低聚物 / film / powder / LC-MS 产物
- 判断是否具有实际塑料降解相关性。
```

第一轮实验目标不是得到最高活性，而是确认：

```text
该新骨架确实具有可检测的催化功能。
```

一旦有弱活性骨架，它就是后续序列优化和能量模型训练的有效起点。

### 4. 骨架内序列生成与 QM/MM-supervised 能量模型

#### 4.1 序列生成对象

在已经实验验证具有弱催化功能的新骨架上，进行固定骨架序列生成：

```text
ProteinMPNN：固定 backbone 生成稳定序列；
LigandMPNN：显式考虑底物、TS-like ligand、中间体或关键水分子。
```

建议生成规模：

```text
5 个 validated scaffolds × 每个 scaffold 200 条序列 = 1000 条设计序列
```

或者：

```text
10 个 scaffold × 每个 scaffold 100 条序列 = 1000 条设计序列
```

序列生成策略：

```text
A. pocket-conservative design
保留一层催化口袋，只改二层和远端背景。

B. pocket-redesign design
允许口袋周围 6–8 Å 内残基改变。

C. ligand-aware design
使用底物或 TS-like ligand 作为设计上下文。

D. transition-state-biased design
以 TS-like state 为设计对象，奖励过渡态稳定化。

E. product-release-aware design
避免产物被过度锁死。
```

#### 4.2 1000 条序列的 QM/MM 标签不是实验活性标签

这里的模型不是直接预测实验活性，而是预测 QM/MM 衍生的能量标签。

准确表述：

```text
QM/MM-supervised reaction-state ensemble energy predictor
```

中文：

```text
基于 QM/MM 标签的反应状态构象集合能量预测模型
```

模型预测的是：

```text
ES、NAC、TS1-like、TI、TS2-like、PS 等状态在不同序列背景下的相对稳定性；
以及由这些状态得到的近似反应势垒、过渡态稳定化和产物困陷倾向。
```

不是直接预测：

```text
kcat
kcat/KM
真实实验转化率
```

#### 4.3 不能使用绝对总能量作为标签

不同序列突变后原子组成、电荷、侧链数量不同，因此绝对总能量不能横向比较。

应该使用相对能量标签：

```text
ΔE‡TS1 = E_TS1-like - E_ES
ΔE_TI = E_TI - E_ES
ΔE‡TS2 = E_TS2-like - E_TI
或 ΔE‡TS2,total = E_TS2-like - E_ES
ΔE_PS = E_PS - E_ES

ΔΔE‡TS1 = ΔE‡TS1(mutant) - ΔE‡TS1(WT)
ΔΔE‡TS2 = ΔE‡TS2(mutant) - ΔE‡TS2(WT)
```

核心标签：

```text
Selective TS stabilization = [G_TS-like - G_ES]_mutant - [G_TS-like - G_ES]_WT
```

#### 4.4 使用 ensemble-based scoring，而不是单结构打分

不把 WT 的单个 ES / TS1 / TS2 结构硬塞进突变体中，而是使用 WT 的反应状态构象集合。

对每条设计序列 `s`，每个状态 `k`，每个构象 `i`：

```text
E(s, k, i)
```

然后用 Boltzmann-like ensemble averaging 得到状态自由能近似：

```text
G(s, k) = -RT ln Σ_i exp[-β E(s, k, i)]
```

如果每个 WT ensemble 构象有采样权重，可以写作：

```text
G(s, k) = -RT ln Σ_i w_i exp[-β E(s, k, i)]
```

这样可以评价：

```text
该序列是否只偶然稳定一个 TS-like pose，
还是稳定整个 TS-like ensemble 中的一批可反应构象。
```

#### 4.5 多保真度标签策略

完整 QM/MM TS 搜索对 1000 条序列成本极高，因此建议采用多保真度策略：

```text
Tier 1：1000 条
TS-like constrained QM/MM single-point 或 constrained optimization labels

Tier 2：100–200 条
局部 QM/MM constrained optimization / reaction-coordinate scan

Tier 3：10–30 条
true TS search、虚频验证、IRC/MEP 或 QM/MM umbrella/metadynamics

Tier 4：top candidates
实验验证
```

命名要严谨：

```text
大规模标签：TS-like ensemble score
高精度标签：true TS energy / free-energy barrier
```

#### 4.6 模型输入特征

建议三类输入：

```text
A. 序列 embedding
- ESM / ProtT5 / MSA-derived per-residue embedding
- pocket residue embedding
- second-shell residue embedding
- global sequence embedding

B. 结构图特征
- pocket residue-ligand graph
- residue-residue contacts
- residue-TS atom distances
- hydrogen-bond geometry
- electrostatic features
- pocket volume and shape

C. 反应状态几何特征
- nucleophile-carbonyl distance
- attack angle
- oxyanion H-bond distances
- proton-transfer distances
- catalytic water position
- electrostatic field along forming/breaking bond
- ligand strain
- product exit geometry
```

#### 4.7 模型类型

数据规模约 1000 条，不适合从零训练大模型。更合理的是：

```text
frozen protein language model embeddings
+
reaction-state geometric descriptors
+
small multi-task regressor / ranking model
```

候选模型：

```text
XGBoost / LightGBM
random forest
ridge / elastic net
Gaussian process
small MLP on frozen embeddings
small GNN for pocket graph
```

训练目标：

```text
multi-task regression + pairwise ranking
```

输出：

```text
y1 = ΔE‡TS1-like
y2 = ΔE_TI
y3 = ΔE‡TS2-like
y4 = ΔE_PS
y5 = geometry penalty
y6 = product trapping penalty
```

更推荐 ranking，而不是只追求绝对能量回归：

```text
目标不是精确预测 ΔE‡ = 13.2 kcal/mol，
而是预测哪些序列应该进入下一轮 QM/MM 或实验验证。
```

### 5. 课题一的最终产出

```text
1. 一个 reaction-state-guided scaffold discovery pipeline；
2. 一批具有不同整体骨架但相似催化机制的新候选塑料降解酶；
3. 至少 1–5 个经实验验证具有弱催化功能的新骨架；
4. 一个 QM/MM-supervised reaction-state ensemble energy predictor；
5. 一批经模型筛选的更优设计序列；
6. 物理解释：哪些序列背景更能稳定 TS-like ensemble、降低 ΔE‡、避免 product trapping。
```

### 6. 课题一的风险和限制

```text
1. 口袋几何相似不等于真实活性。
需要考虑底物进入、产物释放、口袋柔性和 protein expression。

2. WT reaction-state ensemble 只适用于同机制设计。
如果新骨架产生完全不同 proton-transfer pathway 或不同反应中间体，需要重新做 QM/MM path search。

3. 1000 条计算标签不是实验真实活性标签。
模型应定位为 QM/MM teacher 的 surrogate，而不是实验 kcat predictor。

4. 同一 backbone 上的 1000 条序列不代表跨 scaffold 泛化。
需要区分 intra-scaffold optimization 和 cross-scaffold generalization。

5. 最终仍需要实验验证。
计算和模型只能用于优先级排序，不能替代实际活性测定。
```

---

## 课题二：NylC/Nyl10/Nyl50-like 尼龙水解酶对 PA6/PA66 的界面识别、链段进入和催化机制

### 1. 核心科学问题

尼龙酶研究相较 PET 等塑料降解酶仍不充分。NylC 是较经典的尼龙寡聚体内切水解酶；Nyl10 和 Nyl50 是近年报道的 PA66-selective Ntn hydrolase nylon hydrolases。Nyl50 已有结构和 putative substrate-access tunnel 线索，但其 PA66 选择性、产物偏好和底物进入机制尚未被系统解析。

本课题核心问题：

```text
NylC / Nyl10 / Nyl50-like 酶如何识别 PA6/PA66，
如何吸附在聚合物或预处理寡聚体表面，
如何让 PA 链段进入 substrate-access tunnel，
进入后如何进行 N-terminal Thr 介导的酰胺键水解，
以及这些过程如何决定 PA6/PA66 产物分布。
```

### 2. 与课题一的区别

课题二不急于训练通用模型。尼龙体系当前最缺的是清楚的机制状态定义：

```text
1. 什么是 productive NylC/Nyl50-PA66 binding pose？
2. PA66 chain-end 和 internal-loop 哪个更可能进入 tunnel？
3. Nyl50 的 A/D dimer interface 是否参与 PA66 substrate-access tunnel？
4. 为什么 Nyl50 偏 PA66 L2，而 NylC 更偏 PA66 L1？
5. 化学 TS 是否真是限速，还是 chain-entry / accessibility 才是主限速？
```

因此课题二优先做机制模拟和 QM/MM，而不是直接建立大规模活性预测模型。

### 3. 研究对象

建议酶 panel：

```text
1. NylC / NylC-GYAQ：经典 NylC reference；
2. Nyl10：PA66-selective Ntn hydrolase；
3. Nyl50：PA66-selective Ntn hydrolase，已有结构和 tunnel 假说；
4. Nyl12 或其他已报道 Ntn hydrolase：可选对照；
5. Aim 0 或后续搜索得到的新 NylC-like scaffold。
```

需要强调：

```text
Nyl10 和 Nyl50 是天然/筛选得到的 PA66-selective nylon hydrolase scaffold，
不是 NylC 的点突变体。
```

工程化 NylC 变体，例如 NylC-HP、NylC-V3 等，适合作为 benchmark 或后续工程参考，不作为天然骨架搜索种子。

### 4. 底物状态

尼龙底物不应只用单一模型。建议分三层：

```text
A. 固体聚合物表面
- PA6 crystalline slab
- PA6 amorphous slab
- PA66 crystalline slab
- PA66 amorphous slab
- defect-rich / chain-end-rich slab

B. 化学预处理产物
- acid-oligomerized PA6 / PA66
- solvent-dispersed PA6 / PA66
- low-MW oligomer mixture

C. 定义明确的寡聚体
- PA6: Ahx3, Ahx4, Ahx5
- PA66: L1, L2, C1, C2, 1.5-mer / 2.5-mer if available
```

化学预处理的定位：

```text
尼龙不同于 amyloid fibril，可以先用酸、溶剂、热或机械方法松解材料，
再让酶负责选择性水解和产物分布收敛。
```

因此，尼龙酶工程目标不一定是“直接啃完整高结晶 PA66”，也可以是：

```text
高效处理化学预处理后产生的真实 PA6/PA66 寡聚体混合物。
```

### 5. Aim 0：天然新骨架搜索与初筛

如果将骨架发现作为前置工作，可设为 Aim 0。

```text
Aim 0：以 Nyl10 / Nyl50 和 wild-type NylC-like enzymes 为天然 seed，
搜索新的 PA6/PA66-active Ntn hydrolase scaffold。
```

注意：

```text
点突变体不作为 seed；
Nyl10 / Nyl50 作为天然 PA66-selective seed；
NylC 工程化变体只作为 benchmark。
```

流程：

```text
Nyl10 / Nyl50 / natural NylC-like seed set
        ↓
MMseqs2 / HMMER / HHpred sequence-profile search
        ↓
Foldseek / DALI / AlphaFold structure search
        ↓
按 catalytic Ntn motif、自剪切位点、substrate-access tunnel、loop architecture 筛选
        ↓
表达和自剪切成熟验证
        ↓
PA6/PA66 固体和寡聚体底物初筛
        ↓
LC-MS/HPLC 确认 PA6 Ahx1–Ahx5 和 PA66 L1/L2/C1/C2 产物谱
```

### 6. Aim 1：PA6/PA66 材料和预处理底物可及性

目标：

```text
建立 PA6/PA66 不同材料状态下的可及性图谱，
确定哪些链段最可能成为 NylC/Nyl50 的 productive substrate。
```

需要表征：

```text
DSC / WAXS：结晶度和晶型
GPC / SEC：Mn、Mw、链端密度
SEM / AFM：表面形貌和裂纹
BET / particle size：表面积
FTIR / XPS：表面化学和氢键环境
LC-MS blank wash：预存低聚物背景
```

计算模型：

```text
crystalline slab
amorphous slab
defect-rich slab
surface chain-end segment
surface loop segment
internal crystalline segment
```

核心计算量：

```text
amide SASA
PA-PA interchain H-bond number
water accessibility
local chain mobility
interface depth
ΔG_extract：链段从 PA 表面被拉出的自由能
```

### 7. Aim 2：酶-聚合物界面吸附与 productive / nonproductive binding

目标：

```text
解析 NylC/Nyl10/Nyl50/new scaffold 如何吸附到 PA6/PA66 表面，
并区分 productive 和 nonproductive adsorption。
```

模拟体系：

```text
Nyl50 dimer + PA66 slab
Nyl10 model + PA66 slab
NylC tetramer/dimer model + PA6/PA66 slab
new scaffold + PA6/PA66 slab
```

应从多个初始取向开始，而不是只让 active site 面向聚合物。

分析指标：

```text
enzyme-surface contact map
active-site tunnel orientation
surface residence time
productive binding fraction
nonproductive adsorption classes
polymer chain proximity to tunnel entrance
```

productive pose 标准：

```text
1. scissile amide carbonyl C 靠近 N-terminal Thr Oγ；
2. carbonyl O 指向 oxyanion-stabilizing region；
3. leaving amide N 有合理 proton-transfer path；
4. PA chain 在 tunnel 中稳定；
5. 攻击角度接近可反应几何。
```

### 8. Aim 3：PA6/PA66 链段进入 substrate-access tunnel 的自由能

目标：

```text
定量 PA 链段从材料表面或寡聚物状态进入 NylC/Nyl50 active-site tunnel 的自由能代价。
```

关键自由能：

```text
ΔG_extract：链段从 PA 氢键网络中拉出；
ΔG_entry：链段进入酶 substrate-access tunnel；
ΔG_productive：loosely bound → productive catalytic pose。
```

增强采样方法：

```text
umbrella sampling
well-tempered metadynamics
HREMD
path CV / string method
funnel metadynamics
```

关键比较：

```text
NylC vs Nyl50
PA6 vs PA66
chain-end entry vs internal-loop entry
solid surface segment vs acid-oligomerized substrate
L1-producing pose vs L2-producing pose
```

产物选择性假设：

```text
chain-end entry ↑ → 更偏 L1 / Ahx1-Ahx2；
internal-loop entry ↑ → 更偏 L2 / longer oligomers / true endo cleavage；
tunnel 太窄 → 选择性高但转化率低；
tunnel 太宽 → 总活性可能提高但产物分布变宽；
product release 慢 → product trapping 或产物抑制。
```

### 9. Aim 4：QM/MM 催化机制

目标：

```text
在经过 Aim 2–3 确认的 productive complex 上，
解析 N-terminal Thr 介导的 PA6/PA66 酰胺键水解反应机制。
```

不建议直接拿任意 docked oligomer 做 QM/MM。先要有可信的 productive chain-entry complex。

反应步骤：

```text
1. N-terminal Thr Oγ attack on amide carbonyl；
2. tetrahedral intermediate formation；
3. C-N bond cleavage；
4. acyl-enzyme or equivalent intermediate；
5. water-mediated deacylation；
6. product release。
```

比较体系：

```text
NylC + PA6 Ahx4/Ahx5
Nyl50 + PA66 L2-like substrate
Nyl10 + PA66 substrate
new scaffold + preferred substrate
```

核心问题：

```text
Nyl50 的 PA66 选择性来自化学势垒更低，
还是来自 chain-entry / productive binding 更容易？
```

### 10. Aim 5：界面、通道和口袋改造

根据 Aim 1–4 的结果，设计两类改造。

#### 10.1 界面改造

目标：

```text
提高 productive adsorption，减少 nonproductive binding。
```

可改区域：

```text
enzyme surface patches
A/D dimer interface
polymer-facing loops
hydrophobic / aromatic / polar balance
charged residues affecting PA surface binding
```

#### 10.2 口袋 / 通道改造

目标：

```text
降低 chain-entry barrier；
稳定 productive amide geometry；
控制 PA6 Ahx 或 PA66 L1/L2 产物分布；
避免 product trapping。
```

可改区域：

```text
substrate-access tunnel
active-site entrance loop
oxyanion-stabilizing region
subsites along PA chain
product exit channel
```

设计原则：

```text
不要只让底物 binding 更强；
要让 productive state 更稳定、nonproductive state 更不稳定、product release 不被严重阻碍。
```

### 11. 课题二最终产出

```text
1. NylC/Nyl10/Nyl50/new scaffold 的 PA6/PA66 substrate-selectivity 机制；
2. PA66 chain-end vs internal-loop entry 的自由能图谱；
3. Nyl50 A/D dimer tunnel 是否参与 PA66 chain entry 的机制证据；
4. N-terminal Thr 介导 PA6/PA66 酰胺键水解的 QM/MM 反应路径；
5. PA6 Ahx 和 PA66 L1/L2 产物选择性的结构解释；
6. 可用于界面工程和 pocket/tunnel engineering 的突变策略。
```

---

## 两个课题的关系

### 共同思想：多状态设计

两个课题都不应只看一个静态 docking pose。统一思想是：

```text
substrate encounter
→ productive binding
→ chain / pocket entry
→ TS stabilization
→ intermediate handling
→ product release
```

即：

```text
设计的不只是一个能 bind substrate 的酶，
而是一个能沿着正确反应路径稳定关键状态、避免错误状态的酶。
```

### 区别

```text
课题一：
更偏通用方法学。
重点是 reaction-state-guided scaffold discovery 和 QM/MM-supervised energy surrogate。

课题二：
更偏具体体系机制。
重点是 NylC/Nyl10/Nyl50-like nylon hydrolase 的 PA6/PA66 substrate entry、catalysis 和 product selectivity。
```

### 推荐推进顺序

如果资源有限，建议优先推进课题二的机制研究，因为它的体系更具体、文献空白更明确、故事线更容易收敛。

课题一可以作为方法学平台逐步发展：

```text
先在一个已知塑料降解酶骨架上建立 reaction-state ensemble scoring；
再验证新骨架发现；
最后把 QM/MM-supervised energy surrogate 推广到多骨架、多塑料体系。
```

---

## 最简版摘要

```text
课题一：
建立一个以 QM/MM reaction-state ensemble 为核心的塑料降解酶新骨架发现和序列优化平台。先固定催化 motif 搜索不同骨架，实验验证弱活性；再在命中骨架上生成约 1000 条序列，用 ES/TS1/TI/TS2/PS 等 QM/MM 相对能量标签训练 energy surrogate，用于预测不同序列背景对催化能垒和过渡态稳定化的影响。

课题二：
以 NylC/Nyl10/Nyl50-like 尼龙水解酶为对象，解析 PA6/PA66 表面吸附、链段进入 substrate-access tunnel、N-terminal Thr 催化水解和 L1/L2/Ahx 产物选择性机制。通过 MD、enhanced sampling 和 QM/MM 计算，指导酶表面界面改造和 pocket/tunnel 工程，提高尼龙解聚效率并控制产物分布。
```
