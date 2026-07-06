# 课题一：反应状态引导的塑料降解酶新骨架发现与序列优化

_Last updated: 2026-07-06_

## 0. 课题定位

本课题是一个**方法学平台型课题**。它不以某一个具体塑料酶为唯一研究对象，而是希望建立一条可以推广到 PET、聚酯、聚氨酯等塑料水解酶的设计路线。

核心目标是：

```text
从已知塑料降解酶的催化机制和 QM/MM 反应状态出发，
固定最小催化口袋 / 过渡态稳定几何，
寻找或生成整体骨架不同的新酶，
再通过 QM/MM 反应状态能量标签和小型预测模型进行序列优化。
```

一句话概括：

```text
不是继续在已知 PETase/cutinase 骨架上做点突变，
而是从“催化口袋应该如何稳定 ES/TS/中间体”出发，
寻找不同整体 fold 的新骨架，并用 QM/MM 能量标签指导后续序列设计。
```

本课题最终想建立的是：

```text
reaction-state-guided scaffold discovery pipeline
+
QM/MM-supervised reaction-state ensemble energy predictor
```

中文可以称为：

```text
反应状态引导的新骨架发现流程
+
基于 QM/MM 标签的反应状态构象集合能量预测模型
```

注意：这里的模型不是直接预测实验活性，而是预测：

```text
不同序列背景对底物结合态、近攻击构象、过渡态、中间体和产物态的相对稳定性。
```

---

## 1. 背景介绍

塑料酶降解研究在 PET 方向已经发展较快，尤其是 PETase、cutinase、LCC、MHETase 等体系。过去十年主要进展集中在三个方向：

```text
1. 发现天然 PET 降解体系；
2. 解析已知 PET 降解酶结构和催化机制；
3. 在已知骨架上进行定向进化、热稳定化和机器学习辅助突变。
```

这些研究证明了塑料酶降解的可行性，但也暴露出一个明显限制：

```text
大多数研究仍然围绕少数已知 fold 展开，尤其是 cutinase / αβ-hydrolase 类骨架。
```

因此，如果本课题继续做“已知骨架点突变”，创新性不足。更有新意的方向是：

```text
从已知催化机制中抽象出“口袋和过渡态稳定化规则”，
然后用这些规则寻找或生成不同整体骨架的新酶。
```

---

## 2. 目前研究进展：按时间线整理

### 2.1 2016：Yoshida 等人在 Science 发现 PET 降解细菌和 PETase/MHETase 体系

2016 年，Yoshida、Hiraga、Takehana、Taniguchi、Oda 等人在 **Science** 发表文章：

```text
A bacterium that degrades and assimilates poly(ethylene terephthalate)
```

他们从 PET 污染环境中分离出 **Ideonella sakaiensis 201-F6**，证明该菌可以利用 PET 作为碳源和能量来源。

该体系包含两个关键酶：

```text
PETase：将 PET 切成 MHET 和少量 BHET；
MHETase：将 MHET 继续水解为 TPA 和 EG。
```

这一工作的重要性：

```text
1. 第一次系统证明细菌可以利用 PET；
2. 建立了 PETase + MHETase 的生物降解逻辑；
3. 将 PET 酶研究从环境观察推进到分子酶学。
```

但它没有解决：

```text
1. PETase 活性仍然很低；
2. 高结晶 PET 降解效率差；
3. 没有新骨架发现方法；
4. 没有针对过渡态稳定化的系统设计。
```

---

### 2.2 2018：Austin、McGeehan、Beckham 等解析并工程化 PETase

2018 年，Austin、Allen、Donoghue、McGeehan、Beckham 等人在 **PNAS** 发表：

```text
Characterization and engineering of a plastic-degrading aromatic polyesterase
```

他们解析了 PETase 结构，并比较 PETase 与 cutinase-like 酶的活性口袋差异，提出 PETase 的较开放口袋有利于 PET-like 底物进入。

他们做了结构引导突变，证明 PETase 可以通过少数突变改变 PET/PEF 水解性能。

这项工作的重要性：

```text
1. 给出 PETase 结构基础；
2. 证明活性口袋形状和开放程度影响 PET-like 底物水解；
3. 启动了 PETase 结构工程路线。
```

但它仍然局限于：

```text
1. 已知 PETase 骨架内部突变；
2. 主要关注结构和活性差异；
3. 没有系统探索“相同催化口袋、不同整体骨架”；
4. 没有把 QM/MM 多状态反应能量作为可学习标签。
```

---

### 2.3 2020：Tournier 等人在 Nature 工程化 LCC，实现高效 PET 解聚

2020 年，Tournier、Topham、Gilles、Marty 等人在 **Nature** 发表：

```text
An engineered PET depolymerase to break down and recycle plastic bottles
```

该工作以 **leaf-branch compost cutinase, LCC** 为基础，通过工程化提高热稳定性和 PET 解聚效率。在接近 PET 玻璃化转变温度以上的条件下，工程化 LCC 可以高效 depolymerize PET，并将产物用于再聚合。

这项工作的意义：

```text
1. 证明 PET enzymatic recycling 可以接近应用水平；
2. 强调热稳定性对 PET 解聚很重要；
3. 说明真实 PET 材料降解需要考虑温度、结晶度和材料状态；
4. 推动 LCC/cutinase 成为 PET 工业酶路线的重要骨架。
```

但它的不足是：

```text
1. 仍然是在 LCC/cutinase 骨架上优化；
2. 主要通过实验工程获得高性能变体；
3. 没有提出通用新骨架发现策略；
4. 没有建立 reaction-state ensemble 的设计模型。
```

---

### 2.4 2020：Knott、McGeehan、Beckham 等工程化 PETase-MHETase 双酶体系

2020 年，Knott、Erickson、McGeehan、Beckham 等人在 **PNAS** 发表：

```text
Characterization and engineering of a two-enzyme system for plastics depolymerization
```

他们研究 PETase 和 MHETase 的两酶系统，并构建了 PETase-MHETase 连接体或级联设计，目的是减少 MHET 积累，提高 PET 向 TPA/EG 的转化。

这项工作的重要性：

```text
1. 强调 PET 解聚不只是 PETase 单步反应；
2. PETase 释放 MHET 后，需要 MHETase 继续处理；
3. 产物清理和级联反应会影响总体效率；
4. 说明 product release 和 downstream enzyme coupling 是工程重点。
```

但它没有解决：

```text
1. 如何发现新的 PETase-like scaffold；
2. 如何用 QM/MM 反应状态指导新骨架设计；
3. 如何预测不同序列背景对 TS 稳定化和产物释放的影响。
```

---

### 2.5 2022：Lu、Alper 等在 Nature 报道 FAST-PETase

2022 年，Lu、Diaz、Alper 等人在 **Nature** 发表：

```text
Machine learning-aided engineering of hydrolases for PET depolymerization
```

该工作通过机器学习辅助突变设计，获得 **FAST-PETase**，使 PETase 在较温和条件下对多种 PET 底物表现出更高解聚能力。

这项工作的重要性：

```text
1. 证明机器学习可以辅助 PETase 工程；
2. 将突变设计和实验验证结合；
3. 进一步推动 PETase 从基础酶学走向工程应用；
4. 表明少数关键突变可以显著改变稳定性和活性。
```

但对本课题来说，它仍有局限：

```text
1. 它主要还是在 PETase 已知骨架上做突变；
2. 模型主要围绕突变体性能，不是新 scaffold discovery；
3. 没有使用 QM/MM 多反应状态能量作为训练标签；
4. 没有以过渡态 ensemble 为中心训练能量代理模型。
```

---

### 2.6 2022–2024：HotPETase、LCC-ICCG、PET 表面吸附和链段进入模拟

2022 年以后，PET 降解酶研究出现两个趋势：

```text
1. 更稳定、更高温活性的 PET hydrolase 变体；
2. 将 PET 材料表面和酶放在一起做 MD，研究酶如何吸附和形成 productive binding。
```

代表性方向包括：

```text
HotPETase：通过实验或结构引导工程获得更高热稳定性 PETase；
LCC-ICCG：工程化 LCC/cutinase，用于 PET 表面吸附和降解研究；
Sahihi 等 2024 JCIM：LCC-ICCG 在 PET surface 上的全原子 MD 吸附模拟；
Jäckering 等 2024 JACS Au：从 PET bulk 到 hydrolase binding pocket 的链段进入机制模拟。
```

这些研究的重要启发：

```text
1. 酶不是只要有催化中心就能降解塑料；
2. 表面吸附、active-site 朝向、chain entry、productive vs nonproductive binding 都会影响反应；
3. 单个 docking pose 不足以解释真实材料降解；
4. 需要从反应路径和构象 ensemble 角度理解酶-塑料体系。
```

但是，这些研究仍然主要集中于已知 PETase/LCC/PES-H1 骨架，没有系统把：

```text
reaction-state ensemble
+
new scaffold discovery
+
QM/MM energy surrogate
```

结合起来。

---

### 2.7 2023–2025：蛋白生成模型和 fixed-backbone 设计工具成熟

这几年蛋白设计工具快速发展：

```text
RFdiffusion：可在固定 motif 的情况下生成新 backbone；
ProteinMPNN：固定 backbone 生成折叠兼容序列；
LigandMPNN：在显式 ligand / 小分子环境下做序列设计；
AlphaFold / ESMFold：用于验证设计序列是否能回到目标结构；
蛋白语言模型 embedding：用于提取口袋残基和序列背景特征。
```

这些工具让本课题成为可能：

```text
先用 reaction-state motif 寻找或生成新骨架，
再用 MPNN/LigandMPNN 对命中骨架做序列优化，
最后用 QM/MM 能量标签训练小型代理模型进行排序。
```

---

## 3. 现有研究没有做什么

结合上面的时间线，可以明确看到目前缺口：

### 缺口 1：新骨架发现仍然不足

已有 PETase/LCC 工程主要集中在少数骨架内部突变。即使活性很高，仍然不是系统寻找全新 fold。

缺少：

```text
以催化口袋和过渡态稳定几何为核心约束的新 scaffold discovery。
```

### 缺口 2：多数筛选只看单状态或单指标

已有工作常用：

```text
模型底物活性；
PET film 产物总量；
单一 docking pose；
结构稳定性指标。
```

缺少：

```text
ES / NAC / TS / intermediate / product 多状态综合评价。
```

### 缺口 3：QM/MM 结果通常只解释已知酶，而不是作为训练标签

很多 QM/MM 研究用于解释某个已知酶如何催化，但没有进一步把这些结果变成可训练的数据集。

缺少：

```text
用 QM/MM 对大量设计序列生成反应能垒标签，
训练 energy surrogate 来筛更多序列。
```

### 缺口 4：模型底物筛选和真实塑料降解之间存在断层

pNP ester 等模型底物筛选通量高，但不一定对应真实塑料降解。真实 PET film/powder 筛选更可靠，但通量低。

需要建立：

```text
计算打分 → 模型底物快筛 → 真实聚合物 LC-MS 验证 → 材料级验证
```

的闭环。

---

## 4. 我们要做什么：总体研究假设

本课题的核心假设是：

```text
塑料降解酶的催化能力不仅由一个静态活性口袋决定，
而是由底物结合态、近攻击构象、过渡态、中间体和产物态的多状态稳定化共同决定。

如果我们能够从已知 QM/MM 反应路径中提取这些状态的几何和相互作用要求，
就可以寻找整体骨架不同但能够稳定同一反应状态集合的新 scaffold。
```

更具体地说：

```text
1. 口袋几何相似是必要条件；
2. 能够稳定 TS-like ensemble 是更关键条件；
3. 不能过度稳定 ES 或产物态；
4. 底物进入和产物释放必须可行；
5. 新骨架必须能折叠、表达并保持催化 motif。
```

---

## 5. 研究方案

## Aim 1：建立已知塑料降解酶的反应状态构象集合

### 目标

选择一个已有文献报道且 QM/MM 机制清楚的塑料降解酶体系，复现其反应路径，并得到可用于设计的多状态构象集合。

### 为什么要先做这个 Aim

如果只用一个 WT 的 ES 或 TS 结构去筛新骨架，会有严重问题：

```text
1. 酶反应不是单一构象；
2. 过渡态区域也有构象分布；
3. 同一反应路径可能存在多个可行 near-attack pose；
4. 单个 pose 很容易造成假阳性或假阴性。
```

因此需要 reaction-state ensemble。

### 具体做法

#### 1.1 选择体系

候选体系：

```text
PETase + PET-like oligomer；
LCC/cutinase + PET-like oligomer；
MHETase + MHET；
polyesterase + PCL/PLA/PBS oligomer；
PU esterase + urethane / polyester-urethane model compound。
```

选择标准：

```text
1. 有结构；
2. 有反应机制文献；
3. 底物和产物可检测；
4. 可建立 QM/MM 模型；
5. 最好有实验活性数据作为参考。
```

#### 1.2 建立 WT QM/MM 反应路径

需要获得：

```text
ES：底物结合态；
NAC：近攻击构象；
TS1：第一过渡态；
TI：四面体中间体或对应中间体；
TS2：第二过渡态；
PS：产物态。
```

#### 1.3 采样每个状态的 ensemble

方法：

```text
QM/MM umbrella sampling；
QM/MM metadynamics；
constrained reaction-coordinate sampling；
multiple pathway sampling；
near-attack conformation filtering。
```

每个状态保留：

```text
20–100 个代表构象。
```

聚类标准不是全蛋白 RMSD，而是反应相关特征：

```text
亲核原子-羰基碳距离；
攻击角；
oxyanion hole 氢键距离；
general acid/base 几何；
离去基团键长；
关键水分子位置；
底物扭转角；
关键口袋残基 rotamer。
```

### 输出

```text
WT reaction-state ensemble library
```

---

## Aim 2：固定 minimal catalytic motif 搜索或生成新骨架

### 目标

寻找整体 fold 不同，但能承载相同催化 motif 和反应状态几何的新骨架。

### 具体做法

#### 2.1 定义 minimal catalytic motif

固定：

```text
催化残基关键原子；
亲核攻击几何；
oxyanion hole；
general acid/base；
关键水分子；
TS-like ligand 或 intermediate-like ligand；
底物进入方向。
```

不固定：

```text
整个口袋；
所有二层残基；
全蛋白 fold。
```

#### 2.2 搜索或生成骨架

来源：

```text
PDB；
AlphaFold DB；
Foldseek / DALI / TM-align；
RFdiffusion-like motif scaffolding；
未注释结构蛋白。
```

#### 2.3 初步过滤

过滤标准：

```text
1. catalytic motif 能否放入；
2. 反应原子距离角度是否合理；
3. 是否有底物进入通道；
4. 是否有严重 clash；
5. 骨架是否可折叠；
6. 与已知酶 fold 是否不同；
7. 是否适合表达。
```

### 输出

```text
几百个计算候选骨架；
48–96 个进入实验验证的优先骨架。
```

---

## Aim 3：多状态计算打分并实验验证弱活性骨架

### 目标

用多状态 scoring 排除明显不合理骨架，再用实验验证其是否具有基本催化能力。

### 计算打分

对每个候选骨架评价：

```text
ES compatibility；
NAC compatibility；
TS1-like stabilization；
TI compatibility；
TS2-like stabilization；
PS / product release compatibility；
foldability；
access channel。
```

核心指标：

```text
TS-like state 是否相对于 ES 得到选择性稳定。
```

惩罚项：

```text
ES overbinding；
product trapping；
pocket collapse；
blocked access channel；
fold instability。
```

### 实验验证

#### 表达

```text
arrayed gene synthesis；
E. coli BL21(DE3) 或适合表达系统；
pET 系列载体；
96-well expression；
crude lysate 初筛；
top hits 纯化。
```

#### 筛选层级

一级：

```text
pNP ester / fluorescent ester / soluble oligomer
```

二级：

```text
BHET / MHET / PET-like oligomer
```

三级：

```text
真实 PET film / powder / nanoparticles
HPLC/UPLC/LC-MS 检测 TPA/MHET/BHET
```

### 输出

```text
1–5 个实验验证具有弱催化功能的新 scaffold。
```

---

## Aim 4：命中骨架的序列生成和 QM/MM 多状态标签

### 目标

对已验证弱活性的新骨架进行 fixed-backbone sequence design，生成约 1000 条序列，并用 QM/MM 生成多状态能量标签。

### 序列生成

方法：

```text
ProteinMPNN；
LigandMPNN；
带底物 / TS-like ligand / intermediate 的 fixed-backbone design。
```

规模：

```text
5 个 validated scaffolds × 200 条序列 = 1000 条
```

或：

```text
10 个 scaffolds × 100 条序列 = 1000 条
```

### 设计策略

```text
A. pocket-conservative design；
B. pocket-redesign design；
C. ligand-aware design；
D. transition-state-biased design；
E. product-release-aware design。
```

### QM/MM 标签生成

不使用绝对总能量，而使用相对能量：

```text
ΔE‡TS1 = E_TS1-like - E_ES
ΔE_TI = E_TI - E_ES
ΔE‡TS2 = E_TS2-like - E_TI
ΔE_PS = E_PS - E_ES

ΔΔE‡TS1 = ΔE‡TS1(mutant) - ΔE‡TS1(WT)
ΔΔE‡TS2 = ΔE‡TS2(mutant) - ΔE‡TS2(WT)
```

### ensemble 评分

对每条序列 s、状态 k、构象 i：

```text
E(s, k, i)
```

得到：

```text
G(s, k) = -RT ln Σ_i exp[-β E(s, k, i)]
```

### 多保真度策略

```text
1000 条：TS-like constrained QM/MM single-point 或局部约束优化；
100–200 条：局部 QM/MM reaction-coordinate scan；
10–30 条：true TS search / IRC / QM/MM umbrella sampling；
top candidates：实验验证。
```

### 输出

```text
1000 条左右设计序列的 QM/MM 多状态能量数据集。
```

---

## Aim 5：训练 QM/MM energy surrogate 并进行主动学习优化

### 目标

训练一个能快速预测不同序列背景对反应状态能量影响的小型模型，用于优先级排序。

### 输入特征

```text
1. pocket residue embedding；
2. second-shell embedding；
3. sequence background embedding；
4. pocket residue-ligand graph；
5. 反应原子距离/角度；
6. 氢键网络；
7. electrostatic field；
8. pocket volume and shape；
9. product exit geometry。
```

### 输出标签

```text
y1 = ΔE‡TS1-like；
y2 = ΔE_TI；
y3 = ΔE‡TS2-like；
y4 = ΔE_PS；
y5 = geometry penalty；
y6 = product trapping penalty。
```

### 模型类型

```text
XGBoost / LightGBM；
ridge / elastic net；
Gaussian process；
small MLP；
small pocket graph neural network。
```

### 训练方式

```text
multi-task regression + pairwise ranking
```

目标：

```text
不是精确预测绝对能垒，
而是筛选哪些序列值得进入下一轮高精度 QM/MM 或实验。
```

### 主动学习循环

```text
1. 模型预测更多 MPNN/LigandMPNN 设计序列；
2. 选择高分、高不确定性和结构多样的候选；
3. 做更高精度 QM/MM；
4. top hits 做实验；
5. 数据回流更新模型。
```

---

## 6. 高通量筛选体系

### 6.1 一级筛选

底物：

```text
pNP ester；
fluorescent ester；
BHET；
MHET；
soluble polyester oligomer。
```

作用：

```text
排除无表达、无折叠或无基本水解能力的候选。
```

### 6.2 二级筛选

底物：

```text
PET oligomer；
amorphous PET powder；
PET nanoparticles；
low-crystallinity PET film。
```

读出：

```text
HPLC / UPLC / LC-MS 检测 TPA / MHET / BHET。
```

### 6.3 三级验证

```text
PET film / bottle flakes；
GPC / SEC；
SEM / AFM；
质量损失；
产物分布；
热稳定性和残余活性。
```

---

## 7. 预期结果

```text
1. 建立反应状态引导的新骨架发现流程；
2. 获得一批口袋几何相似但整体骨架不同的新候选塑料降解酶；
3. 实验验证 1–5 个弱活性新骨架；
4. 建立约 1000 条序列的 QM/MM 多状态能量数据集；
5. 训练一个 QM/MM-supervised reaction-state ensemble energy predictor；
6. 得到更优设计序列；
7. 解释哪些序列背景更能稳定 TS-like ensemble、降低 ΔE‡、避免 product trapping。
```

---

## 8. 风险和应对

### 风险 1：新骨架没有实验活性

应对：

```text
增加 access channel、product release 和 foldability 过滤；
先用模型底物筛弱活性；
对弱活性骨架再做 MPNN 优化。
```

### 风险 2：1000 条 QM/MM 标签计算量过大

应对：

```text
多保真度标签；
大规模只做 TS-like constrained score；
小规模做真实 TS/free-energy validation。
```

### 风险 3：模型只学会计算打分函数

应对：

```text
明确模型是 QM/MM energy surrogate；
用实验 top hits 校准；
最终不把模型当直接活性预测器。
```

### 风险 4：WT reaction-state ensemble 不适用于新骨架

应对：

```text
只对同机制 scaffold 使用 WT-like ensemble；
对明显不同机制候选重新做 reaction path search。
```

---

## 9. 最终摘要

```text
本课题将从已有塑料降解酶的 QM/MM 反应机制出发，构建 ES/NAC/TS1/TI/TS2/PS 等反应状态构象集合，提取最小催化 motif 和过渡态稳定几何。随后通过结构搜索和生成式 scaffold design，寻找整体 fold 不同但能够承载相同催化口袋的新骨架，并通过多状态打分和实验筛选验证弱活性。在命中骨架上进一步进行 MPNN/LigandMPNN 序列设计，用 QM/MM 为约 1000 条设计序列生成反应状态相对能量标签，并训练反应状态能量代理模型，用于快速筛选更可能降低催化能垒、稳定过渡态且避免产物困陷的设计序列。
```

---

## 10. 关键参考文献

1. Yoshida et al., 2016, **Science**, *A bacterium that degrades and assimilates poly(ethylene terephthalate).*  
2. Austin et al., 2018, **PNAS**, *Characterization and engineering of a plastic-degrading aromatic polyesterase.*  
3. Tournier et al., 2020, **Nature**, *An engineered PET depolymerase to break down and recycle plastic bottles.*  
4. Knott et al., 2020, **PNAS**, *Characterization and engineering of a two-enzyme system for plastics depolymerization.*  
5. Lu et al., 2022, **Nature**, *Machine learning-aided engineering of hydrolases for PET depolymerization.*  
6. Sahihi et al., 2024, **JCIM**, LCC-ICCG cutinase adsorption on PET surface MD.  
7. Jäckering et al., 2024, **JACS Au**, *From Bulk to Binding: Decoding the Entry of PET into Hydrolase Binding Pockets.*  
8. RFdiffusion / motif scaffolding literature for de novo backbone generation.  
9. ProteinMPNN literature for fixed-backbone sequence design.  
10. LigandMPNN literature for ligand-aware sequence design.
