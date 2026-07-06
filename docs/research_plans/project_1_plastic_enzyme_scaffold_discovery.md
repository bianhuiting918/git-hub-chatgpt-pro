# 课题一：反应状态引导的塑料降解酶新骨架发现与序列优化

_Last updated: 2026-07-06_

## 0. 课题定位

本课题是一个**方法学平台型课题**。核心不是继续在某一个已知塑料降解酶骨架上做点突变，而是从已知塑料降解酶的催化口袋和 QM/MM 反应机制出发，寻找或生成**整体骨架不同、但能够承载相同催化几何和反应状态稳定化模式**的新酶骨架。

一句话概括：

```text
用已知塑料降解酶的 QM/MM 反应状态定义“应该具备什么样的催化口袋”，
再寻找或生成不同整体骨架来承载这个口袋，
最后用 QM/MM 反应状态能量作为教师标签训练能量代理模型，指导序列优化。
```

本课题的目标不是一开始就预测实验活性，而是建立一个：

```text
QM/MM-supervised reaction-state ensemble energy predictor
```

中文可称为：

```text
基于 QM/MM 标签的反应状态构象集合能量预测模型
```

它预测的是：

```text
不同序列背景对 ES / NAC / TS1-like / TI / TS2-like / PS 等反应状态的相对稳定性；
不同序列是否能降低近似反应能垒；
不同序列是否会过度稳定底物基态或产物态。
```

它**不是**直接预测：

```text
kcat
kcat/KM
真实实验转化率
```

---

## 1. 背景介绍

塑料酶降解研究已经比较活跃，尤其是 PETase、LCC/cutinase、MHETase、polyesterase、PU esterase 等方向。已有大量研究集中在：

```text
1. 从自然界筛选塑料降解酶；
2. 对已知 PETase / cutinase 骨架进行定向进化；
3. 提高热稳定性；
4. 提高模型底物或真实 PET film 的降解产物释放；
5. 通过结构分析和 QM/MM 解释已知酶的催化机制。
```

这些研究推动了 PET/polyester 降解酶的发展，但也带来一个问题：**大多数工程优化仍然围绕少数已知骨架展开**。例如 PETase、LCC、TfCut2、FAST-PETase、HotPETase 等基本属于 cutinase/αβ-hydrolase 类似空间。继续在这些骨架上做点突变，创新性和可扩展性都会受到限制。

与此同时，蛋白设计领域已经具备了几个关键工具：

```text
1. 结构搜索：PDB / AlphaFold DB / Foldseek / DALI / TM-align；
2. motif scaffolding：在固定催化 motif 的情况下寻找或生成支撑骨架；
3. RFdiffusion/RFdiffusion-like backbone generation；
4. ProteinMPNN / LigandMPNN fixed-backbone sequence design；
5. ESM / ProtT5 等蛋白语言模型 embedding；
6. QM/MM 反应路径和过渡态能量计算。
```

这些工具使一个新的方向变得可能：

```text
不从序列同源性找新酶，
而从“反应几何和过渡态稳定化要求”出发，
寻找整体 fold 不同但功能口袋相似的新骨架。
```

---

## 2. 目前研究进展：与本课题相关的部分

### 2.1 实验筛选层面

塑料降解酶实验筛选通常分三类。

#### 2.1.1 小分子模型底物筛选

常见底物包括：

```text
p-nitrophenyl acetate / butyrate / hexanoate 等 pNP 酯；
荧光酯类底物；
MpNPT 等 MHETase-like 显色底物。
```

检测逻辑：

```text
酶切酯键
→ 释放黄色或荧光产物
→ 酶标仪快速读数
```

优点：

```text
通量高；
便宜；
可做 96/384 孔板；
适合排除完全无水解活性的候选。
```

缺点：

```text
模型底物太小；
容易筛出普通 esterase，而不是 polymer-active enzyme；
不能反映真实塑料表面吸附、链段进入、结晶度、产物释放。
```

因此，小分子模型底物只能作为一级预筛，不能作为最终活性标签。

#### 2.1.2 PET-like 低聚物底物

更接近 PET 的底物包括：

```text
BHET：双羟乙基对苯二甲酸酯；
MHET：单羟乙基对苯二甲酸酯；
PET dimer / trimer / tetramer / hexamer；
HEMT 或其他 PET-like analog。
```

这些底物更适合：

```text
1. 动力学分析；
2. 结构结合研究；
3. docking / MD / QM/MM；
4. 真实 PET 反应的中间层验证。
```

#### 2.1.3 真实塑料底物

最终验证必须用：

```text
PET film；
PET powder；
amorphous PET particles；
PET nanoparticles；
post-consumer PET flakes；
polyester textile particles；
PU dispersion / Impranil；
PCL / PLA / PBS film or powder。
```

检测读出通常是：

```text
HPLC / UPLC / LC-MS 定量 TPA / MHET / BHET 等产物；
GPC / SEC 看聚合物分子量变化；
SEM / AFM 看材料表面变化；
质量损失或力学性能变化作为补充。
```

这一层通量最低，但最能说明真实塑料降解能力。

---

### 2.2 计算设计层面

现有计算研究主要集中于：

```text
1. 已知酶的 docking 和 MD；
2. 已知酶的 QM/MM 反应路径；
3. 单个骨架内的突变筛选；
4. 结合态稳定性分析；
5. 过渡态或四面体中间体稳定化解释。
```

但仍有几个不足：

```text
1. 多数计算仍以单个底物结合构象为中心，缺少 reaction-state ensemble；
2. 多数设计只在已有骨架上做突变，没有系统发现新 scaffold；
3. docking 分数往往奖励底物基态结合，未必奖励过渡态稳定；
4. 很少把 QM/MM 得到的多状态能量作为可学习标签；
5. 很少明确区分“实验活性模型”和“QM/MM 能量代理模型”。
```

---

## 3. 本课题的核心改进点

### 改进点 1：从同源酶搜索转向反应状态约束的新骨架发现

传统路线：

```text
已知酶序列
→ BLAST / HMM 搜索同源序列
→ 表达筛选
```

本课题路线：

```text
已知酶反应机制
→ 提取 minimal catalytic motif 和 reaction-state ensemble
→ 搜索或生成不同整体骨架
→ 筛选能承载该反应几何的新 scaffold
```

目标不是找相似序列，而是找：

```text
口袋几何相似、整体骨架不同的新酶。
```

---

### 改进点 2：从单一 TS 结构转向 reaction-state ensemble

不使用单个野生型 TS1/TS2 结构硬套所有设计序列，而是先在野生型中采样多个反应状态构象集合：

```text
ES ensemble
NAC ensemble
TS1-like ensemble
TI ensemble
TS2-like ensemble
PS ensemble
```

然后对每条设计序列评价其对整个 ensemble 的稳定能力。

这样做的理由是：

```text
酶反应不是单个静态构象；
过渡态附近也存在构象集合；
不同序列背景可能稳定 ensemble 中不同子集；
只看一个结构容易误判。
```

---

### 改进点 3：用 QM/MM 能量标签训练 energy surrogate，而不是直接活性模型

本课题不把 1000 条计算设计序列当作真实实验活性数据。相反：

```text
QM/MM = teacher
ML model = student
```

QM/MM 给每条设计序列生成多状态能量标签；模型学习这些能量标签与口袋序列/结构之间的关系。

模型输出是：

```text
ΔE‡TS1-like
ΔE_TI
ΔE‡TS2-like
ΔE_PS
TS-stabilization score
product-trapping penalty
```

而不是直接输出：

```text
实验 kcat 或真实降解率。
```

---

### 改进点 4：多状态评分避免“只强化底物结合”

酶催化的核心不是把底物基态绑得越紧越好，而是相对于底物基态更好地稳定过渡态。

因此评分应关注：

```text
Selective TS stabilization = [G_TS-like - G_ES]_mutant - [G_TS-like - G_ES]_WT
```

如果这个值为负，说明设计序列相对于野生型更好地稳定 TS-like 状态。

同时需要惩罚：

```text
1. 过度稳定 ES；
2. 产物释放困难；
3. 口袋塌陷；
4. access channel 被堵；
5. 折叠不稳定。
```

---

## 4. 研究方案

### 方向 1：建立野生型反应状态集合

#### 目标

在一个已有文献报道且 QM/MM 机制清楚的塑料降解酶上，复现并扩展反应路径，得到可用于 scaffold 和序列设计的反应状态 ensemble。

#### 具体步骤

```text
1. 选择一个已有 QM/MM 文献体系作为起点；
2. 建立酶-底物复合物；
3. 复现 ES、TS1、TI、TS2、PS 等关键状态；
4. 用 umbrella sampling / metadynamics / constrained reaction-coordinate sampling 采样；
5. 对每个状态进行聚类；
6. 得到每个状态 20–100 个代表构象。
```

#### 需要记录的反应几何

```text
1. 亲核原子到羰基碳距离；
2. 攻击角；
3. oxyanion hole 氢键距离；
4. general acid/base 距离；
5. 离去基团键长；
6. 关键水分子位置；
7. 底物扭转角和应变；
8. 产品释放方向。
```

#### 预期产物

```text
WT reaction-state ensemble library
```

---

### 方向 2：固定 minimal catalytic motif 搜索或生成新骨架

#### 目标

寻找整体 fold 不同、但能够承载相同催化 motif 和反应状态几何的新蛋白骨架。

#### 候选来源

```text
1. PDB 结构数据库；
2. AlphaFold DB 预测结构；
3. Foldseek / DALI / TM-align 结构相似搜索；
4. RFdiffusion-like motif scaffolding；
5. 未注释或远缘功能蛋白。
```

#### 初筛标准

```text
1. 能放置 catalytic motif；
2. 反应关键原子距离和角度合理；
3. 口袋有足够空间容纳底物/过渡态；
4. 底物或低聚物有进入通道；
5. backbone 可折叠；
6. 与已知酶整体 fold 不同。
```

#### 多状态评分

对每个候选骨架评价：

```text
ES compatibility
NAC compatibility
TS1-like stabilization
TI compatibility
TS2-like stabilization
PS / product release compatibility
```

#### 输出

```text
48–96 个优先实验验证的新骨架。
```

---

### 方向 3：实验验证新骨架弱催化功能

#### 筛选层级

```text
一级：小分子模型底物 / 可溶低聚物
二级：PET-like / polyester-like 低聚物
三级：真实聚合物 film / powder / nanoparticles
```

#### 推荐检测方式

```text
1. pNP ester 或荧光酯：快速看是否有基本水解能力；
2. BHET / MHET / PET oligomer：看 PET-like 活性；
3. PET powder / film：HPLC/UPLC/LC-MS 测 TPA/MHET/BHET；
4. GPC / SEM：验证材料层面变化。
```

#### 实验策略

```text
1. 48–96 个候选骨架 arrayed gene synthesis；
2. E. coli BL21(DE3) 或合适表达系统 96-well 表达；
3. crude lysate 初筛；
4. top hits 纯化；
5. 真实底物验证。
```

#### 关键原则

第一轮实验目标不是获得最高活性，而是证明：

```text
该新骨架确实具有可检测的催化功能。
```

---

### 方向 4：命中骨架上的序列生成

#### 目标

对实验验证具有弱活性的新骨架，进行固定骨架序列优化。

#### 方法

```text
ProteinMPNN：生成可折叠序列；
LigandMPNN：在底物、过渡态类似物或中间体存在下生成口袋序列。
```

#### 设计规模

```text
5 个 validated scaffolds × 200 条序列 = 1000 条
或
10 个 scaffolds × 100 条序列 = 1000 条
```

#### 设计类型

```text
A. pocket-conservative design
保留一层催化口袋，只改二层和远端背景。

B. pocket-redesign design
允许口袋周围 6–8 Å 内残基改变。

C. ligand-aware design
显式放入底物或 TS-like ligand。

D. transition-state-biased design
以 TS-like state 为设计上下文。

E. product-release-aware design
避免产物被锁死。
```

---

### 方向 5：QM/MM 多状态标签生成

#### 目标

为 1000 条设计序列生成 QM/MM 衍生能量标签。

#### 标签层级

```text
Tier 1：1000 条
TS-like constrained QM/MM single-point 或 constrained optimization。

Tier 2：100–200 条
局部 QM/MM constrained optimization / reaction-coordinate scan。

Tier 3：10–30 条
true TS search、虚频验证、IRC/MEP 或 QM/MM free-energy sampling。

Tier 4：top candidates
实验验证。
```

#### 标签定义

不使用绝对总能量。使用相对能量：

```text
ΔE‡TS1 = E_TS1-like - E_ES
ΔE_TI = E_TI - E_ES
ΔE‡TS2 = E_TS2-like - E_TI
ΔE_PS = E_PS - E_ES

ΔΔE‡TS1 = ΔE‡TS1(mutant) - ΔE‡TS1(WT)
ΔΔE‡TS2 = ΔE‡TS2(mutant) - ΔE‡TS2(WT)
```

#### ensemble 评分

对每条序列 s、状态 k、构象 i：

```text
E(s, k, i)
```

然后计算：

```text
G(s, k) = -RT ln Σ_i exp[-β E(s, k, i)]
```

或带权重：

```text
G(s, k) = -RT ln Σ_i w_i exp[-β E(s, k, i)]
```

---

### 方向 6：训练 QM/MM 能量代理模型

#### 输入特征

```text
1. pocket residue embedding；
2. second-shell residue embedding；
3. protein language model embedding；
4. pocket residue-ligand graph；
5. 反应原子距离和角度；
6. 氢键网络；
7. electrostatic field；
8. pocket volume / shape；
9. product exit geometry。
```

#### 输出标签

```text
y1 = ΔE‡TS1-like
y2 = ΔE_TI
y3 = ΔE‡TS2-like
y4 = ΔE_PS
y5 = geometry penalty
y6 = product trapping penalty
```

#### 模型类型

```text
XGBoost / LightGBM
ridge / elastic net
Gaussian process
small MLP on frozen embeddings
small graph neural network for pocket graph
```

#### 训练目标

```text
multi-task regression + pairwise ranking
```

目标不是精确预测绝对能垒，而是：

```text
预测哪些序列值得进入下一轮 QM/MM 或实验验证。
```

---

## 5. 高通量筛选体系设计

### 5.1 一级筛选

```text
pNP ester / fluorescent ester
BHET / MHET / PET-like soluble oligomer
```

作用：

```text
排除表达失败、折叠失败或完全无水解活性的候选。
```

### 5.2 二级筛选

```text
真实或半真实聚合物底物：
PET powder / low-crystallinity PET film / PET nanoparticles
```

读出：

```text
HPLC / UPLC / LC-MS 检测 TPA / MHET / BHET。
```

### 5.3 三级验证

```text
GPC / SEC：聚合物分子量变化；
SEM / AFM：材料表面变化；
质量损失：总降解程度；
热稳定性：DSF/nanoDSF 或残余活性；
产物谱：是否有 product trapping 或副产物积累。
```

---

## 6. 预期成果

```text
1. reaction-state-guided scaffold discovery pipeline；
2. 一批不同整体骨架但承载相似催化口袋的新候选塑料降解酶；
3. 至少 1–5 个实验验证具有弱催化功能的新骨架；
4. 1000 条左右设计序列的 QM/MM reaction-state energy dataset；
5. 一个 QM/MM-supervised reaction-state ensemble energy predictor；
6. 一批经模型筛选的更优设计序列；
7. 对“序列背景如何影响 TS-like ensemble 稳定化”的机制解释。
```

---

## 7. 风险和应对

### 风险 1：新骨架口袋像，但没有活性

原因可能是：

```text
底物进不去；
产物出不来；
口袋柔性不对；
催化水/酸碱路径不对；
蛋白不表达或不稳定。
```

应对：

```text
在初筛阶段加入 access channel、product release 和 foldability 过滤；
实验上先确认弱活性，再做大规模序列优化。
```

### 风险 2：1000 条序列标签成本过高

应对：

```text
使用多保真度标签；
1000 条只做 TS-like constrained score；
100–200 条做局部 QM/MM scan；
10–30 条做真实 TS/free-energy validation。
```

### 风险 3：模型只学会计算打分函数，不代表真实活性

应对：

```text
明确模型定位为 QM/MM surrogate；
最终 top candidates 必须做真实底物实验验证；
用实验结果反向校正模型。
```

### 风险 4：同一骨架内数据不能泛化到新骨架

应对：

```text
区分 intra-scaffold model 和 cross-scaffold model；
逐步增加多个 validated scaffolds；
做 scaffold-split 评估。
```

---

## 8. 最终摘要

```text
本课题将以已知塑料降解酶的 QM/MM 反应状态为起点，提取 minimal catalytic motif 和 reaction-state ensemble，搜索或生成能够承载相同催化几何的新骨架。通过多状态打分筛选候选骨架并进行实验验证，获得具有弱催化功能的新 scaffold。随后在命中骨架上进行 MPNN/LigandMPNN 序列生成，并用 QM/MM 计算为约 1000 条设计序列建立 ES/TS1/TI/TS2/PS 相对能量标签，训练反应状态能量代理模型，用于快速预测不同序列背景对过渡态稳定化、反应能垒和产物释放的影响，最终指导实验验证和进一步优化。
```

---

## 9. 关键文献与技术参考

1. PETase / cutinase / LCC 等已知塑料降解酶的结构、突变和 QM/MM 机制文献。  
2. PET 真实底物检测体系：PET film/powder/nanoparticles + HPLC/UPLC/LC-MS 检测 TPA/MHET/BHET。  
3. RFdiffusion / motif scaffolding：用于固定功能 motif 的新骨架生成。  
4. ProteinMPNN：固定骨架序列设计。  
5. LigandMPNN：显式考虑小分子/配体/过渡态类似物环境的序列设计。  
6. Ensemble docking / QM/MM reaction-state ensemble：用于避免单一静态构象打分。
