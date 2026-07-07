# 尼龙酶课题两条并行路线：新骨架搜索与 PA66 分层机制模拟

_Last updated: 2026-07-07_

## 0. 当前研究定位

本文件在原有尼龙酶课题基础上进一步细化最新执行方案。现在课题可以分成两条并行推进的路线：

```text
方向一：新骨架搜索与初筛
从 Nyl10、Nyl50 和天然 NylC-like 序列出发，搜索新的 PA66/PA6-active Ntn hydrolase scaffold，
通过 embedding、结构建模、二聚体预测、稳定性、通道分析和 docking 逐层筛选。

方向二：代表性酶的 PA66 分层机制模拟
选择 NylC、NylC 功能增强突变体、Nyl10、Nyl50，
用 PA66 L2、L4、L8 建立“短产物处理—注册定位—长链段进入”的分层模拟体系。
```

这两条线的关系是：

```text
方向一负责扩展酶空间，找到新的候选 scaffold；
方向二负责解释 NylC/Nyl10/Nyl50 以及候选 scaffold 为什么有不同 PA66 产物偏好，
并为后续界面和通道改造提供机制依据。
```

---

## 1. 方向一：新骨架搜索与初筛

### 1.1 核心目标

目标不是简单找 NylC 的近缘同源物，而是寻找：

```text
1. 保留 Ntn hydrolase 自剪切和 N-terminal Thr 催化架构；
2. 可能具有 Nyl50-like PA66 substrate-access tunnel；
3. 可能形成与 Nyl50 类似的功能性二聚体；
4. 整体序列和结构与已知 NylC/Nyl50 有一定差异；
5. 可能具有不同 PA66 L1/L2 或 PA6 Ahx 产物偏好。
```

最终输出不是几千条序列，而是一个小而有代表性的实验候选 panel：

```text
48–96 条优先候选序列
```

其中包含：

```text
PA66-selective-like candidates
PA6-active-like candidates
broad-active candidates
结构和序列上与 Nyl50/NylC 不同的新 scaffold candidates
```

---

### 1.2 种子序列选择

#### 应作为种子的序列

```text
Nyl10
Nyl50
wild-type / natural NylC-like sequences
其他天然或筛选得到的 Ntn hydrolase nylon hydrolases
```

Nyl10 和 Nyl50 的定位：

```text
Nyl10 / Nyl50 是 PA66-selective Ntn hydrolase scaffold，
不是 NylC 的点突变体，
因此适合作为 PA66-active scaffold discovery 的核心种子。
```

#### 不作为种子的序列

```text
NylC-HP
NylC-V3
NylC-K-TS
其他 NylC directed-evolution 或 rational-design 点突变体
```

这些工程化变体可以作为：

```text
benchmark
positive control
后续工程对照
```

但不应参与定义天然新骨架搜索空间。原因是：

```text
它们代表已知骨架的工程优化路径，
会使搜索偏向已知 NylC 局部突变，而不是发现新 scaffold。
```

---

### 1.3 序列搜索流程

#### Step 1：常规序列同源搜索

工具：

```text
MMseqs2
BLASTP
jackhmmer / HMMER
HHpred / HHblits
```

数据库：

```text
UniRef90 / UniRef100
NCBI nr
MGnify / metagenome proteins
AlphaFold DB sequences
industrial / wastewater / soil / marine metagenome sequences, if available
```

输出：

```text
几千到几万条 Nyl10/Nyl50/NylC-like 候选序列。
```

初步过滤：

```text
1. 序列长度接近 NylC/Nyl50-like Ntn hydrolase；
2. 保留自剪切位点附近 N-terminal Thr/Ser/Cys；
3. 保留催化酸碱网络或等价残基；
4. 去除明显片段化、低复杂度过高、跨膜区过多序列；
5. 去除完全重复和过近序列。
```

---

### 1.4 Pythia / 蛋白语言模型 embedding 过滤

你提出的 `Pythia embedding 0.8 cutoff` 可以作为第二层过滤。这里建议把它明确写成：

```text
对候选序列计算 Pythia 或其他蛋白语言模型 embedding；
以 Nyl10/Nyl50/natural NylC-like active set 为参考，计算 cosine similarity；
用 0.8 作为初筛阈值，保留功能空间上接近 Nyl10/Nyl50/NylC-like 的序列。
```

建议不要把 0.8 作为最终硬标准，而是作为第一版经验阈值：

```text
embedding similarity ≥ 0.8：高可信 NylC/Nyl50-like 区域；
0.65–0.8：远缘探索区域，可少量保留；
< 0.65：除非结构搜索强支持，否则优先排除。
```

这样做的原因：

```text
如果阈值太高，会只留下近缘同源物；
如果阈值太低，会引入大量非 Ntn hydrolase 噪声；
因此 0.8 可以作为第一轮偏保守筛选，但要保留一部分低 identity / 结构相似的远缘候选。
```

embedding 过滤后的目标不是只留下最相似序列，而是得到两类候选：

```text
1. functionally similar candidates：embedding 接近 Nyl10/Nyl50；
2. structurally diverse candidates：embedding 稍远，但保留 Ntn fold 和 tunnel-like features。
```

---

### 1.5 序列 identity 去冗余与骨架多样性选择

如果 embedding 过滤后仍剩下很多序列，需要进一步用 identity 和聚类挑代表。

建议分层：

```text
90% identity cutoff：去除几乎重复序列；
70% identity cutoff：减少过近同源；
50% identity cutoff：得到中等多样性候选；
30–40% identity range：优先挑选可能代表新 scaffold 的远缘候选。
```

实际选择策略：

```text
1. 先按 70% identity 聚类，去冗余；
2. 再按 40–50% identity 选择不同 clade 代表；
3. 每个 sequence cluster 只保留 1–3 条代表；
4. 同时保留 Nyl10-like、Nyl50-like、NylC-like 和远缘 Ntn hydrolase-like 四类。
```

最终得到：

```text
200–500 条结构建模候选。
```

---

### 1.6 结构建模与二聚体装配预测

#### Step 1：单体结构预测

工具：

```text
AlphaFold / ColabFold / ESMFold
```

筛选指标：

```text
1. pLDDT；
2. PAE；
3. Ntn hydrolase fold 是否完整；
4. 自剪切 Thr/Ser/Cys 是否暴露为未来 N-terminal nucleophile；
5. α/β subunit packing 是否合理；
6. catalytic acid/base network 是否空间接近。
```

排除：

```text
1. active-site fold 崩坏；
2. 自剪切位点埋藏异常；
3. 口袋严重塌陷；
4. 预测结构不稳定或低置信。
```

#### Step 2：二聚体建模

因为 Nyl50 的 putative substrate-access tunnel 可能跨二聚体界面，因此必须考察候选酶是否能形成类似二聚体。

方法：

```text
AlphaFold-Multimer / ColabFold-Multimer；
以 Nyl50 A/D dimer 为模板进行结构对齐；
必要时建 A/B-like 和 A/D-like 两种 dimer model。
```

关注指标：

```text
1. ipTM / pTM；
2. interface PAE；
3. buried surface area；
4. interface hydrogen bonds / salt bridges；
5. interface hydrophobic packing；
6. dimer 是否能形成 continuous substrate-access tunnel；
7. catalytic pocket 是否跨界面或受界面稳定。
```

输出：

```text
dimerization-likely candidates
```

---

### 1.7 与 Nyl50 二聚体对齐

你提出“将这些酶 align 成 Nyl50 一样的二聚体”是合理的。

具体做法：

```text
1. 取 Nyl50 A/D dimer 作为 reference；
2. 对候选 dimer 做结构对齐；
3. 比较 active-site Thr、catalytic acid/base、tunnel entrance、subpocket 的空间位置；
4. 判断候选是否具有 Nyl50-like two-subpocket architecture。
```

关键不是 RMSD 越低越好，而是看以下局部结构是否保留：

```text
1. N-terminal nucleophile 位置；
2. catalytic pocket；
3. two larger subpockets；
4. tunnel bottleneck；
5. A/D interface pocket；
6. PA66 chain-entry route。
```

可以定义一个：

```text
Nyl50-like tunnel score
```

包括：

```text
tunnel volume
tunnel bottleneck radius
subpocket distance to catalytic Thr
two-unit PA66 register compatibility
A/D interface contribution
```

---

### 1.8 整体稳定性与二聚体形成趋势

建议同时做快速稳定性评估：

```text
1. FoldX / Rosetta ddG；
2. interface ddG；
3. ProteinMPNN sequence recovery or compatibility；
4. short MD stability check；
5. DSF/nanoDSF 后续实验测 Tm。
```

计算阶段可以建立以下指标：

```text
monomer stability score
dimer interface score
active-site preorganization score
tunnel integrity score
expression risk score
```

---

### 1.9 Docking 与三个方向的优先评分

你提到“docking 选择三个方向分数高的”。建议明确为三类 docking / binding scoring，不要只做一个 docking 分数。

#### 方向 A：PA66 L2 product / substrate docking

目标：

```text
筛选能稳定 PA66 L2 或 L2-producing pose 的候选，
优先发现 Nyl50-like L2-selective scaffold。
```

底物：

```text
PA66 L2
PA66 capped L3/L4 中的 U2–U3 register
```

关注：

```text
L2 是否稳定；
L2 是否容易释放；
U2–U3 键是否能靠近 catalytic Thr。
```

#### 方向 B：PA66 L4 register docking

目标：

```text
比较 L1-producing register 和 L2-producing register。
```

底物：

```text
PA66 L4
capped PA66 L4
```

构象：

```text
Register 1：切 U1–U2，预测 L1；
Register 2：切 U2–U3，预测 L2；
Register 3：切内部更深位置。
```

关注：

```text
哪个 register 的 productive geometry 更稳定；
候选更像 NylC 的 one-unit register，还是 Nyl50 的 two-unit register。
```

#### 方向 C：PA66 L8 / long-chain entry docking

目标：

```text
筛选能够容纳长 PA66 chain segment 的候选，
判断是否有 chain-entry potential。
```

底物：

```text
free PA66 L8
或 tethered PA66 L8 近似模型
```

关注：

```text
长链是否进入 tunnel；
是否出现严重 clash；
是否可以稳定两个 repeat units；
其余链段是否有出口或延伸路径。
```

---

### 1.10 候选排序与最终实验 panel

建议最终综合评分：

```text
Final scaffold priority score =
  embedding similarity score
+ catalytic motif score
+ dimer interface score
+ Nyl50-like tunnel score
+ L2 docking / register score
+ L8 entry compatibility score
+ stability score
+ sequence diversity bonus
- expression risk penalty
- redundancy penalty
```

最终选出：

```text
48–96 条候选
```

构成建议：

```text
20–30 条 Nyl50-like high L2-register candidates；
10–20 条 Nyl10-like PA66-selective candidates；
10–20 条 natural NylC-like candidates；
10–20 条远缘 Ntn hydrolase candidates；
5–8 条 controls：NylC、Nyl10、Nyl50、NylC enhanced mutants。
```

---

## 2. 方向二：代表性酶的 PA66 分层机制模拟

### 2.1 核心目标

用一套明确的底物层级，比较：

```text
NylC
NylC 功能增强突变体 1
NylC 功能增强突变体 2
Nyl10
Nyl50
```

对 PA66 不同长度和不同状态底物的识别、进入和催化差异。

核心假设：

```text
Nyl50 的 L2 产物偏好可能来自 two-unit register；
即 Nyl50 的 dimeric substrate-access tunnel 同时稳定 PA66 链端两个 repeat units，
并将 U2–U3 之间的酰胺键放到 Thr227 附近切割，释放 L2。

NylC 可能偏 one-unit register，
将 U1–U2 键放到 catalytic Thr 附近，释放 L1。
```

---

### 2.2 酶 panel

```text
E1：NylC
E2：NylC 功能增强突变体 1，例如 NylC-HP 或同类突变体
E3：NylC 功能增强突变体 2，例如 NylC-V3 或同类通道突变体
E4：Nyl10
E5：Nyl50
```

结构要求：

```text
Nyl50：优先用 dimer model，尤其 A/D-like dimer；
NylC：根据文献选择 tetramer/dimer model，同时保留成熟 α/β 状态；
Nyl10：若无晶体结构，使用 AlphaFold/AlphaFold-Multimer 结构，并标注结论不如 Nyl50 确定；
NylC mutants：在 NylC 结构上建突变并充分 relax。
```

---

### 2.3 底物分层

#### L2：短产物 / 继续水解模型

用途：

```text
判断 L2 在不同酶中是更像产物，还是更像可继续水解的底物。
```

问题：

```text
Nyl50 产生 L2 后，是否容易释放？
Nyl50 是否不容易继续把 L2 切成 L1？
NylC 是否容易把 L2 继续切成 L1？
Nyl10 是否类似 Nyl50？
```

模拟：

```text
5 个酶 × PA66 L2
```

分析：

```text
L2 binding pose；
L2 是否进入 catalytic pocket；
可切键是否靠近 Thr；
L2 product release tendency；
L2 further-cleavage probability。
```

---

#### L4：register scan 模型

L4 可以写成：

```text
[U1]—[U2]—[U3]—[U4]
```

其中 U 表示一个 PA66 repeat unit。

用途：

```text
测试酶更偏 L1-producing register 还是 L2-producing register。
```

注册方式：

```text
Register A：切 U1–U2，预测释放 L1；
Register B：切 U2–U3，预测释放 L2；
Register C：切 U3–U4，作为更深位置对照。
```

模拟：

```text
5 个酶 × free PA66 L4
5 个酶 × capped PA66 L4
```

为什么要 capped L4：

```text
free L4 两端有自由氨基/羧基，可能引入人工电荷和端基氢键；
capped L4 更接近真实 PA66 内部链段。
```

输出：

```text
每个酶最稳定的 register；
productive pose population；
U1/U2 subpocket occupancy；
L1 vs L2 predicted preference。
```

---

#### L8：长链段初始结合和 chain-entry 模型

L8 可以写成：

```text
[U1]—[U2]—[U3]—[U4]—[U5]—[U6]—[U7]—[U8]
```

用途：

```text
研究较长 PA66 链段如何接触酶表面，
是否能进入 substrate-access tunnel，
以及是否能形成 two-unit register。
```

不建议只使用完全自由 L8。建议做两种：

```text
free L8：模拟可溶长寡聚体；
tethered L8：模拟仍连在 PA66 表面的链段。
```

优先推荐：

```text
tethered PA66 L8
```

因为它更接近真实情况：

```text
PA66 surface — [U1]—[U2]—...—[U8]
```

模拟：

```text
5 个酶 × free L8：短 MD 筛选；
5 个酶 × tethered L8：短 MD 筛选；
从中选 2–3 个代表体系做 enhanced sampling。
```

---

### 2.4 第一轮模拟：口袋兼容性筛选

体系：

```text
5 个酶 × L2
5 个酶 × free L4
5 个酶 × capped L4
```

方法：

```text
docking；
短 MD；
side-chain repacking；
local relaxation；
productive pose filtering。
```

指标：

```text
1. scissile amide carbonyl C 到 catalytic Thr Oγ 距离；
2. 攻击角；
3. carbonyl O 是否被 oxyanion region 稳定；
4. leaving amide N 是否有质子转移路径；
5. PA66 repeat unit 在 subpocket 中的占据时间；
6. L1-register vs L2-register 稳定性。
```

输出：

```text
每个酶的初步 register preference。
```

---

### 2.5 第二轮模拟：L8 初始结合和 chain-entry

体系：

```text
5 个酶 × free L8
5 个酶 × tethered L8
```

方法：

```text
多个初始取向；
短 MD；
adsorption / binding pose clustering；
productive vs nonproductive classification。
```

分析：

```text
enzyme-L8 contact map；
tunnel entrance orientation；
U1/U2 是否靠近 Nyl50-like two subpockets；
chain bending angle；
PA66 self-contact / hydrogen bonding；
nonproductive trapping。
```

输出：

```text
进入增强采样的 2–4 个代表体系。
```

建议优先选择：

```text
Nyl50 + tethered L8；
NylC + tethered L8；
Nyl10 + tethered L8；
最佳 NylC enhanced mutant + tethered L8。
```

---

### 2.6 第三轮模拟：增强采样计算 chain-entry 自由能

方法可选：

```text
umbrella sampling；
well-tempered metadynamics；
HREMD；
path CV；
funnel metadynamics。
```

建议 CV：

```text
CV1：scissile amide carbonyl C 到 catalytic Thr Oγ 距离；
CV2：L8 插入 tunnel 深度；
CV3：U1/U2 two-subpocket occupancy；
CV4：PA66 链段弯曲角；
CV5：PA66-enzyme 接触数；
CV6：PA66-PA66 氢键数量。
```

输出：

```text
ΔG_entry；
ΔG_productive；
nonproductive minima；
L1-register vs L2-register free-energy difference。
```

关键判断：

```text
如果 Nyl50 的 L2-producing register 自由能更低，支持 two-unit register model；
如果 NylC 的 L1-producing register 更低，支持 one-unit register model；
如果 Nyl10 类似 Nyl50，则说明 PA66-selective enzymes 可能共享 register 机制；
如果 Nyl10 不同，则说明 PA66 选择性可以由不同通道/路径实现。
```

---

### 2.7 第四轮模拟：QM/MM 催化机制

不对所有组合做 QM/MM，只做最关键体系。

推荐体系：

```text
1. Nyl50 + PA66 L2
   判断 L2 是否更像产物而不是继续水解底物。

2. NylC + PA66 L2
   判断 NylC 是否容易把 L2 继续切成 L1。

3. Nyl50 + capped L4 的 U2–U3 register
   验证 L2-producing register 的催化可行性。

4. NylC + capped L4 的 U1–U2 register
   验证 L1-producing register 的催化可行性。

5. Nyl10 + predicted preferred register
   如果前面 MD 显示 Nyl10 有代表性。
```

QM/MM 反应步骤：

```text
N-terminal Thr Oγ 攻击酰胺羰基；
四面体中间体形成；
C-N 键断裂；
acyl-enzyme 或等效中间体；
水分子去酰化；
产物释放。
```

分析：

```text
chemical barrier；
TS stabilization；
proton-transfer path；
oxyanion stabilization；
L1/L2 product release tendency。
```

---

## 3. 两条路线如何连接

### 3.1 新骨架搜索为机制模拟提供新对象

方向一筛到的新候选，可以按以下标准进入方向二：

```text
1. docking 中表现出 L2-producing register；
2. 预测有 Nyl50-like dimer tunnel；
3. 二聚体界面稳定；
4. 与 Nyl50/NylC 序列差异较大；
5. AMIDE 或 LC-MS 初筛有 PA66 活性。
```

进入机制模拟的候选数量不宜多：

```text
2–3 个新 scaffold 即可。
```

### 3.2 机制模拟反过来优化新骨架筛选规则

如果方向二发现关键特征，例如：

```text
two-subpocket distance；
U1/U2 occupancy；
A/D interface tunnel radius；
L2-register free-energy；
product release path；
```

这些特征可以回流到方向一，作为新的筛选打分项。

最终形成：

```text
sequence/embedding search
→ structure/dimer/tunnel filtering
→ docking/register scoring
→ experimental screen
→ MD/QM-MM mechanism
→ improved scaffold scoring rules
```

---

## 4. 最小可执行版本

如果先做一个可落地的最小版本，建议：

```text
方向一：
1. Nyl10/Nyl50/NylC seed 搜索；
2. Pythia embedding ≥ 0.8 初筛；
3. identity 聚类去冗余；
4. 选 200–500 条建模；
5. align to Nyl50 dimer；
6. 二聚体和 tunnel 打分；
7. L2/L4/L8 docking；
8. 选 48–96 条实验候选。

方向二：
1. NylC、2 个增强突变体、Nyl10、Nyl50；
2. L2、free L4、capped L4 口袋/register 扫描；
3. free/tethered L8 初始结合；
4. NylC/Nyl50/Nyl10 中选代表做 chain-entry enhanced sampling；
5. 只对关键 register 做 QM/MM。
```

---

## 5. 预期结论类型

本方案最终希望得到的不是单一结论，而是一组机制规则：

```text
1. 哪些 Ntn hydrolase scaffold 具有 Nyl50-like PA66 tunnel；
2. 哪些结构特征支持 L2-producing register；
3. Nyl50 的 L2 偏好是否来自 two-unit register；
4. NylC 的 L1 偏好是否来自 one-unit register 或 L2 further-cleavage；
5. Nyl10 是否和 Nyl50 共享机制，还是代表另一类 PA66-selective route；
6. 哪些界面和通道残基可以改造以提高 PA66 L2 选择性。
```

---

## 6. 一句话总结

```text
本课题现在可以并行推进两条线：
一条从 Nyl10/Nyl50/NylC 出发进行 embedding、identity、结构、二聚体、tunnel 和 docking 多层筛选，寻找新的 PA66-active Ntn hydrolase scaffold；
另一条用 NylC、NylC 增强突变体、Nyl10 和 Nyl50 对 PA66 L2/L4/L8 进行分层模拟，分别解析短产物处理、L1/L2 register 和长链段进入机制。
两条路线最终在“PA66 L2 选择性和 two-unit register”这一机制假设上汇合。
```
