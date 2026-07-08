# 尼龙酶课题两条并行路线：新骨架搜索与 PA66 分层机制模拟

_Last updated: 2026-07-08_

## 0. 当前研究定位

本文件用于组织尼龙酶课题的两条并行路线：

```text
方向一：新骨架搜索与初筛
从 Nyl10、Nyl50 和天然 NylC-like 序列出发，搜索新的 PA66/PA6-active Ntn hydrolase scaffold。
核心筛选逻辑更新为：sequence/embedding 召回 → 结构模型 → 参考寡聚体装配 → interface 几何筛选 → pocket/tunnel 匹配 → 局部 refinement → 少量 docking/MD → 实验候选。

方向二：代表性酶的 PA66 分层机制模拟
选择 NylC、NylC 功能增强突变体、Nyl10、Nyl50，
用 PA66 L2、L4、L8 建立“短产物处理—链端 register—长链段进入”的分层模拟体系。
```

两条线的关系是：

```text
方向一负责扩展酶空间，找到新的候选 scaffold；
方向二负责解释 NylC/Nyl10/Nyl50 以及候选 scaffold 为什么有不同 PA66 产物偏好，
并把机制特征回流到方向一，作为新的筛选打分项。
```

---

## 1. 方向一：新骨架搜索与初筛

### 1.1 核心目标

目标不是简单找 NylC 的近缘同源物，而是寻找：

```text
1. 保留 Ntn hydrolase 自剪切和 N-terminal Thr/Ser/Cys 催化架构；
2. 可能具有 Nyl50-like PA66 substrate-access tunnel；
3. 可能形成与 Nyl50 类似的功能性二聚体，或与 NylC 类似的四聚体/二聚体装配；
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

### 1.2 种子序列和参考结构的分工

#### 用于搜索的 seed sequences

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
因此适合作为 PA66-active scaffold discovery 的核心搜索种子。
```

#### 不作为新骨架搜索 seed 的序列

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
机制模拟对象
```

但不应参与定义天然新骨架搜索空间。原因是：

```text
它们代表已知 NylC 骨架的工程优化路径，
会使搜索偏向 NylC 局部突变，
而不是发现新的天然或远缘 scaffold。
```

#### 用于装配评估的 reference oligomers

新骨架筛选需要单独建立一套 reference oligomer library：

```text
Nyl50 reference：优先整理 monomer、A/D-like dimer、必要时 A/B-like dimer。
NylC reference：整理 mature monomer、可能的 dimer、tetramer。
Nyl10 reference：如果没有实验结构，使用当前可获得的预测结构或同源模型，明确标注 source。
```

每个 reference 必须记录：

```text
reference_id
seed_name: NylC / Nyl50 / Nyl10
oligomer_mode: monomer / dimer / tetramer
chain_pair: A-B / A-D / A-C 等
structure_source: exact / homology proxy / AlphaFold / ESMFold / ColabFold
active_chain
catalytic nucleophile position
interface chains
是否成熟自剪切状态
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

`Pythia embedding 0.8 cutoff` 可以作为第二层过滤。建议明确写成：

```text
对候选序列计算 Pythia 或其他蛋白语言模型 embedding；
以 Nyl10/Nyl50/natural NylC-like active set 为参考，计算 cosine similarity；
用 0.8 作为初筛阈值，保留功能空间上接近 Nyl10/Nyl50/NylC-like 的序列。
```

不要把 0.8 作为最终硬标准，而是作为第一版经验阈值：

```text
embedding similarity ≥ 0.8：高可信 NylC/Nyl50-like 区域；
0.65–0.8：远缘探索区域，可少量保留；
< 0.65：除非结构搜索强支持，否则优先排除。
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

### 1.6 结构来源与第一版不等待新增预测结构

结构预测可以先挂着，不等它。第一版 interface 评估路线只用当前已有结构先跑：

```text
1. 实验结构：优先使用 Nyl50、NylC 或 NylC-like 已有结构；
2. 已有 AlphaFold / ColabFold / ESMFold 结构：作为 homology proxy 或 predicted structure；
3. Nyl10：若没有实验结构，使用当前最可信的预测结构，并明确标注 uncertainty；
4. 后续新增 ESMFold / ColabFold 结构完成后，再追加进入同一 pipeline。
```

每个候选结构必须记录：

```text
candidate_id
sequence_id
structure_source: exact / homology proxy / AlphaFold / ESMFold / ColabFold
model_confidence: pLDDT / pTM / ipTM / PAE if available
global_sequence_identity_to_seed
pocket_embedding_cosine
```

---

### 1.7 结构建模与单体质量筛选

#### 单体结构预测工具

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
6. catalytic acid/base network 是否空间接近；
7. active pocket 是否明显塌陷。
```

排除：

```text
1. active-site fold 崩坏；
2. 自剪切位点埋藏异常；
3. 口袋严重塌陷；
4. 预测结构不稳定或低置信；
5. 多段断裂或无意义 loop 侵入口袋。
```

---

### 1.8 Reference-based oligomer grafting：先固定参考寡聚方式

这里是新骨架筛选的核心更新。

第一轮不要直接让 docking 自由寻找界面。正确问题是：

```text
这个候选如果采用 NylC / Nyl50 / Nyl10 的已知或假定寡聚方式，
是否能形成合理界面、保留活性口袋，并产生 Nyl50-like tunnel？
```

因此第一步应做 reference-based oligomer grafting：

```text
1. 固定 reference oligomer：NylC、Nyl50、Nyl10 的 monomer/dimer/tetramer 模板；
2. 对候选 monomer 与 seed monomer 做结构比对；
3. 读取结构比对产生的旋转/平移矩阵；
4. 将候选 monomer 复制到 reference oligomer 每条 chain 的空间位置；
5. 得到 candidate-as-Nyl50-dimer 或 candidate-as-NylC-tetramer 装配模型；
6. 在这个模型上统计 interface 和 pocket，而不是先自由 docking。
```

结构比对工具：

```text
US-align
TM-align
Foldseek structure alignment
```

输出模型示例：

```text
candidate_001__Nyl50_AD_dimer_grafted.pdb
candidate_001__Nyl50_AB_dimer_grafted.pdb
candidate_001__NylC_ABCD_tetramer_grafted.pdb
candidate_001__Nyl10_dimer_grafted.pdb
```

这个模型不是最终稳定结构，而是回答一个筛选问题：

```text
候选骨架是否兼容参考寡聚体装配？
```

---

### 1.9 第一轮 interface 几何筛选

第一轮只做几何筛选，不急着 docking PA66 L2/L4/L8。

必须统计：

```text
1. interface clash score
2. buried surface area / interface area
3. interface residue contact count
4. interface hydrogen bond count
5. interface salt bridge count
6. interface RMSD to reference
7. monomer RMSD to seed
8. active pocket 是否被界面破坏
9. oligomer assembly 后是否严重重叠
10. catalytic Thr/Ser/Cys 是否仍在可及位置
```

建议的几何定义：

```text
clash：两个不同 chain 的非氢重原子距离 < 2.0 Å 或 < 2.2 Å；
contact：两个不同 chain 的重原子距离 < 4.0–4.5 Å；
hydrogen bond：供体-受体距离和角度同时合理；
salt bridge：带正/负电残基重原子距离 < 4.0 Å；
BSA = SASA(chain A) + SASA(chain B) - SASA(complex)。
```

建议第一轮排除：

```text
1. severe clash 很多，local refinement 难以修复；
2. BSA 过小，几乎没有界面；
3. interface contact count 过低；
4. active pocket 被另一条链直接堵死；
5. catalytic nucleophile 与参考位置偏移过大；
6. 二聚体/四聚体装配后出现大面积 backbone 重叠。
```

这一轮的目标不是判断真实结合自由能，而是淘汰明显不兼容参考寡聚方式的候选。

---

### 1.10 Pocket / tunnel 匹配：在参考装配后立即统计

pocket 形状要在 reference-grafted oligomer 上先统计一次。原因是：

```text
Nyl50-like PA66 tunnel 可能跨二聚体界面；
如果先自由 docking，界面变了，pocket/tunnel 的生物意义也变了；
因此必须先问：候选是否能在 Nyl50-like interface 中形成 tunnel。
```

统计工具：

```text
CAVER / CAVER Analyst
MOLEonline / MOLE
fpocket
P2Rank
自写 Python pocket descriptor
```

重点指标：

```text
tunnel 是否连续；
tunnel bottleneck radius；
tunnel volume；
subpocket number；
two-subpocket geometry；
subpocket-to-catalytic-Thr distance；
tunnel entrance 是否位于合理表面；
A/D interface 是否贡献 pocket wall；
装配后是否有 loop 或 side chain 堵住 tunnel。
```

建议定义：

```text
Nyl50-like tunnel score =
  tunnel continuity
+ bottleneck radius compatibility
+ catalytic pocket accessibility
+ two-subpocket geometry
+ interface contribution
+ PA66 two-unit register compatibility
- steric blockage penalty
```

---

### 1.11 第二轮：constrained local refinement，而不是全自由 docking

第一轮通过的候选，不应该直接做全自由 docking。应该做：

```text
参考约束下的局部 refinement / local docking
```

原则：

```text
不是固定死；
也不是任意重新找界面；
而是在 Nyl50-like / NylC-like 寡聚方式附近做小范围调整。
```

允许：

```text
1. side-chain repacking；
2. 局部 backbone relaxation；
3. monomer 之间小幅平移和旋转；
4. interface loop 小幅调整；
5. 能量最小化。
```

不允许第一轮就发生：

```text
1. 整个 monomer 翻到另一个完全不同界面；
2. tunnel 被重构成与 Nyl50 无关的新口袋；
3. catalytic pocket 远离参考装配位置；
4. 以 docking score 取代 reference interface compatibility。
```

可用工具：

```text
RosettaDock local mode
HADDOCK with interface restraints
LightDock with restraints
OpenMM minimization
GROMACS minimization
FoldX RepairPDB + AnalyseComplex
```

refinement 后重新统计：

```text
clash score
interface area
contact count
interface energy-like score
pocket volume
tunnel bottleneck
active-site geometry
```

如果 refinement 后 interface 变好且 tunnel 保留，则候选进入下一轮。

---

### 1.12 第三轮：少量自由 docking 只作为补充

自由 docking 可以做，但不能作为主筛选。它只用于：

```text
1. 探索 alternative oligomer mode；
2. 处理 reference-grafted 失败但整体结构仍可能二聚化的候选；
3. 与 AlphaFold-Multimer 预测结果交叉验证；
4. 少量 top candidates 的复核。
```

可用工具：

```text
ClusPro
HADDOCK global docking
LightDock global docking
RosettaDock global mode
AlphaFold-Multimer / ColabFold-Multimer
```

自由 docking 的结果必须重新追问：

```text
1. 是否保留 Ntn catalytic pocket？
2. 是否形成 PA66 substrate-access tunnel？
3. 是否能解释 L1/L2 产物偏好？
4. 是否与 Nyl50/NylC 的寡聚机制有关？
5. 是否得到 AlphaFold-Multimer / MD / interface energy 支持？
```

如果只是产生一个低能但无功能意义的界面，不应作为优先候选。

---

### 1.13 Docking 与底物 scoring：只对 interface 合理候选做

PA66 docking 不应该对几万条候选全做。只有满足下面条件的候选才进入底物 docking：

```text
1. monomer fold 合理；
2. 与 seed 单体 RMSD 合理；
3. reference oligomer grafting 后 clash 可接受；
4. interface area/contact count 合理；
5. active pocket 没被破坏；
6. pocket embedding 高；
7. Nyl50-like tunnel 或 NylC-like pocket 保留。
```

#### 方向 A：PA66 L2 product / substrate docking

目标：

```text
筛选能稳定 PA66 L2 或 L2-producing pose 的候选，
优先发现 Nyl50-like L2-selective scaffold。
```

底物：

```text
PA66 L2
PA66 one-end-free L4 中的 U2–U3 register
capped L4 中的 U2–U3 register 作为内部链段对照
```

关注：

```text
L2 是否稳定；
L2 是否容易释放；
U2–U3 键是否能靠近 catalytic Thr；
是否出现 two-unit register。
```

#### 方向 B：PA66 L4 register docking

目标：

```text
比较 L1-producing register 和 L2-producing register。
```

底物：

```text
free-free PA66 L4：可溶寡聚体对照；
one-end-free PA66 L4：链端底物，最重要；
capped PA66 L4：内部链段对照。
```

构象：

```text
Register A：切 U1–U2，预测释放 L1；
Register B：切 U2–U3，预测释放 L2；
Register C：切 U3–U4，作为更深位置对照。
```

关注：

```text
哪个 register 的 productive geometry 更稳定；
候选更像 NylC 的 one-unit register，还是 Nyl50 的 two-unit register；
端基是否造成假阳性结合。
```

#### 方向 C：PA66 L8 / long-chain entry docking

目标：

```text
筛选能够容纳长 PA66 chain segment 的候选，
判断是否有 chain-entry potential。
```

底物：

```text
free PA66 L8：可溶长寡聚体对照；
one-end-free / tethered PA66 L8：表面链端主模型；
fully capped / double-tethered L8：内部 loop 对照。
```

关注：

```text
长链是否进入 tunnel；
是否出现严重 clash；
是否可以稳定两个 repeat units；
其余链段是否有出口或延伸路径。
```

---

### 1.14 候选排序与最终实验 panel

建议最终综合评分：

```text
Final scaffold priority score =
  embedding similarity score
+ catalytic motif score
+ monomer structural quality score
+ reference oligomer compatibility score
+ dimer/tetramer interface geometry score
+ constrained refinement score
+ Nyl50-like tunnel score
+ L2 docking / register score
+ L8 entry compatibility score
+ stability score
+ sequence diversity bonus
- clash penalty
- expression risk penalty
- redundancy penalty
```

每个候选至少输出以下字段：

```text
candidate_id
seed_reference: NylC / Nyl50 / Nyl10
oligomer_mode: monomer / dimer / tetramer
reference_chain_pair: A-B / A-D / A-C / etc.
structure_source: exact / homology proxy / ESMFold / AlphaFold / ColabFold
monomer_RMSD_to_seed
interface_RMSD_to_reference
clash_score
interface_area
contact_count
hydrogen_bond_count
salt_bridge_count
active_pocket_status: intact / partially_blocked / broken
pocket_embedding_cosine
global_sequence_identity
Nyl50_like_tunnel_score
L2_register_score
L8_entry_score
final_priority
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

### 1.15 软件路线

第一版必需：

```text
US-align / TM-align：结构比对和 reference grafting；
Foldseek：结构搜索、结构聚类、快速相似性筛选；
FreeSASA：interface area / BSA；
BioPython / MDAnalysis：clash、contacts、氢键、盐桥统计；
CAVER / MOLE：tunnel 分析；
fpocket / P2Rank：pocket 识别和体积统计；
PyMOL / ChimeraX：人工检查 reference oligomer 和 pocket。
```

可选增强：

```text
Arpeggio：蛋白界面相互作用自动分析；
PDBePISA：少量结构的界面复核；
FoldX：快速 interface energy-like score；
Rosetta InterfaceAnalyzer：界面能量和 packing 复核；
RosettaDock local / HADDOCK restrained / LightDock restrained：局部 refinement；
OpenMM / GROMACS：少量 top hit 的 minimization 和短 MD。
```

PA66 docking / 结合评估：

```text
AutoDock Vina / smina：L2 和部分 L4 初筛；
GNINA：可作为深度学习 scoring 复核；
RosettaLigand：少量复杂体系；
手工 register pose + local minimization：L4/L8 更推荐。
```

---

### 1.16 第一版最小可执行路线：只用当前已有结构

第一版不依赖新增 ESMFold 结果，先用已有结构跑通 pipeline：

```text
Step 1：整理 reference oligomer library
NylC：monomer / dimer / tetramer；
Nyl50：monomer / A-D dimer / A-B dimer；
Nyl10：当前可用预测结构或同源模型。

Step 2：明确每个 reference 的 interface chain pair
例如 Nyl50 A-D、Nyl50 A-B、NylC A-B、NylC A-D、NylC tetramer interfaces。

Step 3：候选 monomer 与 seed monomer 结构比对
US-align / TM-align 输出 RMSD、TM-score、变换矩阵。

Step 4：用 reference transform 生成候选二聚体/四聚体
得到 candidate-as-reference oligomer models。

Step 5：几何统计
clash、BSA、contact count、interface RMSD、active pocket status。

Step 6：pocket/tunnel 统计
CAVER / fpocket / P2Rank，判断 Nyl50-like tunnel 是否保留。

Step 7：保留 top candidates
只有 interface 几何合理、pocket embedding 高、整体 RMSD 合理的候选进入后续 refinement/docking。

Step 8：新增预测结构完成后追加
把新增 ESMFold / ColabFold 结构放入同一输入目录，复用相同 pipeline。
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
E2：NylC-HP 或同类高活性/高稳定突变体
E3：NylC-V3 或同类通道突变体
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

推荐结构：

```text
free zwitterionic PA66 L2
H3N+—[U1]—[U2]—COO-
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

L4 不应只做两端自由或两端 capped。推荐四类：

```text
L4-FF：free-free L4
H3N+—[U1]—[U2]—[U3]—[U4]—COO-
代表完全游离寡聚体。

L4-Nfree：胺端自由，另一端封端
H3N+—[U1]—[U2]—[U3]—[U4]—NHMe
代表 PA66 胺端链端。

L4-Cfree：羧酸端自由，另一端封端
Ac—[U1]—[U2]—[U3]—[U4]—COO-
代表 PA66 羧酸端链端。

L4-capped：两端封端
Ac—[U1]—[U2]—[U3]—[U4]—NHMe
代表 PA66 内部链段。
```

优先级：

```text
第一优先级：L4-Nfree 和 L4-Cfree，一端自由、一端封端，代表链端底物；
第二优先级：L4-capped，用于判断是否可能内部切割；
第三优先级：L4-FF，用于可溶寡聚体对照。
```

输出：

```text
每个酶最稳定的 register；
productive pose population；
U1/U2 subpocket occupancy；
L1 vs L2 predicted preference；
端基是否主导结合。
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

不建议只使用完全自由 L8。建议做：

```text
free L8：模拟可溶长寡聚体；
one-end-free / tethered L8：模拟 PA66 表面链端；
fully capped / double-tethered L8：模拟内部 loop 对照。
```

优先推荐：

```text
one-end-free / tethered PA66 L8
```

因为它更接近真实情况：

```text
PA66 surface — [U8]—...—[U3]—[U2]—[U1]
                                      ↑
                                  自由链端进入酶通道
```

模拟：

```text
5 个酶 × free L8：短 MD 筛选；
5 个酶 × one-end-free/tethered L8：短 MD 筛选；
从中选 2–3 个代表体系做 enhanced sampling。
```

---

### 2.4 第一轮模拟：口袋兼容性筛选

体系：

```text
5 个酶 × L2
5 个酶 × L4-Nfree / L4-Cfree
5 个酶 × L4-capped
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
6. L1-register vs L2-register 稳定性；
7. 端基电荷是否造成人工盐桥。
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
5 个酶 × one-end-free/tethered L8
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

3. Nyl50 + one-end-free L4 的 U2–U3 register
   验证链端 L2-producing register 的催化可行性。

4. NylC + one-end-free L4 的 U1–U2 register
   验证链端 L1-producing register 的催化可行性。

5. Nyl50 + capped L4 的 U2–U3 register
   判断 Nyl50 是否可能内部切割。

6. Nyl10 + predicted preferred register
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
1. reference oligomer grafting 后 interface 几何合理；
2. 预测有 Nyl50-like dimer tunnel 或 NylC-like pocket；
3. constrained refinement 后界面和口袋仍稳定；
4. docking 中表现出 L2-producing register；
5. 与 Nyl50/NylC 序列差异较大；
6. AMIDE 或 LC-MS 初筛有 PA66 活性。
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
端基方向偏好；
chain-end vs internal-chain compatibility。
```

这些特征可以回流到方向一，作为新的筛选打分项。

最终形成：

```text
sequence/embedding search
→ structure filtering
→ reference oligomer grafting
→ interface geometry filtering
→ pocket/tunnel scoring
→ constrained refinement
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
4. 整理当前已有 NylC/Nyl50/Nyl10 参考结构；
5. 明确 reference oligomer chain pairs；
6. 候选 monomer align 到 seed monomer；
7. 按 reference transform 生成 dimer/tetramer grafted models；
8. 统计 clash、BSA、contact count、interface RMSD、active pocket status；
9. 统计 CAVER/fpocket/P2Rank pocket 和 tunnel；
10. 选 interface 几何合理、pocket embedding 高、整体 RMSD 合理的候选；
11. 只对 top candidates 做 constrained refinement 和 L2/L4/L8 docking；
12. 选 48–96 条实验候选。

方向二：
1. NylC、2 个增强突变体、Nyl10、Nyl50；
2. L2、one-end-free L4、capped L4 口袋/register 扫描；
3. free/tethered L8 初始结合；
4. NylC/Nyl50/Nyl10 中选代表做 chain-entry enhanced sampling；
5. 只对关键 register 做 QM/MM。
```

---

## 5. 预期结论类型

本方案最终希望得到的不是单一结论，而是一组机制规则：

```text
1. 哪些 Ntn hydrolase scaffold 具有 Nyl50-like PA66 tunnel；
2. 哪些候选能兼容 Nyl50-like dimer 或 NylC-like tetramer interface；
3. 哪些结构特征支持 L2-producing register；
4. Nyl50 的 L2 偏好是否来自 two-unit register；
5. NylC 的 L1 偏好是否来自 one-unit register 或 L2 further-cleavage；
6. Nyl10 是否和 Nyl50 共享机制，还是代表另一类 PA66-selective route；
7. 哪些界面和通道残基可以改造以提高 PA66 L2 选择性。
```

---

## 6. 一句话总结

```text
本课题现在可以并行推进两条线：
一条从 Nyl10/Nyl50/NylC 出发进行 embedding、identity、结构、reference oligomer grafting、interface 几何、tunnel/pocket 和 docking 多层筛选，寻找新的 PA66-active Ntn hydrolase scaffold；
另一条用 NylC、NylC 增强突变体、Nyl10 和 Nyl50 对 PA66 L2/L4/L8 进行分层模拟，分别解析短产物处理、链端 L1/L2 register 和长链段进入机制。
两条路线最终在“PA66 L2 选择性和 two-unit register”这一机制假设上汇合。
```
