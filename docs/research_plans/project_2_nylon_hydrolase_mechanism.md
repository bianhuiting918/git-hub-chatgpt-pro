# 课题二：NylC/Nyl10/Nyl50-like 尼龙水解酶对 PA6/PA66 的识别、链段进入、催化机制和产物选择性

_Last updated: 2026-07-06_

## 0. 课题定位

本课题是一个**具体体系机制与工程课题**。与课题一的通用方法学不同，本课题聚焦尼龙酶，尤其是 NylC、Nyl10、Nyl50 以及后续可能发现的新 NylC-like / Ntn hydrolase 骨架。

核心问题是：

```text
NylC/Nyl10/Nyl50-like 酶如何识别 PA6/PA66，
如何吸附到聚合物或预处理寡聚体底物上，
如何让 PA 链段进入 substrate-access tunnel，
进入以后如何通过 N-terminal Thr 催化酰胺键水解，
以及这些过程如何决定 PA6/PA66 的产物分布。
```

本课题的重点不是一开始训练大模型，而是先建立清楚的机制图谱：

```text
PA6/PA66 substrate state
→ enzyme-polymer adsorption
→ productive vs nonproductive binding
→ chain-end / internal-loop entry
→ catalytic pose
→ QM/MM hydrolysis
→ product distribution
→ interface / tunnel engineering
```

---

## 1. 背景介绍

尼龙是聚酰胺材料，典型代表为 PA6 和 PA66。

```text
PA6:
[-NH-(CH2)5-CO-]n

PA66:
[-NH-(CH2)6-NH-CO-(CH2)4-CO-]n
```

它们的共同特点是：

```text
1. 主链由酰胺键连接；
2. 酰胺基之间形成强链间氢键；
3. 材料通常是半结晶聚合物；
4. 酶很难接触到内部酰胺键；
5. 反应常受材料可及性和链段进入限制，而不只是化学水解势垒限制。
```

经典尼龙酶体系包括：

```text
NylA：处理 PA6 cyclic dimer，开环；
NylB：处理 PA6 linear dimer / short oligomer，外切；
NylC：处理 DP > 3 的 PA6 linear/cyclic oligomer，内切。
```

其中 NylC 更接近真实 PA6/PA66 解聚问题，因为真实尼龙是高分子链，不是小二聚体。

NylC 属于：

```text
N-terminal nucleophile hydrolase, Ntn hydrolase
```

它先以前体形式表达，经过自剪切生成 α/β 两个亚基，暴露新的 N-terminal Thr。这个 N-terminal Thr 是亲核催化残基，负责攻击酰胺羰基。

Nyl50 也属于同一类 Ntn hydrolase，但它是近年发现的 PA66-selective nylon hydrolase。Nyl50 的结构显示其可能具有跨二聚体界面的 substrate-access tunnel，因此它可能通过不同的通道/界面结构实现 PA66 选择性。

---

## 2. 目前研究进展：与本课题相关的部分

### 2.1 PA6 方向：NylC 已证明可降解，但效率仍低

PA6 酶法解聚已有系统筛选。重要结论是：

```text
NylC-type enzyme 可以水解 PA6 film 或 PA6 powder；
但转化率仍然很低；
加入新酶不一定能重启平台期反应；
主要瓶颈很可能来自 substrate accessibility，而不是酶完全失活。
```

这说明：

```text
只提高催化中心化学反应能力不够；
必须理解 PA6 链段如何暴露、如何进入酶口袋、以及如何形成 productive binding。
```

### 2.2 PA66 方向：Nyl10 / Nyl50 是新的关键起点

近年 95 个 Ntn hydrolase homolog 的筛选发现：

```text
Nyl10 和 Nyl50 对 PA66 有明显选择性；
Nyl50 主要偏向 PA66 linear dimer 类产物；
NylC 更偏向 PA66 linear monomer；
不同 Ntn hydrolase 的产物谱明显不同。
```

这对本课题很重要，因为它说明：

```text
PA66 选择性不是 NylC 一个骨架的简单点突变问题；
天然或筛选得到的 Ntn hydrolase 空间中存在不同产物偏好 scaffold；
Nyl10/Nyl50 可以作为 PA66-selective 机制研究起点。
```

需要特别强调：

```text
Nyl10 和 Nyl50 不是 NylC 点突变体；
它们是独立的 NylC-like / Ntn hydrolase scaffold；
它们适合作为天然 PA66-selective seed 和机制对照。
```

### 2.3 Nyl50 的二聚体和 tunnel 线索

Nyl50 结构分析提示：

```text
1. Nyl50 可能主要以二聚体形式存在；
2. Nyl50 的 putative substrate-access tunnel 可能跨越两个 monomer；
3. Nyl50 的 catalytic residues 与 NylC 类似，但底物/产物选择性不同；
4. 因此选择性很可能来自通道、界面和底物进入路径，而不是催化 Thr 本身。
```

仍然没有解决的问题是：

```text
PA66 链段是否真的从这个 tunnel 进入？
Nyl50 的 A/D dimer interface 是否参与 productive PA66 binding？
Nyl50 为什么偏 PA66 L2，而 NylC 更偏 PA66 L1？
```

这些正是本课题要解决的机制问题。

### 2.4 NylC 工程化进展：从热稳定性到通道工程

已有 NylC 工程化路线包括：

```text
1. 热稳定化变体；
2. directed evolution 小幅提升 PA6 film 活性；
3. substrate-access channel 附近的突变，例如芳香残基引入；
4. 提高 PA6 film 产物释放和温度稳定性。
```

这些结果说明：

```text
NylC 改造不应只盯 catalytic Thr；
口袋入口、通道、loop 和表面界面可能是更重要的改造对象。
```

### 2.5 化学预处理 + 酶精修成为重要路线

尼龙和 amyloid fibril 不同。尼龙是材料，可以在酶处理前进行化学或物理预处理。

可行预处理包括：

```text
酸低聚化；
甲酸处理；
HCl / H2SO4 处理；
溶剂分散；
溶解-反沉淀；
机械粉碎；
热处理和淬冷。
```

预处理的作用：

```text
破坏链间氢键；
降低分子量；
增加链端和水可及性；
生成可溶或半可溶寡聚体；
让 NylC/Nyl50/NylB 更容易发挥选择性水解作用。
```

因此本课题不应局限于“酶直接啃完整高结晶尼龙”。更合理的方向是：

```text
化学预处理负责松解和粗解聚；
NylC/Nyl10/Nyl50-like 酶负责选择性水解和产物分布控制。
```

---

## 3. 当前缺口

### 缺口 1：Nyl10/Nyl50 的 PA66 选择性机制不清楚

已有研究知道 Nyl10/Nyl50 对 PA66 有活性，但还不知道：

```text
1. PA66 链段如何进入 Nyl50 tunnel；
2. Nyl50 是否通过 A/D dimer interface 形成 productive tunnel；
3. Nyl10/Nyl50 与 NylC 的产物偏好差异来自何处；
4. PA66 L1/L2 产物比例由 chain-end entry 还是 internal-loop entry 决定。
```

### 缺口 2：酶-尼龙材料界面模拟不足

现有研究多使用小寡聚体，例如：

```text
PA6 Ahx4/Ahx5；
PA66 L1/L2/C1/C2。
```

但真实底物是：

```text
PA6 film / powder；
PA66 powder / pellet / fiber；
预处理后的低聚物混合物。
```

缺少对以下过程的模拟：

```text
NylC/Nyl50 吸附到 PA6/PA66 表面；
productive vs nonproductive binding；
PA 链段从氢键网络中被拉出；
链段进入 active-site tunnel；
形成可催化构象。
```

### 缺口 3：NylC/Nyl50 的 QM/MM 水解机制不完整

NylB + PA6 linear dimer 的 QM/MM 已有较多研究，但 NylC/Nyl50 对 PA6/PA66 长链或寡聚体的：

```text
N-terminal Thr 攻击；
四面体中间体；
C-N 键断裂；
acyl-enzyme 或等效中间体；
deacylation；
product release；
```

仍缺少系统的 QM/MM 反应路径研究。

### 缺口 4：高通量筛选需要同时看总活性和产物分布

AMIDE/MAFC 可以快速看总胺释放，但不能区分产物。对本课题真正重要的是：

```text
PA6 Ahx1/Ahx2/Ahx3/Ahx4/Ahx5；
PA66 L1/L2/C1/C2；
linear vs cyclic products；
monomer vs dimer vs longer oligomers。
```

因此必须建立：

```text
AMIDE/MAFC 一级筛选 + LC-MS/OPSI-MS 二级产物谱 + 材料级验证
```

的多层筛选体系。

---

## 4. 本课题的核心改进点

### 改进点 1：从“总降解信号”转向“产物选择性机制”

传统筛选常看：

```text
总胺释放；
总产物量；
薄膜质量损失。
```

本课题重点是：

```text
为什么某个酶偏 PA6 或 PA66；
为什么某个酶偏 Ahx2、L1 或 L2；
什么样的口袋/通道决定产物分布。
```

### 改进点 2：从“底物对接口袋”转向“链段进入通道”

尼龙底物不是自由小分子，而是：

```text
高分子链段；
表面 loop；
chain end；
酸低聚化寡聚体；
PA66 L1/L2-like 片段。
```

因此需要研究：

```text
PA 链段如何从材料或低聚物状态进入 enzyme tunnel。
```

### 改进点 3：从“催化位点突变”转向“界面 + 通道 + 口袋协同改造”

改造对象包括：

```text
enzyme-polymer interface；
A/D dimer interface；
substrate-access tunnel；
active-site entrance loop；
oxyanion-stabilizing region；
product exit path。
```

### 改进点 4：把化学预处理产物纳入真实底物体系

本课题不只研究 pristine PA6/PA66，也研究：

```text
acid-oligomerized PA6/PA66；
solvent-dispersed PA6/PA66；
low-MW PA oligomer mixtures。
```

这样更贴近实际 chemo-enzymatic nylon recycling。

---

## 5. 研究方案

## Aim 0：以 Nyl10/Nyl50 为核心的新 Ntn hydrolase 骨架搜索

### 目标

发现新的 PA6/PA66-active Ntn hydrolase scaffold，为后续机制比较和工程改造提供更多骨架。

### 种子序列

```text
Nyl10；
Nyl50；
wild-type / natural NylC-like enzymes；
其他天然 PA6/PA66-active Ntn hydrolase homologs。
```

不作为 seed：

```text
NylC-HP；
NylC-V3；
其他 NylC 点突变体或 directed-evolution 产物。
```

这些工程化变体只作为 benchmark。

### 搜索方法

```text
1. MMseqs2 / BLAST / HMMER 序列搜索；
2. HHpred / profile-profile 搜索远缘同源；
3. Foldseek / DALI / AlphaFold 结构搜索；
4. CAVER / MOLE 分析 substrate-access tunnel；
5. SSN / structure network 聚类。
```

### 候选筛选标准

```text
1. 保留 Ntn hydrolase fold；
2. 保留自剪切位点和 N-terminal Thr/Ser/Cys；
3. 保留 catalytic acid/base network；
4. 有潜在 substrate-access tunnel；
5. 与 NylC/Nyl50 序列有一定距离；
6. 代表不同结构 cluster；
7. 预测可表达、可折叠、可自剪切。
```

### 实验初筛

```text
48–96 个候选基因 arrayed synthesis；
E. coli BL21(DE3) / pET 系统表达；
96-well crude lysate；
PA6/PA66 真实底物和低聚物底物筛选；
AMIDE/MAFC 一级筛；
LC-MS/OPSI-MS 二级确认。
```

### 输出

```text
3–5 个具有不同底物/产物偏好的新 NylC-like scaffold。
```

---

## Aim 1：建立 PA6/PA66 真实底物和预处理底物筛选体系

### 目标

建立与计算模拟和后续工程直接对应的底物体系。

### 底物分层

#### A. 真实固体底物

```text
PA6 film；
washed PA6 powder；
washed PA66 powder；
PA66 pellet / fiber / film；
cryomilled PA6/PA66 powder。
```

#### B. 预处理底物

```text
acid-oligomerized PA6；
acid-oligomerized PA66；
formic-acid-treated PA6/PA66；
HCl / H2SO4 low-MW oligomer mixtures；
solvent-dispersed PA6/PA66 microparticles。
```

#### C. 定义明确的机制底物

```text
PA6: Ahx3, Ahx4, Ahx5；
PA66: L1, L2, C1, C2；
必要时定制合成或从低聚化混合物中分离。
```

### 材料表征

```text
DSC / WAXS：结晶度和晶型；
GPC / SEC：Mn、Mw、链端密度；
SEM / AFM：表面形貌；
BET / particle size：表面积；
FTIR / XPS：氢键和表面化学；
LC-MS blank wash：预存低聚物背景。
```

### 高通量筛选体系

#### 一级：AMIDE/MAFC 总伯胺筛选

```text
PA6 film / PA6 powder / PA66 powder / acid-oligomer mixture
        ↓ 酶反应
释放含 -NH2 的产物
        ↓ MAFC 显色
A494 读数
```

作用：

```text
快速筛总水解能力。
```

限制：

```text
不能区分 Ahx1/Ahx2/Ahx3；
不能区分 PA66 L1/L2；
不能判断 endo/exo；
不能检测 caprolactam。
```

#### 二级：LC-MS / OPSI-MS 产物谱

检测：

```text
PA6:
Ahx1, Ahx2, Ahx3, Ahx4, Ahx5, cyclic oligomers

PA66:
L1, L2, C1, C2, HMD, adipic acid, longer oligomers
```

作用：

```text
判断 product selectivity。
```

#### 三级：材料级验证

```text
GPC/SEC：Mw 是否下降；
SEM/AFM：表面是否侵蚀；
DSC/WAXS：结晶度变化；
mass loss：总失重。
```

### 输出

```text
一套可用于 NylC/Nyl10/Nyl50/new scaffold 的标准化 PA6/PA66 筛选平台。
```

---

## Aim 2：酶-聚合物界面吸附和 productive/nonproductive binding

### 目标

解析 NylC/Nyl10/Nyl50/new scaffold 如何吸附到 PA6/PA66 表面，并区分能否导致催化。

### 模拟体系

```text
Nyl50 dimer + PA66 slab；
Nyl10 model + PA66 slab；
NylC tetramer/dimer + PA6/PA66 slab；
new scaffold + PA6/PA66 slab。
```

### 关键点

不只放一个 active-site-facing 初始构象，而是从多个方向开始：

```text
random orientation adsorption；
active-site-facing pose；
A/D dimer-interface-facing pose；
surface-loop-facing pose；
control nonproductive pose。
```

### 分析指标

```text
enzyme-surface contact map；
polymer-facing residues；
active-site tunnel orientation；
residence time；
productive binding fraction；
nonproductive adsorption classes；
PA chain distance to tunnel entrance。
```

### productive pose 标准

```text
1. scissile amide carbonyl C 靠近 N-terminal Thr Oγ；
2. carbonyl O 指向 oxyanion-stabilizing region；
3. leaving amide N 有质子转移路径；
4. PA chain 在 tunnel 中稳定；
5. 攻击角度合理；
6. 不严重阻碍产物释放。
```

### 输出

```text
每个酶的 productive vs nonproductive adsorption ensemble；
界面残基候选；
适合后续 chain-entry sampling 的初始构象。
```

---

## Aim 3：PA6/PA66 链段进入 substrate-access tunnel 的自由能

### 目标

定量 PA6/PA66 链段从材料表面或低聚物状态进入 NylC/Nyl50 active-site tunnel 的自由能代价。

### 关键自由能

```text
ΔG_extract：链段从 PA 氢键网络中被拉出；
ΔG_entry：链段进入酶 substrate-access tunnel；
ΔG_productive：loosely bound → productive catalytic pose；
ΔG_release：产物释放自由能。
```

### 方法

```text
umbrella sampling；
well-tempered metadynamics；
HREMD；
path CV；
funnel metadynamics；
string method。
```

### 比较对象

```text
NylC vs Nyl50 vs Nyl10；
PA6 vs PA66；
chain-end entry vs internal-loop entry；
solid surface segment vs acid-oligomerized substrate；
L1-producing pose vs L2-producing pose。
```

### 假设

```text
chain-end entry 更可能产生 L1 / Ahx1-Ahx2；
internal-loop entry 更可能产生 L2 / longer oligomers；
Nyl50 的 PA66 选择性可能来自更低的 PA66 chain-entry barrier；
NylC 与 Nyl50 的产物差异可能来自 tunnel geometry 和 product release，而不只是 catalytic Thr。
```

### 输出

```text
PA6/PA66 chain-entry free-energy landscape；
Nyl50 PA66 selectivity 的物理来源；
指导 tunnel / entrance loop 改造的关键位点。
```

---

## Aim 4：QM/MM 催化机制

### 目标

在 Aim 2–3 得到的可信 productive complex 上，研究 N-terminal Thr 介导的 PA6/PA66 酰胺键水解反应。

### 反应步骤

```text
1. N-terminal Thr Oγ 攻击酰胺羰基；
2. 形成四面体中间体；
3. C-N 键断裂；
4. 形成 acyl-enzyme 或等效中间体；
5. 水分子去酰化；
6. 产物释放。
```

### 体系选择

```text
NylC + PA6 Ahx4/Ahx5；
Nyl50 + PA66 L2-like substrate；
Nyl10 + PA66 substrate；
new scaffold + preferred substrate；
必要时比较 capped/internal oligomer 与 free oligomer。
```

### 关键问题

```text
1. 化学 TS 是否是限速步骤？
2. Nyl50 是否比 NylC 更好稳定 PA66 TS？
3. 还是 Nyl50 的优势主要来自 chain-entry 和 productive binding？
4. PA66 L1/L2 产物选择性是否由反应位点位置和 product release 决定？
```

### 输出

```text
NylC/Nyl50-like PA6/PA66 hydrolysis QM/MM energy profile；
关键催化残基和质子转移路径；
化学能垒与链段进入能垒的相对重要性。
```

---

## Aim 5：界面、通道和口袋改造

### 目标

基于 Aim 1–4 的结果，设计更高效、更具产物选择性的 NylC-like 酶。

### 5.1 界面改造

目标：

```text
提高 productive adsorption，降低 nonproductive adsorption。
```

候选区域：

```text
enzyme-polymer surface patches；
A/D dimer interface；
polymer-facing loop；
疏水/芳香/极性残基组合；
影响 PA 表面停留时间和朝向的电荷残基。
```

### 5.2 口袋 / 通道改造

目标：

```text
降低 ΔG_entry；
稳定 productive amide geometry；
控制 L1/L2/Ahx product distribution；
避免 product trapping。
```

候选区域：

```text
substrate-access tunnel；
active-site entrance loop；
oxyanion-stabilizing region；
PA chain subsites；
product exit path。
```

### 5.3 工艺兼容性改造

如果采用酸预处理 + 酶精修路线，改造还应关注：

```text
pH 5.5–8.0 稳定性；
残余 formate / chloride / sulfate tolerance；
残余溶剂 tolerance；
60–75 °C 热稳定性；
高寡聚体浓度下活性；
NylC/Nyl50 与 NylB 级联兼容性。
```

### 实验验证

```text
AMIDE/MAFC：一级总胺筛；
LC-MS/OPSI-MS：产物分布；
GPC/SEC：polymer Mw 变化；
SEM/AFM：表面变化；
DSC/WAXS：材料结晶度变化；
thermal stability assay：Tm 和半衰期。
```

---

## 6. 高通量筛选体系详细设计

### 6.1 表达体系

参考已有 NylC/Nyl10/Nyl50 文献，建议：

```text
E. coli BL21(DE3)；
pET-28a(+) 或 pET21a(+)；
96-well deep-well plate 表达；
20 °C overnight induction；
crude lysate 一级筛；
top hits 再纯化。
```

NylC-like/Nyl50-like 酶还需要确认：

```text
1. 是否可溶表达；
2. 是否发生自剪切成熟；
3. 是否形成正确 dimer/tetramer；
4. crude lysate 信号是否来自目标酶。
```

### 6.2 一级筛：AMIDE/MAFC

底物：

```text
PA6 film disc；
washed PA6 powder；
washed PA66 powder；
acid-oligomerized PA6/PA66 mixture。
```

读出：

```text
A494，总伯胺释放量。
```

适合筛：

```text
hundreds to thousands clones；
通道突变体；
界面突变体；
新 scaffold crude lysate。
```

### 6.3 二级筛：LC-MS / OPSI-MS

检测：

```text
PA6:
Ahx1–Ahx5, cyclic oligomers

PA66:
L1, L2, C1, C2, HMD, adipic acid
```

输出：

```text
total product；
product selectivity；
L2/L1 ratio；
linear/cyclic ratio；
monomer/dimer/longer oligomer ratio。
```

### 6.4 三级验证：材料级

```text
GPC/SEC；
SEM/AFM；
DSC/WAXS；
mass loss；
product accumulation over time；
enzyme re-addition / substrate re-addition platform test。
```

---

## 7. 预期结果

```text
1. 建立 NylC/Nyl10/Nyl50/new scaffold 的 PA6/PA66 substrate-selectivity 图谱；
2. 阐明 Nyl50 PA66 选择性是否来自 A/D dimer tunnel 和 chain-entry advantage；
3. 区分 chain-end cleavage 与 internal-loop cleavage 对 L1/L2 产物的贡献；
4. 建立 N-terminal Thr PA6/PA66 酰胺键水解 QM/MM 机制；
5. 得到可用于界面和通道改造的关键残基；
6. 建立 AMIDE + LC-MS/OPSI-MS + GPC 的尼龙酶筛选闭环。
```

---

## 8. 风险和应对

### 风险 1：真实 PA66 底物信号太低

应对：

```text
使用 acid-oligomerized PA66 mixture；
使用 washed PA66 powder + longer reaction time；
用 AMIDE 一级筛，LC-MS 二级确认；
提高底物表面积或使用 low-MW PA66。
```

### 风险 2：AMIDE 信号高但产物不理想

应对：

```text
所有 hits 必须用 LC-MS/OPSI-MS 看 L1/L2/Ahx 分布；
不能只按总胺选最高信号。
```

### 风险 3：Nyl50 结构中 tunnel 假设不成立

应对：

```text
用 multiple starting orientations MD；
比较 A/B vs A/D dimer models；
用 chain-entry enhanced sampling 验证；
若 tunnel 不参与，则寻找 alternative binding/entry route。
```

### 风险 4：QM/MM 起始结构不可信

应对：

```text
只对 Aim 2–3 确认的 productive complex 做 QM/MM；
不使用任意 docked pose；
比较多个 productive conformations；
做 reaction-coordinate validation。
```

### 风险 5：改造提高 binding 但不提高 catalysis

应对：

```text
区分 productive 和 nonproductive binding；
加入 product release penalty；
LC-MS 检查产物分布；
不要只看总吸附或总胺。
```

---

## 9. 最终摘要

```text
本课题以 NylC、Nyl10、Nyl50 及新筛选 NylC-like 骨架为对象，系统研究其对 PA6/PA66 的界面吸附、链段进入、催化水解和产物选择性机制。通过 AMIDE/MAFC 总胺筛选、LC-MS/OPSI-MS 产物谱、GPC/材料表征建立实验闭环；通过 PA6/PA66 slab MD、productive/nonproductive binding 分析、增强采样和 QM/MM 计算建立机制闭环。最终目标是明确 Nyl50-like PA66 选择性和 L1/L2 产物偏好的分子来源，并指导酶-聚合物界面、A/D dimer tunnel、active-site entrance 和 product exit pathway 的理性改造。
```

---

## 10. 关键文献与技术参考

1. NylA/NylB/NylC 经典 nylon oligomer hydrolase 体系。  
2. NylC 结构和 N-terminal Thr 自剪切机制。  
3. NylC PA6 film 筛选与 NylC_K-TS / NylC 工程化研究。  
4. Nyl10/Nyl50 PA66-selective Ntn hydrolase 筛选与 Nyl50 结构。  
5. AMIDE/MAFC 高通量 nylonase directed evolution 体系。  
6. I.DOT/OPSI-MS PA6/PA66 产物谱筛选体系。  
7. NylC-GYAQ + PA66 L1/L2/C1/C2 MD/SMD 机制研究。  
8. 化学预处理 + NylC/NylB 酶精修的 chemo-enzymatic PA6/PA66 monomerization 研究。
