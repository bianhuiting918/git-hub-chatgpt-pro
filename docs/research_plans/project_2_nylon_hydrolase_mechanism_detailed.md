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

## 2. 目前研究进展：按时间线整理

### 2.1 1990s–2000s：Kinoshita、Negoro 等建立 NylA/NylB/NylC 经典体系

早期尼龙酶研究主要来自能够利用 nylon-6 工业副产物的细菌体系。Kinoshita、Negoro 等研究者解析了 6-aminohexanoate oligomer degradation 的基因和酶学基础。

经典分工是：

```text
NylA：水解 6-aminohexanoate cyclic dimer；
NylB：水解 6-aminohexanoate linear dimer / short oligomers；
NylC：水解较长 linear/cyclic 6-aminohexanoate oligomers。
```

这一阶段完成了：

```text
1. 证明自然界存在尼龙寡聚体降解酶；
2. 建立 NylA/NylB/NylC 基本功能分工；
3. 将 nylon byproduct degradation 与具体酶联系起来。
```

但该阶段主要研究的是：

```text
可溶性 nylon oligomers；
工业副产物；
Ahx dimer / oligomer。
```

没有充分解决：

```text
真实 PA6/PA66 固体材料如何被酶识别和解聚。
```

---

### 2.2 2007–2014：NylC 结构与自剪切机制逐渐明确

NylC 后续被证明属于 N-terminal nucleophile hydrolase, Ntn hydrolase。它先以前体形式表达，随后在 Asn/Thr 附近自剪切，暴露 N-terminal Thr 作为亲核催化残基。

典型结构特征：

```text
1. 一个前体蛋白自剪切成 α 和 β 亚基；
2. β 亚基的新 N 端 Thr 是催化亲核残基；
3. 成熟 protomer 可进一步形成二聚体或四聚体；
4. 寡聚化可能影响稳定性、活性口袋和底物通道。
```

这一阶段完成了：

```text
1. 解释 NylC 为什么需要自剪切；
2. 证明 N-terminal Thr 是催化中心；
3. 确定 NylC 是 Ntn hydrolase 家族成员；
4. 初步说明寡聚化和亚基装配的重要性。
```

但仍没有解决：

```text
1. NylC 如何接触真实 PA6/PA66 固体表面；
2. PA 链段如何进入口袋；
3. PA6/PA66 产物选择性如何形成。
```

---

### 2.3 2014–2015：NylB 的 QM/MM 机制研究较成熟，但对象不是 NylC

Kamiya、Baba、Boero、Matsui、Negoro、Shigeta 等在 2014–2015 年对 NylB 做了 QM/MM 和 metadynamics 研究，研究对象是：

```text
NylB + PA6 linear dimer Ald / Ahx2
```

他们研究了：

```text
1. Ser112/Lys115/Tyr215 催化三联体；
2. Tyr170 对底物构象和四面体中间体崩解的作用；
3. acylation 反应和自由能势垒；
4. 突变体对化学步骤的影响。
```

这项工作的重要性：

```text
证明尼龙寡聚体水解酶可以用 QM/MM 研究化学 TS；
给出 NylB/Ahx2 小分子底物的催化机制。
```

但它不能直接回答本课题问题，因为：

```text
1. NylB 处理的是小线性二聚体；
2. NylB 不是 NylC/Nyl50；
3. 研究的是可溶小底物，不是真实 PA6/PA66 固体；
4. 没有涉及 PA 链段进入材料表面/通道的自由能。
```

因此，本课题的 QM/MM 空白是：

```text
NylC/Nyl50-like enzyme + PA6/PA66 longer chain / oligomer productive complex 的反应机制。
```

---

### 2.4 2023：AMIDE/MAFC 高通量尼龙酶筛选体系建立

Puetz 等在 **ACS Sustainable Chemistry & Engineering** 发表 AMIDE assay，建立了可以筛选 nylon-depolymerizing enzyme 的高通量显色方法。

该体系使用：

```text
MAFC = Meldrum’s acid furfural conjugate
```

检测机制是：

```text
尼龙水解产物中的伯胺 -NH2
        +
MAFC
        ↓
亲核加成 + 呋喃环开环
        ↓
共轭三烯显色产物
        ↓
A494 信号
```

该平台使用真实 PA6 film discs，96-well 表达 NylC 文库，筛选 1700 个克隆，获得 NylC_TS^P27Q/F301L，turnover frequency 提高约 1.9 倍。

这项工作完成了：

```text
1. 证明可用真实 PA6 film 进行 96-well 高通量初筛；
2. 建立总伯胺释放量的快速检测体系；
3. 为 NylC directed evolution 提供工具。
```

但 AMIDE 的局限是：

```text
1. 只能测总伯胺；
2. 不能区分 Ahx1/Ahx2/Ahx3；
3. 不能区分 PA66 L1/L2/C1/C2；
4. 不能判断 endo/exo 机制；
5. 不能直接说明产物分布是否理想。
```

因此本课题需要：

```text
AMIDE 一级筛 + LC-MS/OPSI-MS 二级产物谱 + GPC/材料验证。
```

---

### 2.5 2023–2025：ORNL 建立 I.DOT/OPSI-MS 产物谱筛选体系

ORNL 团队建立了 PAL/I.DOT + OPSI-MS 的高通量质谱方法，用于直接检测 PA6/PA66 水解产物。

该方法可以检测：

```text
PA6:
Ahx1–Ahx5, linear/cyclic oligomers

PA66:
L1, L2, C1, C2, HMD, adipic acid 等
```

这项工作的重要性：

```text
1. 能区分 PA6 和 PA66 具体产物；
2. 能判断产物选择性；
3. 比传统 HPLC 通量更高；
4. 适合筛选几十到几百个 Ntn hydrolase homologs。
```

它解决了 AMIDE 不能区分产物的问题，但也需要质谱平台和标准品。

---

### 2.6 2024：Nature Communications 系统筛选 PA6 活性酶，发现可及性仍是瓶颈

2024 年的 PA6 酶筛选工作系统评估了天然和工程化 nylonase 对 PA6 film/powder 的水解能力。

主要结论：

```text
1. NylC-type enzyme 是目前较有潜力的 PA6 hydrolase；
2. 最好的热稳定 NylC 变体对 PA6 film 的解聚仍然很低；
3. 反应平台期后加新酶不能有效重启；
4. 加新底物可以继续释放产物；
5. 因此限制很可能来自 substrate accessibility，而不是单纯酶失活。
```

这对本课题的启发：

```text
必须研究 PA6/PA66 表面可及性、链段进入和 productive binding；
不能只做催化中心突变。
```

---

### 2.7 2025：Nyl10/Nyl50 被发现为 PA66-selective nylon hydrolases

2025 年 ORNL/Michener/Foster/Bocharova/Chen 等团队筛选了 95 个 Ntn hydrolase homologs，发现 Nyl10 和 Nyl50 对 PA66 有明显选择性。

他们完成了：

```text
1. 证明 Ntn hydrolase 超家族中存在多个 nylon hydrolase；
2. 发现 Nyl10/Nyl50 对 PA66 有更高选择性；
3. 解析 Nyl50 结构；
4. 提出 Nyl50 可能通过 substrate-access tunnel 和 dimer interface 影响 PA66 选择性；
5. 发现不同酶的产物分布不同，例如 Nyl50 更偏 PA66 linear dimer，而 NylC 更偏 PA66 linear monomer。
```

但他们没有做：

```text
1. Nyl50 + PA66 slab 的吸附模拟；
2. PA66 链段进入 Nyl50 tunnel 的自由能；
3. Nyl50 对 PA66 酰胺键水解的 QM/MM；
4. Nyl50 为什么偏 L2 而 NylC 为什么偏 L1 的动态机制。
```

这正是本课题的关键切入点。

---

### 2.8 2025：PA66 高分子酶解机制研究表明真实 PA66 上可能偏 chain-end / surface-limited

2025 年 Polymer Chemistry 的 PA66 研究系统考察了 NylC-GYAQ 对 commercial PA66 的水解。

结论包括：

```text
1. PA66 水解受 molecular weight、surface area、crystallinity 共同影响；
2. commercial PA66 上主要观察到 L1 monomer 增加；
3. 没有明显 L2–L4 中间产物积累；
4. 说明在真实 PA66 高分子上，NylC 可能表现为 chain-end / surface-limited exo-like cleavage，而不是理想 random endo cleavage。
```

他们也做了 NylC-GYAQ + PA66 L1/L2/C1/C2 小底物 MD/SMD，但仍没有模拟真实 PA66 固体表面链段进入。

这对本课题的启发：

```text
要区分 chain-end entry 和 internal-loop entry；
要用链段进入自由能解释 L1/L2 产物偏好；
要比较 NylC 和 Nyl50 的不同进入模式。
```

---

### 2.9 2025：化学预处理 + NylC/NylB 酶精修证明“可及性”是主瓶颈之一

2025 年 PLOS ONE 等工作显示，通过 homogeneous dispersion、acid-induced oligomerization、甲酸/HCl/H2SO4 预处理，可以显著提高 PA6/PA66 对 NylC/NylB 的敏感性。

核心逻辑：

```text
化学预处理：破坏氢键网络、降低分子量、生成可溶/半可溶寡聚体；
酶处理：将寡聚体选择性水解或转化为更集中的产物。
```

这对本课题很重要，因为本课题不应只研究 pristine PA6/PA66，还应研究：

```text
acid-oligomerized PA6/PA66；
low-MW oligomer mixture；
defined oligomers；
真实材料 + 预处理产物。
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

## 4. 我们要做什么：总体研究假设

本课题的核心假设是：

```text
NylC/Nyl10/Nyl50-like 酶对 PA6/PA66 的活性和产物选择性，
主要由材料可及性、酶-聚合物界面吸附、PA 链段进入 substrate-access tunnel、
以及 N-terminal Thr 催化构象共同决定。

Nyl50 的 PA66 选择性可能不是来自更强的化学催化位点，
而是来自更合适的 PA66 chain-entry route、A/D dimer interface tunnel 和产物释放路径。
```

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

## 10. 关键参考文献

1. Kinoshita / Negoro 等早期 NylA/NylB/NylC nylon oligomer hydrolase 研究。  
2. NylC 结构和 N-terminal Thr 自剪切机制相关结构生物学研究。  
3. Kamiya / Boero / Negoro / Shigeta 等，2014–2015，NylB + PA6 dimer QM/MM 机制研究。  
4. Puetz et al., 2023, **ACS Sustainable Chemistry & Engineering**, AMIDE/MAFC 高通量 nylonase directed evolution 体系。  
5. PA6 film enzymatic depolymerization screening, 2024, **Nature Communications**，NylC_K-TS 等 PA6 活性筛选。  
6. Drufva / Michener / Foster / Bocharova / Chen 等，2025，Nyl10/Nyl50 PA66-selective Ntn hydrolase 筛选与 Nyl50 结构。  
7. Polymer Chemistry 2025，NylC-GYAQ 对 commercial PA66 的水解及 MW/surface area/crystallinity 影响。  
8. PLOS ONE 2025，chemical pretreatment + NylC/NylB chemo-enzymatic PA6/PA66 monomerization。  
9. I.DOT/OPSI-MS PA6/PA66 产物谱高通量检测方法。  
10. PET enzyme surface adsorption / chain-entry simulation literature as 方法学类比。
