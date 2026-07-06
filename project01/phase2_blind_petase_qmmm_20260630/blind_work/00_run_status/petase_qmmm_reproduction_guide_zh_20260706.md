# PETase QM/MM 盲复现完整流程中文总说明

更新日期: 2026-07-06

这份 Markdown 是 PETase QM/MM 盲复现工作的中文总说明。它不是只告诉你“下一步做什么”，而是把从结构信息出发、生成 Michaelis complex 假设、建立 Amber QM/MM、寻找 TS 集合、做 committor 验证、再到 umbrella/PMF 的全部路线串起来。读完这份文件，应该能回答四个问题：

- 我们为什么这样设计计算，而不是直接照搬文献坐标；
- 每一轮 attempt 在解决什么科学或数值问题；
- 现在已有结果支持什么结论，哪些结论还不能下；
- 如果要复现、检查或继续推进，应该去服务器上看哪些目录和文件。

注意：GitHub 里只保存小型文本记录、路径、决策逻辑和可复现说明；不保存服务器密码、原始轨迹 `md.nc`、重启文件 `md.rst7` 或大型 `md.out` 日志。

## 0. 总入口和当前结论

服务器项目根目录：

```text
/Dell/Dell14/bianht/petase_blind_qmmm/repo/project01/phase2_blind_petase_qmmm_20260630
```

GitHub 文档入口：

```text
project01/phase2_blind_petase_qmmm_20260630/blind_work/00_run_status/
```

同目录下还有两个英文/表格记录：

```text
petase_blind_qmmm_reproduction_runbook_20260706.md
petase_qmmm_attempt_ledger_20260706.tsv
```

当前总状态：

```text
去酰化 deacylation:
  已完成 38 窗口长 production 和最终 MBAR 分析。
  当前主结论能垒 = 19.061 +/- 0.073 kcal/mol。
  block sensitivity 稳定；没有坏窗口；唯一弱 overlap 在产物侧，不阻塞主能垒。

酰化 acylation:
  cand02 已通过 committor 支持为 TS-core 候选。
  attempts042/043/045 的短 production 与补窗已完成。
  attempt047 正在/已经用于 25 窗口长 production，完成后才应给最终酰化 PMF。
```

## 1. 这项复现的核心原则

用户的目标不是“把文献的结果复制一遍”，而是假设一开始没有那篇文章，只能根据 PETase 结构和通用酶催化化学，从第一性原理出发解析机制。文献只能在方法学层面提供启发，例如：可以用 QM/MM、committor、umbrella/PMF；不能用它的 TS 坐标、窗口、反应坐标、能垒或机理结论作为输入。

允许使用：

```text
PETase 结构信息
通用丝氨酸水解酶化学
通用活性中心几何判断
对接和 pose 筛选
Amber/Sander QM/MM 方法
ATESA / committor / umbrella / MBAR 这些计算思想
```

在验证完成前不允许使用：

```text
文献 TS 坐标
文献 Michaelis complex 坐标
文献 umbrella windows
文献具体 reaction coordinate 公式
文献 trajectory / restart
文献能垒、速控步结论、关键残基运动结论
```

整个流程用的内部 grill-me 标准是：

```text
1. 现有数据真的支持这个结论吗？
2. 这个结论最可能在哪一点误导我们？
3. 下一步计算是否在减少最大的剩余不确定性？
```

## 2. 从零出发的完整研究路线

如果从 PETase 结构出发，完整路线不是直接 umbrella sampling，而是分成下面几层。

### 2.1 建立反应假设

先根据丝氨酸水解酶的一般化学确定可能的反应类型：

```text
酰化 acylation:
  SerO 进攻底物酯羰基 C；
  SerH 转移给 His；
  C-Oleaving 键断裂；
  形成酰酶中间体。

去酰化 deacylation:
  水分子进攻酰酶羰基 C；
  水 H 转移给 His；
  SerO-C 酰酶键断裂；
  释放水解产物并再生活性 Ser。
```

这一步只用通用化学，不使用文献的具体反应态。

### 2.2 构造 Michaelis complex 或酰酶初始构象

起点不是 PMF，而是合理的反应前几何。要通过结构、对接、pose 过滤和短程动力学筛掉不合理构象。

对酰化，需要底物酯羰基满足：

```text
SerO 靠近羰基 C，方向接近亲核进攻；
SerH 与 His NE2 有可转移几何；
离去 O 与羰基 C 保持成键但可以被拉长；
底物整体在结合口袋中不严重冲突。
```

对去酰化，需要酰酶中间体和水分子满足：

```text
SerO-C 酰酶键已经存在；
反应水靠近羰基 C；
水的一个 H 指向 His NE2；
羰基 O 能被 oxyanion hole 稳定；
水不是随机 bulk water，而是能形成可反应氢键网络的 active-site water。
```

### 2.3 建立 Amber/Sander QM/MM 模型

最终正式路线选择 Amber/Sander，而不是早先 GMX/CP2K 的势能测试路线。原因是这里的主目标是动态 TS 验证和 PMF，需要与 ATESA/Amber restraint/trajectory 分析更直接对接。

正式生产设置：

```text
QM/MM 程序: Amber/Sander
QM 方法: DFTB3
温度: 310 K
生产 timestep: dt = 0.0001 ps
Thermostat: Langevin, gamma_ln = 5.0
轨迹输出: ntwx = 20
```

这里的 `dt = 0.0001 ps` 比文献式 1 fs 更保守。原因是盲探索阶段的反应性 restraint 容易导致数值不稳定，小步长减少了 vlimit、SCC 和爆炸轨迹造成的假失败。

### 2.4 先找 TS-core，再做 PMF

关键逻辑：不能一上来就把某条坐标当成反应坐标做 PMF。需要先证明某些构象确实在 reactant/product 分界附近。

因此路线是：

```text
粗略路径/扫描 -> TS-like seed -> ATESA/aimless shooting -> focused committor -> TS-core 候选 -> umbrella windows -> MBAR/PMF
```

committor 的作用是回答：从这个结构出发，随机化速度后，向反应物和产物两边走的概率是否接近 0.5。如果接近，它才是 TS 集合附近的候选，而不只是某个高能构象或人为拉出来的中间态。

## 2A. 关键技术概念详细解释

这一节专门解释后文反复出现的技术词。读者如果不熟悉 QM/MM、ATESA、aimless shooting、committor、umbrella sampling、PMF、MBAR，应该先读这一节，再看后面的 attempt 复盘。

### 2A.1 QM/MM 是什么

QM/MM 是把体系分成两层来算：反应中心用量子力学 QM，周围蛋白和水用分子力学 MM。

为什么需要 QM：

```text
化学反应涉及成键、断键、质子转移和电荷重排。
普通经典力场不能自然描述键断裂和新键生成。
因此 SerO-C、C-Oleaving、水 O-C、H-His 等反应键必须放在 QM 区或至少由 QM 描述。
```

为什么还需要 MM：

```text
整个 PETase + 底物 + 水盒太大，全部用 QM 不现实。
蛋白远端主要提供构象、电场、疏水环境和 steric constraint，可以用 MM 近似。
```

在本项目中的含义：

```text
QM 区负责真实反应化学：亲核进攻、酰酶键断裂/形成、质子转移。
MM 区负责酶环境：口袋形状、氢键网络、蛋白电场、溶剂约束。
```

常见误解：

```text
QM/MM 不等于“比经典 MD 更高级所以一定正确”。
结果仍依赖 QM 区选择、半经验方法精度、初始构象、采样长度和 reaction coordinate。
所以我们还需要 committor、overlap、block sensitivity 等验证。
```

### 2A.2 Amber/Sander、DFTB3、SCC-DFTB 分别是什么

Amber 是分子模拟生态；Sander 是 Amber 里可以跑 QM/MM MD 的程序之一。

DFTB3 是一种半经验量子方法，可以看作比常规 DFT 便宜很多的近似 QM。它适合做较长时间的 QM/MM 动力学，但精度通常低于高等级 DFT。

SCC-DFTB 中的 SCC 是 self-consistent charge，自洽电荷。每一步 QM 计算需要迭代得到稳定电荷分布。如果某一步 SCC 没收敛，Amber 会警告：这一小步的能量和力不够准确。

在本项目中怎么判断 SCC warning：

```text
偶发 1-2 次 SCC warning：
  不一定代表整条 trajectory 报废，但该窗口必须标记，完成后重点复查 CV 和能量连续性。

连续大量 SCC warning：
  说明 QM 电子结构或几何有严重问题，通常需要低 dt 重跑、调整初始构象或补窗。

伴随 vlimit、NaN、SANDER BOMB：
  视为坏窗口或失败窗口。
```

当前例子：

```text
acylation attempt047 的 w007_ts208 在 NSTEP=18800 附近出现 SCC-DFTB convergence warning。
它目前没有 vlimit 或硬错误，所以先让它跑完；完成后 CV/error summary 会把它作为重点窗口复查。
```

### 2A.3 timestep、production、equilibration、restart 是什么

`timestep` 是 MD 每一步前进的时间长度。本项目正式 production 用：

```text
dt = 0.0001 ps = 0.1 fs
```

为什么这么小：

```text
反应性 QM/MM + umbrella restraint 会让某些键和质子运动很快。
如果 timestep 太大，容易出现 vlimit、能量漂移、结构爆炸或 SCC 不收敛。
小 timestep 更慢，但对盲探索更稳。
```

`production` 指正式用于统计的采样，不是测试跑。比如：

```text
attempt040 deacylation: 38 windows x 40000 steps
attempt047 acylation: 25 windows x 40000 steps
```

`restart` 或 `md.rst7` 是某个 MD 窗口最后的坐标和速度/盒子信息，可用于：

```text
继续延长 production；
作为补窗 seed；
从坏窗口前一阶段重启；
构造下一轮 umbrella window。
```

### 2A.4 Michaelis complex 是什么

Michaelis complex 是酶和底物结合后、反应真正发生前的复合物。对酰化来说，它应该是“底物已结合，Ser 还没进攻或刚准备进攻”的状态。

不是任意 docked pose 都是 Michaelis complex。合理的 Michaelis complex 至少要满足：

```text
底物酯羰基 C 朝向 SerO；
SerO-C 距离和进攻角允许亲核进攻；
SerH 与 His NE2 能形成质子转移路径；
离去 O 仍与羰基 C 成键；
底物在口袋里没有明显冲突；
短 MD 后这些几何不会立即崩掉。
```

对去酰化，起点更接近酰酶中间体，而不是普通 Michaelis complex：

```text
SerO-C 酰酶键已形成；
反应水在活性位点；
水 O 指向羰基 C；
水 H 指向 His；
羰基 O 有 oxyanion-hole 稳定。
```


### 2A.4b 如何诱导形成 Michaelis complex

Michaelis complex 不是“随便 dock 一个底物进去”就完成了。真正可用于反应计算的 Michaelis complex 要同时满足结合合理性和反应几何合理性。我们的盲法思路是：先用结构和通用化学产生多个候选，再通过几何过滤、短动力学、必要的软约束诱导，把体系带到一个不靠文献坐标也能自洽的反应前态。

完整诱导流程如下。

**第一步：准备蛋白和底物**

目的：得到没有明显结构错误、质子化状态合理、能进入对接/MD 的体系。

输入：

```text
PETase 晶体结构或 AlphaFold/修复后的结构；
底物或模型底物结构；
催化三联体 Ser-His-Asp；
活性口袋附近关键残基和可能反应水。
```

动作：

```text
1. 检查缺失残基、缺失原子、异常构象和晶体水。
2. 判断催化 His 的质子化方式，保证 NE2 能作为接收 SerH 或水 H 的碱。
3. 保留或重新评估活性位点附近结构水。
4. 为底物生成 3D 构象和力场参数。
```

输出：

```text
可对接的蛋白结构；
可参数化的底物结构；
明确的催化原子编号候选。
```

失败标准：

```text
His 质子化导致无法接受质子；
SerO 被错误质子化或方向明显不合理；
底物参数化失败；
活性位点缺失关键原子或构象严重冲突。
```

**第二步：生成多个 substrate poses**

目的：不要假设唯一结合方式。先产生多个可能的底物结合构象，再筛选。

动作：

```text
1. 以活性 Ser、His 和口袋为中心定义 docking box。
2. 允许底物酯键附近有构象自由度。
3. 保留多个排名靠前 pose，而不是只拿 docking score 第一名。
4. 对每个 pose 记录底物羰基 C、离去 O、酯 O、芳香环/链段相对口袋的位置。
```

为什么 docking score 不能直接决定 Michaelis complex：

```text
Docking score 主要衡量结合，不保证反应几何。
一个结合很紧的 pose 可能 SerO 完全无法进攻羰基 C。
所以 docking 只用于生成候选，不是最终判据。
```

**第三步：几何过滤**

目的：用通用丝氨酸水解酶化学筛掉不能反应的 pose。

酰化 Michaelis complex 的过滤指标：

```text
SerO_C 距离：SerO 到底物羰基 C 应在可诱导进攻范围内，通常不应远到无法靠近。
attack angle：SerO -> C=O 的进攻方向应接近 Bürgi-Dunitz 亲核进攻方向，而不是从羰基平面错误方向进入。
SerH_His 距离：SerH 到 His NE2 应能形成质子转移前体。
C_Oleave：离去 C-O 键应仍是成键状态，不能一开始就是产物。
oxyanion stabilization：羰基 O 应朝向可能的 oxyanion hole 氢键供体。
clash：底物不能和蛋白严重重叠。
```

去酰化酰酶/水构象的过滤指标：

```text
SerO_C：酰酶键应存在或可保持。
OW_C：候选水 O 应能靠近酰基 C。
Hwater_His：水的一个 H 应能朝向 His NE2。
WatO-C-Ocarbonyl angle：水进攻方向应接近羰基亲核进攻方向。
water identity：反应水最好是活性位点稳定水或短 MD 中稳定进入口袋的水，而非瞬时 bulk water。
```

输出：

```text
通过几何过滤的 pose/rank 列表；
被拒绝 pose 的原因；
每个候选 pose 的关键距离和角度。
```

**第四步：短 classical MD 或 restrained relaxation**

目的：判断 pose 是否只是静态对接偶然结果，还是能在蛋白环境中保持反应前几何。

动作：

```text
1. 对体系加水、加离子，做能量最小化。
2. 对蛋白骨架和底物重原子做轻约束，逐步升温到 310 K。
3. 在短 MD 中观察底物是否逃出口袋，羰基是否仍面向 SerO/水。
4. 对关键反应几何可以使用很弱的 flat-bottom restraint，防止一开始就漂走，但不能直接拉成产物。
```

为什么可以用弱约束：

```text
Michaelis complex 是反应前态，不是平衡 docking pose。
如果完全不约束，短 MD 可能因为初始 pose 还没放松就漂走；
如果约束太强，又会人为制造反应态。
所以只允许用来保持合理接触，不允许把旧键拉断或新键强行拉成。
```

通过标准：

```text
底物仍在活性口袋；
SerO_C / OW_C / H-His 等关键距离没有完全跑坏；
没有明显 steric clash 或水盒/蛋白结构异常；
His-Asp-Ser 氢键网络保持合理；
羰基 O 仍有被 oxyanion hole 稳定的可能。
```

失败标准：

```text
底物逃逸；
羰基翻转到 SerO 无法进攻方向；
His 与 SerH/水 H 断开太远；
口袋塌陷或底物严重冲突；
反应水离开且没有替代水进入。
```

**第五步：QM/MM 预松弛或短 restrained QM/MM**

目的：经典力场不能准确描述即将反应的键和电荷。进入正式 TS/PMF 前，需要确认在 QM/MM 力场下反应前几何仍稳定。

动作：

```text
1. 选定 QM 区：反应 Ser/His/Asp 相关原子、底物反应片段、必要水分子。
2. 用 DFTB3 QM/MM 对候选构象做短最小化或短 MD。
3. 保持弱几何 restraint，只用于防止候选立即漂走。
4. 检查 SCC、vlimit、温度和关键 CV。
```

通过标准：

```text
QM/MM 下结构不爆；
SCC warning 不连续出现；
关键 CV 保持反应前态而不是直接变成产物；
反应中心没有明显不合理电荷/键长。
```

**第六步：选出可以进入 TS search 的 Michaelis-like seed**

最终 seed 不是“看起来最像文献”的结构，而是满足以下条件的结构：

```text
1. 从结构出发独立产生；
2. docking/pose 几何合理；
3. 短 MD 后仍保持 active-site binding；
4. QM/MM 下数值稳定；
5. 反应原子距离和角度足以通过温和 restraint 或 shooting 接近 TS region；
6. 没有已经偏到产物态。
```

这些 seed 才能进入后续 scan、aimless shooting、committor 或 umbrella window 构建。

### 2A.4c 本项目中 Michaelis complex/酰酶构象的记录应怎么看

相关 GitHub/服务器记录分布在：

```text
blind_work/01_system_setup/
blind_work/02_classical_md/
blind_work/03_mechanism_tree/
blind_work/04_qmmm_exploration/
```

应该重点看：

```text
stage1_ligand_and_protonation_execution_protocol.md
stage1_pose_geometry_filter_protocol.md
ligand_model_manifest.tsv
pose_generation_queue.tsv
gs_pose_manifest.tsv
rejected_pose_manifest.tsv
productive_conformer_manifest.tsv
candidate_cv_sets.tsv
path_screening_table.tsv
ts_like_guess_manifest.tsv
```

每个文件的意义：

```text
gs_pose_manifest.tsv:
  记录通过几何筛选的 ground-state / Michaelis-like pose。

rejected_pose_manifest.tsv:
  记录为什么某些 pose 被拒绝，例如距离、角度、冲突或质子转移几何不合理。

productive_conformer_manifest.tsv:
  记录短 MD 或构象筛选后仍能保持反应几何的 conformer。

candidate_cv_sets.tsv:
  记录我们认为可能描述反应进程的 CV 组合，不等于最终 reaction coordinate。

path_screening_table.tsv:
  记录早期 QM/MM path 或 scan 的筛选结果。

ts_like_guess_manifest.tsv:
  记录几何上可能接近 TS 的 guess，后续还必须用 shooting/committor 验证。
```

如果要判断 Michaelis complex 是否可信，不能只看一个 PDB。至少要同时看：

```text
1. 初始 pose 的几何指标；
2. 短 MD 后这些指标是否保持；
3. 是否有 rejected pose 作为反例；
4. 进入 QM/MM 后是否仍稳定；
5. 后续 TS search 是否能从该构象找到 boundary behavior。
```

### 2A.5 TS、TS-like、TS-core、TS ensemble 的区别

TS 是 transition state，过渡态。严格地说，过渡态是反应自由能面上分隔 reactant basin 和 product basin 的分界面附近状态。

在真实酶 QM/MM 动力学里，我们通常不会只得到一个“完美 TS 点”，而是得到一组接近分界面的构象。

几个词的区别：

```text
TS-like:
  几何上看起来像过渡态，例如一个键半断半连。
  但它不一定真的有 50% 概率去反应物/产物两侧。

TS-core:
  经过 committor 或 shooting 证据支持，确实接近反应分界面的核心候选。

TS ensemble:
  一组 TS-core 附近构象，而不是单个结构。
  酶反应在有限温度下本来就是构象集合，不应过度迷信一个坐标。
```

本项目为什么强调 TS-core：

```text
如果只根据某个 C-O 键长等于 1.9 Å 就说它是 TS，很容易错。
正确做法是看从这个点出发的短动力学会不会有一部分到 reactant、一部分到 product。
```

### 2A.6 Basin 和 endpoint classification 是什么

`basin` 是稳定区域。简单说：

```text
reactant basin：反应物区域。
product basin：产物区域。
undecided：短时间内没有清楚到任一侧，或几何混乱。
```

endpoint classification 是把一条短 trajectory 的终点按几何规则分类。分类规则不能来自文献结论，而应来自化学定义。

去酰化 endpoint 例子：

```text
product-like:
  水 O-C 新键形成，SerO-C 旧酰酶键断裂，水 H 已转移或接近 His。

reactant-like:
  SerO-C 酰酶键仍在，水 O-C 未形成，水 H 仍在水附近。

undecided:
  部分键变化但没有完整反应，或关键几何混乱。
```

酰化 endpoint 例子：

```text
product-like:
  SerO-C 键形成，C-Oleaving 键断裂，SerH 转给 His。

reactant-like:
  SerO-C 未形成，C-Oleaving 仍在，SerH 仍在 SerO。
```

### 2A.7 Committor 和 pB 是什么

committor 是判断一个构象是不是 TS 附近的核心工具。

做法：

```text
拿同一个候选结构；
随机化多组速度；
从这个结构出发跑很多条短 trajectory；
看每条最后到 reactant basin 还是 product basin；
product 的比例叫 pB。
```

解释：

```text
pB ≈ 0：基本都会回反应物，说明太 reactant-like。
pB ≈ 1：基本都会去产物，说明太 product-like。
pB ≈ 0.5：两边概率接近，说明接近反应分界面，也就是 TS ensemble 附近。
```

本项目例子：

```text
acylation cand02:
  forward/product-like = 18
  reactant-like = 13
  other = 1
  pB_committed = 0.580645
  pB_all = 0.5625
```

这说明 cand02 不是纯 product-like，也不是纯 reactant-like，而是接近 TS-core。

常见误解：

```text
pB 不需要精确等于 0.500000。
有限 shots 下会有统计波动。
关键是它不能明显偏向 0 或 1，并且 endpoint 分类要化学合理。
```

### 2A.8 Aimless shooting 是什么

aimless shooting 是 transition path sampling 的一种做法，用来从候选分界点附近自动寻找反应路径和 TS ensemble。

核心思想：

```text
选一个可能接近 TS 的结构作为 shooting point；
随机给它一组速度；
从这个点向前跑一条 trajectory；
把速度反向或重新随机化，再向后跑一条 trajectory；
看这两端是不是分别落到 reactant 和 product basin。
```

如果一对 forward/backward trajectories 分别到达两个不同 basin，这个 shooting move 就有可能被接受，因为它说明 shooting point 靠近连接反应物和产物的过渡路径。

为什么叫 aimless：

```text
它不需要预先精确知道最佳 reaction coordinate。
它通过不断随机速度和移动 shooting point，让系统自己告诉我们哪些结构在反应通道附近。
```

在本项目中：

```text
deacyl_attempt003 和 deacyl_attempt004 用 aimless shooting/ATESA 找去酰化边界 seed。
这些 seed 后来不是直接当最终结果，而是进入 focused committor 和 PMF 建路。
```

常见误解：

```text
aimless shooting 不是 umbrella sampling。
它不是为了直接给自由能曲线，而是为了找真实动力学边界附近的 seed/TS ensemble。
```

### 2A.9 ATESA 是什么

ATESA 是一个自动化 enhanced sampling/transition path sampling 工作流工具。它可以帮忙管理 aimless shooting、endpoint 分类、accepted move、候选 reaction coordinate 选择等任务。

本项目里 ATESA 的作用：

```text
1. 从多个 TS-like seed 开始自动发起 shooting trajectories。
2. 根据 endpoint 分类判断 shooting move 是否连接 reactant/product。
3. 收集 accepted 或 near-boundary seeds。
4. 给后续 focused committor 和 PMF 提供候选起点。
```

ATESA 不是魔法黑箱。它仍然依赖：

```text
合理的初始 seed；
合理的 endpoint classification；
足够长的 commitment trajectory；
正确的 QM/MM 稳定性；
后续人工/脚本复查。
```

为什么 attempt003/004 之后还要 focused committor：

```text
ATESA accepted seed 说明它有边界行为，但不等于已经精确定义最终 TS ensemble。
所以我们继续用 focused committor 检查 seed 1、seed 4 等候选的 pB。
```

### 2A.10 LMax、reaction coordinate 和 CV 是什么

CV 是 collective variable，集体变量。它是用少数几何量描述复杂反应进程的方式，比如距离、角度、距离差。

reaction coordinate 是真正能区分反应进程的坐标。一个 CV 不一定就是好的 reaction coordinate。

LMax 通常指在候选变量组合中，通过似然或分类模型寻找最能区分 reactant/product endpoint 的 reaction-coordinate 形式。

在本项目中实际含义：

```text
我们会先列出可能重要的距离/角度：SerO-C、C-Oleaving、H-His、水 O-C 等。
然后通过 shooting/committor/endpoint 数据判断哪些变量真的区分两侧 basin。
```

常见误解：

```text
eta = 一个旧键距离 - 一个新键距离。
eta 很有用，但不一定足够。
如果质子转移、进攻角和成键/断键不同步，只用 eta 会产生假峰或错误 overlap 判断。
```

### 2A.11 Umbrella sampling 是什么

umbrella sampling 是为了计算自由能曲线而做的 biased sampling。

问题背景：

```text
普通 MD 很难自发跨越反应能垒。
如果只等自然反应，可能跑很久都看不到一次反应。
```

umbrella 的做法：

```text
沿反应路径选很多 window；
每个 window 用 restraint 把体系拉在某个目标几何附近；
每个 window 单独跑 MD；
最后把所有 biased distributions 合并，去掉 bias，得到 unbiased PMF。
```

在本项目中一个 window 是什么：

```text
一个特定目标几何，例如 SerO-C、C-Oleaving、SerH-His 等距离的组合；
一个对应的 umbrella.RST restraint 文件；
一条 md.nc trajectory；
一个 md.out 日志；
一个完成后的 md.rst7 restart。
```

为什么需要 bridge/refine/gap windows：

```text
相邻 window 的采样分布必须有 overlap。
如果两个 window 之间没有重叠，MBAR 无法可靠判断两者自由能差。
所以 attempt013/014/015/016/038/043/045 都是在修补 overlap 或几何 gap。
```

### 2A.12 PMF 是什么

PMF 是 potential of mean force，可以理解为沿某个反应坐标或状态序列的自由能剖面。

PMF 不是单点电子能，也不是单个结构能量。它包含：

```text
反应中心化学能变化；
蛋白环境响应；
溶剂和构象熵；
被选择的 CV/path 上的统计平均。
```

能垒通常从 PMF 上读：

```text
barrier = 最高自由能状态 - reactant reference 自由能
```

本项目去酰化用：

```text
reactant reference = r150
barrier ≈ br208_mid_from_br206 / br208_mid_from_r210 - r150
barrier = 19.061 +/- 0.073 kcal/mol
```

常见误解：

```text
PMF 结果不只看图上最高点。
要同时看窗口健康、overlap、block sensitivity、bias 公式是否正确。
```

### 2A.13 WHAM、MBAR 是什么

WHAM 和 MBAR 都是把多个 umbrella windows 合并成无偏自由能的方法。

WHAM：

```text
更传统，常用于一维 reaction coordinate 的直方图合并。
依赖 binning；bin 选得不好会影响结果。
```

MBAR：

```text
更通用，直接使用各个样本在各个 bias potential 下的 reduced potential。
可以处理多状态、多窗口，通常比 WHAM 更适合复杂 path/state 分析。
```

本项目为什么用 MBAR：

```text
酰化和去酰化都不是简单 1D 反应。
多维 restraint 和 path/state free energies 比单纯 eta 直方图更可靠。
```

MBAR 的前提：

```text
相邻状态要有足够 overlap；
每个窗口 bias potential 要算对；
温度 beta 要一致；
坏窗口不能混进去。
```

### 2A.14 Amber `umbrella.RST` 和 piecewise bias 是什么

Amber 的 `&rst` restraint 不是简单的无限 harmonic restraint。它由 `r1/r2/r3/r4` 和 `rk2/rk3` 定义一个分段势能：

```text
r2 到 r3：平台区，能量为 0。
r1 到 r2：左侧 harmonic 区。
r3 到 r4：右侧 harmonic 区。
r1 外侧、r4 外侧：线性延伸区。
```

为什么这很重要：

```text
如果分析时把所有 restraint 都当成简单 harmonic，会算错 bias energy。
错误的 bias energy 会直接导致错误 PMF，甚至假高能垒。
```

本项目早期去酰化就踩过这个坑，所以 attempt037 开始使用 exact Amber piecewise bias 重新分析。

### 2A.15 Overlap 是什么

overlap 是相邻 umbrella windows 采样分布的重叠程度。

直观理解：

```text
window A 采样范围和 window B 采样范围有交集：MBAR 能比较两者自由能。
两个 window 完全不重叠：中间缺信息，自由能差会不可靠。
```

本项目常用两个诊断：

```text
neighbor_overlap_mbar_matrix.tsv:
  MBAR 意义上的 overlap，更直接关系自由能可靠性。

neighbor_overlap_s_hist.tsv:
  沿路径投影的 histogram overlap，是几何诊断。
```

判断：

```text
overlap < 0.03 且在 barrier region：通常要补窗。
overlap < 0.03 但在远产物侧，且主 barrier 稳定：可暂时不补。
```

当前去酰化：

```text
h119 -> h104 = 0.0198
这是产物侧，不在主 barrier，所以不阻塞去酰化主结论。
```

### 2A.16 Block sensitivity 是什么

block sensitivity 是把 trajectory 分成不同时间块，分别重算 barrier，看结果是否稳定。

本项目最常用：

```text
first_half: 每个窗口前半段帧
second_half: 每个窗口后半段帧
```

如果两半结果接近，说明采样至少在这个尺度上比较稳定。

去酰化结果：

```text
first_half  = 19.136 kcal/mol
second_half = 18.965 kcal/mol
差约 0.17 kcal/mol
```

解释：主 barrier 稳定，没有明显只靠某一半 trajectory 支撑。

### 2A.17 Hysteresis 是什么

hysteresis 是路径依赖或迟滞。比如从 reactant 方向 seed 的窗口和从 product 方向 seed 的窗口，即使 restraint target 接近，采样到的构象也可能不同。

出现 hysteresis 说明：

```text
该区域构象转换慢；
单个窗口采样可能不充分；
需要双向 seed、bridge window、更长 production 或更合适 CV。
```

attempt038 中 midpoint bridge 有 mild hysteresis，但合并 attempt039 后 barrier-region overlap 被修复，所以进入全路径长 production。

### 2A.18 vlimit exceeded、SANDER BOMB、hard error 是什么

`vlimit exceeded` 是 Amber 中速度/位移异常的信号，常见于结构被 restraint 拉坏、timestep 太大、碰撞严重或 QM/MM 力异常。

`SANDER BOMB` 是 Sander 明确停止的严重错误。

本项目 bad window 判据：

```text
NSTEP 未达到目标；
vlimit exceeded > 0；
SCC warning 大量连续出现；
SANDER BOMB / fatal / segmentation / forrtl / NaN；
md.rst7 或 md.nc 缺失；
CV 明显跑偏，不能代表目标 window。
```

注意：

```text
TESTING RELATIVE ERROR
Ewald error estimate
```

这类通常是 Amber 常规诊断，不是 hard error。

### 2A.19 `md.out`、`md.nc`、`md.rst7` 分别看什么

```text
md.out:
  文字日志。看 NSTEP、温度、SCC warning、vlimit、SANDER BOMB、是否正常结束。

md.nc:
  NetCDF trajectory。用于计算 CV distribution、PMF、MBAR。

md.rst7:
  最终 restart。用于继续跑、补窗、延长 production。

umbrella.RST:
  restraint 定义。PMF 分析必须解析它来计算 bias energy。

window_manifest.tsv:
  窗口来源和目标几何。用于追踪每个 window 来自哪个 attempt、哪个 seed。
```

### 2A.20 为什么要保留失败尝试

失败 attempt 不是噪音。它们说明我们为什么没有走某条路。

例子：

```text
attempt014/015/016 说明单纯补 1D eta gap 不能完全解决去酰化 PMF；
attempt026 的 r150 vlimit 说明 reactant reference 必须低 dt 重跑；
attempt037 说明 harmonic-bias 简化会误导 PMF；
attempt046 说明酰化短 production 还不足以给最终 barrier。
```

这些记录能防止后来的人重复同样错误，也能解释为什么最终路线看起来复杂。


## 2B. 从头到尾每一步应该怎么做、看什么、怎么判断

这一节把整个项目按“操作步骤”重写一遍。后面的 attempt 章节是历史复盘；本节是接手者复现时的执行框架。

### Step 1：定义问题和盲法边界

目的：避免把文献答案带入输入。

输入：

```text
PETase 结构；
通用丝氨酸水解酶机理知识；
可用计算资源；
允许参考的方法学。
```

动作：

```text
明确哪些信息不能用：文献 TS、文献窗口、文献坐标、文献能垒。
明确最终要比较什么：我们独立得到的 TS ensemble 和 PMF，与文献结果做后验比较。
```

输出：

```text
blind rule；
反应分支：酰化和去酰化；
初始机制假设。
```

通过标准：

```text
任何用于 seed/window/RC 的具体数值都能从我们自己的结构、化学或计算中解释，而不是来自文献结果。
```

### Step 2：准备蛋白、底物和质子化状态

目的：得到可模拟且化学合理的 PETase-底物体系。

输入：

```text
蛋白结构文件；
底物结构或 SMILES；
催化残基 Ser/His/Asp；
可能活性水。
```

动作：

```text
1. 补全蛋白缺失原子/氢。
2. 判断 His 质子化状态，使其能作为一般碱。
3. 构造底物 3D 结构并参数化。
4. 检查 SerO、HisNE2、Asp 与底物反应原子的编号。
5. 保留候选结构水，但不默认任何水一定是反应水。
```

输出文件应看：

```text
ligand_smiles.tsv
ligand_model_manifest.tsv
stage1_ligand_and_protonation_execution_protocol.md
```

失败时怎么处理：

```text
参数化失败：换底物片段模型或检查电荷/键级。
His 方向不对：重新质子化或旋转 His tautomer。
结构缺失严重：回到结构修复，而不是进入 QM/MM。
```

### Step 3：生成 Michaelis complex 候选

目的：从多个可能结合姿势中找出反应几何合理的候选。

输入：

```text
准备好的蛋白；
底物 3D 构象；
活性位点定义；
几何过滤规则。
```

动作：

```text
1. 以 catalytic Ser/His 附近为 docking box。
2. 生成多个 pose，而不是只保留 docking score 第一。
3. 对每个 pose 测 SerO_C、attack angle、SerH_His、C_Oleave、羰基 O 指向。
4. 对不满足反应几何的 pose 写入 rejected_pose_manifest。
5. 对通过 pose 写入 gs_pose_manifest。
```

输出文件应看：

```text
pose_generation_queue.tsv
gs_pose_manifest.tsv
rejected_pose_manifest.tsv
stage1_pose_geometry_filter_protocol.md
```

通过标准：

```text
SerO 能从合理方向接近羰基 C；
His 能接受 SerH；
离去键仍在；
羰基 O 有 oxyanion-hole 稳定可能；
底物无严重冲突。
```

失败时怎么处理：

```text
所有 pose 都不合理：重新定义 docking box 或底物构象。
pose 结合好但反应几何差：不能作为 Michaelis complex，只能作为非反应结合态。
```

### Step 4：用短 MD/约束松弛诱导稳定反应前态

目的：把静态 docked pose 转为能在蛋白环境中稳定存在的 Michaelis-like conformer。

输入：

```text
gs_pose_manifest 中的候选 pose；
Amber topology；
水盒和离子；
轻约束方案。
```

动作：

```text
1. 最小化水和侧链冲突。
2. 逐步升温到 310 K。
3. 对底物和关键反应距离使用弱 flat-bottom restraint。
4. 运行短 classical MD。
5. 每隔一定帧测反应几何和底物 RMSD/口袋稳定性。
```

输出文件应看：

```text
md_replicate_queue.tsv
productive_conformer_manifest.tsv
stage2_classical_md_protocol.md
```

通过标准：

```text
底物仍在口袋；
反应距离没有完全漂走；
His/Ser/Asp 网络仍合理；
没有蛋白局部结构崩坏；
至少有一个 conformer 能进入 QM/MM。
```

失败时怎么处理：

```text
底物漂走：回到 pose/docking。
关键距离漂远但结合仍稳定：可能需要换反应构象或温和诱导。
口袋结构崩坏：检查参数、质子化或初始 clash。
```

### Step 5：建立机制假设树和候选 CV

目的：不要只假设一条反应坐标；先列出可能机理和变量。

输入：

```text
productive conformers；
通用 serine hydrolase chemistry；
反应原子编号。
```

动作：

```text
1. 为酰化列出 SerO 进攻、SerH->His、C-Oleaving 断裂等事件。
2. 为去酰化列出水进攻、水 H->His、SerO-C 断裂等事件。
3. 为每个事件定义距离、角度、距离差等候选 CV。
4. 把候选 CV 写入 candidate_cv_sets。
```

输出文件应看：

```text
mechanism_hypotheses.yaml
candidate_cv_sets.tsv
stage3_stage4_mechanism_qmmm_protocol.md
```

通过标准：

```text
每个 CV 都有明确化学含义；
CV 覆盖成键、断键、质子转移和进攻角；
没有把单个文献 RC 当成唯一答案。
```

### Step 6：低成本 QM/MM 探索和 TS-like guess

目的：用短 QM/MM scan 或 restrained dynamics 判断哪些路径有可能反应。

输入：

```text
productive conformers；
candidate CV sets；
Amber/Sander QM/MM 输入。
```

动作：

```text
1. 沿候选 CV 做短 restrained scan 或 path screening。
2. 监控是否出现合理成键/断键/质子转移趋势。
3. 排除一拉就崩或直接进入错误产物的路径。
4. 记录 TS-like guess，但不把它直接当 TS。
```

输出文件应看：

```text
path_screening_table.tsv
ts_like_guess_manifest.tsv
generate_stage4_amber_qmmm_inputs.py
```

通过标准：

```text
路径能在 QM/MM 下稳定推进；
关键几何变化符合化学预期；
没有大量 SCC/vlimit/硬错误；
存在可用于 shooting 的 TS-like seed。
```

失败时怎么处理：

```text
一开始就是产物：seed 太后期，回到更早构象。
始终不反应：可能 pose 不对或 CV 不对。
数值爆炸：降低 dt、减弱 restraint 或重新松弛结构。
```

### Step 7：ATESA / aimless shooting 找边界 seed

目的：验证 TS-like guess 是否真的靠近 reactant/product 分界。

输入：

```text
TS-like guess restart；
endpoint classification rules；
Sander QM/MM input；
ATESA settings。
```

动作：

```text
1. 从 seed 随机化速度。
2. 向前/向后或多方向发起短 trajectory。
3. 按几何规则分类 endpoint。
4. 接受能连接 reactant/product 两侧的 shooting move。
5. 收集 accepted/near-boundary seeds。
```

输出目录例子：

```text
blind_work/05_atesa_deacylation/formal_as_attempt003_midseed_16x2
blind_work/05_atesa_deacylation/formal_as_attempt004_from_accepted3_3x8
```

通过标准：

```text
存在 accepted shooting moves；
endpoint 分类不是全部 product 或全部 reactant；
accepted seed 的几何符合化学定义；
没有严重数值失败。
```

失败时怎么处理：

```text
全部到 product：seed 太后期。
全部到 reactant：seed 太早期。
endpoint 多为 undecided：trajectory 太短或分类规则/seed 不合适。
```

### Step 8：focused committor 验证 TS-core

目的：对少数候选 seed 做更直接的 pB 测试。

输入：

```text
ATESA accepted/near-boundary seeds；
随机速度集合；
endpoint classification rules。
```

动作：

```text
1. 对每个 seed 发起多条短 commitment trajectory。
2. 统计 product/reactant/undecided 数量。
3. 计算 pB_committed 和 pB_all。
4. 选择 pB 接近 0.5 的 seed 作为 TS-core 候选。
```

输出目录例子：

```text
blind_work/07_focused_committor/deacylation_lmax_seed12_8shot_attempt006_shortpaths
blind_work/07_focused_committor/deacylation_seed01_seed04_longcommit_attempt007_2x16
```

通过标准：

```text
pB 不接近 0 或 1；
product/reactant 两侧都有足够样本；
undecided 不占主导；
endpoint 几何分类可信。
```

失败时怎么处理：

```text
pB 偏 1：从更早结构找 seed。
pB 偏 0：从更晚结构找 seed。
undecided 多：延长 commitment 或改 endpoint 判据。
```

### Step 9：从 TS-core 构建 umbrella path

目的：把 TS-core 附近向 reactant 和 product 两侧扩展成一串可采样窗口。

输入：

```text
TS-core restart；
reactant/product 方向 short trajectories；
关键 CV；
umbrella restraint 模板。
```

动作：

```text
1. 按 CV 选取从 reactant 到 product 的 state/window。
2. 每个 window 设置 target distances/angles。
3. 从相邻或对应 state 的 rst7 作为 seed。
4. 复制或生成 umbrella.RST。
5. 先短 production/smoke test，再进入长 production。
```

输出文件应看：

```text
window_manifest.tsv
windows/w*/umbrella.RST
windows/w*/md.out
windows/w*/md.nc
windows/w*/md.rst7
```

通过标准：

```text
每个 window 能跑；
CV 分布围绕 target；
相邻 window 有 overlap；
无 SCC/vlimit/硬错误。
```

### Step 10：补窗、bridge、refine

目的：修复 overlap 或几何断层，而不是盲目加长所有窗口。

输入：

```text
初步 MBAR overlap；
CV distribution；
bad window summary；
相邻窗口 final rst7。
```

动作：

```text
1. 找出弱 overlap pair。
2. 判断弱 overlap 是否在 barrier region。
3. 如果在 barrier region，取两侧 CV target 中点或双 seed 建 bridge window。
4. 如果只是远产物侧且 barrier 稳定，暂时不补。
5. 对数值坏窗，优先低 dt 重跑该窗口。
```

本项目例子：

```text
attempt038 修补 deacylation br206->r210 barrier-region weak overlap。
attempt043/045 修补 acylation attempt042 的 bridge/refine gaps。
```

通过标准：

```text
补窗后 weak overlap 消失或不再影响 barrier；
block sensitivity 稳定；
没有新增坏窗口。
```

### Step 11：长 production

目的：用统一长度的窗口采样得到可报告 PMF。

输入：

```text
已修补 overlap 的 window set；
每个 window final rst7；
对应 umbrella.RST；
正式 production input。
```

动作：

```text
1. 每个 window 跑相同步数，例如 40000 steps。
2. 控制总并发 <=64。
3. 定期扫描 NSTEP、TEMP、SCC、vlimit、hard errors。
4. 不在 production 未完成前报告最终 barrier。
```

本项目例子：

```text
deacyl_attempt040: 38 windows x 40000 steps。
acyl_attempt047: 25 windows x 40000 steps。
```

通过标准：

```text
全部窗口达到目标 NSTEP；
md.nc 和 md.rst7 齐全；
SCC/vlimit/hard error 在可接受范围；
没有明显掉队或中断窗口。
```

### Step 12：CV/error summary

目的：先判断数据能不能进 MBAR，不要把坏窗口直接合并。

输入：

```text
complex.prmtop；
windows/w*/md.nc；
windows/w*/md.out；
windows/w*/umbrella.RST。
```

动作：

```text
1. 从 trajectory 计算每帧 CV。
2. 从 md.out 解析 NSTEP、温度、SCC、vlimit、hard errors。
3. 标记 bad windows。
4. 输出 cv_error_summary.tsv。
```

通过标准：

```text
每个窗口帧数合理；
CV 没有跑飞；
没有未完成窗口；
坏窗数为 0，或坏窗有清晰处理方案。
```

### Step 13：MBAR/PMF 分析

目的：去除 umbrella bias，得到自由能。

输入：

```text
CV time series；
umbrella.RST 解析结果；
温度 310 K；
window/state 列表。
```

动作：

```text
1. 用 exact Amber piecewise bias 计算每帧在每个 window 下的 bias energy。
2. 转成 reduced potential。
3. 运行 MBAR。
4. 输出 state_free_energies.tsv、neighbor_overlap_mbar_matrix.tsv、path/eta PMF。
5. 做 block sensitivity。
```

通过标准：

```text
barrier region overlap 合格；
state barrier 和 block sensitivity 稳定；
path-bin PMF 没有由空 bin 支撑的假峰；
不确定性可解释。
```

### Step 14：最终机理判断和文献后验比较

目的：把计算结果变成可解释的机理论断。

输入：

```text
acylation final PMF；
deacylation final PMF；
TS/committor evidence；
窗口健康和 overlap 诊断；
文献结果，只作为后验比较。
```

动作：

```text
1. 分别报告 acylation 和 deacylation barrier。
2. 比较两者差值和可能速控步。
3. 说明每条结论的数据证据。
4. 再和文献比较相同点/差异点。
5. 对残余风险单独列出：DFTB3 精度、QM 区大小、采样长度、RC 投影、SCC warning 等。
```

通过标准：

```text
结论可追溯到具体文件和具体 attempt；
没有把 preliminary PMF 当最终 PMF；
没有把文献结果当输入；
所有弱点都有明确说明。
```

## 3. 原子编号和主要 CV

Amber 原子编号是 1-based。

### 3.1 去酰化 atom/CV

```text
SerO        1911
His NE2     2992
Acyl C      3843
Acyl O      3844
Water O     3863
Water H2    3865

C_SerO  = distance(3843, 1911)
OW_C    = distance(3863, 3843)
H2_His  = distance(3865, 2992)
eta     = C_SerO - OW_C
angle   = angle(3863, 3843, 3844)
```

解释：

```text
C_SerO  大：SerO-C 酰酶键断裂；小：酰酶键保留。
OW_C    小：水 O 接近羰基 C；大：水未进攻。
H2_His  小：水 H 转给 His；大：质子仍在水上。
eta     用于描述旧键/新键相对进程，但不能单独代表完整反应。
angle   判断水进攻羰基的几何是否合理。
```

### 3.2 酰化 atom/CV

```text
SerO        1911
SerH        1912
His NE2     2992
Carbonyl C  3843
Leaving O   3845

SerO_C      = distance(1911, 3843)
C_Oleave    = distance(3843, 3845)
SerH_SerO   = distance(1912, 1911)
SerH_His    = distance(1912, 2992)
eta         = SerO_C - C_Oleave
```

解释：

```text
SerO_C      小：SerO 正在/已经进攻羰基 C。
C_Oleave    大：离去 C-O 键正在/已经断裂。
SerH_His    小：SerH 已转移给 His。
SerH_SerO   大：SerH 离开 SerO。
eta         描述亲核成键与离去键断裂的相对进程。
```

## 4. 总目录结构应该怎么读

服务器上最关键的目录层级：

```text
blind_work/01_system_setup/          结构、配体、pose 和早期 setup
blind_work/02_classical_md/          经典 MD 和 productive conformer gate
blind_work/03_mechanism_tree/        机理假设树和候选 CV 集
blind_work/04_qmmm_exploration/      低成本 QM/MM 探索与 path screening
blind_work/05_atesa_deacylation/     去酰化 ATESA / aimless shooting
blind_work/07_focused_committor/     focused committor 和 long commitment
blind_work/09_umbrella_pmf/          umbrella/PMF 正式主线
```

GitHub 上不会保存原始大文件，只保存说明、表格、脚本和路径。服务器目录才有 `md.nc`、`md.rst7`、`md.out`。

## 5. setup 阶段：从结构到可计算体系

这一阶段对应 ledger 里的：

```text
stage1_structure_pose
amber_protocol_choice
```

目的：不依赖文献坐标，自己从 PETase 结构和通用反应几何出发构造候选复合物。

主要工作：

```text
1. 准备 PETase 结构和底物模型。
2. 生成或筛选底物 pose。
3. 用几何规则过滤 active-site pose。
4. 做短 MD 或构象筛选，寻找能保持反应几何的 productive conformer。
5. 建 Amber topology 和 QM/MM 输入。
```

GitHub 可看：

```text
blind_work/01_system_setup/stage1_remote_execution_instructions.md
blind_work/01_system_setup/stage1_ligand_and_protonation_execution_protocol.md
blind_work/01_system_setup/stage1_pose_geometry_filter_protocol.md
blind_work/02_classical_md/stage2_classical_md_protocol.md
```

这些文件解释早期结构、配体、pose 和 classical MD gate 的规则。它们不是最终能垒结果，但决定后面 QM/MM 是否从一个合理体系开始。

## 6. 去酰化路线完整复盘

去酰化主线已经走到最终 PMF，是目前最完整的一条线。

### 6.1 ATESA 和 seed 寻找：attempt003/004

```text
deacyl_attempt003:
  formal ATESA midseed aimless shooting。
  目的：测试中间 seed 是否能产生 accepted shooting move。
  结果：产生了一些 accepted/near-boundary seeds，但还不足以直接作为 TS ensemble。

deacyl_attempt004:
  从 attempt003 accepted seeds 继续 focused ATESA。
  目的：提高 shooting 效率。
  结果：seed pool 改善，但仍不能直接进入最终 PMF。
```

这一步的科学意义：我们不是把某个 C-O 键长强行设成 TS，而是先看它在动力学上是否能走向两侧 basin。

### 6.2 focused committor：attempt006/007

```text
deacyl_attempt006:
  focused committor short paths。
  目的：估计 seed 候选的 pB，找 TS-core。
  结果：seed 1 和 seed 4 成为重要边界候选。

deacyl_attempt007:
  对 seed 1 和 seed 4 做 long commitment。
  目的：检查 short committor 是否因为路径太短而误分 basin。
  结果：支持 deacylation TS-core region，但个别 trajectory 很慢或偏后期。
```

决策：long commitment 支持可用 TS-core 后，不再无限等待病态慢轨迹，而是转向 umbrella/PMF 建路。

### 6.3 初始 PMF 中心/桥接尝试：attempt012-016

```text
deacyl_attempt012:
  4 个 central PMF production windows，5000 步。
  目的：测试中央 umbrella windows 是否数值稳定。
  结果：能跑，但路径覆盖不足。

deacyl_attempt013:
  w03-w04 bridge windows。
  目的：修补 central windows 之间的 overlap。
  结果：改善一部分，但 eta 0.70-0.90 仍有 gap。

deacyl_attempt014:
  eta 0.75 和 0.85 gap windows。
  目的：针对低 overlap 区补窗。
  结果：仍不足以稳定 projected PMF。

deacyl_attempt015:
  flank windows。
  目的：一次多跑相邻窗口，提高效率。
  结果：覆盖改善但中心 gap 仍脆弱。

deacyl_attempt016:
  eta 0.70-0.90 strong seeded gap windows。
  目的：测试更强 restraint 是否能修复假峰/低权重。
  结果：帮助诊断出 1D eta 不够，需要多维 path treatment。
```

这一段的经验：不是所有 PMF 问题都能靠单个 eta 坐标补窗解决。反应涉及亲核进攻、旧键断裂、质子转移和角度，低维投影可能制造假峰。

### 6.4 扩展正式路径：attempt026/027/029/031/033

```text
deacyl_attempt026:
  20-window 310 K A23/A25 production。
  目的：建立更宽的 formal 310 K path。
  结果：大部分可用，但 r150 窗口有严重 vlimit。

deacyl_attempt027:
  r150 backup low-dt 1 window，40000 步。
  目的：替换坏的 r150 reactant-side reference。
  结果：成功，无 vlimit；作为 reactant reference。

deacyl_attempt029:
  refined bridge 8 windows。
  目的：修复 bridge region 的几何/overlap。
  结果：8/8 可用。

deacyl_attempt031:
  H-transfer bridge 8 windows。
  目的：覆盖亲核进攻后的质子转移分支。
  结果：8/8 可用。

deacyl_attempt033:
  h139/hb130 one-window patch。
  目的：测试一个具体 H-transfer gap。
  结果：无错误完成。
```

关键点：r150 是去酰化 PMF 的 reactant reference。reference 窗口坏了会污染整个相对自由能，所以先用低 dt 重跑它是必要的。

### 6.5 发现并修正 PMF 分析错误：attempt037

```text
deacyl_attempt037:
  Amber-piecewise 310 K path MBAR。
  目的：用 exact Amber piecewise restraint energy 重新分析。
  输入：attempt026/027/029/031/033。
  结果：能垒约 19.27 kcal/mol，但 br206->r210 overlap 仍弱。
```

这里有一个重要方法学修正：早期 PMF 用了过度简化的 harmonic bias 和旧温度假设，可能产生假高能垒。正式分析必须解析每个 `umbrella.RST`，使用 Amber `&rst` 的 piecewise 势能。

Amber piecewise bias 逻辑：

```python
def rst_energy(x, b):
    r1, r2, r3, r4 = b["r1"], b["r2"], b["r3"], b["r4"]
    rk2, rk3 = b["rk2"], b["rk3"]
    e = np.zeros_like(x, dtype=float)
    lo = x < r1
    lq = (x >= r1) & (x < r2)
    uq = (x > r3) & (x <= r4)
    hi = x > r4
    e[lq] = rk2 * (x[lq] - r2) ** 2
    e[uq] = rk3 * (x[uq] - r3) ** 2
    e[lo] = rk2 * (r1 - r2) ** 2 + 2 * rk2 * (r1 - r2) * (x[lo] - r1)
    e[hi] = rk3 * (r4 - r3) ** 2 + 2 * rk3 * (r4 - r3) * (x[hi] - r4)
    return e
```

### 6.6 修补 barrier-region overlap：attempt038/039

```text
deacyl_attempt038:
  br206-r210 midpoint bridge two seeds。
  目的：修复 attempt037 发现的 barrier-region weak overlap。
  结果：两个 midpoint windows 完成，有 mild hysteresis 但可用。

deacyl_attempt039:
  合并 attempt037 + attempt038 preliminary MBAR。
  目的：确认 midpoint bridge 是否修复 br206->r210。
  结果：barrier 约 19.57 kcal/mol；br206->r210 弱连接被修复；无 <0.03 弱连接。
```

决策：因为 barrier-region overlap 被修复，下一步不是继续小补窗，而是全路径统一长 production。

### 6.7 最终长 production 和分析：attempt040/041

最终 production：

```text
/Dell/Dell14/bianht/petase_blind_qmmm/repo/project01/phase2_blind_petase_qmmm_20260630/blind_work/09_umbrella_pmf/deacylation_seed1_attempt040_merged037038_38win_longprod_40000
```

最终分析：

```text
/Dell/Dell14/bianht/petase_blind_qmmm/repo/project01/phase2_blind_petase_qmmm_20260630/blind_work/09_umbrella_pmf/deacylation_seed1_attempt041_from040_longprod_pymbar
```

结果：

```text
38/38 windows completed
76000 frames total, 2000 per window
Bad windows: 0
SCC warnings: 0
vlimit exceeded: 0
Hard errors: 0
State barrier from r150: 19.061 +/- 0.073 kcal/mol
Barrier state: br208_mid_from_br206 / br208_mid_from_r210
Block sensitivity:
  first half  = 19.136 +/- 0.103 kcal/mol
  second half = 18.965 +/- 0.103 kcal/mol
Weak MBAR overlap: h119 -> h104 = 0.0198
```

最终判断：

```text
去酰化主能垒稳定在约 19.1 kcal/mol。
唯一弱 overlap 在产物侧 h119->h104，不在 barrier region。
因此不优先补该窗口；主线转向酰化长 production。
```

## 7. 去酰化结果应该看哪些文件

进入最终分析目录：

```bash
cd /Dell/Dell14/bianht/petase_blind_qmmm/repo/project01/phase2_blind_petase_qmmm_20260630/blind_work/09_umbrella_pmf/deacylation_seed1_attempt041_from040_longprod_pymbar
```

### 7.1 `summary.md`

这是最终摘要。看这些字段：

```text
state_barrier_from_r150_kcal
state_barrier_uncertainty_kcal
state_barrier_label
bad_windows_count
mbar_neighbor_weak_lt0p03
weak_mbar_links
```

当前关键值：

```text
state_barrier_from_r150_kcal = 19.061
state_barrier_uncertainty_kcal = 0.073
state_barrier_label = br208_mid_from_br206
bad_windows_count = 0
weak_mbar_links = h119->h104:0.0198
```

### 7.2 `state_free_energies.tsv`

这是主能垒判断文件。关键列：

```text
state
label
C_target
OW_target
H_target
F_rel_r150_kcal
dF_rel_r150_kcal
```

看法：按 `F_rel_r150_kcal` 从大到小排序，最高点就是 state free energy 给出的 barrier。

### 7.3 `block_sensitivity_state_barrier.tsv`

用于判断采样稳定性。当前：

```text
first_half  = 19.136 kcal/mol
second_half = 18.965 kcal/mol
```

两半差约 0.17 kcal/mol，主峰位置一致，所以去酰化 PMF 主能垒稳定。

### 7.4 `neighbor_overlap_mbar_matrix.tsv` 和 `neighbor_overlap_s_hist.tsv`

用于判断是否需要补窗。规则：

```text
弱 overlap 在 barrier region: 必须补窗。
弱 overlap 在产物侧且主 barrier 稳定: 不优先补。
```

当前最弱是 `h119 -> h104 = 0.0198`，属于产物侧，不阻塞主结论。

### 7.5 `attempt040_cv_error_summary.tsv`

用于判断每个 window 是否坏。关键列：

```text
N_frames
last_nstep
scc_warn
vlimit
hard_errors
ended
C_mean
OW_mean
H_mean
eta_mean
angle_mean
```

坏窗判断：

```text
last_nstep < 40000
或 scc_warn > 0
或 vlimit > 0
或 hard_errors > 0
或 ended = False
```

当前坏窗数为 0。

## 8. 酰化路线完整复盘

酰化还没有最终 PMF，但 TS-core 和短 PMF/补窗路线已经完成，当前主线是 attempt047 长 production。

### 8.1 TS-core 候选 cand02

```text
acyl_committor_cand02:
  目的：验证 acylation TS candidate 是否接近 pB 0.5。
  结果：32 shots 中 forward/product-like 18，reactant-like 13，other 1。
  pB_committed = 0.580645
  pB_all       = 0.5625
```

解释：cand02 不是因为它“长得像文献 TS”被选中，而是因为 committor 支持它位于 reactant/product 分界附近。

### 8.2 初始 4D path：attempt040/041/042

```text
acyl_attempt040:
  7-window 4D path smoke test，5000 步。
  目的：确认 4D path 数值上能跑。
  结果：可行。

acyl_attempt041:
  13-window 4D path stage，10000 步。
  目的：准备 production seeds。
  结果：生成 final restarts。

acyl_attempt042:
  13-window 4D path production，20000 步。
  目的：第一版主要酰化 PMF。
  结果：13/13 完成，但有若干 weak overlap。
```

这里用 4D path，因为酰化不只是一个键断裂/成键，还包括 SerH 向 His 的质子转移。

### 8.3 补桥和 refine：attempt043/044/045/046

```text
acyl_attempt043:
  7 bridge windows，20000 步。
  目的：填补 attempt042 的 weak overlap gaps。
  结果：7/7 完成。

acyl_attempt044:
  合并 042/043 preliminary MBAR。
  目的：检查 bridge 是否修复 PMF。
  结果：中段和产物侧仍需 refine。

acyl_attempt045:
  5 refine windows，20000 步。
  目的：补 m197-mid185 和 m172-bm054 等 gap。
  结果：5/5 完成。

acyl_attempt046:
  合并 042/043/045 preliminary MBAR/overlap。
  目的：诊断长 production 前的窗口质量。
  结果：仍有 reactant-side 和 product-side 弱 overlap；eta PMF 只能作 preliminary。
```

决策：此时最大不确定性不是再补一个小窗口，而是酰化整体只有 20000-step 短 production，和去酰化 40000-step 长 production 不对等。所以启动 attempt047。

### 8.4 当前长 production：attempt047

目录：

```text
/Dell/Dell14/bianht/petase_blind_qmmm/repo/project01/phase2_blind_petase_qmmm_20260630/blind_work/09_umbrella_pmf/acylation_cand02_attempt047_merged042043045_25win_longprod_40000
```

设计：

```text
25 windows
40000 steps/window
MAX_WORKERS = 25
Seeds 来自 attempts042/043/045 final restarts
每个窗口复制其 source umbrella.RST
输入文件: input/acyl_4d_path_longprod_40000.in
```

关键输入文件：

```text
complex.prmtop
window_manifest.tsv
input/acyl_4d_path_longprod_40000.in
windows/w*/umbrella.RST
windows/w*/md.nc
windows/w*/md.rst7
```

完成后必须做：

```text
1. CV/error summary。
2. 解析每个 umbrella.RST。
3. 用 exact Amber piecewise bias 计算 reduced potential。
4. 在 310 K 下跑 MBAR。
5. 输出 state free energies、eta/path diagnostic PMF、neighbor overlap。
6. 判断 barrier estimate 和 uncertainty。
7. 与 deacylation attempt041 做双反应对比。
```

## 9. 酰化实时/完成后应该看哪些文件

进入 attempt047 目录：

```bash
cd /Dell/Dell14/bianht/petase_blind_qmmm/repo/project01/phase2_blind_petase_qmmm_20260630/blind_work/09_umbrella_pmf/acylation_cand02_attempt047_merged042043045_25win_longprod_40000
```

### 9.1 进程和完成情况

```bash
ps -u "$USER" -o pid,ppid,stat,etime,cmd | grep -E 'run_acyl_25win_longprod_40000.py|sander' | grep -v grep
wc -l umbrella_run_progress.tsv
find windows -name md.out | wc -l
find windows -name md.nc | wc -l
find windows -name md.rst7 | wc -l
```

运行中通常看到：

```text
25 个 sander
1 个 run_acyl_25win_longprod_40000.py driver
md.out = 25
md.nc = 25
md.rst7 可能尚未生成，完成后应为 25
```

### 9.2 `umbrella_run_progress.tsv`

窗口完成后会写入。关键列：

```text
state
label
status
returncode
seconds
job_dir
out
rst7
nc
```

判断：

```text
status = done 且 returncode = 0：窗口正常完成。
status = failed 或 returncode != 0：检查该窗口 md.out/run.log。
```

### 9.3 `windows/w*/md.out`

看这些：

```text
NSTEP 是否到 40000
TEMP(K) 是否在 310 K 附近波动
SCC convergence warning 是否出现
vlimit exceeded 是否出现
SANDER BOMB / segmentation / forrtl 等硬错误是否出现
```

注意不要误判：

```text
TESTING RELATIVE ERROR
Ewald error estimate
```

这类是 Amber 常规数值检查，不等于硬错误。

### 9.4 `window_manifest.tsv`

关键列：

```text
state
label
source_attempt
source_window_index
SerO_C_target
C_Oleave_target
SerH_SerO_target
SerH_His_target
force_bond
force_h
seed_rst7
rst_src
```

用途：追踪每个 long-production window 来源于哪个短 production 或 bridge/refine window。

## 10. 如何从服务器复查去酰化最终结论

```bash
cd /Dell/Dell14/bianht/petase_blind_qmmm/repo/project01/phase2_blind_petase_qmmm_20260630/blind_work/09_umbrella_pmf/deacylation_seed1_attempt041_from040_longprod_pymbar
cat summary.md
head state_free_energies.tsv
cat block_sensitivity_state_barrier.tsv
sort -k7,7g neighbor_overlap_mbar_matrix.tsv | head
```

你应该看到的核心结论：

```text
barrier around 19.061 kcal/mol
uncertainty around 0.073 kcal/mol
bad_windows_count = 0
only weak link h119->h104 on product side
```

如果和这个不一致，优先检查：

```text
是不是读错了 analysis 目录；
是不是用了旧的 harmonic-bias PMF；
是不是没有合并 attempt038 midpoint bridge；
是不是温度或 beta 设置不是 310 K；
是不是漏掉 attempt027 r150 backup。
```

## 11. 如何复现 attempt ledger 的逻辑

GitHub 表格：

```text
project01/phase2_blind_petase_qmmm_20260630/blind_work/00_run_status/petase_qmmm_attempt_ledger_20260706.tsv
```

关键列：

```text
attempt
branch
purpose
why_this_was_done
input_or_seed
output_dir_or_file
outcome
next_decision
```

读法：

```text
只看 setup:      过滤 branch = setup
只看去酰化:     过滤 branch = deacylation
只看酰化:       过滤 branch = acylation
找最终去酰化:   deacyl_attempt040 和 deacyl_attempt041
找当前酰化:     acyl_committor_cand02 到 acyl_attempt047
找失败原因:     看 outcome 和 next_decision
```

这张表的作用不是替代本说明，而是提供每一步的最简审计 trail。

## 12. 为什么不能只看 projected path-bin PMF

在窗口覆盖不均匀、某些 path bins 为空或很稀疏时，投影 PMF 容易被采样空洞和 binning 误导。因此当前正式判断优先级是：

```text
1. state free energies 和不确定性；
2. block sensitivity；
3. neighbor MBAR overlap；
4. CV/error summary；
5. path-bin 或 eta PMF 作为诊断图，不单独作为最终 barrier。
```

如果 state barrier 稳定但 path-bin PMF 有奇怪尖峰，先检查 occupied bins、overlap 和 bias 公式，不要直接相信尖峰。

## 13. 什么情况下需要补窗

需要补窗：

```text
barrier region 的相邻窗口 MBAR overlap < 0.03；
block sensitivity 显示前半/后半 barrier 差异大；
某窗口 NSTEP 未完成或有 vlimit/SCC/hard error；
CV distribution 跑偏，窗口没有采到 restraint 附近；
反应关键几何断裂，例如质子转移窗口中 H 没有保持在合理通道。
```

暂时不补窗：

```text
弱 overlap 在远离 barrier 的产物侧；
state barrier 和 block sensitivity 已稳定；
补窗不能减少当前最大不确定性；
另一条反应线还缺少同等级别 production。
```

去酰化 `h119 -> h104` 就属于暂时不补窗。

## 14. 现在还不能下的结论

在 acylation attempt047 完成并分析前，不应该说：

```text
酰化最终能垒是多少；
哪一步一定是速控步；
我们的 acylation/deacylation 差值已经与文献一致或不一致；
酰化 PMF 的 product-side/reactant-side weak overlap 已经完全解决。
```

可以说：

```text
去酰化当前长 production 结果约 19.1 kcal/mol；
酰化 cand02 是经过 committor 支持的 TS-core 候选；
酰化短 PMF 已完成多轮 bridge/refine，但最终能垒等待 attempt047。
```

## 15. 与文献比较应该放在最后

文献可作为方法学参考，但最终比较只能在我们自己的 TS/PMF 完成后进行。比较时应分开写：

```text
方法是否相似：Amber/Sander QM/MM、DFTB3、310 K、committor、PMF。
输入是否独立：没有使用文献 TS、文献窗口或文献坐标作为 seed。
结果是否接近：比较 barrier 和 acylation/deacylation 差值。
差异原因：QM region、timestep、初始构象、水分子选择、CV/path、采样长度、bias 形式。
```

不要把“结果接近文献”写成复现成功的唯一证据。更重要的是每一步是否从独立结构和动力学验证出发。

## 16. 最小复现检查清单

如果别人接手，只想确认这套结果是否站得住，按这个顺序看：

```text
1. GitHub README 和本中文说明，理解主线。
2. petase_qmmm_attempt_ledger_20260706.tsv，确认每个 attempt 的目的和转向。
3. 服务器 deacylation attempt041 summary.md，确认 19.061 kcal/mol。
4. 服务器 deacylation state_free_energies.tsv，确认最高 state。
5. 服务器 deacylation block_sensitivity_state_barrier.tsv，确认前后半稳定。
6. 服务器 deacylation neighbor_overlap_mbar_matrix.tsv，确认弱 overlap 不在 barrier。
7. 服务器 acylation attempt047，确认是否全部 25 窗口完成。
8. attempt047 完成后运行 acylation CV/error + MBAR。
9. 最后再写 acylation/deacylation 双反应对比。
```

## 17. 一句话总结

当前我们已经用盲法路线完成了去酰化 TS-core 到 PMF 的闭环，得到稳定长 production 能垒约 `19.061 +/- 0.073 kcal/mol`；酰化已经有 committor 支持的 cand02 TS-core 和 25 窗口长 production 主线，最终酰化 barrier 需要等 attempt047 完成后用同一套 exact Amber piecewise MBAR 流程重算。
