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
