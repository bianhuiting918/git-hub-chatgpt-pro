# PETase QM/MM 复现记录中文说明

更新日期: 2026-07-06

这份文件是给后续接手者看的中文入口。它解释现在 PETase QM/MM 复现到哪里了、每类数据应该看哪个文件、每个文件的关键列是什么意思、以及下一步该根据什么标准推进。

注意: GitHub 里只保存小型文本记录和可复现路径，不保存服务器密码、原始轨迹 `md.nc`、重启文件 `md.rst7` 或大型日志。

## 1. 当前总状态

服务器项目根目录:

```text
/Dell/Dell14/bianht/petase_blind_qmmm/repo/project01/phase2_blind_petase_qmmm_20260630
```

GitHub 文档入口:

```text
project01/phase2_blind_petase_qmmm_20260630/blind_work/00_run_status/
```

最重要的三个 GitHub 记录文件:

```text
petase_blind_qmmm_reproduction_runbook_20260706.md
petase_qmmm_attempt_ledger_20260706.tsv
petase_qmmm_reproduction_guide_zh_20260706.md
```

当前结论:

```text
去酰化 deacylation:
  已完成长 production 和 MBAR 分析。
  当前能垒 = 19.061 +/- 0.073 kcal/mol。
  结果比较稳定。

酰化 acylation:
  正在跑 attempt047 长 production。
  25 个窗口, 每个 40000 步。
  完成后需要重新做 CV/error 汇总和 MBAR/PMF 分析。
```

## 2. 如果只想快速知道做到哪了

看 GitHub:

```text
project01/phase2_blind_petase_qmmm_20260630/blind_work/00_run_status/petase_blind_qmmm_reproduction_runbook_20260706.md
```

重点看这些小节:

```text
Deacylation Final State
Acylation Current State
Current Barrier Summary
How to Reproduce the Current Results
Decision Logic Used
```

这个文件回答:

- 用的是什么服务器和软件环境;
- 原子编号和 CV 是什么;
- 去酰化最终结果是什么;
- 酰化现在跑到哪个 attempt;
- 为什么不继续补去酰化某些窗口;
- 为什么现在转向酰化长 production。

## 3. 如果想知道每一次 attempt 为什么做

看 GitHub:

```text
project01/phase2_blind_petase_qmmm_20260630/blind_work/00_run_status/petase_qmmm_attempt_ledger_20260706.tsv
```

这是一个 TSV 表。每一行是一个关键 attempt。

关键列:

```text
attempt
  attempt 名称, 例如 deacyl_attempt040, acyl_attempt047。

branch
  属于 setup / deacylation / acylation 哪条线。

purpose
  这一步具体想解决什么问题。

why_this_was_done
  为什么当时要做这一步, 也就是决策依据。

input_or_seed
  输入来自哪里, 例如哪个 attempt 的 final restart 或 accepted seed。

output_dir_or_file
  服务器上的输出目录或文件名。

outcome
  这一步结果如何, 成功、失败、部分成功或发现了什么问题。

next_decision
  基于这个结果, 下一步为什么这么做。
```

怎么读:

- 想复盘去酰化路线, 过滤 `branch = deacylation`。
- 想复盘酰化路线, 过滤 `branch = acylation`。
- 想找最终结果, 看 `deacyl_attempt040`, `deacyl_attempt041`, `acyl_attempt047`。
- 想知道失败/返工原因, 看 `outcome` 和 `next_decision`。

## 4. 去酰化最终结果看哪些文件

服务器最终 production 目录:

```text
/Dell/Dell14/bianht/petase_blind_qmmm/repo/project01/phase2_blind_petase_qmmm_20260630/blind_work/09_umbrella_pmf/deacylation_seed1_attempt040_merged037038_38win_longprod_40000
```

服务器最终分析目录:

```text
/Dell/Dell14/bianht/petase_blind_qmmm/repo/project01/phase2_blind_petase_qmmm_20260630/blind_work/09_umbrella_pmf/deacylation_seed1_attempt041_from040_longprod_pymbar
```

### 4.1 总结文件

```text
summary.md
```

看什么:

```text
state_barrier_from_r150_kcal
state_barrier_uncertainty_kcal
state_barrier_label
bad_windows_count
mbar_neighbor_weak_lt0p03
weak_mbar_links
```

当前关键值:

```text
state_barrier_from_r150_kcal = 19.061
state_barrier_uncertainty_kcal = 0.073
state_barrier_label = br208_mid_from_br206
bad_windows_count = 0
weak_mbar_links = h119->h104:0.0198
```

解释:

- `state_barrier_from_r150_kcal` 是主要能垒结论。
- `bad_windows_count = 0` 说明 38 个 production 窗口都可用。
- `weak_mbar_links = h119->h104` 是产物侧弱 overlap, 不在主能垒附近, 所以目前不阻塞主线。

### 4.2 各窗口自由能

```text
state_free_energies.tsv
```

关键列:

```text
state
  窗口编号。

label
  窗口标签, 例如 r150_backup, br206, br208_mid_from_br206。

C_target, OW_target, H_target
  该窗口 restraint 的目标几何。

F_rel_r150_kcal
  相对 r150 reactant reference 的自由能。

dF_rel_r150_kcal
  MBAR 不确定性估计。
```

怎么看能垒:

```text
按 F_rel_r150_kcal 从大到小排序。
最高点就是 state free energy 给出的 barrier。
```

当前最高点:

```text
br208_mid_from_br206 / br208_mid_from_r210
F_rel_r150_kcal = 19.060878
```

### 4.3 前半段/后半段稳定性

```text
block_sensitivity_state_barrier.tsv
```

关键列:

```text
block
  first_half 或 second_half。

barrier_label
  这一半数据得到的最高点窗口。

barrier_kcal
  这一半数据得到的能垒。
```

当前结果:

```text
first_half  = 19.136 kcal/mol
second_half = 18.965 kcal/mol
```

解释:

两半差约 0.17 kcal/mol, 主峰位置一致, 说明去酰化 PMF 的主 barrier 比较稳。

### 4.4 窗口 overlap

```text
neighbor_overlap_mbar_matrix.tsv
neighbor_overlap_s_hist.tsv
```

关键列:

```text
label_i, label_j
  相邻两个窗口。

mbar_overlap_sym
  MBAR overlap 对称值。越低说明相邻窗口重叠越差。

hist_overlap
  沿路径投影的直方图 overlap, 作为几何诊断。
```

当前最弱 MBAR overlap:

```text
h119 -> h104 = 0.0198
```

判断:

- 如果弱 overlap 在 barrier 附近, 需要补窗。
- 如果弱 overlap 只在产物侧, 且 barrier/block sensitivity 稳定, 不优先补。
- 当前属于第二种, 所以没有继续补 h119 -> h104。

### 4.5 每个窗口是否坏了

```text
attempt040_cv_error_summary.tsv
```

关键列:

```text
N_frames
  每个窗口帧数。当前应为 2000。

C_mean, OW_mean, H_mean, eta_mean, angle_mean
  该窗口采样到的平均几何。

last_nstep
  最后步数。当前应为 40000。

scc_warn
  SCC 收敛 warning 数。

vlimit
  vlimit exceeded 次数。

hard_errors
  fatal/segmentation/forrtl 等硬错误数。

ended
  md.out 是否正常结束。
```

判断坏窗:

```text
last_nstep < 40000
或 scc_warn > 0
或 vlimit > 0
或 hard_errors > 0
或 ended = False
```

当前:

```text
bad windows = 0
```

## 5. 酰化现在看哪些文件

当前正在跑的服务器目录:

```text
/Dell/Dell14/bianht/petase_blind_qmmm/repo/project01/phase2_blind_petase_qmmm_20260630/blind_work/09_umbrella_pmf/acylation_cand02_attempt047_merged042043045_25win_longprod_40000
```

### 5.1 运行进度

```text
umbrella_run_progress.tsv
```

如果正在运行, 这个文件可能暂时只有表头或为空。窗口完成后才会逐行写入。

关键列:

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

判断:

```text
status = done, returncode = 0  表示窗口正常完成。
status = failed 或 returncode != 0 表示需要检查该窗口。
```

### 5.2 每个窗口实时状态

窗口目录:

```text
windows/w000_r260/
windows/w001_r247/
...
windows/w024_prod138/
```

每个窗口主要看:

```text
md.out
md.nc
md.rst7
umbrella.RST
run.log
```

文件含义:

```text
md.out
  Sander 输出。看 NSTEP、温度、SCC warning、vlimit、fatal error。

md.nc
  轨迹。完成后用于计算 CV 和 MBAR。

md.rst7
  最终 restart。完成后用于后续补窗或延长 production。

umbrella.RST
  该窗口的 restraint 定义。分析 PMF 时必须解析它, 不能只用简单 harmonic 近似。

run.log
  启动命令和标准输出/错误。
```

### 5.3 输入设置

```text
input/acyl_4d_path_longprod_40000.in
window_manifest.tsv
```

`window_manifest.tsv` 关键列:

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

解释:

- `SerO_C_target` 和 `C_Oleave_target` 描述成键/断键方向。
- `SerH_SerO_target` 和 `SerH_His_target` 描述质子转移方向。
- `seed_rst7` 是该窗口从哪个已有 final restart 开始。
- `rst_src` 是该窗口 restraint 从哪里复制来的。

### 5.4 attempt047 完成后应生成的新分析

建议输出目录:

```text
/Dell/Dell14/bianht/petase_blind_qmmm/repo/project01/phase2_blind_petase_qmmm_20260630/blind_work/09_umbrella_pmf/acylation_cand02_attempt048_from047_longprod_mbar
```

应生成这些文件:

```text
attempt047_cv_error_summary.tsv
state_free_energies.tsv
neighbor_overlap_mbar_matrix.tsv
neighbor_overlap_eta_hist.tsv
eta_pmf.tsv
block_sensitivity_state_barrier.tsv
summary.md
```

分析逻辑:

1. 用 `md.nc` 和 `complex.prmtop` 计算 4 个 CV:
   - SerO_C
   - C_Oleave
   - SerH_SerO
   - SerH_His
2. 从 `umbrella.RST` 解析每个窗口真实 Amber piecewise bias。
3. 用 310 K 的 beta 做 MBAR。
4. 先看 `state_free_energies.tsv` 得到主 barrier。
5. 再看 overlap 判断是否需要补窗。
6. eta PMF 只作为诊断, 如果有空 bin 或弱 overlap, 不直接当最终 barrier。

## 6. 如何判断下一步是否补窗

补窗不是看到弱 overlap 就立刻做。按下面顺序判断:

### 情况 A: 必须补窗

```text
弱 overlap 在 barrier 附近;
或 block sensitivity 前后半差异很大;
或最高点位置在相邻窗口之间跳动;
或窗口有 SCC/vlimit/hard error;
或某些窗口未到目标 NSTEP。
```

动作:

```text
在弱 overlap 两侧目标 CV 的中点增加 1-3 个 bridge windows。
从两侧 final rst7 分别 seed。
必要时降低 dt 或增强 restraint。
```

### 情况 B: 暂时不补窗

```text
弱 overlap 在产物侧或反应物远端;
barrier 区 overlap 合格;
block sensitivity 稳定;
主 barrier 位置不变。
```

动作:

```text
先汇总主结论, 不消耗 CPU 去补不影响 barrier 的窗口。
```

当前去酰化就是情况 B。

### 情况 C: 需要改 RC

```text
窗口 overlap 看似可以, 但 path-bin PMF 出现巨大假峰;
或 eta 投影有空 bin;
或不同几何变量显示明显滞后/迟滞。
```

动作:

```text
不要只加密 1D eta 窗口。
改用 state free energies 或 3D/4D path coordinate。
检查 Amber piecewise bias 是否正确。
```

早期 deacylation 就遇到过这个问题, 后来改用 exact Amber piecewise bias 和 state free energies 作为主判据。

## 7. 监控命令模板

在服务器上运行:

```bash
cd /Dell/Dell14/bianht/petase_blind_qmmm/repo/project01/phase2_blind_petase_qmmm_20260630/blind_work/09_umbrella_pmf/acylation_cand02_attempt047_merged042043045_25win_longprod_40000

date
ps -u "$USER" -o pid,ppid,stat,etime,cmd | grep -E 'run_acyl_25win_longprod_40000.py|sander' | grep -v grep
wc -l umbrella_run_progress.tsv
find windows -name md.out | wc -l
find windows -name md.nc | wc -l
find windows -name md.rst7 | wc -l
```

检查错误:

```bash
grep -R -i "SCC\\|vlimit exceeded\\|fatal\\|segmentation\\|forrtl" windows/*/md.out
```

注意:

```text
Ewald error estimate
TESTING RELATIVE ERROR
```

这两类是 Amber 常见诊断行, 不要当成硬错误。

## 8. 最终应该如何写结果

不要只写一个 barrier 数字。最终结果至少应包含:

```text
1. 结构和反应路径来源:
   是否来自 blind docking/MD/QM-MM, 是否用 paper 坐标。

2. TS/committor 证据:
   哪些 seed 接近 pB 0.5, 哪些被拒绝。

3. PMF 协议:
   QM 区、温度、dt、窗口数、步数、CV、restraint 形式。

4. 窗口健康:
   bad windows, SCC, vlimit, hard errors。

5. overlap:
   barrier 区是否有弱 overlap。

6. barrier:
   acylation barrier, deacylation barrier, uncertainty, block sensitivity。

7. 残余风险:
   例如短采样、QM 区偏小、DFTB3 精度、RC 投影假峰等。
```

当前已可写:

```text
Deacylation barrier = 19.061 +/- 0.073 kcal/mol
```

当前还不能最终写:

```text
Acylation final barrier
Rate-limiting step
Final comparison with paper
```

这些要等 `attempt047` 完成并分析后再写。

