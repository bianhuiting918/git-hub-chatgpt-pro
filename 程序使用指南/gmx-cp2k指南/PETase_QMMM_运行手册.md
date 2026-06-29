# PETase GMX-CP2K/QMMM 运行手册

更新时间：2026-06-29

本文档记录在 CPU 服务器上运行 PETase GMX-CP2K/QMMM 计算的实际流程。目标是让新的 Codex 窗口可以直接恢复、监控、扩展和分析当前计算。

## 1. 服务器和关键路径

CPU 服务器：

```bash
ssh bianht@210.73.40.29
```

常用工作根目录：

```bash
/Dell/Dell14/bianht/gromacs_cp2k
/Dell/Dell14/bianht/enzyme/orbmol-comput/all_pdb_with_kcat
```

GMX-CP2K wrapper：

```bash
/Dell/Dell14/bianht/gmx_cp2k_patched.sh
```

当前 PETase WT 无水 QMMM 目录：

```bash
/Dell/Dell14/bianht/gromacs_cp2k/PETase_WT_QMMM_vacuum_total
```

当前 PETase 突变体批量 QMMM 目录：

```bash
/Dell/Dell14/bianht/gromacs_cp2k/PETase_QMMM_vacuum_total_batch
```

源结构目录：

```bash
/Dell/Dell14/bianht/enzyme/orbmol-comput/all_pdb_with_kcat/aligned_complexes/PETase/step2_minimization
/Dell/Dell14/bianht/enzyme/orbmol-comput/all_pdb_with_kcat/aligned_complexes/PETase/step2_ts_guess
```

## 2. GMX-CP2K 调用方式

不要直接调用宿主机上的 `/data/gromacs-2021.1-plumed-qmmm/build/bin/gmx_mpi_d`，因为它实际在 bwrap/singularity 风格环境里运行。统一使用 wrapper：

```bash
/Dell/Dell14/bianht/gmx_cp2k_patched.sh <gmx-subcommand> [args...]
```

示例：

```bash
/Dell/Dell14/bianht/gmx_cp2k_patched.sh grompp \
  -f em-qmmm-gs-vacuum.mdp \
  -c gs_complex.gro \
  -p topol_gs_vac.top \
  -n index_vac.ndx \
  -o gs_lg3_qmmm_vacuum_total.tpr \
  -maxwarn 8
```

运行 QMMM 单点：

```bash
OMP_NUM_THREADS=4 /Dell/Dell14/bianht/gmx_cp2k_patched.sh mdrun \
  -s gs_lg3_qmmm_vacuum_total.tpr \
  -deffnm gs_lg3_qmmm_vacuum_total \
  -ntomp 4 \
  -nsteps 1
```

## 3. 当前 QMMM 科学口径

当前用于和 OrbMol 比较的 QMMM 口径是：

- 全蛋白 MM。
- 配体 QM。
- 不加水，不加离子。
- GS 配体名：`LG3`。
- TS 配体名：`LG4`。
- QM 区域：配体 23 个原子。
- 最终能量来自 GROMACS 日志里的 `Potential Energy`。

计算势能差：

```text
DeltaE_kJ_mol = Potential_TS - Potential_GS
DeltaE_kcal_mol = DeltaE_kJ_mol / 4.184
DeltaE_eV = DeltaE_kJ_mol / 96.48533212
```

重要：不要用 CP2K 输出中的

```text
ENERGY| Total FORCE_EVAL ( QMMM ) energy [a.u.]
```

作为最终势垒。该值用于判断 CP2K 是否完成，但最终比较应使用 GROMACS/QMMM 的整体 `Potential Energy`。

## 4. WT 单体系结果和提取方法

WT 目录：

```bash
cd /Dell/Dell14/bianht/gromacs_cp2k/PETase_WT_QMMM_vacuum_total
```

检查完成状态：

```bash
cat wt_gs_lg3_qmmm_vacuum_total.exit
cat wt_ts_lg4_qmmm_vacuum_total.exit
grep -E 'Potential Energy|Finished mdrun|Fatal error|Error' \
  wt_gs_lg3_qmmm_vacuum_total.log \
  wt_ts_lg4_qmmm_vacuum_total.log
```

已得到的 WT 结果：

```text
GS Potential Energy = -336953.241152172 kJ/mol
TS Potential Energy = -337726.990758486 kJ/mol
TS-GS = -773.749606314 kJ/mol
TS-GS = -184.930594243 kcal/mol
TS-GS = -8.019349567 eV
```

## 5. 突变体批量队列

批量目录：

```bash
cd /Dell/Dell14/bianht/gromacs_cp2k/PETase_QMMM_vacuum_total_batch
```

主要脚本：

```bash
petase_qmmm_vacuum_batch.py
petase_qmmm_global2_worker.py
monitor_qmmm_progress.sh
```

主要输出：

```bash
qmmm_vacuum_total_results.csv
qmmm_vs_orbmol_merged.csv
qmmm_vs_orbmol_summary.txt
qmmm_vs_orbmol_allatom_step2_merged.csv
qmmm_vs_orbmol_allatom_step2_summary.txt
```

当前调度器 PID：

```bash
cat batch_worker.pid
ps -p "$(cat batch_worker.pid)" -o pid,stat,etime,pcpu,pmem,rss,cmd
```

当前定时监控 PID：

```bash
cat monitor_qmmm_progress.pid
ps -p "$(cat monitor_qmmm_progress.pid)" -o pid,stat,etime,pcpu,pmem,rss,cmd
```

监控日志：

```bash
tail -120 qmmm_progress_timer.log
```

当前推荐调度策略：

- 全局并发 2 个 QMMM 单点起步。
- 确认稳定后可升到 4 个 QMMM 单点。
- 不建议一开始超过 8 个单点。

## 6. 单突变体准备流程

脚本会自动执行下列工作：

1. 从 `step2_minimization` 读取 GS PDB。
2. 从 `step2_ts_guess` 读取 TS PDB。
3. 去掉水，只保留蛋白和配体。
4. 对蛋白运行 `pdb2gmx`，生成 GROMOS54A7 蛋白拓扑。
5. 添加配体 placeholder topology。
6. 生成无水 `gs_complex.gro` 和 `ts_complex.gro`。
7. 动态生成 `index_vac.ndx`，其中 QMatoms 是配体 23 个原子。
8. 动态生成 CP2K input，更新 `MM_INDEX`。
9. 运行 `grompp` 生成 `.tpr`。
10. 运行 `mdrun -nsteps 1` 得到 QMMM 单点势能。

手动准备一个突变体：

```bash
cd /Dell/Dell14/bianht/gromacs_cp2k/PETase_QMMM_vacuum_total_batch
./petase_qmmm_vacuum_batch.py prepare --key IsPETase_M1-A248P
```

手动运行一个突变体：

```bash
./petase_qmmm_vacuum_batch.py run-key IsPETase_M1-A248P
```

手动汇总：

```bash
./petase_qmmm_vacuum_batch.py collect
```

## 7. 原子范围

当前 QMMM 包含：

- MM：全蛋白。
- QM：配体 23 个原子。
- 不含水。

示例：

```text
WT:
  MM atoms = 1..2493
  QM atoms = 2494..2516
  total = 2516

IsPETase_M1-A248P:
  MM atoms = 1..2500
  QM atoms = 2501..2523
  total = 2523

IsPETase_M1-D186H:
  MM atoms = 1..2504
  QM atoms = 2505..2527
  total = 2527
```

QM 配体原子名：

```text
O10 O12 O14 O18
C21 C22 C23 C24 C25 C26 C27 C28 C29 C30
H42 H43 H44 H45 H46 H47 H48 H49 H50
```

## 8. 结果提取

从 GROMACS 日志提取：

```bash
grep 'Potential Energy' gs_lg3_qmmm_vacuum_total.log
grep 'Potential Energy' ts_lg4_qmmm_vacuum_total.log
```

批量结果表：

```bash
cd /Dell/Dell14/bianht/gromacs_cp2k/PETase_QMMM_vacuum_total_batch
column -s, -t qmmm_vacuum_total_results.csv | less -S
```

状态计数：

```bash
awk -F, 'NR>1{c[$2]++} END{for(k in c) print k,c[k]}' qmmm_vacuum_total_results.csv
```

失败检查：

```bash
find runs -name FAILED -print -exec cat {} \;
```

## 9. 和 OrbMol 全原子结果比较

之前的 OrbMol 全原子 step2 结果在：

```bash
/Dell/Dell14/bianht/enzyme/orbmol-comput/all_pdb_with_kcat/results_orb_v3_all_residues/step2_gs_energies.csv
/Dell/Dell14/bianht/enzyme/orbmol-comput/all_pdb_with_kcat/results_orb_v3_all_residues/step2_ts_energies.csv
```

OrbMol 全原子对比口径：

```text
OrbMol all-atom DeltaE = E_total(TS) - E_total(GS)
```

注意两者不是完全同构：

- OrbMol 全原子结果包含全蛋白、配体和 3 个 WAT。
- 当前 QMMM 结果包含全蛋白和配体 QM，不含水。

当前已生成的合并表：

```bash
/Dell/Dell14/bianht/gromacs_cp2k/PETase_QMMM_vacuum_total_batch/qmmm_vs_orbmol_allatom_step2_merged.csv
/Dell/Dell14/bianht/gromacs_cp2k/PETase_QMMM_vacuum_total_batch/qmmm_vs_orbmol_allatom_step2_summary.txt
```

## 10. 资源估算和并发

服务器实时基线示例：

```text
total memory = 2.0 TiB
available memory = about 1.7 TiB
CPU cores = 192
single QMMM RSS = about 24-26 GiB
```

建议：

- 稳妥：2 个 QMMM 单点并发。
- 推荐尝试：4 个 QMMM 单点并发。
- 最高谨慎：6-8 个 QMMM 单点并发。
- 不建议超过 8 个，容易出现 CPU 争用、I/O 压力和 SCF 变慢。

每个突变体需要 GS 和 TS 两个单点。实测：

```text
WT: about 3-3.5 h
IsPETase_M1-A248P: about 4 h 26 min
IsPETase_M1-D186H TS: about 3 h 19 min
```

保守估计一个突变体完整 barrier 需要 4-5 小时；并发数越高，总 wall time 近似按并发缩短，但单点可能因 CPU 争用变慢。

## 11. 常见问题

### 11.1 新窗口找不到 `gmx_mpi_d`

使用 wrapper，不要直接找宿主机二进制：

```bash
/Dell/Dell14/bianht/gmx_cp2k_patched.sh --version
```

### 11.2 `grompp` 报 `No default Proper Dih. types`

当前批处理脚本会自动处理：

- 只删除 `ERROR ... No default Proper Dih. types` 对应 topology 行。
- 保留 `[ molecules ]` 中的配体行。
- 每次删除前会写 topology 备份。

不要手动大范围删除 topology 行。

### 11.3 不能用 CP2K 中间能量做最终比较

如果只看到：

```text
ENERGY| Total FORCE_EVAL ( QMMM ) energy [a.u.]
```

说明 CP2K 已经给出一次 QMMM 能量，但最终对比仍要等 GROMACS 日志写出：

```text
Potential Energy = ...
Finished mdrun
```

### 11.4 旧 worker 和新 worker 冲突

检查：

```bash
ps -ef | grep -E 'petase_qmmm_.*worker|gmx_mpi_d mdrun' | grep -v grep
cat batch_worker.pid
```

只保留一个调度器父进程。正在跑的 `gmx_mpi_d mdrun` 子进程不要随便杀，除非明确要停止计算。

## 12. 快速恢复命令

```bash
cd /Dell/Dell14/bianht/gromacs_cp2k/PETase_QMMM_vacuum_total_batch

date
free -h
uptime

ps -p "$(cat batch_worker.pid)" -o pid,stat,etime,pcpu,pmem,rss,cmd
ps -eo pid,ppid,stat,etime,pcpu,pmem,rss,cmd --sort=-rss | \
  egrep 'gmx_mpi_d|petase_qmmm|global2|mdrun' | grep -v egrep | head -30

tail -120 batch_worker_global2.nohup.log
tail -120 qmmm_progress_timer.log

awk -F, 'NR>1{c[$2]++} END{for(k in c) print k,c[k]}' qmmm_vacuum_total_results.csv 2>/dev/null || true
find runs -name FAILED -print -exec cat {} \;
```

