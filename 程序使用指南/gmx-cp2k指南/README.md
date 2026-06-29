# GMX-CP2K 指南

本目录记录 CPU 服务器上 GROMACS 调 CP2K 进行 PETase QMMM 计算的实际可运行方法。

## 文档入口

- [PETase_QMMM_运行手册.md](./PETase_QMMM_运行手册.md)
  - 环境和路径
  - `gmx_cp2k_patched.sh` 调用方式
  - WT 单体系 QMMM 运行
  - 突变体批量 QMMM 队列
  - 结果提取和 OrbMol 全原子结果比较
  - 并发数和内存估算
  - 常见错误处理

## 当前计算口径

当前 PETase QMMM 批量计算采用：

- 无水体系，不加入溶剂盒。
- 全蛋白作为 MM 区域。
- 配体 `LG3` 或 `LG4` 的 23 个原子作为 QM 区域。
- 最终势能差使用 GROMACS/QMMM 日志里的整体 `Potential Energy`：

```text
DeltaE = Potential Energy(TS) - Potential Energy(GS)
```

不要用 CP2K `ENERGY| Total FORCE_EVAL` 作为最终对比值；该值不是用户当前要求的完整 MM+QM 总势能口径。

