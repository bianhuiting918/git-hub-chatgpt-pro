# GitHub bridge workspace

This repository is intended as a shared workspace between Codex and ChatGPT Pro web.

Use it for files that should be visible to both tools. Temporary downloads, installers, model checkpoints, raw QM/MM trajectories, large PLACER outputs, and scratch files should stay out of version control.

## Active research projects

This workspace currently contains two linked enzyme-design research projects prepared for ChatGPT Pro + Codex collaboration.

The current design is a **true-TS teacher → predicted-TS student** workflow:

1. [`projects/01-specialized-ts-aware-scorer`](projects/01-specialized-ts-aware-scorer/README.md)  
   **课题一 / true-TS teacher**：针对具体酶或具体酶家族，例如丝氨酸水解酶、PETase、CALB 等，在已知或可由 QM/MM、DFT、约束优化确认的真实 TS（过渡态）和 GS（基态）条件下，结合 PLACER 构象筛选、QM/MM 构象采集和势垒计算，训练 TS-aware barrier scorer，预测 `ΔG‡`、`ΔΔG‡` 和突变体/设计体排序。

2. [`projects/02-general-enzyme-prediction`](projects/02-general-enzyme-prediction/README.md)  
   **课题二 / predicted-TS student**：推理时不使用真实 TS，而是使用 TS 构象预测模型产生的 TS geometry / TS embedding / TS prior，再经 PLACER 筛选构象，学习逼近课题一 true-TS teacher 的势垒预测和机制分解。第一阶段目标是某一催化类型内的 computed barrier proxy / catalytic potential ranking，而不是直接预测实验 `kcat` 或 `kcat/KM`。

## Boundary between the two projects

```text
Project 01:
  GS + true/refined TS + QM/MM or DFT labels
  → high-confidence true-TS barrier scorer

Project 02:
  GS + reaction prior + predicted TS embedding/geometry + PLACER screening
  → practical predicted-TS catalytic-potential predictor
```

Project 01 is the high-confidence teacher and upper-bound model. Project 02 is the practical inference model for cases where a true TS is not available.

## Codex working convention

- Treat each project directory as an independent work package.
- Read the project `README.md` before editing code or adding files.
- Follow the project `CODEX_TASKS.md` files for implementation steps.
- Do not commit temporary downloads, model checkpoints, large raw datasets, PLACER ensemble dumps, QM/MM trajectories, installers, or scratch notebooks unless a project document explicitly requests a small synthetic example artifact.
