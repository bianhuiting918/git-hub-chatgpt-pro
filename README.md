# GitHub bridge workspace

This repository is intended as a shared workspace between Codex and ChatGPT Pro web.

Use it for files that should be visible to both tools. Temporary downloads, installers, and scratch files should stay out of version control.

## Active research projects

This workspace currently contains two enzyme-design research projects prepared for ChatGPT Pro + Codex collaboration:

1. [`projects/01-specialized-ts-aware-scorer`](projects/01-specialized-ts-aware-scorer/README.md)  
   专用酶项目：在已知 TS（过渡态）/GS（基态）的固定反应体系内，学习酶突变体或生成酶对 TS 的稳定能力，并输出 TS-GS 能垒或活性排序。

2. [`projects/02-general-enzyme-prediction`](projects/02-general-enzyme-prediction/README.md)  
   通用酶预测项目：在没有人工指定 TS 作为必需输入的情况下，利用蛋白 embedding（嵌入表示）、构象集合和反应先验，预测不同酶-底物体系的反应活性潜力。

## Codex working convention

- Treat each project directory as an independent work package.
- Read the project `README.md` before editing code or adding files.
- Follow the project `CODEX_TASKS.md` files for implementation steps.
- Do not commit temporary downloads, model checkpoints, large raw datasets, installers, or scratch notebooks unless a project document explicitly requests a small example artifact.
