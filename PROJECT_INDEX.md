# Enzyme reaction learning project index

Created: 2026-06-28

This repository contains two related but separate research projects. They should be developed in parallel but not mixed into one code path too early.

## Project 01 — Specialized TS-aware enzyme scorer

Path: `projects/01-specialized-ts-aware-scorer/`

Purpose: build a high-precision scorer for a fixed reaction/enzyme family where TS（transition state，过渡态）and GS（ground state，基态）are known or can be constructed consistently.

Core use case:

- PETase（聚对苯二甲酸乙二醇酯水解酶）
- PUase / polyurethane hydrolase（聚氨酯水解酶）
- nylonase（尼龙酶）
- laccase（漆酶，后期作为氧化还原体系）
- CALB（Candida antarctica lipase B，南极假丝酵母脂肪酶B，可作为水解酶基准）

Main question:

> Given a known GS/TS reaction model, can we judge whether a generated enzyme or mutant stabilizes TS relative to GS?

This project may use TS as an input because the target reaction is fixed and the goal is mechanism-conditioned evaluation and reranking.

## Project 02 — General enzyme prediction model

Path: `projects/02-general-enzyme-prediction/`

Purpose: build a more general predictor that uses protein embedding（蛋白嵌入表示）, ligand/complex conformation（配体/复合物构象）, and reaction prior（反应先验）to estimate activity without requiring a manually supplied TS complex at inference time.

Main question:

> Can we transfer mechanism-aware signals learned from Project 01 into a broader enzyme activity predictor?

This project should treat explicit TS information as optional training supervision or a generated prior, not as a mandatory inference input.

## Relationship between the two projects

Project 01 is the controlled, high-precision mechanism scorer.

Project 02 is the generalization layer.

The intended path is:

```text
Project 01: known TS/GS → mechanism-aware labels and scorer
        ↓ distillation / transfer
Project 02: GS + protein/ligand/reaction prior → general activity prediction
```

## Shared principle

Do not train a black-box `embedding → activity` model as the primary claim. The central claim is:

> Decompose enzyme activity into reaction-state geometry, electrostatic TS stabilization, conformational ensemble accessibility, and energetic penalties; then learn how these components explain TS-GS barrier changes.
