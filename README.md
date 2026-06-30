# GitHub bridge workspace

This repository is intended as a shared workspace between Codex and ChatGPT Pro web.

Use it for files that should be visible to both tools. Temporary downloads, installers, model checkpoints, raw QM/MM trajectories, large PLACER outputs, and scratch files should stay out of version control.

## Active research projects

This workspace currently contains two linked enzyme-design research projects prepared for ChatGPT Pro + Codex collaboration.

The current design is a **true-TS teacher 鈫?predicted-TS student** workflow:

1. [`projects/01-specialized-ts-aware-scorer`](projects/01-specialized-ts-aware-scorer/README.md)  
   **璇鹃涓€ / true-TS teacher**锛氶拡瀵瑰叿浣撻叾鎴栧叿浣撻叾瀹舵棌锛屼緥濡備笣姘ㄩ吀姘磋В閰躲€丳ETase銆丆ALB 绛夛紝鍦ㄥ凡鐭ユ垨鍙敱 QM/MM銆丏FT銆佺害鏉熶紭鍖栫‘璁ょ殑鐪熷疄 TS锛堣繃娓℃€侊級鍜?GS锛堝熀鎬侊級鏉′欢涓嬶紝缁撳悎 PLACER 鏋勮薄绛涢€夈€丵M/MM 鏋勮薄閲囬泦鍜屽娍鍨掕绠楋紝璁粌 TS-aware barrier scorer锛岄娴?`螖G鈥銆乣螖螖G鈥 鍜岀獊鍙樹綋/璁捐浣撴帓搴忋€?

2. [`projects/02-general-enzyme-prediction`](projects/02-general-enzyme-prediction/README.md)  
   **璇鹃浜?/ predicted-TS student**锛氭帹鐞嗘椂涓嶄娇鐢ㄧ湡瀹?TS锛岃€屾槸浣跨敤 TS 鏋勮薄棰勬祴妯″瀷浜х敓鐨?TS geometry / TS embedding / TS prior锛屽啀缁?PLACER 绛涢€夋瀯璞★紝瀛︿範閫艰繎璇鹃涓€ true-TS teacher 鐨勫娍鍨掗娴嬪拰鏈哄埗鍒嗚В銆傜涓€闃舵鐩爣鏄煇涓€鍌寲绫诲瀷鍐呯殑 computed barrier proxy / catalytic potential ranking锛岃€屼笉鏄洿鎺ラ娴嬪疄楠?`kcat` 鎴?`kcat/KM`銆?

## Boundary between the two projects

```text
Project 01:
  GS + true/refined TS + QM/MM or DFT labels
  鈫?high-confidence true-TS barrier scorer

Project 02:
  GS + reaction prior + predicted TS embedding/geometry + PLACER screening
  鈫?practical predicted-TS catalytic-potential predictor
```

Project 01 is the high-confidence teacher and upper-bound model. Project 02 is the practical inference model for cases where a true TS is not available.

## Codex working convention

- Treat each project directory as an independent work package.
- Read the project `README.md` before editing code or adding files.
- Follow the project `CODEX_TASKS.md` files for implementation steps.
- Do not commit temporary downloads, model checkpoints, large raw datasets, PLACER ensemble dumps, QM/MM trajectories, installers, or scratch notebooks unless a project document explicitly requests a small synthetic example artifact.

## Shared Entries

- [Sugon SCNet GROMACS+CP2K QM/MM entry](docs/scnet-gmx-cp2k.md)
