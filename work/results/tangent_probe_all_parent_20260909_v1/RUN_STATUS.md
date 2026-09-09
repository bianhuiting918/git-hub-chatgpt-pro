# Full-parent tangent-probe run, 2026-09-09

Status at submission: RUNNING. User requested all surface ester/amide accessibility before plotting.

Remote project: /data/bht2/polymer_material_reference_20260827/simulation_slabs/interface_surface_robustness_20260908_v1

CPU process PID at launch: 516427.
Existing Python: /data/bht2/polymer_material_reference_20260827/simulation_slabs/unbiased_100chain_v1/env/bin/python
Existing script SHA256: 1762eda47dce4a5d461143397c9ca2df94a381dde0ccf6afe0b7d64f25a555ea

Command (run from remote project):
```
/data/bht2/polymer_material_reference_20260827/simulation_slabs/unbiased_100chain_v1/env/bin/python -u tangent_probe_v1.py --output tangent_probe_all_parent_20260909_v1 --pilot 8400 --directions 1024
```

The reused CLI calls its requested sample count --pilot. Here 8400 equals/exceeds every material catalogue and selects ALL 8400 PET, 6800 PA6, 6800 PA66 model bonds. PET list identity and uniqueness were verified after launch. Verify PA6/PA66 lists when created. The script's generic PILOT status strings must not be interpreted as a random subsample for this run.

No density/support filter; includes interior bonds so no surface candidate is lost. These 22000 model bonds are NOT the surface denominator. Caps and loose tails still require explicit labels before dense-body interface fitting. Do not call no sampled access a proof of buried location.

Unchanged radii 0.5,1,1.5,2,3,4,5,6,8,10 A, 1024 directions, all atoms retained, exact swept-segment checks, 120 A straight-path cap, XY periodicity, unresolved and right-censoring flags. Six analytic tests passed again before launch. The 64-direction screen changes ordering only in this all-catalogue run; no record is excluded by it.

Output: tangent_probe_all_parent_20260909_v1/
Log: tangent_probe_all_parent_20260909_v1.stdout.log
Existing script appends RUN_LOG.jsonl and RUNBOOK.md. Per-site NPZ files persist before final summaries. Do not overwrite/restart this output directory; inspect process and records before recovery.

Completion gates: full unique coverage per material, ten radius records per site, masks/counts/monotonicity, source identity, prior pilot numerical agreement, unresolved and cap/tail scope. Plot only after these gates. Planned figures merge top/bottom in a carbonyl-centered oriented frame, stratify NEW geometric exposure rather than legacy support or water exposure. No new triad/donor scan is authorized by this phase. No figures or full-run numerical conclusions are claimed in this submission note.
