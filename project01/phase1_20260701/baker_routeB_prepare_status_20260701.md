# Baker Route B Prepare Status

Date: 2026-07-01 00:30 CST

## Current Result

The Baker-style serine-hydrolase input package has been prepared on the GPU filesystem:

```text
/data/bht/project01_baker_serhyd_routeB_20260701
```

The route is labeled:

```text
baker_theozyme_new_backbone
```

This is intentionally separate from the older route:

```text
/data/bht/project01_phase1_reset_gpu/denovo_scaffold_routeB/serine_hydrolase
```

The older route has useful scripts and diagnostic evidence, but it is not counted as completion for this Baker-theozyme reset because it used the old `denovo_SER_hydrolase_full_input.pdb` / fixed-pocket style entrance route.

## Remote Evidence

Prepared input manifest:

```text
/data/bht/project01_baker_serhyd_routeB_20260701/manifests/baker_serhyd_routeB_input_manifest.tsv
```

Prepare status:

```text
/data/bht/project01_baker_serhyd_routeB_20260701/manifests/baker_serhyd_routeB_prepare_status.json
```

Source repository commit:

```text
2442801d2c4a37994197faa8f20b459381e5dbd6
```

Staged Baker files include:

- `mu1.params`
- `1LNS_mu1.cst`
- `simple_theozyme.pdb`
- `bn1.params`
- `theozyme.cst`
- `theozyme.pdb`
- `super_af2_bu2.pdb`

## GPU Status

At prepare time:

```text
gpu_util_percent = 100
gpu_memory_used_mib = 4617
gpu_memory_total_mib = 81920
```

The active GPU process was a `pxdesign` Python process. The Baker smoke was not launched at that moment to avoid competing with a fully utilized GPU.

## Next Action

Run the checked-in script when GPU utilization drops below the threshold:

```bash
bash project01/phase1_20260701/scripts/launch_baker_theozyme_smoke.sh /data/bht/project01_baker_serhyd_routeB_20260701
```

Default behavior: if GPU utilization is above `MAX_GPU_UTIL=40`, the script writes `BLOCKED_GPU_BUSY_OR_UNAVAILABLE` instead of launching.

This remains a compute-resource block, not a biological or design failure.

## Waiting Launcher

The checked-in waiting launcher can be used for long GPU waits:

```bash
nohup bash project01/phase1_20260701/scripts/wait_and_launch_baker_theozyme_smoke.sh \
  /data/bht/project01_baker_serhyd_routeB_20260701 \
  > /data/bht/project01_baker_serhyd_routeB_20260701/logs/baker_theozyme_smoke_monitor.nohup.log 2>&1 &
```

Defaults:

```text
MAX_GPU_UTIL=40
INTERVAL_SECONDS=300
MAX_WAIT_MINUTES=360
```

It writes only lightweight status/log files and invokes `launch_baker_theozyme_smoke.sh` once GPU utilization drops below the threshold.

Current monitor:

```text
remote_pid = 4100072
remote_status = /data/bht/project01_baker_serhyd_routeB_20260701/manifests/baker_theozyme_smoke_monitor_status.json
remote_log = /data/bht/project01_baker_serhyd_routeB_20260701/logs/baker_theozyme_smoke_monitor.log
status_at_2026-07-01T00:40:59+08:00 = WAITING_GPU_BUSY
gpu_util_percent = 99
gpu_memory_used_mib = 4619
```

This monitor does not count any candidate as evaluated until the launcher actually starts and RFAA writes output.

## Launcher Check

The launcher was pulled and syntax-checked on the GPU node from commit `087a779`.

Result:

```text
status = BLOCKED_GPU_BUSY_OR_UNAVAILABLE
gpu_util_percent = 100
gpu_memory_used_mib = 4617
max_gpu_util_percent = 40
```

Remote status file:

```text
/data/bht/project01_baker_serhyd_routeB_20260701/manifests/baker_theozyme_smoke_status.json
```
