#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import time
import traceback
from pathlib import Path


ROOT = Path("/data/bht/project01_baker_serhyd_routeB_20260701")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--selected-tsv", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--status-tsv", required=True)
    parser.add_argument("--summary-json", required=True)
    parser.add_argument("--memory-fraction", type=float, default=0.20)
    parser.add_argument("--chunk-size", type=int, default=128)
    parser.add_argument("--max-total", type=int, default=0, help="0 means no limit")
    parser.add_argument("--limit-per-bin", type=int, default=0, help="0 means no per-bin limit")
    parser.add_argument("--bins", nargs="*", default=["90", "80", "70", "60", "50"])
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def select_rows(rows: list[dict[str, str]], bins: list[str], limit_per_bin: int, max_total: int) -> list[dict[str, str]]:
    selected: list[dict[str, str]] = []
    by_bin: dict[str, int] = {bin_id: 0 for bin_id in bins}
    for row in rows:
        bin_id = row.get("bin", "")
        if bin_id not in by_bin:
            continue
        if limit_per_bin and by_bin[bin_id] >= limit_per_bin:
            continue
        selected.append(row)
        by_bin[bin_id] += 1
        if max_total and len(selected) >= max_total:
            break
    return selected


def main() -> int:
    args = parse_args()
    selected_tsv = Path(args.selected_tsv)
    out_dir = Path(args.out_dir)
    status_tsv = Path(args.status_tsv)
    summary_json = Path(args.summary_json)
    out_dir.mkdir(parents=True, exist_ok=True)
    status_tsv.parent.mkdir(parents=True, exist_ok=True)
    summary_json.parent.mkdir(parents=True, exist_ok=True)

    all_rows = read_tsv(selected_tsv)
    rows = select_rows(all_rows, args.bins, args.limit_per_bin, args.max_total)

    import torch
    import esm
    import pytorch_lightning.utilities.seed as pl_seed

    if not hasattr(pl_seed, "seed_everything"):
        from lightning_fabric.utilities.seed import seed_everything

        pl_seed.seed_everything = seed_everything

    torch.cuda.set_per_process_memory_fraction(args.memory_fraction, 0)
    model = esm.pretrained.esmfold_v1().eval().cuda()
    model.set_chunk_size(args.chunk_size)

    fields = [
        "sample_id",
        "parent_sample_id",
        "bin",
        "identity",
        "mutation_count",
        "status",
        "elapsed_s",
        "bytes",
        "pdb",
        "error",
    ]
    done: list[dict[str, object]] = []
    started = time.strftime("%Y-%m-%dT%H:%M:%S")
    print(
        json.dumps(
            {
                "event": "start",
                "started": started,
                "selected_tsv": str(selected_tsv),
                "input_rows": len(all_rows),
                "run_rows": len(rows),
                "memory_fraction": args.memory_fraction,
                "chunk_size": args.chunk_size,
            },
            sort_keys=True,
        ),
        flush=True,
    )

    for index, row in enumerate(rows, 1):
        sample_id = row["sample_id"]
        out = out_dir / f"{sample_id}.esmfold.pdb"
        if out.exists() and out.stat().st_size > 1000:
            rec = {
                "sample_id": sample_id,
                "parent_sample_id": row.get("parent_sample_id", ""),
                "bin": row.get("bin", ""),
                "identity": row.get("identity", ""),
                "mutation_count": row.get("mutation_count", ""),
                "status": "SKIP_EXISTS",
                "elapsed_s": 0,
                "bytes": out.stat().st_size,
                "pdb": str(out),
                "error": "",
            }
            done.append(rec)
            print(f"[{index}/{len(rows)}] {sample_id} SKIP_EXISTS", flush=True)
            continue
        t0 = time.time()
        try:
            with torch.no_grad():
                pdb = model.infer_pdb(row["sequence"])
            out.write_text(pdb)
            rec = {
                "sample_id": sample_id,
                "parent_sample_id": row.get("parent_sample_id", ""),
                "bin": row.get("bin", ""),
                "identity": row.get("identity", ""),
                "mutation_count": row.get("mutation_count", ""),
                "status": "OK",
                "elapsed_s": round(time.time() - t0, 2),
                "bytes": out.stat().st_size,
                "pdb": str(out),
                "error": "",
            }
            print(f"[{index}/{len(rows)}] {sample_id} OK bytes={rec['bytes']} elapsed_s={rec['elapsed_s']}", flush=True)
        except Exception as exc:
            rec = {
                "sample_id": sample_id,
                "parent_sample_id": row.get("parent_sample_id", ""),
                "bin": row.get("bin", ""),
                "identity": row.get("identity", ""),
                "mutation_count": row.get("mutation_count", ""),
                "status": "FAIL",
                "elapsed_s": round(time.time() - t0, 2),
                "bytes": 0,
                "pdb": str(out),
                "error": repr(exc),
            }
            print(f"[{index}/{len(rows)}] {sample_id} FAIL {exc!r}", flush=True)
            traceback.print_exc()
        done.append(rec)
        with status_tsv.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(done)

    counts: dict[str, int] = {}
    counts_by_bin: dict[str, dict[str, int]] = {}
    for row in done:
        status = str(row["status"])
        bin_id = str(row["bin"])
        counts[status] = counts.get(status, 0) + 1
        counts_by_bin.setdefault(bin_id, {})
        counts_by_bin[bin_id][status] = counts_by_bin[bin_id].get(status, 0) + 1
    summary = {
        "status": "DONE",
        "started": started,
        "finished": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "selected_tsv": str(selected_tsv),
        "out_dir": str(out_dir),
        "status_tsv": str(status_tsv),
        "input_rows": len(all_rows),
        "run_rows": len(rows),
        "counts": counts,
        "counts_by_bin": counts_by_bin,
        "memory_fraction": args.memory_fraction,
        "chunk_size": args.chunk_size,
        "evaluated_universe": "Rows selected for this ESMFold run; unselected manifest rows are NOT_EVALUATED.",
    }
    summary_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
