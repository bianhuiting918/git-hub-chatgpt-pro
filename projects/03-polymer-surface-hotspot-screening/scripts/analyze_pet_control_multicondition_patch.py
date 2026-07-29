#!/usr/bin/env python3
import argparse
import csv
import datetime
import hashlib
import json
import math
import pathlib
import random
import statistics
import tempfile
import os
from collections import defaultdict

ROOT = pathlib.Path("/work/home/acshdt1dks/polymer_surface_hotspot_screen_20260725")
PATCH_METRICS = ROOT / "results/experimental_control31_patch_relative_stickiness_20260729_v1/patch_relative_stickiness_metrics.tsv"
ACTIVITY = pathlib.Path("/work/home/acshdt1dks/petase_orbmol_lg1_lg4_layer8343_20260721/inputs/activity/activity_energy_long_authority.tsv")
OUT = ROOT / "results/experimental_control31_multicondition_patch_analysis_20260729_v1"
PRIMARY_METRICS = (
    "r14_relative_external_sticky_fraction",
    "r14_composite_catalytic_candidate_percentile",
    "r14_mean_catalytic_candidate_percentile_ACOA",
)
ALL_AGGREGATE_METRICS = tuple(
    f"r{radius}_{suffix}"
    for radius in (6, 10, 14)
    for suffix in (
        "relative_external_sticky_fraction",
        "composite_catalytic_candidate_percentile",
        "mean_catalytic_candidate_percentile_ACOA",
    )
)
MIN_CONDITION_N = 10
N_PERM = 20000
N_BOOT = 5000
SEED = 20260729


def read_tsv(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def write_tsv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", newline="", dir=path.parent, delete=False
    ) as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
        tmp = pathlib.Path(handle.name)
    os.replace(tmp, path)


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", dir=path.parent, delete=False
    ) as handle:
        json.dump(value, handle, indent=2, sort_keys=True, ensure_ascii=False, allow_nan=False)
        handle.write("\n")
        tmp = pathlib.Path(handle.name)
    os.replace(tmp, path)


def rankdata(values):
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(values):
        j = i + 1
        while j < len(values) and values[order[j]] == values[order[i]]:
            j += 1
        average_rank = (i + j - 1) / 2 + 1
        for k in range(i, j):
            ranks[order[k]] = average_rank
        i = j
    return ranks


def pearson(x, y):
    if len(x) != len(y) or len(x) < 2:
        return float("nan")
    mean_x = statistics.fmean(x)
    mean_y = statistics.fmean(y)
    dx = [value - mean_x for value in x]
    dy = [value - mean_y for value in y]
    denominator = math.sqrt(sum(v * v for v in dx) * sum(v * v for v in dy))
    if denominator == 0:
        return float("nan")
    return sum(a * b for a, b in zip(dx, dy)) / denominator


def spearman(x, y):
    return pearson(rankdata(x), rankdata(y))


def spearman_stats(x, y, seed):
    rho = spearman(x, y)
    if not math.isfinite(rho):
        return rho, float("nan"), float("nan"), float("nan"), 0
    rng = random.Random(seed)
    rank_x = rankdata(x)
    rank_y = rankdata(y)
    hits = 0
    for _ in range(N_PERM):
        permuted = rank_y[:]
        rng.shuffle(permuted)
        if abs(pearson(rank_x, permuted)) >= abs(rho) - 1e-15:
            hits += 1
    boots = []
    n = len(x)
    for _ in range(N_BOOT):
        indices = [rng.randrange(n) for _ in range(n)]
        value = spearman([x[i] for i in indices], [y[i] for i in indices])
        if math.isfinite(value):
            boots.append(value)
    boots.sort()
    if not boots:
        return rho, (hits + 1) / (N_PERM + 1), float("nan"), float("nan"), 0

    def quantile(fraction):
        position = fraction * (len(boots) - 1)
        lower = int(position)
        upper = min(lower + 1, len(boots) - 1)
        weight = position - lower
        return boots[lower] * (1 - weight) + boots[upper] * weight

    return (
        rho,
        (hits + 1) / (N_PERM + 1),
        quantile(0.025),
        quantile(0.975),
        len(boots),
    )


def bh(rows, pkey, qkey):
    finite = sorted(
        [(index, float(row[pkey])) for index, row in enumerate(rows)
         if row.get(pkey) not in ("", None) and math.isfinite(float(row[pkey]))],
        key=lambda item: item[1],
    )
    last = 1.0
    total = len(finite)
    for rank in range(total, 0, -1):
        index, pvalue = finite[rank - 1]
        last = min(last, pvalue * total / rank)
        rows[index][qkey] = last
    for row in rows:
        row.setdefault(qkey, "")


def parse_condition(value):
    parts = value.split("|")
    if len(parts) != 4:
        raise ValueError(f"unexpected condition: {value}")
    return tuple(parts)


def condition_percentiles(rows):
    ordered = sorted(rows, key=lambda row: row["sequence_md5"])
    values = [float(row["activity_value"]) for row in ordered]
    ranks = rankdata(values)
    denominator = len(values) - 1
    if denominator == 0:
        scaled = [0.5]
    else:
        scaled = [(rank - 1) / denominator for rank in ranks]
    return {row["sequence_md5"]: value for row, value in zip(ordered, scaled)}


def safe_number(value):
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {key: safe_number(item) for key, item in value.items()}
    if isinstance(value, list):
        return [safe_number(item) for item in value]
    return value


def self_test():
    assert rankdata([3, 1, 1, 2]) == [4.0, 1.5, 1.5, 3.0]
    assert abs(spearman([1, 2, 3], [3, 2, 1]) + 1.0) < 1e-12
    assert parse_condition("Table D3|PET_unspecified|50.0C|H7.5") == (
        "Table D3", "PET_unspecified", "50.0C", "H7.5"
    )
    percentiles = condition_percentiles([
        {"sequence_md5": "a", "activity_value": 0.0},
        {"sequence_md5": "b", "activity_value": 0.0},
        {"sequence_md5": "c", "activity_value": 2.0},
    ])
    assert percentiles == {"a": 0.25, "b": 0.25, "c": 1.0}
    rows = [{"p": 0.01}, {"p": 0.04}, {"p": 0.03}]
    bh(rows, "p", "q")
    assert [round(row["q"], 12) for row in rows] == [0.03, 0.04, 0.04]
    print("SELF_TEST_PASS")


def load_features():
    rows = read_tsv(PATCH_METRICS)
    if len(rows) != 30 or len({row["sequence_md5"] for row in rows}) != 30:
        raise SystemExit("patch feature denominator is not 30 unique MD5")
    required = {"sequence_md5", "protein_names", *ALL_AGGREGATE_METRICS}
    if not required.issubset(rows[0]):
        raise SystemExit("patch feature columns missing")
    output = {}
    for row in rows:
        record = {
            "sequence_md5": row["sequence_md5"],
            "protein_names": row["protein_names"],
        }
        for metric in ALL_AGGREGATE_METRICS:
            record[metric] = float(row[metric])
        output[row["sequence_md5"]] = record
    return output


def load_activity(feature_md5s):
    source = read_tsv(ACTIVITY)
    filtered = []
    seen = set()
    for row in source:
        if (
            row["dataset"] != "NATURAL2022"
            or row["activity_metric"] != "sum aromatic products"
            or row["activity_unit"] != "mg/L"
            or row["sequence_md5"] not in feature_md5s
        ):
            continue
        key = (row["sequence_md5"], row["condition"])
        if key in seen:
            raise SystemExit(f"duplicate activity key: {key}")
        seen.add(key)
        table, substrate, temperature, buffer_name = parse_condition(row["condition"])
        filtered.append({
            "sequence_md5": row["sequence_md5"],
            "protein_names": row["protein_names"],
            "condition": row["condition"],
            "table": table,
            "substrate": substrate,
            "temperature": temperature,
            "buffer": buffer_name,
            "activity_value": float(row["activity_value"]),
            "activity_unit": row["activity_unit"],
        })
    if not filtered:
        raise SystemExit("no joined NATURAL2022 activity rows")
    return filtered


def joined_rows(activity_rows, features, table):
    output = []
    for row in activity_rows:
        if row["table"] != table:
            continue
        merged = dict(row)
        merged.update({
            metric: features[row["sequence_md5"]][metric]
            for metric in ALL_AGGREGATE_METRICS
        })
        output.append(merged)
    return sorted(output, key=lambda row: (row["condition"], row["sequence_md5"]))


def condition_associations(rows):
    grouped = defaultdict(list)
    for row in rows:
        grouped[row["condition"]].append(row)
    output = []
    for condition_index, condition in enumerate(sorted(grouped)):
        records = grouped[condition]
        if len(records) < MIN_CONDITION_N:
            continue
        table, substrate, temperature, buffer_name = parse_condition(condition)
        y = [row["activity_value"] for row in records]
        for metric_index, metric in enumerate(PRIMARY_METRICS):
            x = [row[metric] for row in records]
            rho, pvalue, low, high, boot_n = spearman_stats(
                x, y, SEED + 1009 * condition_index + 97 * metric_index
            )
            output.append({
                "condition": condition,
                "table": table,
                "substrate": substrate,
                "temperature": temperature,
                "buffer": buffer_name,
                "metric": metric,
                "n": len(records),
                "missing_from_patch_panel_n": 30 - len(records),
                "zero_activity_n": sum(row["activity_value"] == 0 for row in records),
                "spearman_rho": rho,
                "p_permutation": pvalue,
                "ci95_boot_low": low,
                "ci95_boot_high": high,
                "bootstrap_valid": boot_n,
                "analysis_role": "D3_CONDITION_REPLICATION",
            })
    bh(output, "p_permutation", "q_bh_d3_condition_metric_family")
    return sorted(output, key=lambda row: (
        row["temperature"], row["buffer"], row["metric"]
    ))


def aggregate_condition_ranks(rows, features, group_label, group_value):
    grouped = defaultdict(list)
    for row in rows:
        grouped[row["condition"]].append(row)
    values = defaultdict(list)
    conditions_used = []
    for condition in sorted(grouped):
        records = grouped[condition]
        if len(records) < MIN_CONDITION_N:
            continue
        conditions_used.append(condition)
        percentiles = condition_percentiles(records)
        for md5, percentile in percentiles.items():
            values[md5].append(percentile)
    output = []
    for md5 in sorted(values):
        record = {
            "group_type": group_label,
            "group_value": group_value,
            "sequence_md5": md5,
            "protein_names": features[md5]["protein_names"],
            "mean_condition_activity_percentile": statistics.fmean(values[md5]),
            "median_condition_activity_percentile": statistics.median(values[md5]),
            "n_conditions": len(values[md5]),
            "eligible_condition_count": len(conditions_used),
        }
        record.update({metric: features[md5][metric] for metric in ALL_AGGREGATE_METRICS})
        output.append(record)
    return output, conditions_used


def aggregate_associations(records, metrics, role, seed_offset):
    output = []
    y = [row["mean_condition_activity_percentile"] for row in records]
    for index, metric in enumerate(metrics):
        x = [row[metric] for row in records]
        rho, pvalue, low, high, boot_n = spearman_stats(
            x, y, SEED + seed_offset + 131 * index
        )
        output.append({
            "group_type": records[0]["group_type"],
            "group_value": records[0]["group_value"],
            "metric": metric,
            "n": len(records),
            "eligible_condition_count": records[0]["eligible_condition_count"],
            "min_conditions_per_protein": min(row["n_conditions"] for row in records),
            "max_conditions_per_protein": max(row["n_conditions"] for row in records),
            "spearman_rho": rho,
            "p_permutation": pvalue,
            "ci95_boot_low": low,
            "ci95_boot_high": high,
            "bootstrap_valid": boot_n,
            "analysis_role": role if metric in PRIMARY_METRICS else "RADIUS_SENSITIVITY",
        })
    return output


def direction_consistency(rows):
    output = []
    for metric in PRIMARY_METRICS:
        values = [
            float(row["spearman_rho"]) for row in rows
            if row["metric"] == metric and math.isfinite(float(row["spearman_rho"]))
        ]
        expected = "POSITIVE" if "external" in metric else "NEGATIVE"
        count = sum(value > 0 for value in values) if expected == "POSITIVE" else sum(
            value < 0 for value in values
        )
        output.append({
            "metric": metric,
            "expected_direction": expected,
            "finite_condition_n": len(values),
            "expected_direction_n": count,
            "expected_direction_fraction": count / len(values) if values else "",
            "median_condition_rho": statistics.median(values) if values else "",
            "min_condition_rho": min(values) if values else "",
            "max_condition_rho": max(values) if values else "",
            "interpretation": "DESCRIPTIVE_CORRELATED_ENDPOINTS",
        })
    return output


def main():
    features = load_features()
    activity_rows = load_activity(set(features))
    d3_rows = joined_rows(activity_rows, features, "Table D3")
    d6_rows = joined_rows(activity_rows, features, "Table D6")
    if not d3_rows:
        raise SystemExit("no Table D3 joined rows")
    OUT.mkdir(parents=True, exist_ok=True)

    d3_join_path = OUT / "d3_joined_records.tsv"
    joined_fields = [
        "sequence_md5", "protein_names", "condition", "table", "substrate",
        "temperature", "buffer", "activity_value", "activity_unit",
        *PRIMARY_METRICS,
    ]
    write_tsv(
        d3_join_path,
        [{field: row[field] for field in joined_fields} for row in d3_rows],
        joined_fields,
    )

    d3_condition = condition_associations(d3_rows)
    d3_condition_path = OUT / "d3_condition_associations.tsv"
    write_tsv(d3_condition_path, d3_condition, list(d3_condition[0]))

    d3_aggregate, d3_conditions_used = aggregate_condition_ranks(
        d3_rows, features, "TABLE", "Table D3"
    )
    d3_aggregate_path = OUT / "d3_crosscondition_aggregate.tsv"
    write_tsv(d3_aggregate_path, d3_aggregate, list(d3_aggregate[0]))

    d3_aggregate_assoc = aggregate_associations(
        d3_aggregate, ALL_AGGREGATE_METRICS, "D3_CROSSCONDITION_PRIMARY", 300000
    )
    bh(d3_aggregate_assoc, "p_permutation", "q_bh_d3_aggregate_metric_family")
    d3_aggregate_assoc = sorted(
        d3_aggregate_assoc,
        key=lambda row: (0 if row["metric"] in PRIMARY_METRICS else 1,
                         row["p_permutation"], row["metric"]),
    )
    d3_aggregate_assoc_path = OUT / "d3_crosscondition_associations.tsv"
    write_tsv(d3_aggregate_assoc_path, d3_aggregate_assoc, list(d3_aggregate_assoc[0]))

    directions = direction_consistency(d3_condition)
    direction_path = OUT / "d3_direction_consistency.tsv"
    write_tsv(direction_path, directions, list(directions[0]))

    d6_by_form = defaultdict(list)
    for row in d6_rows:
        d6_by_form[row["substrate"]].append(row)
    d6_aggregate_all = []
    d6_associations = []
    d6_condition_counts = {}
    for form_index, form in enumerate(sorted(d6_by_form)):
        aggregate, conditions_used = aggregate_condition_ranks(
            d6_by_form[form], features, "TABLE_D6_SUBSTRATE_FORM", form
        )
        if len(aggregate) < MIN_CONDITION_N:
            continue
        d6_condition_counts[form] = len(conditions_used)
        d6_aggregate_all.extend(aggregate)
        d6_associations.extend(
            aggregate_associations(
                aggregate, PRIMARY_METRICS, "D6_LOW_POWER_EXPLORATORY",
                500000 + 10000 * form_index,
            )
        )
    if not d6_aggregate_all or not d6_associations:
        raise SystemExit("no eligible Table D6 substrate-form aggregate")
    bh(d6_associations, "p_permutation", "q_bh_d6_form_metric_family")
    d6_associations = sorted(
        d6_associations, key=lambda row: (row["group_value"], row["p_permutation"], row["metric"])
    )
    d6_aggregate_path = OUT / "d6_crosscondition_aggregate.tsv"
    write_tsv(d6_aggregate_path, d6_aggregate_all, list(d6_aggregate_all[0]))
    d6_association_path = OUT / "d6_form_associations.tsv"
    write_tsv(d6_association_path, d6_associations, list(d6_associations[0]))

    output_paths = {
        "d3_joined_records": d3_join_path,
        "d3_condition_associations": d3_condition_path,
        "d3_crosscondition_aggregate": d3_aggregate_path,
        "d3_crosscondition_associations": d3_aggregate_assoc_path,
        "d3_direction_consistency": direction_path,
        "d6_crosscondition_aggregate": d6_aggregate_path,
        "d6_form_associations": d6_association_path,
    }
    primary_d3 = [
        row for row in d3_aggregate_assoc if row["metric"] in PRIMARY_METRICS
    ]
    audit = {
        "status": "MULTICONDITION_PATCH_ANALYSIS_PASS",
        "created_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "evaluated_universe": "Nature 2022 PET controls with exact-canonical valid PATCH_PASS",
        "patch_exact_canonical_n": len(features),
        "protein_202_excluded": True,
        "protein_202_status": "NONCANONICAL_MAPPING_REQUIRED",
        "material_family": "PET",
        "pet_nylon_denominators_separate": True,
        "activity_dataset": "NATURAL2022",
        "activity_metric": "sum aromatic products",
        "activity_unit": "mg/L",
        "activity_endpoint_caveat": "Endpoint sum of aromatic products; not a single-product readout, kinetics, or barrier.",
        "patch_recomputation": False,
        "table_d3_joined_rows": len(d3_rows),
        "table_d3_eligible_conditions": len(d3_conditions_used),
        "table_d3_crosscondition_protein_n": len(d3_aggregate),
        "table_d6_joined_rows": len(d6_rows),
        "table_d6_eligible_conditions_by_form": d6_condition_counts,
        "table_d6_low_power_exploratory": True,
        "table_d4_excluded": True,
        "table_d4_exclusion_reason": "Only 3-4 records per condition.",
        "minimum_condition_n": MIN_CONDITION_N,
        "permutation_replicates": N_PERM,
        "bootstrap_replicates": N_BOOT,
        "random_seed": SEED,
        "primary_metrics": list(PRIMARY_METRICS),
        "d3_crosscondition_primary_associations": primary_d3,
        "d3_direction_consistency": directions,
        "input_patch_metrics": {"path": str(PATCH_METRICS), "sha256": sha256(PATCH_METRICS)},
        "input_activity_authority": {"path": str(ACTIVITY), "sha256": sha256(ACTIVITY)},
        "outputs": {
            name: {"path": str(path), "sha256": sha256(path)}
            for name, path in output_paths.items()
        },
    }
    audit_path = OUT / "ANALYSIS_PASS.json"
    write_json(audit_path, safe_number(audit))
    print(json.dumps(safe_number({
        "status": audit["status"],
        "audit_path": str(audit_path),
        "audit_sha256": sha256(audit_path),
        "d3_primary": primary_d3,
        "d3_direction_consistency": directions,
        "d6_associations": d6_associations,
    }), indent=2, ensure_ascii=False, allow_nan=False))
    print("MULTICONDITION_ANALYSIS_PASS")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
    else:
        main()
