#!/usr/bin/env python3
"""Render immutable NylC M1 replica MDPs with a frozen seed allow-list."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


APPROVED_SEEDS = {"26711", "26723", "26737"}
KEY_RE = re.compile(r"^\s*([^;=]+?)\s*=\s*([^;]*?)(\s*;.*)?$")


def parse_mdp(text: str) -> dict[str, str]:
    values = {}
    for line_number, raw in enumerate(text.splitlines(), 1):
        body = raw.split(";", 1)[0].strip()
        if not body or "=" not in body:
            continue
        key, value = body.split("=", 1)
        normalized = key.strip().lower()
        if normalized in values:
            raise ValueError(f"duplicate MDP key {normalized} at line {line_number}")
        values[normalized] = value.strip()
    return values


def render_mdp(template: str, replacements: dict[str, str]) -> str:
    values = parse_mdp(template)
    normalized_replacements = {
        str(key).strip().lower(): str(value).strip()
        for key, value in replacements.items()
    }
    missing = sorted(set(normalized_replacements) - set(values))
    if missing:
        raise ValueError(f"MDP key is not present in template: {', '.join(missing)}")
    if (
        "gen-seed" in normalized_replacements
        and normalized_replacements["gen-seed"] not in APPROVED_SEEDS
    ):
        raise ValueError(
            f"gen-seed must be one of the approved values: {sorted(APPROVED_SEEDS)}"
        )
    if "define" in normalized_replacements and "define" not in values:
        raise ValueError("cannot inject define into a fully unrestrained MDP")

    output = []
    replaced = set()
    for raw in template.splitlines():
        match = KEY_RE.match(raw)
        if match is None:
            output.append(raw)
            continue
        key = match.group(1).strip().lower()
        if key not in normalized_replacements:
            output.append(raw)
            continue
        comment = match.group(3) or ""
        output.append(
            f"{match.group(1).rstrip()} = {normalized_replacements[key]}{comment}"
        )
        replaced.add(key)
    if replaced != set(normalized_replacements):
        raise AssertionError("not every requested MDP key was replaced exactly once")
    return "\n".join(output) + "\n"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--template", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--set", action="append", default=[])
    args = parser.parse_args()

    if args.output.exists():
        print(f"refusing to overwrite existing MDP: {args.output}", file=sys.stderr)
        return 2
    replacements = {}
    for assignment in args.set:
        if "=" not in assignment:
            print(f"invalid --set assignment: {assignment}", file=sys.stderr)
            return 2
        key, value = assignment.split("=", 1)
        normalized = key.strip().lower()
        if normalized in replacements:
            print(f"duplicate --set key: {normalized}", file=sys.stderr)
            return 2
        replacements[normalized] = value.strip()
    try:
        rendered = render_mdp(args.template.read_text(), replacements)
    except ValueError as error:
        print(str(error), file=sys.stderr)
        return 2
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
