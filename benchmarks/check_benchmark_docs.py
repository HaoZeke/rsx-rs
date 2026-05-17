#!/usr/bin/env python3
"""Validate benchmark documentation against the tracked benchmark CSV."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


REQUIRED_TASKS = (
    "pixi run -e benchmark generate-data",
    "pixi run -e benchmark run-benchmarks",
    "pixi run -e benchmark plot-benchmarks",
    "pixi run -e benchmark test-package",
    "pixi run -e benchmark package-figshare",
)

STALE_COMMANDS = (
    "python3 benchmarks/generate_data.py",
    "bash benchmarks/run_benchmarks.sh",
    "python3 benchmarks/plot_benchmarks.py",
)


def load_csv(path: Path) -> dict[tuple[str, str], dict[str, float]]:
    rows: dict[tuple[str, str], dict[str, float]] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"scale", "command", "impl", "time_secs"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} is missing required columns: {', '.join(sorted(missing))}")
        for row in reader:
            key = (row["scale"], row["command"])
            rows.setdefault(key, {})[row["impl"]] = float(row["time_secs"])
    return rows


def extract_raw_table(path: Path) -> dict[tuple[str, str], tuple[str, str, str]]:
    in_table = False
    rows: dict[tuple[str, str], tuple[str, str, str]] = {}

    for line in path.read_text().splitlines():
        if line.strip() == "** Raw timing data (seconds, median of 3 runs)":
            in_table = True
            continue
        if in_table and line.startswith("** "):
            break
        if not in_table or not line.startswith("|"):
            continue
        if line.startswith("|-") or "Scale" in line and "Command" in line:
            continue
        cells = [cell.strip() for cell in line.strip().strip("|").split("|")]
        if len(cells) != 5:
            raise ValueError(f"unexpected raw timing row in {path}: {line}")
        scale, command, cpp, rust, speedup = cells
        rows[(scale, command)] = (cpp, rust, speedup)

    if not rows:
        raise ValueError(f"could not find raw timing table in {path}")
    return rows


def expected_doc_rows(csv_rows: dict[tuple[str, str], dict[str, float]]) -> dict[tuple[str, str], tuple[str, str, str]]:
    expected: dict[tuple[str, str], tuple[str, str, str]] = {}
    for key, impls in csv_rows.items():
        if "cpp" not in impls or "rust" not in impls:
            continue
        cpp = impls["cpp"]
        rust = impls["rust"]
        expected[key] = (f"{cpp:.3f}", f"{rust:.3f}", f"{cpp / rust:.1f}x")
    return expected


def compare_rows(
    expected: dict[tuple[str, str], tuple[str, str, str]],
    actual: dict[tuple[str, str], tuple[str, str, str]],
) -> list[str]:
    errors: list[str] = []
    for key in sorted(set(expected).difference(actual)):
        errors.append(f"missing documentation row for {key[0]} {key[1]}")
    for key in sorted(set(actual).difference(expected)):
        errors.append(f"unexpected documentation row for {key[0]} {key[1]}")
    for key in sorted(set(expected).intersection(actual)):
        if expected[key] != actual[key]:
            exp = " | ".join(expected[key])
            got = " | ".join(actual[key])
            errors.append(f"{key[0]} {key[1]}: expected {exp}, found {got}")
    return errors


def check_reproduction_commands(path: Path) -> list[str]:
    text = path.read_text()
    errors: list[str] = []
    for task in REQUIRED_TASKS:
        if task not in text:
            errors.append(f"missing reproduction task: {task}")
    for command in STALE_COMMANDS:
        if command in text:
            errors.append(f"stale direct reproduction command remains: {command}")
    return errors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", default="benchmarks/results/benchmark_results.csv", type=Path)
    parser.add_argument("--docs", default="docs/orgmode/reference/benchmarks.org", type=Path)
    args = parser.parse_args()

    expected = expected_doc_rows(load_csv(args.csv))
    actual = extract_raw_table(args.docs)
    errors = compare_rows(expected, actual)
    errors.extend(check_reproduction_commands(args.docs))

    if errors:
        for error in errors:
            print(error, file=sys.stderr)
        return 1

    print(f"{args.docs} matches {args.csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
