#!/usr/bin/env python3
"""Merge per-dataset literature benchmark CSV files."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def merge_csv_files(inputs: list[Path], output: Path, replace_by: str = "dataset") -> None:
    if not inputs:
        raise ValueError("at least one input CSV is required")
    output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames: list[str] | None = None
    merged_rows: list[dict[str, str]] = []
    latest_groups: dict[str, tuple[tuple[int, str], list[dict[str, str]]]] = {}
    for path in sorted(inputs):
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames is None:
                continue
            if fieldnames is None:
                fieldnames = reader.fieldnames
            elif reader.fieldnames != fieldnames:
                raise ValueError(f"{path} has incompatible columns")
            rows = list(reader)
        if replace_by and fieldnames and replace_by in fieldnames:
            stamp = (path.stat().st_mtime_ns, str(path))
            grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
            for row in rows:
                grouped[row.get(replace_by, "")].append(row)
            for key, group_rows in grouped.items():
                if key not in latest_groups or stamp >= latest_groups[key][0]:
                    latest_groups[key] = (stamp, group_rows)
        else:
            merged_rows.extend(rows)
    if fieldnames is None:
        raise ValueError("input CSV files did not contain a header")
    if latest_groups:
        merged_rows = []
        for key in sorted(latest_groups):
            merged_rows.extend(latest_groups[key][1])
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(merged_rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", default=Path("benchmarks/results/slurm"), type=Path)
    parser.add_argument("--pattern", default="literature_speed_comparison_*.csv")
    parser.add_argument("--output", default=Path("benchmarks/results/literature_speed_comparison.csv"), type=Path)
    parser.add_argument("--replace-by", default="dataset", help="Keep only the newest shard for each value in this column")
    args = parser.parse_args()

    inputs = sorted(args.input_dir.glob(args.pattern))
    merge_csv_files(inputs, args.output, replace_by=args.replace_by)
    print(f"Merged {len(inputs)} files into {args.output}")


if __name__ == "__main__":
    main()
