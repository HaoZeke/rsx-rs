#!/usr/bin/env python3
"""Merge per-dataset literature benchmark CSV files."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def merge_csv_files(inputs: list[Path], output: Path) -> None:
    if not inputs:
        raise ValueError("at least one input CSV is required")
    output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames: list[str] | None = None
    rows: list[dict[str, str]] = []
    for path in sorted(inputs):
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames is None:
                continue
            if fieldnames is None:
                fieldnames = reader.fieldnames
            elif reader.fieldnames != fieldnames:
                raise ValueError(f"{path} has incompatible columns")
            rows.extend(reader)
    if fieldnames is None:
        raise ValueError("input CSV files did not contain a header")
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", default=Path("benchmarks/results/slurm"), type=Path)
    parser.add_argument("--pattern", default="literature_speed_comparison_*.csv")
    parser.add_argument("--output", default=Path("benchmarks/results/literature_speed_comparison.csv"), type=Path)
    args = parser.parse_args()

    inputs = sorted(args.input_dir.glob(args.pattern))
    merge_csv_files(inputs, args.output)
    print(f"Merged {len(inputs)} files into {args.output}")


if __name__ == "__main__":
    main()
