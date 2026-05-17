#!/usr/bin/env python3
"""Compare rsx and C++ RADSex on downloaded literature FASTQ datasets."""

from __future__ import annotations

import argparse
import csv
import shutil
import sys
from pathlib import Path

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from benchmarks.run_literature_benchmarks import (
    command_output,
    relpath,
    resolve_binary,
    run_command,
    summarize_output,
)


RESULT_COLUMNS = [
    "dataset",
    "impl",
    "command",
    "min_depth",
    "elapsed_seconds",
    "markers",
    "rows",
    "significant_markers",
    "output_bytes",
    "output_path",
]


def parse_min_depths(value: str) -> list[int]:
    depths = [int(item) for item in value.split(",") if item.strip()]
    if not depths:
        raise argparse.ArgumentTypeError("at least one minimum depth is required")
    return depths


def discover_datasets(workdir: Path, selected: list[str] | None) -> list[str]:
    if selected:
        return selected
    return sorted(path.name for path in workdir.iterdir() if (path / "samples").exists() and (path / "popmap.tsv").exists())


def prune_results(path: Path, dataset_names: set[str]) -> None:
    if not path.exists() or not dataset_names:
        return
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        rows = [row for row in reader if row.get("dataset") not in dataset_names]
        fieldnames = reader.fieldnames or RESULT_COLUMNS
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    tmp.replace(path)


class ResultWriter:
    def __init__(self, path: Path) -> None:
        self.path = path
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.handle = path.open("a", newline="")
        self.writer = csv.DictWriter(self.handle, fieldnames=RESULT_COLUMNS)
        if path.stat().st_size == 0:
            self.writer.writeheader()
            self.handle.flush()

    def close(self) -> None:
        self.handle.close()

    def write(self, row: dict[str, str]) -> None:
        self.writer.writerow({column: row.get(column, "") for column in RESULT_COLUMNS})
        self.handle.flush()


def command_path(root: Path, command: str, min_depth: int | None = None) -> Path:
    if command == "process":
        return root / "markers_table.tsv"
    if command == "depth":
        return root / "depth.tsv"
    suffix = ".fa" if command == "signif" else ".tsv"
    return root / f"{command}_{min_depth}{suffix}"


def run_impl(
    dataset: str,
    dataset_dir: Path,
    impl: str,
    binary: str,
    threads: int,
    groups: str,
    min_depths: list[int],
    force: bool,
    writer: ResultWriter,
) -> None:
    impl_dir = dataset_dir / "comparison" / impl
    logs_dir = impl_dir / "logs"
    impl_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)
    markers_table = command_path(impl_dir, "process")
    if force or not markers_table.exists():
        elapsed = run_command(
            [
                binary,
                "process",
                "-i",
                str(dataset_dir / "samples"),
                "-o",
                str(markers_table),
                "-T",
                str(threads),
                "-d",
                "1",
            ],
            logs_dir / "process.log",
        )
    else:
        elapsed = 0.0
    writer.write(
        {
            "dataset": dataset,
            "impl": impl,
            "command": "process",
            "min_depth": "1",
            "elapsed_seconds": f"{elapsed:.6f}",
        }
        | summarize_output("process", markers_table)
    )

    depth_path = command_path(impl_dir, "depth")
    if force or not depth_path.exists():
        elapsed = run_command(
            [
                binary,
                "depth",
                "-t",
                str(markers_table),
                "-p",
                str(dataset_dir / "popmap.tsv"),
                "-o",
                str(depth_path),
            ],
            logs_dir / "depth.log",
        )
    else:
        elapsed = 0.0
    writer.write(
        {
            "dataset": dataset,
            "impl": impl,
            "command": "depth",
            "elapsed_seconds": f"{elapsed:.6f}",
        }
        | summarize_output("depth", depth_path)
    )

    for min_depth in min_depths:
        for command in ("freq", "distrib", "signif"):
            output = command_path(impl_dir, command, min_depth)
            args = [
                binary,
                command,
                "-t",
                str(markers_table),
                "-o",
                str(output),
                "-d",
                str(min_depth),
            ]
            if command in {"distrib", "signif"}:
                args.extend(["-p", str(dataset_dir / "popmap.tsv"), "-G", groups])
            if command == "signif":
                args.append("-a")
            if force or not output.exists():
                elapsed = run_command(args, logs_dir / f"{command}_{min_depth}.log")
            else:
                elapsed = 0.0
            writer.write(
                {
                    "dataset": dataset,
                    "impl": impl,
                    "command": command,
                    "min_depth": str(min_depth),
                    "elapsed_seconds": f"{elapsed:.6f}",
                }
                | summarize_output(command, output)
            )


def require_binary(label: str, binary: str) -> str:
    resolved = resolve_binary(binary)
    if shutil.which(resolved) or Path(resolved).exists():
        return resolved
    raise FileNotFoundError(f"{label} binary not found: {binary}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workdir", default=Path("benchmarks/literature-workdir"), type=Path)
    parser.add_argument("--results", default=Path("benchmarks/results/literature_speed_comparison.csv"), type=Path)
    parser.add_argument("--dataset", action="append", help="Dataset name to compare; repeat for multiple datasets")
    parser.add_argument("--rsx", default=str(Path("target/release/rsx")))
    parser.add_argument("--radsex", default="radsex")
    parser.add_argument("--threads", default=4, type=int)
    parser.add_argument("--groups", default="male,female")
    parser.add_argument("--min-depths", default="1,2,5,10", type=parse_min_depths)
    parser.add_argument("--append", action="store_true")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    datasets = discover_datasets(args.workdir, args.dataset)
    if not args.append:
        prune_results(args.results, set(datasets))
    writer = ResultWriter(args.results)
    try:
        rsx = require_binary("rsx", args.rsx)
        radsex = require_binary("RADSex", args.radsex)
        for dataset in datasets:
            dataset_dir = args.workdir / dataset
            run_impl(dataset, dataset_dir, "rust", rsx, args.threads, args.groups, args.min_depths, args.force, writer)
            run_impl(dataset, dataset_dir, "cpp", radsex, args.threads, args.groups, args.min_depths, args.force, writer)
    finally:
        writer.close()
    print(f"Wrote {args.results}")


if __name__ == "__main__":
    main()
