#!/usr/bin/env python3
"""Combine per-dataset low-depth SLURM shards into the paper sweep CSV."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd


SHARD_RE = re.compile(r"^lowdepth_(?P<dataset>.+)_d(?P<depth>[0-9]+)\.csv$")


def collect_shards(shard_dir: Path) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for path in sorted(shard_dir.glob("lowdepth_*_d*.csv")):
        if not SHARD_RE.match(path.name):
            continue
        frames.append(pd.read_csv(path))
    if not frames:
        return pd.DataFrame()
    combined = pd.concat(frames, ignore_index=True)
    return combined.sort_values(["dataset", "min_depth", "mode"]).reset_index(drop=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--shard-dir",
        type=Path,
        default=Path("benchmarks/results"),
        help="Directory holding lowdepth_<dataset>_d<depth>.csv shards",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("benchmarks/results/literature_mode_effects_sweep.csv"),
        help="Combined sweep CSV path",
    )
    parser.add_argument(
        "--expected-shards",
        type=int,
        default=16,
        help="Expected shard count; set to 0 to disable the check",
    )
    args = parser.parse_args()

    shards = [p for p in args.shard_dir.glob("lowdepth_*_d*.csv") if SHARD_RE.match(p.name)]
    if args.expected_shards and len(shards) != args.expected_shards:
        print(
            f"expected {args.expected_shards} low-depth shards under {args.shard_dir}, found {len(shards)}",
            file=sys.stderr,
        )
        return 2

    combined = collect_shards(args.shard_dir)
    if combined.empty:
        print(f"no low-depth shards found under {args.shard_dir}", file=sys.stderr)
        return 2

    args.output.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(args.output, index=False)
    print(f"Wrote {args.output} ({len(combined)} rows from {len(shards)} shards)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
