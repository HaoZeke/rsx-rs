#!/usr/bin/env python3
"""Combine per-dataset low-depth SLURM shards into the paper sweep CSV."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd


SHARD_RE = re.compile(r"^lowdepth_(?P<dataset>.+)_d(?P<depth>[0-9]+)\.csv$")


def iter_shards(shard_dir: Path) -> list[Path]:
    def key(path: Path) -> tuple[str, int]:
        match = SHARD_RE.match(path.name)
        if not match:
            return (path.name, -1)
        meta = match.groupdict()
        return (meta["dataset"], int(meta["depth"]))

    return sorted(
        (p for p in shard_dir.glob("lowdepth_*_d*.csv") if SHARD_RE.match(p.name)),
        key=key,
    )


def collect_shards(shard_dir: Path) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for path in iter_shards(shard_dir):
        frames.append(pd.read_csv(path, dtype=str))
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


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

    shards = iter_shards(args.shard_dir)
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
