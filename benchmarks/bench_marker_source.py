#!/usr/bin/env python3
"""Benchmark the three MarkerTableSource code paths against each literature
dataset:

1. =MarkerTable.from_path= -- the legacy mmap'd TSV streamer used by the CLI.
2. =MarkerTable.from_dataframe= with =RSX_FORCE_SPILL=0= -- the in-memory
   Arrow path (RecordBatches stay in RAM, the Beissinger estimator agrees).
3. =MarkerTable.from_dataframe= with =RSX_FORCE_SPILL=1= -- the spilled
   Parquet path (the estimator forced to disk).

For each (dataset, source) combination we run a small fixed suite of commands
(=freq=, =depth=, =triage=) and record wall time + peak RSS in megabytes.
The output CSV is read by the paper docs and the figure (depth-stability
SVG counterpart) lands in =docs/figures/=.
"""

from __future__ import annotations

import argparse
import gc
import os
import resource
import sys
import time
from pathlib import Path

import pandas as pd

# Lazy import: only needed for the actual benchmark, not the CLI parse.
def _import_pyrsx():
    import pyrsx  # noqa: F401

    return pyrsx


def peak_rss_mb() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


def reset_peak_rss() -> float:
    # ru_maxrss is monotonic per-process; we just record a baseline and
    # report deltas. RSS does not strictly drop after Python frees objects,
    # so use this as a high-water-mark probe rather than a strict measure.
    return peak_rss_mb()


def time_call(fn) -> tuple[float, float, float]:
    gc.collect()
    rss_before = reset_peak_rss()
    t0 = time.perf_counter()
    fn()
    elapsed = time.perf_counter() - t0
    rss_after = peak_rss_mb()
    return elapsed, rss_before, rss_after


def load_popmap(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", header=None, names=["individual", "group"])


def bench_one(
    dataset: str,
    markers_path: Path,
    popmap_path: Path,
    source: str,
    spill: str,
) -> list[dict]:
    """Run the suite of commands and return a row per command."""
    pyrsx = _import_pyrsx()
    from pyrsx import MarkerTable, TriageParams

    rows: list[dict] = []
    # Source construction is also measured because it includes the
    # Arrow IPC serialisation for the DataFrame paths.
    os.environ["RSX_FORCE_SPILL"] = spill

    if source == "path":
        construct = lambda: MarkerTable.from_path(str(markers_path))
        popmap_arg = str(popmap_path)
    else:
        markers_df = pd.read_csv(markers_path, sep="\t", comment="#")
        popmap_df = load_popmap(popmap_path)
        construct = lambda: MarkerTable.from_dataframe(markers_df)
        popmap_arg = popmap_df

    # Measure construction.
    elapsed, rss_b, rss_a = time_call(construct)
    rows.append(_row(dataset, source, spill, "construct", elapsed, rss_b, rss_a))

    mt = construct()  # rebuild once for the actual command runs

    def _freq():
        return mt.freq(min_depth=10)

    def _depth():
        return mt.depth(popmap=popmap_arg, min_frequency=0.1)

    def _triage():
        params = TriageParams(
            min_depth=10,
            posterior_threshold=0.9,
            prior=0.01,
            linked_prob=0.9,
            group1="male",
            group2="female",
        )
        return mt.triage(popmap=popmap_arg, params=params)

    for label, fn in [("freq", _freq), ("depth", _depth), ("triage", _triage)]:
        try:
            elapsed, rss_b, rss_a = time_call(fn)
            rows.append(_row(dataset, source, spill, label, elapsed, rss_b, rss_a))
        except Exception as exc:  # noqa: BLE001
            rows.append(_row(dataset, source, spill, label, float("nan"), float("nan"), float("nan"), error=str(exc)))

    return rows


def _row(
    dataset: str,
    source: str,
    spill: str,
    command: str,
    elapsed: float,
    rss_before: float,
    rss_after: float,
    error: str = "",
) -> dict:
    return {
        "dataset": dataset,
        "source": source,
        "force_spill": spill,
        "command": command,
        "elapsed_seconds": elapsed,
        "rss_before_mb": rss_before,
        "rss_after_mb": rss_after,
        "rss_delta_mb": rss_after - rss_before if rss_before == rss_before else float("nan"),
        "error": error,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--workdir",
        type=Path,
        default=Path("benchmarks/literature-workdir"),
        help="literature-workdir root with per-dataset folders",
    )
    parser.add_argument(
        "--datasets",
        nargs="+",
        default=[
            "danio_albolineatus",
            "notothenia_rossii",
            "plecoglossus_altivelis",
            "tinca_tinca",
        ],
    )
    parser.add_argument(
        "--results",
        type=Path,
        default=Path("benchmarks/results/marker_source_paths.csv"),
    )
    args = parser.parse_args()

    all_rows: list[dict] = []
    for ds in args.datasets:
        markers = args.workdir / ds / "markers_table.tsv"
        popmap = args.workdir / ds / "popmap.tsv"
        if not (markers.exists() and popmap.exists()):
            print(f"skipping {ds}: missing markers or popmap", file=sys.stderr)
            continue
        for source, spill in [("path", "0"), ("dataframe", "0"), ("dataframe", "1")]:
            print(f"=== {ds} source={source} force_spill={spill} ===", flush=True)
            try:
                rows = bench_one(ds, markers, popmap, source, spill)
            except Exception as exc:  # noqa: BLE001
                print(f"    {ds}/{source}/{spill} failed: {exc}", file=sys.stderr)
                rows = [_row(ds, source, spill, "construct", float("nan"), float("nan"), float("nan"), error=str(exc))]
            all_rows.extend(rows)
            for row in rows:
                print(
                    f"    {row['command']:<10} {row['elapsed_seconds']:>8.3f}s  "
                    f"RSS {row['rss_before_mb']:>8.1f} -> {row['rss_after_mb']:>8.1f} MB  "
                    f"{row['error']}".rstrip()
                )

    args.results.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(all_rows).to_csv(args.results, index=False)
    print(f"Wrote {args.results} ({len(all_rows)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
