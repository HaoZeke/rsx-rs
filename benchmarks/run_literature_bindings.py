#!/usr/bin/env python3
"""Run Python-binding feature checks on downloaded literature datasets."""

from __future__ import annotations

import argparse
import csv
import sys
import time
from pathlib import Path
from typing import Callable


RESULT_COLUMNS = [
    "dataset",
    "feature",
    "api_call",
    "min_depth",
    "elapsed_seconds",
    "input_path",
    "output_path",
    "output_rows",
    "output_bytes",
    "summary",
]


def resolve_markers_table(dataset_dir: Path) -> Path:
    candidates = [
        dataset_dir / "markers_table.tsv",
        dataset_dir / "comparison" / "rust" / "markers_table.tsv",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise FileNotFoundError(f"markers table not found under {dataset_dir}")


def count_data_rows(path: Path) -> int:
    if path.is_dir():
        return sum(count_data_rows(child) for child in path.iterdir() if child.is_file())
    rows = 0
    saw_header = False
    with path.open() as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            if not saw_header and ("\t" in stripped) and not stripped.split("\t", 1)[0].isdigit():
                saw_header = True
                continue
            saw_header = True
            rows += 1
    return rows


def output_bytes(path: Path) -> int:
    if path.is_dir():
        return sum(child.stat().st_size for child in path.iterdir() if child.is_file())
    return path.stat().st_size if path.exists() else 0


def timed_call(action: Callable[[], None]) -> float:
    start = time.perf_counter()
    action()
    return time.perf_counter() - start


def write_rows(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=RESULT_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)


def feature_row(
    dataset: str,
    feature: str,
    api_call: str,
    min_depth: int,
    elapsed: float,
    input_path: Path,
    output_path: Path,
    summary: str,
) -> dict[str, str]:
    return {
        "dataset": dataset,
        "feature": feature,
        "api_call": api_call,
        "min_depth": str(min_depth),
        "elapsed_seconds": f"{elapsed:.6f}",
        "input_path": str(input_path),
        "output_path": str(output_path),
        "output_rows": str(count_data_rows(output_path)),
        "output_bytes": str(output_bytes(output_path)),
        "summary": summary,
    }


def run_dataset(dataset_dir: Path, dataset: str, output_dir: Path, min_depth: int, n_components: int, bayes: bool) -> list[dict[str, str]]:
    import pyrsx

    markers_table = resolve_markers_table(dataset_dir)
    popmap = dataset_dir / "popmap.tsv"
    if not popmap.exists():
        raise FileNotFoundError(f"popmap not found: {popmap}")

    dataset_output = output_dir / dataset
    dataset_output.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, str]] = []

    freq_out = dataset_output / f"freq_d{min_depth}.tsv"
    elapsed = timed_call(lambda: pyrsx.freq(str(markers_table), str(freq_out), min_depth=min_depth))
    rows.append(
        feature_row(
            dataset,
            "Python bindings",
            "pyrsx.freq",
            min_depth,
            elapsed,
            markers_table,
            freq_out,
            "Depth-frequency distribution from a published RADSex workflow marker table.",
        )
    )

    depth_out = dataset_output / "depth.tsv"
    elapsed = timed_call(lambda: pyrsx.depth(str(markers_table), str(popmap), str(depth_out)))
    rows.append(
        feature_row(
            dataset,
            "Python bindings",
            "pyrsx.depth",
            min_depth,
            elapsed,
            markers_table,
            depth_out,
            "Per-sample depth QC from the same literature marker table and popmap.",
        )
    )

    distrib_out = dataset_output / f"distrib_fisher_fdr_d{min_depth}.tsv"
    elapsed = timed_call(
        lambda: pyrsx.distrib(
            str(markers_table),
            str(popmap),
            str(distrib_out),
            min_depth=min_depth,
            test="fisher",
            correction="fdr",
        )
    )
    rows.append(
        feature_row(
            dataset,
            "Statistical tests",
            "pyrsx.distrib(test='fisher', correction='fdr')",
            min_depth,
            elapsed,
            markers_table,
            distrib_out,
            "Fisher exact test and Benjamini-Hochberg correction exposed through Python.",
        )
    )

    pca_dir = dataset_output / f"pca_d{min_depth}"
    elapsed = timed_call(lambda: pyrsx.pca(str(markers_table), str(pca_dir), min_depth=min_depth, n_components=n_components))
    rows.append(
        feature_row(
            dataset,
            "Streaming PCA",
            "pyrsx.pca",
            min_depth,
            elapsed,
            markers_table,
            pca_dir,
            "Depth-matrix PCA using the O(n_individuals^2) Gram workflow.",
        )
    )

    if bayes:
        signif_out = dataset_output / f"signif_bayes_none_d{min_depth}.tsv"
        elapsed = timed_call(
            lambda: pyrsx.signif(
                str(markers_table),
                str(popmap),
                str(signif_out),
                min_depth=min_depth,
                correction="none",
                test="chisq",
                bayes=True,
            )
        )
        rows.append(
            feature_row(
                dataset,
                "Bayesian marker evidence",
                "pyrsx.signif(bayes=True)",
                min_depth,
                elapsed,
                markers_table,
                signif_out,
                "Bayes factor and posterior columns produced through Python bindings.",
            )
        )

    return rows


def discover_datasets(workdir: Path, selected: list[str] | None) -> list[str]:
    if selected:
        return selected
    return sorted(path.name for path in workdir.iterdir() if (path / "popmap.tsv").exists())


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workdir", default=Path("benchmarks/literature-workdir"), type=Path)
    parser.add_argument("--results", default=Path("benchmarks/results/literature_binding_results.csv"), type=Path)
    parser.add_argument("--output-dir", default=Path("benchmarks/literature-workdir/bindings"), type=Path)
    parser.add_argument("--dataset", action="append", help="Dataset name to analyze; repeat for multiple datasets")
    parser.add_argument("--min-depth", default=10, type=int)
    parser.add_argument("--n-components", default=4, type=int)
    parser.add_argument("--no-bayes", action="store_true")
    args = parser.parse_args()

    datasets = discover_datasets(args.workdir, args.dataset)
    rows: list[dict[str, str]] = []
    for dataset in datasets:
        rows.extend(
            run_dataset(
                args.workdir / dataset,
                dataset,
                args.output_dir,
                args.min_depth,
                args.n_components,
                bayes=not args.no_bayes,
            )
        )
    write_rows(args.results, rows)
    print(f"Wrote {args.results}")


if __name__ == "__main__":
    if __package__ in {None, ""}:
        sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    main()
