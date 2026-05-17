#!/usr/bin/env python3
"""Run rsx analysis modes on downloaded literature marker tables."""

from __future__ import annotations

import argparse
import csv
import math
import sys
import time
from collections import defaultdict
from pathlib import Path
from statistics import mean
from typing import Callable


RESULT_COLUMNS = [
    "dataset",
    "min_depth",
    "mode",
    "api_call",
    "elapsed_seconds",
    "output_path",
    "tested_markers",
    "output_rows",
    "significant_markers",
    "posterior_gt_0_5",
    "posterior_gt_0_9",
    "bayes_factor_gt_10",
    "pc1_variance_fraction",
    "pc2_variance_fraction",
    "sex_loading_delta_pc1",
    "sex_loading_delta_pc2",
    "singleton_fraction",
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


def timed_call(action: Callable[[], None]) -> float:
    start = time.perf_counter()
    action()
    return time.perf_counter() - start


def write_rows(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=RESULT_COLUMNS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def parse_metadata(line: str) -> dict[str, str]:
    if not line.startswith("#"):
        return {}
    metadata: dict[str, str] = {}
    for item in line[1:].strip().split(";"):
        if ":" not in item:
            continue
        key, value = item.split(":", 1)
        metadata[key.strip()] = value.strip()
    return metadata


def as_int(value: str | None) -> int:
    if value in {None, ""}:
        return 0
    return int(float(value))


def as_float(value: str | None) -> float:
    if value in {None, ""}:
        return 0.0
    return float(value)


def read_table(path: Path) -> tuple[dict[str, str], list[str], list[dict[str, str]]]:
    metadata: dict[str, str] = {}
    header: list[str] = []
    rows: list[dict[str, str]] = []
    with path.open(newline="") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                metadata |= parse_metadata(line)
                continue
            if not header:
                header = line.split("\t")
                continue
            values = line.split("\t")
            rows.append(dict(zip(header, values, strict=False)))
    return metadata, header, rows


def summarize_freq(path: Path) -> dict[str, str]:
    _, _, rows = read_table(path)
    retained_markers = sum(as_int(row.get("Count")) for row in rows)
    singleton_markers = sum(as_int(row.get("Count")) for row in rows if row.get("Frequency") == "1")
    singleton_fraction = singleton_markers / retained_markers if retained_markers else 0.0
    return {
        "tested_markers": str(retained_markers),
        "output_rows": str(len(rows)),
        "singleton_fraction": f"{singleton_fraction:.6g}",
        "summary": f"{retained_markers} retained marker-frequency observations; {singleton_fraction:.2%} are singletons.",
    }


def summarize_depth(path: Path) -> dict[str, str]:
    _, _, rows = read_table(path)
    reads = [as_int(row.get("Reads")) for row in rows]
    retained = [as_int(row.get("Retained")) for row in rows]
    group_counts: dict[str, int] = defaultdict(int)
    for row in rows:
        group_counts[row.get("Group", "")] += 1
    mean_reads = mean(reads) if reads else 0.0
    mean_retained = mean(retained) if retained else 0.0
    groups = ", ".join(f"{group}:{count}" for group, count in sorted(group_counts.items()) if group)
    return {
        "output_rows": str(len(rows)),
        "summary": f"{len(rows)} samples ({groups}); mean reads {mean_reads:.0f}; mean retained markers {mean_retained:.0f}.",
    }


def summarize_distrib(path: Path) -> dict[str, str]:
    metadata, _, rows = read_table(path)
    significant_markers = 0
    significant_cells = 0
    tested_markers = as_int(metadata.get("n_markers"))
    for row in rows:
        markers = as_int(row.get("Markers"))
        if row.get("Signif") == "True" and markers > 0:
            significant_cells += 1
            significant_markers += markers
    return {
        "tested_markers": str(tested_markers),
        "output_rows": str(len(rows)),
        "significant_markers": str(significant_markers),
        "summary": f"{significant_cells} significant distribution cells covering {significant_markers} markers.",
    }


def summarize_signif(path: Path) -> dict[str, str]:
    metadata, header, rows = read_table(path)
    has_bayes = "Posterior_SexLinked" in header and "Bayes_Factor" in header
    posterior_gt_0_5 = 0
    posterior_gt_0_9 = 0
    bayes_factor_gt_10 = 0
    if has_bayes:
        for row in rows:
            posterior = as_float(row.get("Posterior_SexLinked"))
            bayes_factor = as_float(row.get("Bayes_Factor"))
            posterior_gt_0_5 += int(posterior > 0.5)
            posterior_gt_0_9 += int(posterior > 0.9)
            bayes_factor_gt_10 += int(bayes_factor > 10.0)
    summary = {
        "tested_markers": str(as_int(metadata.get("n_markers"))),
        "output_rows": str(len(rows)),
        "significant_markers": str(len(rows)),
    }
    if has_bayes:
        summary |= {
            "posterior_gt_0_5": str(posterior_gt_0_5),
            "posterior_gt_0_9": str(posterior_gt_0_9),
            "bayes_factor_gt_10": str(bayes_factor_gt_10),
            "summary": (
                f"{len(rows)} marker rows; {posterior_gt_0_9} rows with posterior > 0.9; "
                f"{bayes_factor_gt_10} rows with Bayes factor > 10."
            ),
        }
    else:
        summary["summary"] = f"{len(rows)} marker rows."
    return summary


def load_popmap(path: Path) -> dict[str, str]:
    groups: dict[str, str] = {}
    with path.open() as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            sample, group = stripped.split("\t")[:2]
            groups[sample] = group
    return groups


def component_mean_delta(rows: list[dict[str, str]], groups: dict[str, str], component: str, group1: str, group2: str) -> float:
    values: dict[str, list[float]] = {group1: [], group2: []}
    for row in rows:
        group = groups.get(row["individual"])
        if group in values:
            values[group].append(as_float(row.get(component)))
    if not values[group1] or not values[group2]:
        return 0.0
    return abs(mean(values[group1]) - mean(values[group2]))


def summarize_pca(path: Path, popmap: Path, group1: str, group2: str) -> dict[str, str]:
    _, _, eigen_rows = read_table(path / "eigenvalues.tsv")
    _, _, loading_rows = read_table(path / "loadings.tsv")
    groups = load_popmap(popmap)
    pc1_variance = as_float(eigen_rows[0].get("variance_fraction")) if eigen_rows else 0.0
    pc2_variance = as_float(eigen_rows[1].get("variance_fraction")) if len(eigen_rows) > 1 else 0.0
    pc1_delta = component_mean_delta(loading_rows, groups, "PC1", group1, group2)
    pc2_delta = component_mean_delta(loading_rows, groups, "PC2", group1, group2)
    return {
        "output_rows": str(len(eigen_rows) + len(loading_rows)),
        "pc1_variance_fraction": f"{pc1_variance:.6g}",
        "pc2_variance_fraction": f"{pc2_variance:.6g}",
        "sex_loading_delta_pc1": f"{pc1_delta:.6g}",
        "sex_loading_delta_pc2": f"{pc2_delta:.6g}",
        "summary": (
            f"PC1 explains {pc1_variance:.2%} and PC2 explains {pc2_variance:.2%}; "
            f"|mean {group1}-{group2}| loading deltas are {pc1_delta:.4g} and {pc2_delta:.4g}."
        ),
    }


def base_row(
    dataset: str,
    min_depth: int,
    mode: str,
    api_call: str,
    elapsed: float,
    output_path: Path,
    summary: dict[str, str],
) -> dict[str, str]:
    row = {column: "" for column in RESULT_COLUMNS}
    row.update(
        {
            "dataset": dataset,
            "min_depth": str(min_depth),
            "mode": mode,
            "api_call": api_call,
            "elapsed_seconds": f"{elapsed:.6f}",
            "output_path": str(output_path),
        }
    )
    row.update(summary)
    return row


def run_dataset(
    dataset_dir: Path,
    dataset: str,
    output_dir: Path,
    min_depth: int,
    group1: str,
    group2: str,
    n_components: int,
) -> list[dict[str, str]]:
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
    rows.append(base_row(dataset, min_depth, "frequency_qc", "pyrsx.freq", elapsed, freq_out, summarize_freq(freq_out)))

    depth_out = dataset_output / "depth.tsv"
    elapsed = timed_call(lambda: pyrsx.depth(str(markers_table), str(popmap), str(depth_out)))
    rows.append(base_row(dataset, min_depth, "depth_qc", "pyrsx.depth", elapsed, depth_out, summarize_depth(depth_out)))

    distrib_modes = [
        ("distribution_chisq_bonferroni", "chisq", "bonferroni"),
        ("distribution_fisher_fdr", "fisher", "fdr"),
        ("distribution_gtest_fdr", "gtest", "fdr"),
    ]
    for mode, test, correction in distrib_modes:
        distrib_out = dataset_output / f"{mode}_d{min_depth}.tsv"
        elapsed = timed_call(
            lambda test=test, correction=correction, distrib_out=distrib_out: pyrsx.distrib(
                str(markers_table),
                str(popmap),
                str(distrib_out),
                min_depth=min_depth,
                group1=group1,
                group2=group2,
                test=test,
                correction=correction,
            )
        )
        rows.append(
            base_row(
                dataset,
                min_depth,
                mode,
                f"pyrsx.distrib(test='{test}', correction='{correction}')",
                elapsed,
                distrib_out,
                summarize_distrib(distrib_out),
            )
        )

    strict_out = dataset_output / f"signif_chisq_bonferroni_d{min_depth}.tsv"
    elapsed = timed_call(
        lambda: pyrsx.signif(
            str(markers_table),
            str(popmap),
            str(strict_out),
            min_depth=min_depth,
            group1=group1,
            group2=group2,
            correction="bonferroni",
            test="chisq",
        )
    )
    rows.append(
        base_row(
            dataset,
            min_depth,
            "marker_extraction_chisq_bonferroni",
            "pyrsx.signif(test='chisq', correction='bonferroni')",
            elapsed,
            strict_out,
            summarize_signif(strict_out),
        )
    )

    bayes_out = dataset_output / f"signif_bayes_chisq_none_d{min_depth}.tsv"
    elapsed = timed_call(
        lambda: pyrsx.signif(
            str(markers_table),
            str(popmap),
            str(bayes_out),
            min_depth=min_depth,
            group1=group1,
            group2=group2,
            correction="none",
            test="chisq",
            bayes=True,
        )
    )
    rows.append(
        base_row(
            dataset,
            min_depth,
            "bayesian_marker_table",
            "pyrsx.signif(test='chisq', correction='none', bayes=True)",
            elapsed,
            bayes_out,
            summarize_signif(bayes_out),
        )
    )

    pca_dir = dataset_output / f"pca_d{min_depth}"
    elapsed = timed_call(lambda: pyrsx.pca(str(markers_table), str(pca_dir), min_depth=min_depth, n_components=n_components))
    rows.append(base_row(dataset, min_depth, "streaming_pca", "pyrsx.pca", elapsed, pca_dir, summarize_pca(pca_dir, popmap, group1, group2)))

    return rows


def discover_datasets(workdir: Path, selected: list[str] | None) -> list[str]:
    if selected:
        return selected
    datasets: list[str] = []
    for path in workdir.iterdir():
        if not path.is_dir() or not (path / "popmap.tsv").exists():
            continue
        try:
            resolve_markers_table(path)
        except FileNotFoundError:
            continue
        datasets.append(path.name)
    return sorted(datasets)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workdir", default=Path("benchmarks/literature-workdir"), type=Path)
    parser.add_argument("--results", default=Path("benchmarks/results/literature_mode_effects.csv"), type=Path)
    parser.add_argument("--output-dir", default=Path("benchmarks/literature-workdir/modes"), type=Path)
    parser.add_argument("--dataset", action="append", help="Dataset name to analyze; repeat for multiple datasets")
    parser.add_argument("--min-depth", default=10, type=int)
    parser.add_argument("--group1", default="male")
    parser.add_argument("--group2", default="female")
    parser.add_argument("--n-components", default=4, type=int)
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
                args.group1,
                args.group2,
                args.n_components,
            )
        )
    write_rows(args.results, rows)
    print(f"Wrote {args.results}")


if __name__ == "__main__":
    if __package__ in {None, ""}:
        sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    main()
