#!/usr/bin/env python3
"""Summarize Bayesian sex-linkage evidence from literature distrib outputs."""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


RESULT_COLUMNS = [
    "dataset",
    "min_depth",
    "group1",
    "group2",
    "total_group1",
    "total_group2",
    "cells",
    "markers",
    "markers_bf_gt_10",
    "markers_bf_gt_100",
    "markers_posterior_gt_0_5",
    "markers_posterior_gt_0_9",
    "posterior_marker_mass",
    "mean_posterior",
    "top_cell_group1",
    "top_cell_group2",
    "top_cell_markers",
    "top_cell_bayes_factor",
    "top_cell_posterior",
    "input_path",
]


@dataclass(frozen=True)
class DistribCell:
    group1: int
    group2: int
    markers: int


@dataclass(frozen=True)
class PrevalencePrior:
    """Fixed or Beta prevalence family matching ``rsxcore::PrevalencePrior``."""

    family: str
    probability: float | None = None
    alpha: float | None = None
    beta_shape: float | None = None

    @classmethod
    def fixed(cls, probability: float) -> "PrevalencePrior":
        if not math.isfinite(probability) or not 0.0 < probability < 1.0:
            raise ValueError("fixed prevalence probability must be strictly between zero and one")
        return cls(family="fixed", probability=probability)

    @classmethod
    def beta(cls, alpha: float, beta: float) -> "PrevalencePrior":
        if not math.isfinite(alpha) or not math.isfinite(beta) or alpha <= 0.0 or beta <= 0.0:
            raise ValueError("Beta prevalence shapes must be finite and greater than zero")
        return cls(family="beta", alpha=alpha, beta_shape=beta)


def parse_min_depths(value: str) -> list[int]:
    depths = [int(item) for item in value.split(",") if item.strip()]
    if not depths:
        raise argparse.ArgumentTypeError("at least one minimum depth is required")
    return depths


def read_popmap(path: Path) -> dict[str, int]:
    groups: dict[str, int] = {}
    with path.open() as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            fields = stripped.split()
            if len(fields) < 2:
                continue
            groups[fields[1]] = groups.get(fields[1], 0) + 1
    if len(groups) < 2:
        raise ValueError(f"{path} must contain at least two groups")
    return groups


def iter_data_lines(path: Path) -> Iterable[list[str]]:
    with path.open() as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            yield stripped.split("\t")


def read_distrib(path: Path) -> tuple[str, str, list[DistribCell]]:
    iterator = iter_data_lines(path)
    try:
        header = next(iterator)
    except StopIteration as error:
        raise ValueError(f"{path} does not contain a distribution table") from error
    if len(header) < 3 or header[2] != "Markers":
        raise ValueError(f"{path} has an unexpected header: {header}")
    group1, group2 = header[0], header[1]
    cells: list[DistribCell] = []
    for fields in iterator:
        if len(fields) < 3:
            continue
        cells.append(DistribCell(int(fields[0]), int(fields[1]), int(fields[2])))
    return group1, group2, cells


def log_beta_binom(k: int, n: int) -> float:
    return math.lgamma(k + 1) + math.lgamma(n - k + 1) - math.lgamma(n + 2)


def bayes_factor_2x2(n_g1: int, n_g2: int, total_g1: int, total_g2: int) -> float:
    log_h1 = log_beta_binom(n_g1, total_g1) + log_beta_binom(n_g2, total_g2)
    log_h0 = log_beta_binom(n_g1 + n_g2, total_g1 + total_g2)
    return math.exp(log_h1 - log_h0)


def logsumexp2(a: float, b: float) -> float:
    max_value = max(a, b)
    if math.isinf(max_value):
        return max_value
    return max_value + math.log(math.exp(a - max_value) + math.exp(b - max_value))


def log_prevalence_marginal(k: int, n: int, prior: PrevalencePrior) -> float:
    if prior.family == "fixed":
        probability = prior.probability
        if probability is None:
            raise ValueError("fixed prevalence family requires probability")
        return k * math.log(probability) + (n - k) * math.log1p(-probability)
    if prior.family == "beta":
        if prior.alpha is None or prior.beta_shape is None:
            raise ValueError("Beta prevalence family requires alpha and beta")
        return (
            math.lgamma(k + prior.alpha)
            + math.lgamma(n - k + prior.beta_shape)
            - math.lgamma(n + prior.alpha + prior.beta_shape)
            - math.lgamma(prior.alpha)
            - math.lgamma(prior.beta_shape)
            + math.lgamma(prior.alpha + prior.beta_shape)
        )
    raise ValueError(f"unsupported prevalence family: {prior.family}")


def posterior_sex_linked_with_model(
    n_g1: int,
    n_g2: int,
    total_g1: int,
    total_g2: int,
    *,
    linkage_prior: float,
    group1_linked_weight: float,
    linked: PrevalencePrior,
    null: PrevalencePrior,
) -> float:
    """Evaluate the directional model used by the Rust implementation."""

    if not 0 <= n_g1 <= total_g1 or not 0 <= n_g2 <= total_g2:
        raise ValueError("marker counts must lie within their group totals")
    if not math.isfinite(linkage_prior) or not 0.0 < linkage_prior < 1.0:
        raise ValueError("linkage_prior must be strictly between zero and one")
    if not math.isfinite(group1_linked_weight) or not 0.0 < group1_linked_weight < 1.0:
        raise ValueError("group1_linked_weight must be strictly between zero and one")

    linked_total = total_g1 + total_g2
    ll_g1_linked = log_prevalence_marginal(
        n_g1 + (total_g2 - n_g2), linked_total, linked
    )
    ll_g2_linked = log_prevalence_marginal(
        (total_g1 - n_g1) + n_g2, linked_total, linked
    )
    ll_linked = logsumexp2(
        ll_g1_linked + math.log(group1_linked_weight),
        ll_g2_linked + math.log1p(-group1_linked_weight),
    )
    ll_null = log_prevalence_marginal(n_g1 + n_g2, linked_total, null)
    log_odds = ll_linked - ll_null + math.log(linkage_prior / (1.0 - linkage_prior))
    if log_odds > 20.0:
        return 1.0
    if log_odds < -20.0:
        return 0.0
    return 1.0 / (1.0 + math.exp(-log_odds))


def posterior_sex_linked(
    n_g1: int,
    n_g2: int,
    total_g1: int,
    total_g2: int,
    pi: float,
    p_sex: float,
) -> float:
    return posterior_sex_linked_with_model(
        n_g1,
        n_g2,
        total_g1,
        total_g2,
        linkage_prior=pi,
        group1_linked_weight=0.5,
        linked=PrevalencePrior.fixed(p_sex),
        null=PrevalencePrior.fixed(0.5),
    )


def analyze_distrib(
    dataset: str,
    min_depth: int,
    distrib_path: Path,
    popmap_path: Path,
    pi: float,
    p_sex: float,
) -> dict[str, str]:
    group1, group2, cells = read_distrib(distrib_path)
    group_totals = read_popmap(popmap_path)
    total_g1 = group_totals[group1]
    total_g2 = group_totals[group2]
    markers = sum(cell.markers for cell in cells)
    bf10 = 0
    bf100 = 0
    post50 = 0
    post90 = 0
    posterior_mass = 0.0
    top_cell: tuple[DistribCell, float, float] | None = None

    for cell in cells:
        bf = bayes_factor_2x2(cell.group1, cell.group2, total_g1, total_g2)
        posterior = posterior_sex_linked(cell.group1, cell.group2, total_g1, total_g2, pi, p_sex)
        if bf > 10.0:
            bf10 += cell.markers
        if bf > 100.0:
            bf100 += cell.markers
        if posterior > 0.5:
            post50 += cell.markers
        if posterior > 0.9:
            post90 += cell.markers
        posterior_mass += posterior * cell.markers
        candidate = (cell, bf, posterior)
        if top_cell is None or (posterior, bf, cell.markers) > (top_cell[2], top_cell[1], top_cell[0].markers):
            top_cell = candidate

    top = top_cell or (DistribCell(0, 0, 0), 0.0, 0.0)
    mean_posterior = posterior_mass / markers if markers else 0.0
    return {
        "dataset": dataset,
        "min_depth": str(min_depth),
        "group1": group1,
        "group2": group2,
        "total_group1": str(total_g1),
        "total_group2": str(total_g2),
        "cells": str(len(cells)),
        "markers": str(markers),
        "markers_bf_gt_10": str(bf10),
        "markers_bf_gt_100": str(bf100),
        "markers_posterior_gt_0_5": str(post50),
        "markers_posterior_gt_0_9": str(post90),
        "posterior_marker_mass": f"{posterior_mass:.6f}",
        "mean_posterior": f"{mean_posterior:.9f}",
        "top_cell_group1": str(top[0].group1),
        "top_cell_group2": str(top[0].group2),
        "top_cell_markers": str(top[0].markers),
        "top_cell_bayes_factor": f"{top[1]:.9g}",
        "top_cell_posterior": f"{top[2]:.9f}",
        "input_path": str(distrib_path),
    }


def discover_datasets(workdir: Path, selected: list[str] | None) -> list[str]:
    if selected:
        return selected
    return sorted(path.name for path in workdir.iterdir() if (path / "popmap.tsv").exists())


def analyze_workdir(
    workdir: Path,
    datasets: list[str],
    min_depths: list[int],
    pi: float,
    p_sex: float,
    allow_missing: bool,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    missing: list[Path] = []
    for dataset in datasets:
        dataset_dir = workdir / dataset
        popmap_path = dataset_dir / "popmap.tsv"
        if not popmap_path.exists():
            missing.append(popmap_path)
            continue
        for min_depth in min_depths:
            distrib_path = dataset_dir / f"distrib_{min_depth}.tsv"
            if not distrib_path.exists():
                missing.append(distrib_path)
                continue
            rows.append(analyze_distrib(dataset, min_depth, distrib_path, popmap_path, pi, p_sex))
    if missing and not allow_missing:
        raise FileNotFoundError("missing Bayesian evidence inputs:\n" + "\n".join(str(path) for path in missing))
    return rows


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=RESULT_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workdir", default=Path("benchmarks/literature-workdir"), type=Path)
    parser.add_argument("--results", default=Path("benchmarks/results/literature_bayesian_evidence.csv"), type=Path)
    parser.add_argument("--dataset", action="append", help="Dataset name to analyze; repeat for multiple datasets")
    parser.add_argument("--min-depths", default="1,2,5,10", type=parse_min_depths)
    parser.add_argument("--pi", default=0.01, type=float, help="Prior marker sex-linkage probability")
    parser.add_argument("--p-sex", default=0.9, type=float, help="Expected presence probability in the linked sex")
    parser.add_argument("--allow-missing", action="store_true", help="Skip missing dataset outputs")
    args = parser.parse_args()

    datasets = discover_datasets(args.workdir, args.dataset)
    rows = analyze_workdir(args.workdir, datasets, args.min_depths, args.pi, args.p_sex, args.allow_missing)
    write_csv(args.results, rows)
    print(f"Wrote {args.results}")


if __name__ == "__main__":
    main()
