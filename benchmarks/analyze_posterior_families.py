#!/usr/bin/env python3
"""Evaluate configured posterior families on literature distribution tables."""

from __future__ import annotations

import argparse
import csv
import math
import tomllib
from dataclasses import dataclass
from pathlib import Path

try:
    from benchmarks.analyze_bayesian_evidence import (
        PrevalencePrior,
        read_distrib,
        read_popmap,
        posterior_sex_linked_with_model,
    )
except ModuleNotFoundError:
    from analyze_bayesian_evidence import (
        PrevalencePrior,
        read_distrib,
        read_popmap,
        posterior_sex_linked_with_model,
    )


RESULT_COLUMNS = [
    "dataset",
    "min_depth",
    "profile",
    "linkage_prior",
    "group1_linked_weight",
    "linked_family",
    "linked_probability",
    "linked_alpha",
    "linked_beta",
    "null_family",
    "null_probability",
    "null_alpha",
    "null_beta",
    "threshold",
    "markers",
    "markers_posterior_gt_threshold",
    "posterior_marker_mass",
    "mean_posterior",
    "source_table",
]


@dataclass(frozen=True)
class PosteriorProfile:
    name: str
    linkage_prior: float
    group1_linked_weight: float
    linked: PrevalencePrior
    null: PrevalencePrior


def _check_keys(mapping: dict, allowed: set[str], location: str) -> None:
    unexpected = sorted(set(mapping) - allowed)
    if unexpected:
        raise ValueError(f"unexpected {location} keys: {', '.join(unexpected)}")


def _read_prevalence(value: object, location: str) -> PrevalencePrior:
    if not isinstance(value, dict):
        raise ValueError(f"{location} must be a TOML table")
    family = value.get("family")
    if family == "fixed":
        _check_keys(value, {"family", "probability"}, location)
        if "probability" not in value:
            raise ValueError(f"{location}.probability is required for fixed family")
        return PrevalencePrior.fixed(float(value["probability"]))
    if family == "beta":
        _check_keys(value, {"family", "alpha", "beta"}, location)
        if "alpha" not in value or "beta" not in value:
            raise ValueError(f"{location}.alpha and {location}.beta are required for beta family")
        return PrevalencePrior.beta(float(value["alpha"]), float(value["beta"]))
    raise ValueError(f"{location}.family must be 'fixed' or 'beta'")


def read_profiles(path: Path) -> tuple[float, list[PosteriorProfile]]:
    with path.open("rb") as handle:
        document = tomllib.load(handle)
    _check_keys(document, {"threshold", "profile"}, "top-level")
    threshold = float(document.get("threshold", 0.9))
    if not math.isfinite(threshold) or not 0.0 < threshold < 1.0:
        raise ValueError("threshold must be strictly between zero and one")
    raw_profiles = document.get("profile")
    if not isinstance(raw_profiles, list) or not raw_profiles:
        raise ValueError("at least one [[profile]] table is required")

    profiles: list[PosteriorProfile] = []
    names: set[str] = set()
    for index, raw in enumerate(raw_profiles):
        if not isinstance(raw, dict):
            raise ValueError(f"profile[{index}] must be a TOML table")
        _check_keys(
            raw,
            {"name", "linkage_prior", "group1_linked_weight", "linked", "null"},
            f"profile[{index}]",
        )
        missing = [
            key
            for key in ("name", "linkage_prior", "group1_linked_weight", "linked", "null")
            if key not in raw
        ]
        if missing:
            raise ValueError(f"profile[{index}] is missing: {', '.join(missing)}")
        name = str(raw["name"])
        if not name or name in names:
            raise ValueError(f"profile names must be non-empty and unique: {name!r}")
        names.add(name)
        linkage_prior = float(raw["linkage_prior"])
        group1_linked_weight = float(raw["group1_linked_weight"])
        if not 0.0 < linkage_prior < 1.0:
            raise ValueError(f"profile[{index}].linkage_prior must lie between zero and one")
        if not 0.0 < group1_linked_weight < 1.0:
            raise ValueError(f"profile[{index}].group1_linked_weight must lie between zero and one")
        profiles.append(
            PosteriorProfile(
                name=name,
                linkage_prior=linkage_prior,
                group1_linked_weight=group1_linked_weight,
                linked=_read_prevalence(raw["linked"], f"profile[{index}].linked"),
                null=_read_prevalence(raw["null"], f"profile[{index}].null"),
            )
        )
    return threshold, profiles


def _prior_columns(prior: PrevalencePrior, prefix: str) -> dict[str, str | float]:
    return {
        f"{prefix}_family": prior.family,
        f"{prefix}_probability": "" if prior.probability is None else prior.probability,
        f"{prefix}_alpha": "" if prior.alpha is None else prior.alpha,
        f"{prefix}_beta": "" if prior.beta_shape is None else prior.beta_shape,
    }


def analyze_profiles(
    workdir: Path,
    datasets: list[str],
    min_depths: list[int],
    threshold: float,
    profiles: list[PosteriorProfile],
) -> list[dict[str, str | int | float]]:
    rows: list[dict[str, str | int | float]] = []
    for dataset in datasets:
        dataset_dir = workdir / dataset
        group_totals = read_popmap(dataset_dir / "popmap.tsv")
        for min_depth in min_depths:
            table_path = dataset_dir / f"distrib_{min_depth}.tsv"
            group1, group2, cells = read_distrib(table_path)
            total_g1 = group_totals[group1]
            total_g2 = group_totals[group2]
            marker_total = sum(cell.markers for cell in cells)
            for profile in profiles:
                qualifying = 0
                posterior_mass = 0.0
                for cell in cells:
                    posterior = posterior_sex_linked_with_model(
                        cell.group1,
                        cell.group2,
                        total_g1,
                        total_g2,
                        linkage_prior=profile.linkage_prior,
                        group1_linked_weight=profile.group1_linked_weight,
                        linked=profile.linked,
                        null=profile.null,
                    )
                    posterior_mass += posterior * cell.markers
                    if posterior > threshold:
                        qualifying += cell.markers
                rows.append(
                    {
                        "dataset": dataset,
                        "min_depth": min_depth,
                        "profile": profile.name,
                        "linkage_prior": profile.linkage_prior,
                        "group1_linked_weight": profile.group1_linked_weight,
                        **_prior_columns(profile.linked, "linked"),
                        **_prior_columns(profile.null, "null"),
                        "threshold": threshold,
                        "markers": marker_total,
                        "markers_posterior_gt_threshold": qualifying,
                        "posterior_marker_mass": f"{posterior_mass:.9f}",
                        "mean_posterior": f"{posterior_mass / marker_total if marker_total else 0.0:.12f}",
                        "source_table": str(table_path.relative_to(workdir)),
                    }
                )
    return rows


def parse_depths(value: str) -> list[int]:
    depths = [int(part) for part in value.split(",") if part.strip()]
    if not depths or any(depth < 1 for depth in depths):
        raise argparse.ArgumentTypeError("minimum depths must be positive integers")
    return depths


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workdir", type=Path, default=Path("benchmarks/literature-workdir"))
    parser.add_argument("--profiles", type=Path, default=Path("benchmarks/posterior_family_profiles.toml"))
    parser.add_argument(
        "--results",
        type=Path,
        default=Path("benchmarks/results/literature_posterior_family_sensitivity.csv"),
    )
    parser.add_argument("--dataset", action="append", dest="datasets", required=True)
    parser.add_argument("--min-depths", type=parse_depths, default=[10])
    args = parser.parse_args()

    threshold, profiles = read_profiles(args.profiles)
    rows = analyze_profiles(args.workdir, args.datasets, args.min_depths, threshold, profiles)
    args.results.parent.mkdir(parents=True, exist_ok=True)
    with args.results.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=RESULT_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} rows to {args.results}")


if __name__ == "__main__":
    main()
