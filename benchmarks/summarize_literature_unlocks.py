#!/usr/bin/env python3
"""Summarize biological gains unlocked by rsx mode outputs."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from benchmarks.analyze_literature_modes import as_float, as_int, read_table


SUMMARY_COLUMNS = [
    "dataset",
    "min_depth",
    "tested_markers",
    "strict_candidates",
    "posterior_gt_0_9",
    "strict_and_posterior",
    "strict_only",
    "posterior_only",
    "bayes_factor_gt_10",
    "bayes_factor_only",
    "singleton_fraction",
    "pc1_variance_fraction",
    "sex_loading_delta_pc1",
    "unlock_class",
    "biological_interpretation",
]

CANDIDATE_COLUMNS = [
    "dataset",
    "min_depth",
    "rank",
    "id",
    "sequence",
    "strict_call",
    "posterior_gt_0_9",
    "bayes_factor_gt_10",
    "posterior_sex_linked",
    "bayes_factor",
    "group1",
    "group1_present",
    "group1_total",
    "group1_penetrance",
    "group2",
    "group2_present",
    "group2_total",
    "group2_penetrance",
    "bias_direction",
    "candidate_class",
]


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


def marker_key(row: dict[str, str]) -> str:
    return f"{row.get('id', '')}\t{row.get('sequence', '')}"


def format_float(value: float) -> str:
    if value == 0.0:
        return "0"
    return f"{value:.6g}"


def parse_candidate_table(
    path: Path,
    popmap: dict[str, str],
    min_depth: int,
    group1: str,
    group2: str,
    strict_call: bool,
) -> dict[str, dict[str, object]]:
    _, header, rows = read_table(path)
    sample_columns = [
        column
        for column in header
        if column not in {"id", "sequence", "Bayes_Factor", "Posterior_SexLinked"}
    ]
    group1_total = sum(1 for sample in sample_columns if popmap.get(sample) == group1)
    group2_total = sum(1 for sample in sample_columns if popmap.get(sample) == group2)

    candidates: dict[str, dict[str, object]] = {}
    for row in rows:
        group1_present = sum(
            1 for sample in sample_columns if popmap.get(sample) == group1 and as_float(row.get(sample)) >= min_depth
        )
        group2_present = sum(
            1 for sample in sample_columns if popmap.get(sample) == group2 and as_float(row.get(sample)) >= min_depth
        )
        if group1_present > group2_present:
            direction = f"{group1}-biased"
        elif group2_present > group1_present:
            direction = f"{group2}-biased"
        else:
            direction = "balanced"
        posterior = as_float(row.get("Posterior_SexLinked"))
        bayes_factor = as_float(row.get("Bayes_Factor"))
        key = marker_key(row)
        candidates[key] = {
            "id": row.get("id", ""),
            "sequence": row.get("sequence", ""),
            "strict_call": strict_call,
            "posterior_sex_linked": posterior,
            "bayes_factor": bayes_factor,
            "group1": group1,
            "group1_present": group1_present,
            "group1_total": group1_total,
            "group1_penetrance": group1_present / group1_total if group1_total else 0.0,
            "group2": group2,
            "group2_present": group2_present,
            "group2_total": group2_total,
            "group2_penetrance": group2_present / group2_total if group2_total else 0.0,
            "bias_direction": direction,
        }
    return candidates


def merge_candidate_tables(
    strict_candidates: dict[str, dict[str, object]],
    bayes_candidates: dict[str, dict[str, object]],
) -> dict[str, dict[str, object]]:
    merged = {key: dict(value) for key, value in bayes_candidates.items()}
    for key, value in strict_candidates.items():
        if key in merged:
            merged[key]["strict_call"] = True
        else:
            merged[key] = dict(value)
    return merged


def classify_candidate(candidate: dict[str, object]) -> str:
    strict_call = bool(candidate.get("strict_call"))
    posterior = float(candidate.get("posterior_sex_linked") or 0.0)
    bayes_factor = float(candidate.get("bayes_factor") or 0.0)
    if strict_call and posterior > 0.9:
        return "strict+posterior"
    if strict_call:
        return "strict_only"
    if posterior > 0.9:
        return "posterior_only"
    if bayes_factor > 10.0:
        return "bayes_factor_only"
    return "exploratory"


def candidate_sort_key(candidate: dict[str, object]) -> tuple[object, ...]:
    class_order = {
        "strict+posterior": 0,
        "posterior_only": 1,
        "strict_only": 2,
        "bayes_factor_only": 3,
        "exploratory": 4,
    }
    candidate_class = classify_candidate(candidate)
    penetrance_delta = abs(float(candidate["group1_penetrance"]) - float(candidate["group2_penetrance"]))
    marker_id = str(candidate.get("id", ""))
    numeric_id = int(marker_id) if marker_id.isdigit() else math.inf
    return (
        class_order[candidate_class],
        -float(candidate.get("posterior_sex_linked") or 0.0),
        -float(candidate.get("bayes_factor") or 0.0),
        -penetrance_delta,
        numeric_id,
        marker_id,
    )


def summarize_unlock_class(
    strict_count: int,
    posterior_count: int,
    strict_and_posterior: int,
    posterior_only: int,
    bayes_factor_only: int,
) -> str:
    if strict_and_posterior > 0 and posterior_only > 0:
        return "strict recovery plus posterior expansion"
    if strict_and_posterior > 0:
        return "strict recovery with posterior support"
    if strict_count > 0:
        return "strict recovery with candidate ranking"
    if posterior_count > 0:
        return "posterior triage without strict reclassification"
    if bayes_factor_only > 0:
        return "QC restraint on Bayes-factor-only signal"
    return "strict null with QC context"


def biological_interpretation(summary: dict[str, str]) -> str:
    strict_count = as_int(summary["strict_candidates"])
    posterior_only = as_int(summary["posterior_only"])
    bayes_factor_only = as_int(summary["bayes_factor_only"])
    singleton_fraction = as_float(summary.get("singleton_fraction"))
    pc1 = as_float(summary.get("pc1_variance_fraction"))
    pieces: list[str] = []
    if strict_count:
        pieces.append(f"strict extraction recovers {strict_count} marker(s)")
    else:
        pieces.append("strict extraction is empty")
    if posterior_only:
        pieces.append(f"posterior ranking adds {posterior_only} follow-up marker(s)")
    elif as_int(summary["posterior_gt_0_9"]):
        pieces.append("posterior ranking supports the strict marker set")
    else:
        pieces.append("posterior ranking does not add high-posterior markers")
    if bayes_factor_only:
        pieces.append(f"{bayes_factor_only} Bayes-factor-only row(s) are held below the posterior threshold")
    if singleton_fraction:
        pieces.append(f"singleton fraction is {singleton_fraction:.1%}")
    if pc1:
        pieces.append(f"PC1 depth variance is {pc1:.1%}")
    return "; ".join(pieces) + "."


def summarize_dataset(
    workdir: Path,
    mode_dir: Path,
    dataset: str,
    min_depth: int,
    group1: str,
    group2: str,
    mode_effects: dict[str, str],
    top_candidates: int,
) -> tuple[dict[str, str], list[dict[str, str]]]:
    popmap = load_popmap(workdir / dataset / "popmap.tsv")
    dataset_dir = mode_dir / dataset
    strict_table = dataset_dir / f"signif_chisq_bonferroni_d{min_depth}.tsv"
    bayes_table = dataset_dir / f"signif_bayes_chisq_none_d{min_depth}.tsv"
    if not strict_table.exists():
        raise FileNotFoundError(f"strict marker table not found: {strict_table}")
    if not bayes_table.exists():
        raise FileNotFoundError(f"Bayesian marker table not found: {bayes_table}")

    strict_candidates = parse_candidate_table(strict_table, popmap, min_depth, group1, group2, strict_call=True)
    bayes_candidates = parse_candidate_table(bayes_table, popmap, min_depth, group1, group2, strict_call=False)
    candidates = merge_candidate_tables(strict_candidates, bayes_candidates)

    strict_keys = {key for key, value in candidates.items() if bool(value.get("strict_call"))}
    posterior_keys = {
        key for key, value in candidates.items() if float(value.get("posterior_sex_linked") or 0.0) > 0.9
    }
    bayes_factor_keys = {key for key, value in candidates.items() if float(value.get("bayes_factor") or 0.0) > 10.0}
    strict_and_posterior = strict_keys & posterior_keys
    posterior_only = posterior_keys - strict_keys
    bayes_factor_only = {
        key
        for key in bayes_factor_keys
        if key not in strict_keys and float(candidates[key].get("posterior_sex_linked") or 0.0) <= 0.9
    }

    summary = {
        "dataset": dataset,
        "min_depth": str(min_depth),
        "tested_markers": str(as_int(mode_effects.get("tested_markers"))),
        "strict_candidates": str(len(strict_keys)),
        "posterior_gt_0_9": str(len(posterior_keys)),
        "strict_and_posterior": str(len(strict_and_posterior)),
        "strict_only": str(len(strict_keys - posterior_keys)),
        "posterior_only": str(len(posterior_only)),
        "bayes_factor_gt_10": str(len(bayes_factor_keys)),
        "bayes_factor_only": str(len(bayes_factor_only)),
        "singleton_fraction": format_float(as_float(mode_effects.get("singleton_fraction"))),
        "pc1_variance_fraction": format_float(as_float(mode_effects.get("pc1_variance_fraction"))),
        "sex_loading_delta_pc1": format_float(as_float(mode_effects.get("sex_loading_delta_pc1"))),
    }
    summary["unlock_class"] = summarize_unlock_class(
        len(strict_keys),
        len(posterior_keys),
        len(strict_and_posterior),
        len(posterior_only),
        len(bayes_factor_only),
    )
    summary["biological_interpretation"] = biological_interpretation(summary)

    candidate_rows: list[dict[str, str]] = []
    for rank, candidate in enumerate(sorted(candidates.values(), key=candidate_sort_key)[:top_candidates], start=1):
        candidate_class = classify_candidate(candidate)
        posterior = float(candidate.get("posterior_sex_linked") or 0.0)
        bayes_factor = float(candidate.get("bayes_factor") or 0.0)
        row = {
            "dataset": dataset,
            "min_depth": str(min_depth),
            "rank": str(rank),
            "id": str(candidate.get("id", "")),
            "sequence": str(candidate.get("sequence", "")),
            "strict_call": str(bool(candidate.get("strict_call"))),
            "posterior_gt_0_9": str(posterior > 0.9),
            "bayes_factor_gt_10": str(bayes_factor > 10.0),
            "posterior_sex_linked": format_float(posterior),
            "bayes_factor": format_float(bayes_factor),
            "group1": str(candidate["group1"]),
            "group1_present": str(candidate["group1_present"]),
            "group1_total": str(candidate["group1_total"]),
            "group1_penetrance": format_float(float(candidate["group1_penetrance"])),
            "group2": str(candidate["group2"]),
            "group2_present": str(candidate["group2_present"]),
            "group2_total": str(candidate["group2_total"]),
            "group2_penetrance": format_float(float(candidate["group2_penetrance"])),
            "bias_direction": str(candidate["bias_direction"]),
            "candidate_class": candidate_class,
        }
        row[f"{candidate['group1']}_penetrance"] = row["group1_penetrance"]
        row[f"{candidate['group2']}_penetrance"] = row["group2_penetrance"]
        candidate_rows.append(row)
    return summary, candidate_rows


def load_mode_effects(path: Path) -> dict[tuple[str, int], dict[str, str]]:
    if not path.exists():
        return {}
    result: dict[tuple[str, int], dict[str, str]] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            key = (row["dataset"], as_int(row.get("min_depth")))
            target = result.setdefault(key, {})
            mode = row.get("mode")
            if mode == "bayesian_marker_table":
                target["tested_markers"] = row.get("tested_markers", "")
            elif mode == "frequency_qc":
                target["singleton_fraction"] = row.get("singleton_fraction", "")
            elif mode == "streaming_pca":
                target["pc1_variance_fraction"] = row.get("pc1_variance_fraction", "")
                target["sex_loading_delta_pc1"] = row.get("sex_loading_delta_pc1", "")
    return result


def discover_datasets(mode_dir: Path, selected: list[str] | None) -> list[str]:
    if selected:
        return selected
    return sorted(path.name for path in mode_dir.iterdir() if path.is_dir())


def write_rows(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    extras = sorted({key for row in rows for key in row if key not in columns})
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns + extras, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workdir", default=Path("benchmarks/literature-workdir"), type=Path)
    parser.add_argument("--mode-dir", default=Path("benchmarks/literature-workdir/modes"), type=Path)
    parser.add_argument("--mode-effects", default=Path("benchmarks/results/literature_mode_effects.csv"), type=Path)
    parser.add_argument("--summary", default=Path("benchmarks/results/literature_bio_unlocks.csv"), type=Path)
    parser.add_argument("--candidates", default=Path("benchmarks/results/literature_candidate_triage.csv"), type=Path)
    parser.add_argument("--dataset", action="append", help="Dataset name to summarize; repeat for multiple datasets")
    parser.add_argument("--min-depth", default=10, type=int)
    parser.add_argument("--group1", default="male")
    parser.add_argument("--group2", default="female")
    parser.add_argument("--top-candidates", default=30, type=int)
    args = parser.parse_args()

    mode_effects = load_mode_effects(args.mode_effects)
    summaries: list[dict[str, str]] = []
    candidates: list[dict[str, str]] = []
    for dataset in discover_datasets(args.mode_dir, args.dataset):
        summary, candidate_rows = summarize_dataset(
            args.workdir,
            args.mode_dir,
            dataset,
            args.min_depth,
            args.group1,
            args.group2,
            mode_effects.get((dataset, args.min_depth), {}),
            args.top_candidates,
        )
        summaries.append(summary)
        candidates.extend(candidate_rows)
    write_rows(args.summary, summaries, SUMMARY_COLUMNS)
    write_rows(args.candidates, candidates, CANDIDATE_COLUMNS)
    print(f"Wrote {args.summary} and {args.candidates}")


if __name__ == "__main__":
    if __package__ in {None, ""}:
        sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    main()
