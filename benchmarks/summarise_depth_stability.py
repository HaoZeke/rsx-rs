#!/usr/bin/env python3
"""Collapse the (dataset, depth) low-depth sweep into a stability summary.

Reads benchmarks/results/literature_mode_effects_sweep.csv (128 rows, the
output of running analyze_literature_modes.py at min_depth in {3, 5, 8, 10}
across the four literature datasets) and emits, per (dataset, depth):

- significant_markers (strict chisq/Bonferroni)
- posterior_gt_0_9    (Bayesian)
- bayes_factor_gt_10  (Bayesian)

Writes:
- benchmarks/results/literature_depth_stability.csv
- docs/figures/literature_depth_stability.{svg,pdf} (line plot per dataset)
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("benchmarks/results/literature_mode_effects_sweep.csv"),
    )
    parser.add_argument(
        "--summary",
        type=Path,
        default=Path("benchmarks/results/literature_depth_stability.csv"),
    )
    parser.add_argument(
        "--figure-dir",
        type=Path,
        default=Path("docs/figures"),
    )
    args = parser.parse_args()

    df = pd.read_csv(args.input)

    # Keep only the rows that carry call counts.
    call_modes = {
        "marker_extraction_chisq_bonferroni": "significant_markers",
        "bayesian_marker_table": "bayesian_marker_table",
    }
    sub = df[df["mode"].isin(call_modes)].copy()

    def _as_int(value) -> int:
        if value is None or pd.isna(value):
            return 0
        return int(value)

    rows = []
    for (dataset, min_depth), group in sub.groupby(["dataset", "min_depth"]):
        record = {"dataset": dataset, "min_depth": int(min_depth)}
        for _, r in group.iterrows():
            if r["mode"] == "marker_extraction_chisq_bonferroni":
                record["significant_markers"] = _as_int(r.get("significant_markers"))
            if r["mode"] == "bayesian_marker_table":
                record["posterior_gt_0_9"] = _as_int(r.get("posterior_gt_0_9"))
                record["bayes_factor_gt_10"] = _as_int(r.get("bayes_factor_gt_10"))
        rows.append(record)

    summary = pd.DataFrame(rows).sort_values(["dataset", "min_depth"]).reset_index(drop=True)
    summary = summary[
        [
            "dataset",
            "min_depth",
            "significant_markers",
            "posterior_gt_0_9",
            "bayes_factor_gt_10",
        ]
    ]
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.summary, index=False)
    print(f"Wrote {args.summary} ({len(summary)} rows)")
    print(summary.to_string(index=False))

    try:
        from plotnine import (
            aes,
            element_blank,
            element_rect,
            element_text,
            facet_wrap,
            geom_line,
            geom_point,
            ggplot,
            labs,
            scale_y_log10,
            theme,
            theme_minimal,
        )
    except ImportError:
        print("plotnine not available; skipping figure")
        return 0

    long = summary.melt(
        id_vars=["dataset", "min_depth"],
        value_vars=["significant_markers", "posterior_gt_0_9", "bayes_factor_gt_10"],
        var_name="call",
        value_name="count",
    )
    long["count_plot"] = long["count"].clip(lower=0.5)

    palette = {
        "significant_markers": "#004D40",
        "posterior_gt_0_9": "#FF655D",
        "bayes_factor_gt_10": "#F1DB4B",
    }
    label_map = {
        "significant_markers": "Strict (chisq/Bonferroni)",
        "posterior_gt_0_9": "Posterior > 0.9",
        "bayes_factor_gt_10": "Bayes factor > 10",
    }
    long["call_label"] = long["call"].map(label_map)

    p = (
        ggplot(long, aes(x="min_depth", y="count_plot", color="call_label", group="call_label"))
        + geom_line(size=1.0)
        + geom_point(size=2.0)
        + facet_wrap("~dataset", scales="free_y", ncol=2)
        + scale_y_log10()
        + labs(
            x="min_depth",
            y="marker count (log scale, zero shown as 0.5)",
            color="call",
            title="Depth stability of the call layers",
            subtitle="Per-dataset marker counts across min_depth in {3, 5, 8, 10}",
        )
        + theme_minimal()
        + theme(
            figure_size=(8.4, 6.0),
            panel_grid_minor=element_blank(),
            panel_background=element_rect(fill="white"),
            strip_text=element_text(face="bold"),
        )
    )

    args.figure_dir.mkdir(parents=True, exist_ok=True)
    svg = args.figure_dir / "literature_depth_stability.svg"
    pdf = args.figure_dir / "literature_depth_stability.pdf"
    p.save(str(svg), verbose=False)
    p.save(str(pdf), verbose=False)
    print(f"Wrote {svg}")
    print(f"Wrote {pdf}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
