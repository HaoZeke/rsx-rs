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
import os
import sys
from pathlib import Path

import pandas as pd


FALLBACK_RUHI_COLORS = {
    "coral": "#FF655D",
    "sunshine": "#F1DB4B",
    "teal": "#004D40",
    "sky": "#1E88E5",
    "magenta": "#D81B60",
}


def load_ruhi_colors() -> dict[str, str]:
    """Match plot_literature_benchmarks.py: prefer the shared chemparseplot
    palette, fall back to the inline RUHI colours otherwise."""
    try:
        from chemparseplot.plot.theme import RUHI_COLORS

        return dict(RUHI_COLORS)
    except ImportError:
        root = os.environ.get("CHEMPARSEPLOT_ROOT")
        if root:
            sys.path.insert(0, root)
            try:
                from chemparseplot.plot.theme import RUHI_COLORS

                return dict(RUHI_COLORS)
            except ImportError:
                pass
    return dict(FALLBACK_RUHI_COLORS)


def base_theme():
    """Shared paper theme, identical to base_theme() in
    plot_literature_benchmarks.py."""
    from plotnine import element_blank, element_line, element_rect, element_text, theme, theme_bw

    return theme_bw(base_size=10) + theme(
        figure_size=(7.2, 4.6),
        panel_grid_major=element_line(color="#E6E6E6", size=0.35),
        panel_grid_minor=element_blank(),
        legend_title=element_text(size=9),
        legend_text=element_text(size=8),
        axis_text_x=element_text(rotation=0),
        axis_text_y=element_text(size=8),
        strip_background=element_rect(fill="#F7F7F7", colour="#CCCCCC"),
    )


def dataset_label(key: str) -> str:
    """Format dataset keys exactly as the sibling plots do, e.g.
    danio_albolineatus -> "Danio Albolineatus"."""
    return key.replace("_", " ").title()


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
            facet_wrap,
            geom_line,
            geom_point,
            ggplot,
            labs,
            scale_color_manual,
            scale_y_log10,
        )
    except ImportError:
        print("plotnine not available; skipping figure")
        return 0

    colors = load_ruhi_colors()

    long = summary.melt(
        id_vars=["dataset", "min_depth"],
        value_vars=["significant_markers", "posterior_gt_0_9", "bayes_factor_gt_10"],
        var_name="call",
        value_name="count",
    )
    long["count_plot"] = long["count"].clip(lower=0.5)
    long["dataset_label"] = long["dataset"].map(dataset_label)

    label_map = {
        "significant_markers": "Strict (chi-square / Bonferroni)",
        "posterior_gt_0_9": "Posterior > 0.9",
        "bayes_factor_gt_10": "Bayes factor > 10",
    }
    call_order = ["Strict (chi-square / Bonferroni)", "Posterior > 0.9", "Bayes factor > 10"]
    long["call_label"] = pd.Categorical(
        long["call"].map(label_map), categories=call_order, ordered=True
    )
    palette = {
        "Strict (chi-square / Bonferroni)": colors["teal"],
        "Posterior > 0.9": colors["magenta"],
        "Bayes factor > 10": colors["sky"],
    }

    p = (
        ggplot(long, aes(x="min_depth", y="count_plot", color="call_label", group="call_label"))
        + geom_line(size=0.7)
        + geom_point(size=1.9)
        + facet_wrap("~ dataset_label", scales="free_y", ncol=2)
        + scale_color_manual(values=palette)
        + scale_y_log10()
        + labs(
            x="Minimum read depth",
            y="Markers (log10, zero shown as 0.5)",
            color="",
            title="Depth stability of the call layers",
            subtitle="Per-dataset marker counts across minimum read depth in {3, 5, 8, 10}",
        )
        + base_theme()
    )

    args.figure_dir.mkdir(parents=True, exist_ok=True)
    width, height = 7.2, 5.2
    svg = args.figure_dir / "literature_depth_stability.svg"
    pdf = args.figure_dir / "literature_depth_stability.pdf"
    p.save(str(svg), width=width, height=height, units="in", dpi=300, verbose=False)
    _strip_trailing_whitespace(svg)
    p.save(str(pdf), width=width, height=height, units="in", dpi=300, verbose=False)
    print(f"Wrote {svg}")
    print(f"Wrote {pdf}")
    return 0


def _strip_trailing_whitespace(path: Path) -> None:
    lines = path.read_text().splitlines()
    path.write_text("\n".join(line.rstrip() for line in lines) + "\n")


if __name__ == "__main__":
    raise SystemExit(main())
