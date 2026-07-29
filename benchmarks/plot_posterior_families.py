#!/usr/bin/env python3
"""Render posterior-family sensitivity across literature datasets."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from plotnine import (
    aes,
    element_blank,
    element_rect,
    element_text,
    geom_text,
    geom_tile,
    ggplot,
    labs,
    scale_fill_gradient,
    theme,
    theme_minimal,
)


DATASET_LABELS = {
    "danio_albolineatus": "D. albolineatus",
    "notothenia_rossii": "N. rossii",
    "plecoglossus_altivelis": "P. altivelis",
    "tinca_tinca": "T. tinca",
}


def prepare_plot_rows(rows: list[dict] | pd.DataFrame) -> pd.DataFrame:
    frame = pd.DataFrame(rows).copy()
    required = {"dataset", "profile", "linked_family", "markers_posterior_gt_threshold"}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"posterior-family results are missing: {', '.join(missing)}")
    frame["markers_posterior_gt_threshold"] = pd.to_numeric(
        frame["markers_posterior_gt_threshold"], errors="raise"
    ).astype(int)
    frame["plot_count"] = frame["markers_posterior_gt_threshold"] + 1
    frame["count_label"] = frame["markers_posterior_gt_threshold"].astype(str)
    frame["dataset_label"] = frame["dataset"].map(DATASET_LABELS).fillna(frame["dataset"])
    frame["family_label"] = frame["linked_family"].map(
        {"fixed": "Fixed prevalence", "beta": "Beta integrated"}
    )
    if frame["family_label"].isna().any():
        raise ValueError("posterior-family results contain an unsupported linked family")
    return frame


def render_plot(frame: pd.DataFrame, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    plot = (
        ggplot(frame, aes("profile", "dataset_label", fill="plot_count"))
        + geom_tile(color="white", size=0.6)
        + geom_text(aes(label="count_label"), size=8)
        + scale_fill_gradient(low="#F7FBFF", high="#006D77", trans="log10", name="Markers + 1")
        + labs(
            x="Configured posterior profile",
            y="Literature dataset",
            title="Posterior-family sensitivity on real RAD-seq marker tables",
            subtitle="Cell labels give markers with posterior probability above 0.9 at depth 10",
        )
        + theme_minimal()
        + theme(
            figure_size=(9.0, 4.5),
            panel_grid=element_blank(),
            panel_background=element_rect(fill="white"),
            axis_text_x=element_text(rotation=25, ha="right"),
        )
    )
    plot.save(str(output_dir / "literature_posterior_family_sensitivity.svg"), verbose=False)
    plot.save(str(output_dir / "literature_posterior_family_sensitivity.pdf"), verbose=False)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("benchmarks/results/literature_posterior_family_sensitivity.csv"),
    )
    parser.add_argument("--output", type=Path, default=Path("docs/figures"))
    args = parser.parse_args()

    frame = prepare_plot_rows(pd.read_csv(args.input))
    render_plot(frame, args.output)
    print(f"Wrote posterior-family figures under {args.output}")


if __name__ == "__main__":
    main()
