#!/usr/bin/env python3
"""Aggregate repeated CUDA benchmarks and render speedup and timing plots."""

from __future__ import annotations

import argparse
import csv
import io
from pathlib import Path

import pandas as pd
from plotnine import (
    aes,
    element_blank,
    element_rect,
    element_text,
    geom_hline,
    geom_line,
    geom_point,
    ggplot,
    labs,
    scale_color_manual,
    scale_x_log10,
    scale_y_log10,
    theme,
    theme_minimal,
)


NUMERIC_COLUMNS = [
    "markers",
    "cpu_total_s",
    "cuda_setup_s",
    "cuda_h2d_s",
    "cuda_kernel_s",
    "cuda_d2h_s",
    "cuda_total_s",
    "h2d_bytes",
    "d2h_bytes",
    "h2d_gb_s",
    "d2h_gb_s",
    "kernel_speedup",
    "total_speedup",
    "output_buffer_reused",
    "max_abs_error",
]

COLORS = {
    "CPU total": "#D55E00",
    "CUDA total": "#0072B2",
    "CUDA kernel": "#009E73",
    "Kernel": "#009E73",
    "End-to-end": "#0072B2",
}


def load_benchmark(path: Path) -> pd.DataFrame:
    lines = []
    header_seen = False
    with path.open() as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            if line.startswith("markers,"):
                if header_seen:
                    continue
                header_seen = True
            lines.append(line)
    if not header_seen:
        raise ValueError(f"{path} does not contain a CUDA benchmark header")
    rows = list(csv.DictReader(io.StringIO("".join(lines))))
    frame = pd.DataFrame(rows)
    for column in NUMERIC_COLUMNS:
        frame[column] = pd.to_numeric(frame[column], errors="raise")
    return frame


def summarize_benchmark(frame: pd.DataFrame) -> pd.DataFrame:
    if frame.empty:
        raise ValueError("CUDA benchmark contains no measurements")
    if float(frame["max_abs_error"].max()) > 2.0e-15:
        raise ValueError("CUDA benchmark failed CPU/GPU numerical agreement")
    aggregate = {
        column: "median"
        for column in NUMERIC_COLUMNS
        if column not in {"markers", "max_abs_error"}
    }
    aggregate["max_abs_error"] = "max"
    summary = frame.groupby(["device", "markers"], as_index=False).agg(aggregate)
    counts = frame.groupby(["device", "markers"]).size().rename("repetitions").reset_index()
    return summary.merge(counts, on=["device", "markers"]).sort_values(["device", "markers"])


def render_plots(summary: pd.DataFrame, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    speedup = summary.melt(
        id_vars=["device", "markers"],
        value_vars=["kernel_speedup", "total_speedup"],
        var_name="metric",
        value_name="speedup",
    )
    speedup["metric"] = speedup["metric"].map(
        {"kernel_speedup": "Kernel", "total_speedup": "End-to-end"}
    )
    speedup_plot = (
        ggplot(speedup, aes("markers", "speedup", color="metric"))
        + geom_hline(yintercept=10.0, linetype="dashed", color="#666666")
        + geom_line(size=1.0)
        + geom_point(size=2.4)
        + scale_x_log10()
        + scale_y_log10()
        + scale_color_manual(values=COLORS)
        + labs(
            x="Markers evaluated",
            y="Speedup over CPU (log scale)",
            color="Measurement boundary",
            title="CUDA acceleration of chi-squared significance evaluation",
            subtitle="Medians across repeated runs; dashed line marks 10-fold speedup",
        )
        + theme_minimal()
        + theme(
            figure_size=(7.4, 4.8),
            panel_grid_minor=element_blank(),
            panel_background=element_rect(fill="white"),
            legend_title=element_text(face="bold"),
        )
    )

    timing = summary.melt(
        id_vars=["device", "markers"],
        value_vars=["cpu_total_s", "cuda_total_s", "cuda_kernel_s"],
        var_name="phase",
        value_name="seconds",
    )
    timing["phase"] = timing["phase"].map(
        {"cpu_total_s": "CPU total", "cuda_total_s": "CUDA total", "cuda_kernel_s": "CUDA kernel"}
    )
    timing_plot = (
        ggplot(timing, aes("markers", "seconds", color="phase"))
        + geom_line(size=1.0)
        + geom_point(size=2.4)
        + scale_x_log10()
        + scale_y_log10()
        + scale_color_manual(values=COLORS)
        + labs(
            x="Markers evaluated",
            y="Elapsed seconds (log scale)",
            color="Measurement",
            title="CPU and CUDA significance timing boundaries",
            subtitle="CUDA total includes setup, transfers, kernel execution, and result retrieval",
        )
        + theme_minimal()
        + theme(
            figure_size=(7.4, 4.8),
            panel_grid_minor=element_blank(),
            panel_background=element_rect(fill="white"),
            legend_title=element_text(face="bold"),
        )
    )

    for stem, plot in [("cuda_speedup", speedup_plot), ("cuda_timing", timing_plot)]:
        plot.save(str(output_dir / f"{stem}.svg"), verbose=False)
        plot.save(str(output_dir / f"{stem}.pdf"), verbose=False)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument(
        "--summary",
        type=Path,
        default=Path("benchmarks/results/cuda_benchmark_summary.csv"),
    )
    parser.add_argument("--output", type=Path, default=Path("docs/figures"))
    args = parser.parse_args()

    summary = summarize_benchmark(load_benchmark(args.input))
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.summary, index=False)
    render_plots(summary, args.output)
    print(f"Wrote {args.summary} and CUDA figures under {args.output}")


if __name__ == "__main__":
    main()
