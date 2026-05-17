#!/usr/bin/env python3
"""Generate paper figures from literature benchmark result CSVs."""

from __future__ import annotations

import argparse
import csv
import os
import sys
from collections import defaultdict
from pathlib import Path


FALLBACK_RUHI_COLORS = {
    "coral": "#FF655D",
    "sunshine": "#F1DB4B",
    "teal": "#004D40",
    "sky": "#1E88E5",
    "magenta": "#D81B60",
}


NUMERIC_COLUMNS = {
    "elapsed_seconds",
    "samples",
    "total_spots",
    "total_bases",
    "total_sra_bytes",
    "markers",
    "rows",
    "significant_markers",
    "output_bytes",
}


def load_ruhi_colors() -> dict[str, str]:
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


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def as_float(value: str | None) -> float:
    if value in {None, ""}:
        return 0.0
    return float(value)


def as_int(value: str | None) -> int:
    return int(as_float(value))


def summarize_dataset_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row["dataset"]].append(row)
    summary_rows: list[dict[str, str]] = []
    for dataset, dataset_rows in sorted(grouped.items()):
        metadata = next((row for row in dataset_rows if row["command"] == "metadata"), dataset_rows[0])
        process = next((row for row in dataset_rows if row["command"] == "process"), {})
        download_seconds = sum(as_float(row.get("elapsed_seconds")) for row in dataset_rows if row["command"] == "download")
        process_seconds = sum(as_float(row.get("elapsed_seconds")) for row in dataset_rows if row["command"] == "process")
        analysis_seconds = sum(
            as_float(row.get("elapsed_seconds"))
            for row in dataset_rows
            if row["command"] not in {"metadata", "download", "process"}
        )
        total_seconds = download_seconds + process_seconds + analysis_seconds
        markers = as_int(process.get("markers"))
        total_bases = as_int(metadata.get("total_bases"))
        significant_markers = sum(as_int(row.get("significant_markers")) for row in dataset_rows if row["command"] == "signif")
        summary_rows.append(
            {
                "dataset": dataset,
                "samples": str(as_int(metadata.get("samples"))),
                "spots": str(as_int(metadata.get("total_spots"))),
                "bases": str(total_bases),
                "markers": str(markers),
                "download_seconds": f"{download_seconds:.3f}",
                "process_seconds": f"{process_seconds:.3f}",
                "analysis_seconds": f"{analysis_seconds:.3f}",
                "total_seconds": f"{total_seconds:.3f}",
                "markers_per_second": f"{markers / process_seconds:.3f}" if process_seconds else "",
                "mbases_per_second": f"{total_bases / process_seconds / 1_000_000:.3f}" if process_seconds else "",
                "significant_markers": str(significant_markers),
            }
        )
    return summary_rows


def write_summary(path: Path, rows: list[dict[str, str]]) -> None:
    columns = [
        "dataset",
        "samples",
        "spots",
        "bases",
        "markers",
        "download_seconds",
        "process_seconds",
        "analysis_seconds",
        "total_seconds",
        "markers_per_second",
        "mbases_per_second",
        "significant_markers",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def prepare_results_frame(input_path: Path):
    import pandas as pd
    from siuba import _, group_by, summarize

    frame = pd.read_csv(input_path)
    for column in NUMERIC_COLUMNS.intersection(frame.columns):
        frame[column] = pd.to_numeric(frame[column], errors="coerce").fillna(0.0)
    frame["dataset_label"] = frame["dataset"].str.replace("_", " ").str.title()
    frame["min_depth_label"] = frame["min_depth"].fillna("").astype(str).str.replace(".0", "", regex=False)

    runtime = frame[frame["command"].isin(["download", "process", "depth", "freq", "distrib", "signif"])].copy()
    runtime["phase"] = runtime["command"].map(
        {
            "download": "Download",
            "process": "Process",
            "depth": "Analysis",
            "freq": "Analysis",
            "distrib": "Analysis",
            "signif": "Analysis",
        }
    )
    phase_runtime = runtime >> group_by(_.dataset, _.dataset_label, _.phase) >> summarize(elapsed_seconds=_.elapsed_seconds.sum())
    summary = (
        runtime
        >> group_by(_.dataset, _.dataset_label)
        >> summarize(total_seconds=_.elapsed_seconds.sum(), samples=_.samples.max(), total_bases=_.total_bases.max())
    )
    return frame, phase_runtime, summary


def order_by_total(frame, phase_runtime, summary):
    order = list(summary.sort_values("total_seconds")["dataset_label"])
    frame["dataset_label"] = frame["dataset_label"].astype("category").cat.set_categories(order, ordered=True)
    phase_runtime["dataset_label"] = phase_runtime["dataset_label"].astype("category").cat.set_categories(order, ordered=True)
    summary["dataset_label"] = summary["dataset_label"].astype("category").cat.set_categories(order, ordered=True)
    return frame, phase_runtime, summary


def base_theme():
    from plotnine import element_rect, element_text, theme, theme_bw

    return theme_bw(base_size=10) + theme(
        figure_size=(7.2, 4.6),
        panel_grid_minor=element_rect(fill="white", colour="white"),
        legend_title=element_text(size=9),
        legend_text=element_text(size=8),
        axis_text_x=element_text(rotation=0),
        axis_text_y=element_text(size=8),
        strip_background=element_rect(fill="#F7F7F7", colour="#CCCCCC"),
    )


def save_plot(plot, output_dir: Path, stem: str, width: float = 7.2, height: float = 4.6) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    for suffix in ("svg", "pdf"):
        plot.save(output_dir / f"{stem}.{suffix}", width=width, height=height, units="in", dpi=300, verbose=False)


def plot_runtime_breakdown(phase_runtime, output_dir: Path, colors: dict[str, str]) -> None:
    from plotnine import aes, coord_flip, geom_col, ggplot, labs, scale_fill_manual

    plot_data = phase_runtime.copy()
    plot_data["elapsed_minutes"] = plot_data["elapsed_seconds"] / 60.0
    palette = {
        "Download": colors["sky"],
        "Process": colors["teal"],
        "Analysis": colors["coral"],
    }
    plot = (
        ggplot(plot_data, aes(x="dataset_label", y="elapsed_minutes", fill="phase"))
        + geom_col(width=0.72)
        + coord_flip()
        + scale_fill_manual(values=palette)
        + labs(x="", y="Elapsed minutes", fill="")
        + base_theme()
    )
    save_plot(plot, output_dir, "literature_runtime_breakdown")


def plot_marker_throughput(frame, output_dir: Path, colors: dict[str, str]) -> None:
    from plotnine import aes, coord_flip, geom_col, geom_point, ggplot, labs, scale_color_manual, scale_fill_manual

    process = frame[frame["command"] == "process"].copy()
    process["markers_millions"] = process["markers"] / 1_000_000.0
    process["markers_per_second"] = process["markers"] / process["elapsed_seconds"].where(process["elapsed_seconds"] > 0)
    process["markers_per_second_millions"] = process["markers_per_second"] / 1_000_000.0
    plot = (
        ggplot(process, aes(x="dataset_label", y="markers_millions"))
        + geom_col(fill=colors["teal"], width=0.7)
        + geom_point(aes(y="markers_per_second_millions", color="'Throughput'"), size=2.5)
        + coord_flip()
        + scale_fill_manual(values=[colors["teal"]])
        + scale_color_manual(values={"Throughput": colors["magenta"]})
        + labs(x="", y="Markers or markers/s (millions)", color="")
        + base_theme()
    )
    save_plot(plot, output_dir, "literature_marker_throughput")


def plot_downstream_times(frame, output_dir: Path, colors: dict[str, str]) -> None:
    from plotnine import aes, facet_wrap, geom_line, geom_point, ggplot, labs, scale_color_manual, scale_y_log10

    downstream = frame[frame["command"].isin(["freq", "distrib", "signif"])].copy()
    if downstream.empty:
        return
    downstream["min_depth_numeric"] = downstream["min_depth"].astype(float).astype(int)
    palette = {
        label: color
        for label, color in zip(
            downstream["dataset_label"].cat.categories,
            [colors["teal"], colors["sky"], colors["coral"], colors["magenta"], colors["sunshine"]],
            strict=False,
        )
    }
    plot = (
        ggplot(downstream, aes(x="min_depth_numeric", y="elapsed_seconds", color="dataset_label", group="dataset_label"))
        + geom_line(linewidth=0.7)
        + geom_point(size=1.9)
        + facet_wrap("~ command", scales="free_y")
        + scale_color_manual(values=palette)
        + scale_y_log10()
        + labs(x="Minimum depth", y="Elapsed seconds (log10)", color="")
        + base_theme()
    )
    save_plot(plot, output_dir, "literature_downstream_times")


def plot_bayesian_evidence(path: Path, output_dir: Path, colors: dict[str, str]) -> None:
    if not path.exists():
        return
    import pandas as pd
    from plotnine import aes, facet_wrap, geom_col, ggplot, labs, scale_fill_manual

    evidence = pd.read_csv(path)
    if evidence.empty:
        return
    evidence["dataset_label"] = evidence["dataset"].str.replace("_", " ").str.title()
    long = evidence.melt(
        id_vars=["dataset_label", "min_depth"],
        value_vars=["markers_bf_gt_10", "markers_posterior_gt_0_9"],
        var_name="metric",
        value_name="markers",
    )
    long["metric"] = long["metric"].map(
        {
            "markers_bf_gt_10": "Bayes factor > 10",
            "markers_posterior_gt_0_9": "Posterior > 0.9",
        }
    )
    long["markers"] = long["markers"] / 1_000_000.0
    plot = (
        ggplot(long, aes(x="factor(min_depth)", y="markers", fill="metric"))
        + geom_col(position="dodge", width=0.72)
        + facet_wrap("~ dataset_label", scales="free_y")
        + scale_fill_manual(values={"Bayes factor > 10": colors["sky"], "Posterior > 0.9": colors["magenta"]})
        + labs(x="Minimum depth", y="Markers (millions)", fill="")
        + base_theme()
    )
    save_plot(plot, output_dir, "literature_bayesian_evidence", width=7.2, height=5.2)


def plot_speed_comparison(path: Path, output_dir: Path, colors: dict[str, str]) -> None:
    if not path.exists():
        return
    import pandas as pd
    from plotnine import aes, coord_flip, geom_col, ggplot, labs, scale_fill_manual

    comparison = pd.read_csv(path)
    if comparison.empty:
        return
    pivot = comparison.pivot_table(
        index=["dataset", "command", "min_depth"],
        columns="impl",
        values="elapsed_seconds",
        aggfunc="sum",
    ).reset_index()
    if "cpp" not in pivot or "rust" not in pivot:
        return
    pivot = pivot[(pivot["rust"] > 0) & (pivot["cpp"] > 0)].copy()
    pivot["speedup"] = pivot["cpp"] / pivot["rust"]
    pivot["dataset_label"] = pivot["dataset"].str.replace("_", " ").str.title()
    pivot["command_depth"] = pivot["command"] + " d" + pivot["min_depth"].fillna("").astype(str).str.replace(".0", "", regex=False)
    plot = (
        ggplot(pivot, aes(x="command_depth", y="speedup", fill="dataset_label"))
        + geom_col(width=0.72)
        + coord_flip()
        + scale_fill_manual(values=[colors["teal"], colors["sky"], colors["coral"], colors["magenta"], colors["sunshine"]])
        + labs(x="", y="C++ RADSex / rsx elapsed time", fill="")
        + base_theme()
    )
    save_plot(plot, output_dir, "literature_radsex_speedups", width=7.2, height=5.2)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default=Path("benchmarks/results/literature_benchmark_results.csv"), type=Path)
    parser.add_argument("--bayesian", default=Path("benchmarks/results/literature_bayesian_evidence.csv"), type=Path)
    parser.add_argument("--comparison", default=Path("benchmarks/results/literature_speed_comparison.csv"), type=Path)
    parser.add_argument("--output", default=Path("docs/figures"), type=Path)
    args = parser.parse_args()

    rows = read_csv_rows(args.input)
    write_summary(args.output / "literature_benchmark_summary.tsv", summarize_dataset_rows(rows))

    colors = load_ruhi_colors()
    frame, phase_runtime, summary = prepare_results_frame(args.input)
    frame, phase_runtime, summary = order_by_total(frame, phase_runtime, summary)
    plot_runtime_breakdown(phase_runtime, args.output, colors)
    plot_marker_throughput(frame, args.output, colors)
    plot_downstream_times(frame, args.output, colors)
    plot_bayesian_evidence(args.bayesian, args.output, colors)
    plot_speed_comparison(args.comparison, args.output, colors)
    print(f"Wrote figures to {args.output}")


if __name__ == "__main__":
    main()
