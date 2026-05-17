#!/usr/bin/env python3
"""Render a per-dataset (prior x linked-probability) heatmap from the
detailed triage TSVs produced by the prior-sensitivity SLURM sweep.

The 96 input files live under benchmarks/results/slurm/ and follow the
pattern triage_<dataset>_pi<prior>_psex<linked_prob>.tsv. For each
combination we count how many markers cleared the posterior > 0.9
threshold (the standard Bayesian-call line used elsewhere in the paper).

Outputs:
- docs/figures/literature_prior_sensitivity_heatmap.{svg,pdf}
- benchmarks/results/prior_sensitivity_from_triage.csv
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd
from plotnine import (
    aes,
    element_blank,
    element_rect,
    element_text,
    facet_wrap,
    geom_text,
    geom_tile,
    ggplot,
    labs,
    scale_fill_gradient,
    theme,
    theme_minimal,
)

FILE_RE = re.compile(
    r"^triage_(?P<dataset>.+?)_pi(?P<pi>[0-9.]+)_psex(?P<psex>[0-9.]+)\.tsv$"
)


def summarise_one(path: Path) -> dict[str, float | int | str]:
    df = pd.read_csv(path, sep="\t", comment="#")
    n_posterior_gt_0_9 = int((df["Posterior_SexLinked"] > 0.9).sum())
    n_bf_gt_10 = int((df["Bayes_Factor"] > 10).sum())
    n_strict = int(df["Strict_Call"].astype(str).str.lower().eq("true").sum())
    n_total = int(len(df))
    return {
        "n_posterior_gt_0_9": n_posterior_gt_0_9,
        "n_bayes_factor_gt_10": n_bf_gt_10,
        "n_strict_call": n_strict,
        "n_candidate_rows": n_total,
    }


def collect(slurm_dir: Path) -> pd.DataFrame:
    rows = []
    for tsv in sorted(slurm_dir.glob("triage_*_pi*_psex*.tsv")):
        m = FILE_RE.match(tsv.name)
        if not m:
            continue
        meta = m.groupdict()
        summary = summarise_one(tsv)
        rows.append(
            {
                "dataset": meta["dataset"],
                "prior": float(meta["pi"]),
                "linked_prob": float(meta["psex"]),
                **summary,
            }
        )
    return pd.DataFrame(rows).sort_values(["dataset", "prior", "linked_prob"])


def plot_heatmap(summary: pd.DataFrame, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    # Treat prior and linked_prob as ordered categoricals so the heatmap
    # cells sit on a regular grid (plotnine cannot tile floats reliably).
    prior_levels = sorted(summary["prior"].unique())
    psex_levels = sorted(summary["linked_prob"].unique())
    plot_df = summary.copy()
    plot_df["prior"] = pd.Categorical(
        plot_df["prior"].astype(str), categories=[str(v) for v in prior_levels], ordered=True
    )
    plot_df["linked_prob"] = pd.Categorical(
        plot_df["linked_prob"].astype(str), categories=[str(v) for v in psex_levels], ordered=True
    )

    p = (
        ggplot(plot_df, aes(x="prior", y="linked_prob", fill="n_posterior_gt_0_9"))
        + geom_tile(color="white", size=0.4)
        + geom_text(aes(label="n_posterior_gt_0_9"), size=7, color="#212121")
        + facet_wrap("~dataset", ncol=2)
        + scale_fill_gradient(low="#FFFFFF", high="#004D40", name="markers w/ posterior > 0.9")
        + labs(
            x="prior probability sex-linked (pi)",
            y="linked-sex prevalence (p_sex)",
            title="Prior sensitivity of the Bayesian triage layer",
            subtitle="Counts of qualifying markers per dataset across the 6x4 (pi, p_sex) grid",
        )
        + theme_minimal()
        + theme(
            figure_size=(8.4, 7.2),
            panel_grid=element_blank(),
            panel_background=element_rect(fill="white"),
            strip_text=element_text(face="bold"),
        )
    )

    svg_path = output_dir / "literature_prior_sensitivity_heatmap.svg"
    pdf_path = output_dir / "literature_prior_sensitivity_heatmap.pdf"
    p.save(str(svg_path), verbose=False)
    p.save(str(pdf_path), verbose=False)
    print(f"Wrote {svg_path}")
    print(f"Wrote {pdf_path}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--slurm-dir",
        type=Path,
        default=Path("benchmarks/results/slurm"),
        help="Directory holding triage_<dataset>_pi<pi>_psex<psex>.tsv files",
    )
    parser.add_argument(
        "--summary-csv",
        type=Path,
        default=Path("benchmarks/results/prior_sensitivity_from_triage.csv"),
        help="Where to write the (dataset, prior, p_sex) -> counts summary",
    )
    parser.add_argument(
        "--figure-dir",
        type=Path,
        default=Path("docs/figures"),
        help="Directory for the heatmap SVG/PDF",
    )
    args = parser.parse_args()

    if not args.slurm_dir.is_dir():
        print(f"slurm directory does not exist: {args.slurm_dir}", file=sys.stderr)
        return 2

    summary = collect(args.slurm_dir)
    if summary.empty:
        print(
            f"no triage TSVs found under {args.slurm_dir}; nothing to do",
            file=sys.stderr,
        )
        return 2

    args.summary_csv.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.summary_csv, index=False)
    print(f"Wrote {args.summary_csv} ({len(summary)} rows)")

    plot_heatmap(summary, args.figure_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
