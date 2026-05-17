"""Plotting helpers for pyrsx high-level API.

Defaults to plotnine (ggplot2) + the user's preferred "ruhi" colorscheme.
"""

from __future__ import annotations

from typing import Any

# Ruhi colorscheme (from chemparseplot, user's personal palette)
RUHI_COLORS = {
    "coral": "#FF655D",
    "sunshine": "#F1DB4B",
    "teal": "#004D40",
    "sky": "#1E88E5",
    "magenta": "#D81B60",
}

# Order that looks good for categorical evidence classes
RUHI_PALETTE = [
    RUHI_COLORS["teal"],
    RUHI_COLORS["sky"],
    RUHI_COLORS["magenta"],
    RUHI_COLORS["coral"],
    RUHI_COLORS["sunshine"],
]


def ruhi_theme():
    """Return a plotnine theme using the ruhi colors and Atkinson's Hyperlegible font."""
    from plotnine import theme, element_text, element_line, element_rect

    return (
        theme(
            text=element_text(family="Atkinson Hyperlegible", size=11),
            plot_background=element_rect(fill="white"),
            panel_background=element_rect(fill="white"),
            panel_grid_major=element_line(color="#f0f0f0"),
            axis_text=element_text(color="black"),
            axis_title=element_text(color="black"),
            legend_background=element_rect(fill="white"),
        )
        + theme(figure_size=(8, 5))
    )


def plot_evidence(
    df: Any,
    *,
    x: str = "dataset",
    fill: str = "candidate_class",
    theme: Any | None = None,
    **kwargs: Any,
) -> Any:
    """
    Evidence class breakdown plot using plotnine + the ruhi colorscheme.

    This is the default plotting experience for TriageResult.
    """
    from plotnine import (
        ggplot,
        aes,
        geom_bar,
        labs,
        scale_fill_manual,
        position_dodge,
        theme as p9_theme,
    )

    if theme is None:
        theme = ruhi_theme()

    p = (
        ggplot(df, aes(x=x, fill=fill))
        + geom_bar(position=position_dodge, stat="count", **kwargs)
        + scale_fill_manual(values=RUHI_PALETTE)
        + labs(
            x="",
            y="Number of markers",
            fill="Evidence class",
        )
        + theme
    )
    return p
