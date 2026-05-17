"""pyrsx: Python bindings for rsx (RAD-seq sex determination toolkit).

This package exposes both the low-level (fast, thin) Rust bindings and a
higher-level, idiomatic Python API designed for notebooks, workflows,
and data analysis (the recommended interface for most users).

Example
-------
>>> import pyrsx as rsx
>>> table = rsx.MarkerTable.from_path("markers.tsv")
>>> result = table.triage(popmap="popmap.tsv", min_depth=10)
>>> result.plot_evidence()
"""

# Low-level direct bindings (still available for power users / legacy code)
from .pyrsx import (  # noqa: F401
    depth,
    distrib,
    freq,
    merge,
    pca,
    process,
    signif,
    triage,
)

# High-level idiomatic API (recommended)
from .api.markers import MarkerTable  # noqa: F401
from .api.params import TriageParams  # noqa: F401
from .api.results import TriageResult  # noqa: F401

__all__ = [
    # Low-level (for compatibility)
    "process",
    "distrib",
    "signif",
    "triage",
    "freq",
    "depth",
    "merge",
    "pca",
    # High-level (new recommended API)
    "MarkerTable",
    "TriageResult",
    "TriageParams",
]
