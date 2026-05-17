"""High-level MarkerTable class for the pyrsx Python API (Option B)."""

from __future__ import annotations

import tempfile
from pathlib import Path
from typing import Any, Literal

import narwhals as nw

from pyrsx._adapters import from_narwhals, is_dataframe_like, to_narwhals
from pyrsx._rsx import process as _process_lowlevel  # type: ignore[attr-defined]
# We will gradually add more low-level imports as we implement methods


class MarkerTable:
    """
    High-level representation of a RAD-seq marker depth table.

    This is the central object in the pyrsx high-level API. It can be
    constructed from files or any narwhals-compatible DataFrame and
    provides a fluent, Pythonic interface to the rsx analysis commands.

    Examples
    --------
    >>> table = MarkerTable.from_path("markers.tsv")
    >>> triage = table.triage(popmap="popmap.tsv", min_depth=10)
    >>> df = triage.to_pandas()
    """

    def __init__(
        self,
        data: Any | None = None,
        *,
        path: str | Path | None = None,
        backend: Literal["pandas", "polars", "pyarrow", "auto"] = "auto",
    ) -> None:
        """
        Internal constructor. Prefer the classmethods `from_path` or
        `from_dataframe`.
        """
        if (data is None) == (path is None):
            raise ValueError("Exactly one of `data` or `path` must be provided.")

        self._path: Path | None = Path(path) if path is not None else None
        self._df: nw.DataFrame | None = None
        self._backend = backend

        if data is not None:
            if not is_dataframe_like(data):
                raise TypeError(
                    f"Expected a DataFrame-like object, got {type(data)}"
                )
            self._df = to_narwhals(data)

    # ------------------------------------------------------------------ #
    # Constructors
    # ------------------------------------------------------------------ #

    @classmethod
    def from_path(
        cls,
        path: str | Path,
        *,
        backend: Literal["pandas", "polars", "pyarrow", "auto"] = "auto",
    ) -> MarkerTable:
        """Load a marker table from a TSV/CSV/Parquet file on disk."""
        return cls(path=path, backend=backend)

    @classmethod
    def from_dataframe(
        cls,
        df: Any,
        *,
        backend: Literal["pandas", "polars", "pyarrow", "auto"] = "auto",
    ) -> MarkerTable:
        """Create a MarkerTable directly from a DataFrame."""
        return cls(data=df, backend=backend)

    # ------------------------------------------------------------------ #
    # Properties & introspection
    # ------------------------------------------------------------------ #

    @property
    def n_markers(self) -> int:
        if self._df is not None:
            return len(self._df)
        # Fall back to counting lines if we only have a path (cheap for now)
        # In a real implementation we would cache a header count.
        with open(self._path) as f:  # type: ignore[arg-type]
            return sum(1 for _ in f) - 1  # subtract header

    @property
    def n_individuals(self) -> int:
        if self._df is not None:
            # Assume first two columns are id + sequence, rest are samples
            return len(self._df.columns) - 2
        # Parse header from file
        with open(self._path) as f:  # type: ignore[arg-type]
            header = next(f).strip().split("\t")
            return len(header) - 2

    def __repr__(self) -> str:
        return (
            f"MarkerTable(n_markers={self.n_markers}, "
            f"n_individuals={self.n_individuals})"
        )

    def summary(self) -> str:
        """Return a human-readable one-line summary."""
        return (
            f"MarkerTable with {self.n_markers} markers across "
            f"{self.n_individuals} individuals"
        )

    # ------------------------------------------------------------------ #
    # Core analysis methods (stubs that will grow)
    # ------------------------------------------------------------------ #

    def triage(
        self,
        *,
        popmap: Any,
        min_depth: int = 10,
        posterior_threshold: float = 0.9,
        bayes_factor_threshold: float = 10.0,
        prior: float = 0.01,
        linked_prob: float = 0.9,
        group1: str = "M",
        group2: str = "F",
        **kwargs: Any,
    ) -> "TriageResult":
        """
        Run the hybrid strict + Bayesian triage (the core biological output of rsx).

        This is the first fully working method of the high-level "right" API.
        It accepts DataFrames for both the marker table and popmap, calls the
        fast Rust implementation (via temporary files for the first version),
        and returns a rich `TriageResult` object.
        """
        from .results import TriageResult
        import tempfile
        from pathlib import Path
        import pandas as pd

        # --- Resolve inputs to paths (first version uses temp files) ---
        if self._df is not None:
            markers_path = Path(tempfile.NamedTemporaryFile(suffix=".tsv", delete=False).name)
            pd_df = from_narwhals(self._df, backend="pandas")
            pd_df.to_csv(markers_path, sep="\t", index=False)
        else:
            markers_path = self._path  # type: ignore[assignment]

        if is_dataframe_like(popmap):
            popmap_path = Path(tempfile.NamedTemporaryFile(suffix=".tsv", delete=False).name)
            pd_pop = from_narwhals(to_narwhals(popmap), backend="pandas")
            pd_pop.to_csv(popmap_path, sep="\t", index=False)
        else:
            popmap_path = Path(popmap)

        output_path = Path(tempfile.NamedTemporaryFile(suffix="_triage.tsv", delete=False).name)

        # Call the existing fast low-level binding (exposed at package top level)
        import pyrsx as _pyrsx
        _lowlevel_triage = _pyrsx.triage

        _lowlevel_triage(
            str(markers_path),
            str(popmap_path),
            str(output_path),
            min_depth=min_depth,
            posterior_threshold=posterior_threshold,
            bayes_factor_threshold=bayes_factor_threshold,
            prior=prior,
            linked_prob=linked_prob,
            group1=group1,
            group2=group2,
        )

        # Read result back and wrap
        result_df = pd.read_csv(output_path, sep="\t", comment="#")
        nw_df = to_narwhals(result_df)

        params = {
            "min_depth": min_depth,
            "posterior_threshold": posterior_threshold,
            "bayes_factor_threshold": bayes_factor_threshold,
            "prior": prior,
            "linked_prob": linked_prob,
            "group1": group1,
            "group2": group2,
        }

        return TriageResult(
            _df=nw_df,
            params=params,
            _input_backend=self._backend,
        )

    # Future methods: .pca(), .depth_stats(), .distrib(), etc.

    # ------------------------------------------------------------------ #
    # Export / conversion
    # ------------------------------------------------------------------ #

    def to_dataframe(
        self,
        *,
        backend: Literal["pandas", "polars", "pyarrow", "auto"] = "auto",
    ) -> Any:
        """Return the underlying data as a native DataFrame."""
        if self._df is not None:
            return from_narwhals(self._df, backend=backend)
        # If we only have a path, the caller probably wants to keep using the path
        raise RuntimeError(
            "MarkerTable was constructed from a path and has no in-memory "
            "DataFrame. Use `.to_path()` or load it with `from_dataframe`."
        )
