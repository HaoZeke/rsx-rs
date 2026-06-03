"""High-level result objects for pyrsx (Option B design).

All results hold their data as a **narwhals DataFrame** (see `_df` / `.df`).
This makes them backend-agnostic by design:

- The bindings (and internal helpers like `_read_core_tsv`) produce or
  consume pyarrow Tables where possible for efficiency.
- Everything is then wrapped via `to_narwhals(...)` (or comes from
  `*_from_arrow` / Arrow paths that already return pyarrow).
- Users get the full narwhals API plus easy conversions:
  `.to_pandas()`, `.to_polars()`, `to_dataframe(backend=...)`.
- Pandas (or polars, etc.) is **never** a required internal dependency
  for reading rsx core outputs or running analysis — it is only a
  user-requested output backend (the "standard narwhals approach").

Designed to be excellent Python citizens:
- Frozen dataclasses where it makes sense (immutability + hashing)
- Properties instead of public _underscored fields
- Delegation to the underlying DataFrame so siuba / plotnine users have a great experience
- Rich __repr__ and summary methods
- `__arrow_c_stream__` support on most results for modern Arrow consumers.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Literal

import narwhals as nw

from pyrsx._adapters import from_narwhals, to_narwhals


@dataclass(frozen=True, kw_only=True)
class TriageResult:
    """
    Rich, immutable result object returned by `MarkerTable.triage(...)`.

    This is a first-class Python object, not just a bag of data.
    The `_df` is **always** a narwhals DataFrame (backend agnostic):
    - For in-memory `MarkerTable` (Arrow path): comes from `*_to_arrow_from_arrow`
      → pyarrow Table → `to_narwhals(...)`.
    - For path-backed `MarkerTable`: the low-level CLI binding writes a TSV;
      we read it with the pure-pyarrow `_read_core_tsv` helper → `to_narwhals(...)`.

    This follows the standard narwhals approach everywhere. The concrete
    backend is pyarrow by default (lightweight, no forced pandas/polars for
    internal I/O of rsx artifacts), but you can convert on demand.

    It supports:
    - Direct attribute delegation to the underlying narwhals DataFrame
      (so you can do `result.posterior_sex_linked`, `result >> siuba...` etc.)
    - `.df` property returning a narwhals DataFrame (great for siuba/plotnine)
    - Explicit `.to_pandas()`, `.to_polars()`, `.to_dataframe(backend=...)`
    - Plotting methods using your preferred plotnine + ruhi theme
    - Full provenance via the `params` dataclass
    - `__arrow_c_stream__` for zero-copy consumption by Polars, DuckDB, etc.
    """

    _df: nw.DataFrame
    params: dict[str, Any]
    _input_backend: Literal["pandas", "polars", "pyarrow", "auto"] = "auto"

    def __post_init__(self) -> None:
        if not isinstance(self._df, nw.DataFrame):
            object.__setattr__(self, "_df", to_narwhals(self._df))

    # ------------------------------------------------------------------ #
    # Data access (the important part for siuba/plotnine users)
    # ------------------------------------------------------------------ #

    @property
    def df(self) -> nw.DataFrame:
        """Narwhals DataFrame — the best object to feed into siuba or plotnine."""
        return self._df

    def __getattr__(self, name: str) -> Any:
        """Delegate unknown attributes to the underlying narwhals DataFrame.

        This is the magic that makes the object feel like a DataFrame for
        siuba users: `result.posterior_sex_linked`, `result.query(...)`, etc.
        """
        if name.startswith("_"):
            raise AttributeError(name)
        try:
            return getattr(self._df, name)
        except AttributeError:
            raise AttributeError(
                f"{type(self).__name__!r} object has no attribute {name!r}"
            ) from None

    # ------------------------------------------------------------------ #
    # Export
    # ------------------------------------------------------------------ #

    def to_dataframe(self, *, backend: str | None = None) -> Any:
        b = backend or self._input_backend
        return from_narwhals(self._df, backend=b)  # type: ignore[arg-type]

    def to_pandas(self) -> Any:
        return self.to_dataframe(backend="pandas")

    def to_polars(self) -> Any:
        return self.to_dataframe(backend="polars")

    def __arrow_c_stream__(self, requested_schema=None):
        """Support the Arrow C Data Interface / stream protocol for zero-copy consumers.

        Allows Polars, DuckDB, PyArrow, etc. to consume the result without
        going through pandas or full materialization in many cases.
        """
        try:
            import pyarrow as pa
            table = self.to_dataframe(backend="pyarrow")
            if hasattr(table, "__arrow_c_stream__"):
                return table.__arrow_c_stream__(requested_schema)
            return pa.table(table).__arrow_c_stream__(requested_schema)
        except Exception as e:
            raise NotImplementedError(
                f"Arrow C stream not supported for this result in current backend: {e}"
            ) from e

    # ------------------------------------------------------------------ #
    # Introspection
    # ------------------------------------------------------------------ #

    def summary(self) -> str:
        n = len(self._df)
        post_col = "posterior_sex_linked"
        post = (self._df[post_col] > 0.9).sum() if post_col in self._df.columns else 0
        return f"TriageResult(n_rows={n}, posterior>0.9≈{int(post)})"

    def __repr__(self) -> str:
        return f"<{self.summary()}>"

    # ------------------------------------------------------------------ #
    # Plotting (plotnine + ruhi by default — your preferred stack)
    # ------------------------------------------------------------------ #

    def plot_evidence(self, **kwargs: Any) -> Any:
        """Evidence class breakdown using plotnine + ruhi colors."""
        from pyrsx.plot import plot_evidence as _plot_evidence
        pdf = self.to_pandas()
        return _plot_evidence(pdf, **kwargs)


@dataclass(frozen=True, kw_only=True)
class PcaResult:
    """Result of `MarkerTable.pca(...)`."""

    _df: nw.DataFrame
    params: dict[str, Any] = field(default_factory=dict)
    _input_backend: Literal["pandas", "polars", "pyarrow", "auto"] = "auto"

    @property
    def df(self) -> nw.DataFrame:
        return self._df

    def __getattr__(self, name: str) -> Any:
        if name.startswith("_"):
            raise AttributeError(name)
        return getattr(self._df, name)

    def to_dataframe(self, *, backend: str | None = None) -> Any:
        b = backend or self._input_backend
        return from_narwhals(self._df, backend=b)

    def to_pandas(self) -> Any:
        return self.to_dataframe(backend="pandas")

    def to_polars(self) -> Any:
        return self.to_dataframe(backend="polars")

    def __arrow_c_stream__(self, requested_schema=None):
        """Support the Arrow C Data Interface / stream protocol for zero-copy consumers.

        Allows Polars, DuckDB, PyArrow, etc. to consume the result without
        going through pandas or full materialization in many cases.
        """
        try:
            import pyarrow as pa
            table = self.to_dataframe(backend="pyarrow")
            if hasattr(table, "__arrow_c_stream__"):
                return table.__arrow_c_stream__(requested_schema)
            return pa.table(table).__arrow_c_stream__(requested_schema)
        except Exception as e:
            raise NotImplementedError(
                f"Arrow C stream not supported for this result in current backend: {e}"
            ) from e

    def summary(self) -> str:
        return f"PcaResult(n_rows={len(self._df)})"

    def __repr__(self) -> str:
        return f"<{self.summary()}>"

    def plot(self, **kwargs: Any) -> Any:
        """Quick PCA scatter using plotnine + ruhi theme."""
        from pyrsx.plot import ruhi_theme
        from plotnine import ggplot, aes, geom_point, labs

        pdf = self.to_pandas()
        cols = [c for c in pdf.columns if str(c).upper().startswith("PC")]
        if len(cols) >= 2:
            x, y = cols[0], cols[1]
        else:
            x, y = pdf.columns[0], pdf.columns[1]

        p = (
            ggplot(pdf, aes(x=x, y=y))
            + geom_point(alpha=0.6, color="#004D40")
            + labs(title="PCA (streaming Gram)")
            + ruhi_theme()
        )
        return p


# ------------------------------------------------------------------ #
# ------------------------------------------------------------------ #
# Generic result for simple table-returning commands
# (freq, depth, distrib, signif, subset, etc.)
# ------------------------------------------------------------------ #

@dataclass(frozen=True, kw_only=True)
class TableResult:
    """
    Generic, clean result object for commands that mainly return a table
    (freq, depth, distrib, signif, etc.).

    The `_df` is a narwhals DataFrame produced via the standard backend-
    agnostic path: pyarrow Table (from `*_from_arrow` or `_read_core_tsv`)
    → `to_narwhals(...)`. See the module docstring and `MarkerTable` for
    the full story on why this design keeps pandas/polars as pure output
    backends.

    This class provides one representation for simple table-returning commands.
    """
    _df: nw.DataFrame
    command: str                              # e.g. "freq", "depth", "distrib", "signif"
    params: dict[str, Any] = field(default_factory=dict)

    @property
    def df(self) -> nw.DataFrame:
        return self._df

    def __getattr__(self, name: str) -> Any:
        if name.startswith("_"):
            raise AttributeError(name)
        return getattr(self._df, name)

    def to_dataframe(self, *, backend: str = "pandas") -> Any:
        return from_narwhals(self._df, backend=backend)

    def to_pandas(self) -> Any:
        return self.to_dataframe(backend="pandas")

    def to_polars(self) -> Any:
        return self.to_dataframe(backend="polars")

    def __arrow_c_stream__(self, requested_schema=None):
        """Support the Arrow C Data Interface / stream protocol for zero-copy consumers.

        Allows Polars, DuckDB, PyArrow, etc. to consume the result without
        going through pandas or full materialization in many cases.
        """
        try:
            import pyarrow as pa
            table = self.to_dataframe(backend="pyarrow")
            if hasattr(table, "__arrow_c_stream__"):
                return table.__arrow_c_stream__(requested_schema)
            return pa.table(table).__arrow_c_stream__(requested_schema)
        except Exception as e:
            raise NotImplementedError(
                f"Arrow C stream not supported for this result in current backend: {e}"
            ) from e

    def __repr__(self) -> str:
        return f"<TableResult command={self.command!r} shape={len(self._df)}x{len(self._df.columns)}>"
