"""High-level result objects for pyrsx (Option B design).

Designed to be excellent Python citizens:
- Frozen dataclasses where it makes sense (immutability + hashing)
- Properties instead of public _underscored fields
- Delegation to the underlying DataFrame so siuba / plotnine users have a great experience
- Rich __repr__ and summary methods
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
    It supports:
    - Direct attribute delegation to the underlying narwhals DataFrame
      (so you can do `result.posterior_sex_linked`, `result >> siuba...` etc.)
    - `.df` property returning a narwhals DataFrame (great for siuba/plotnine)
    - Explicit `.to_pandas()`, `.to_polars()`, `.to_dataframe(backend=...)`
    - Plotting methods using your preferred plotnine + ruhi theme
    - Full provenance via the `params` dataclass
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
# Lightweight result wrappers for other commands
# ------------------------------------------------------------------ #

@dataclass(frozen=True, kw_only=True)
class FreqResult:
    """Result of `MarkerTable.freq(...)`."""
    _df: nw.DataFrame
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
        return self.to_dataframe("pandas")


@dataclass(frozen=True, kw_only=True)
class DepthResult:
    """Result of `MarkerTable.depth(...)`."""
    _df: nw.DataFrame
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
        return self.to_dataframe("pandas")


@dataclass(frozen=True, kw_only=True)
class DistribResult:
    """Result of `MarkerTable.distrib(...)`."""
    _df: nw.DataFrame
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
        return self.to_dataframe("pandas")


@dataclass(frozen=True, kw_only=True)
class SignifResult:
    """Result of `MarkerTable.signif(...)`."""
    _df: nw.DataFrame
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
        return self.to_dataframe("pandas")
