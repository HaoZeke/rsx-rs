"""High-level result objects for pyrsx (Option B design)."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Literal

import narwhals as nw

from pyrsx._adapters import from_narwhals, to_narwhals


@dataclass
class TriageResult:
    """
    Rich result object returned by `MarkerTable.triage(...)`.

    Contains the triage output as a DataFrame (in the flavor the user
    prefers) plus the exact parameters that were used (for full
    reproducibility).
    """

    _df: nw.DataFrame
    params: dict[str, Any] = field(default_factory=dict)
    _input_backend: Literal["pandas", "polars", "pyarrow", "auto"] = "auto"

    def __post_init__(self) -> None:
        if not isinstance(self._df, nw.DataFrame):
            self._df = to_narwhals(self._df)

    @property
    def df(self) -> nw.DataFrame:
        """The raw narwhals DataFrame (for advanced users)."""
        return self._df

    def to_dataframe(self, *, backend: str | None = None) -> Any:
        """Return the triage results as a native DataFrame."""
        b = backend or self._input_backend
        return from_narwhals(self._df, backend=b)  # type: ignore[arg-type]

    def to_pandas(self) -> Any:
        return self.to_dataframe(backend="pandas")

    def to_polars(self) -> Any:
        return self.to_dataframe(backend="polars")

    def summary(self) -> str:
        """Human-readable one-line summary."""
        n = len(self._df)
        post = (self._df["posterior_sex_linked"] > 0.9).sum() if "posterior_sex_linked" in self._df.columns else 0
        return f"TriageResult(n_rows={n}, posterior>0.9={post})"

    def __repr__(self) -> str:
        return f"<{self.summary()}>"

    def plot_evidence(self, **kwargs: Any) -> Any:
        """Plot the evidence class breakdown using plotnine + the ruhi colorscheme."""
        from pyrsx.plot import plot_evidence as _plot_evidence

        # Convert to pandas for plotnine (plotnine works great with pandas)
        pdf = self.to_pandas()
        return _plot_evidence(pdf, **kwargs)
