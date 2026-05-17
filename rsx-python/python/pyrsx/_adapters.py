"""DataFrame adapters for pyrsx high-level API.

Uses narwhals to provide a lightweight, zero-dependency (beyond narwhals)
way to accept and return pandas, polars, pyarrow, and any other
narwhals-compatible DataFrame without forcing a hard dependency on any
particular library.

This follows modern Python scientific ecosystem best practices for
"write once, run anywhere" DataFrame code.
"""

from __future__ import annotations

from typing import Any, Literal, Protocol, TypeVar

import narwhals as nw
try:
    # narwhals < 2.x
    from narwhals.typing import IntoDataFrame, DataFrame as NWDataFrame  # type: ignore[attr-defined]
except ImportError:
    # narwhals >= 2.x renamed DataFrame -> DataFrameT
    from narwhals.typing import IntoDataFrame, DataFrameT as NWDataFrame  # type: ignore[attr-defined]

DF = TypeVar("DF", bound=IntoDataFrame)


class DataFrameAdapter(Protocol):
    """Protocol for objects that can be converted to/from narwhals DataFrames."""

    def __narwhals_dataframe__(self) -> Any: ...


def to_narwhals(df: IntoDataFrame | NWDataFrame) -> NWDataFrame:
    """Convert any supported DataFrame-like to a narwhals DataFrame."""
    return nw.from_native(df, strict=False)  # type: ignore[return-value]


def from_narwhals(
    df: NWDataFrame,
    *,
    backend: Literal["pandas", "polars", "pyarrow", "auto"] = "auto",
) -> Any:
    """
    Convert a narwhals DataFrame back to the desired native backend.

    If backend="auto", try to return the same flavor as the original input
    when possible (best-effort). Falls back to pandas.
    """
    if backend == "auto":
        # Narwhals can sometimes preserve the original implementation
        try:
            return df.to_native()
        except Exception:
            backend = "pandas"

    if backend == "pandas":
        return df.to_pandas()
    if backend == "polars":
        return df.to_polars()
    if backend == "pyarrow":
        return df.to_arrow()

    # Fallback
    return df.to_pandas()


def is_dataframe_like(obj: Any) -> bool:
    """Check whether an object can be treated as a DataFrame by our adapters."""
    try:
        nw.from_native(obj, strict=False)
        return True
    except Exception:
        return False
