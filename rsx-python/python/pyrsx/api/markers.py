"""High-level MarkerTable class for the pyrsx Python API (Option B)."""

from __future__ import annotations

import tempfile
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal

import narwhals as nw

from pyrsx._adapters import from_narwhals, is_dataframe_like, to_narwhals

if TYPE_CHECKING:
    from .params import TriageParams
    from .results import TriageResult
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
        if (data is None) == (path is None):
            raise ValueError("Exactly one of `data` or `path` must be provided.")

        self._path: Path | None = Path(path) if path is not None else None
        self._df: nw.DataFrame | None = None
        self._backend: Literal[...] = backend

        if data is not None:
            if not is_dataframe_like(data):
                raise TypeError(f"Expected DataFrame-like, got {type(data)}")
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
        with open(self._path) as f:  # type: ignore[arg-type]
            return sum(1 for _ in f) - 1

    @property
    def n_individuals(self) -> int:
        if self._df is not None:
            return len(self._df.columns) - 2
        with open(self._path) as f:  # type: ignore[arg-type]
            header = next(f).strip().split("\t")
            return len(header) - 2

    @property
    def data(self) -> nw.DataFrame | None:
        """The underlying data as a narwhals DataFrame (preferred over .df for consistency)."""
        return self._df

    @property
    def df(self) -> nw.DataFrame | None:
        """Alias for .data — convenient for siuba/plotnine users."""
        return self._df

    def __getattr__(self, name: str) -> Any:
        if name.startswith("_"):
            raise AttributeError(name)
        if self._df is not None:
            return getattr(self._df, name)
        raise AttributeError(
            f"MarkerTable (path-backed) has no attribute '{name}'. "
            "Load it with from_dataframe() or access via .df after materializing."
        )

    def __len__(self) -> int:
        return self.n_markers

    def __repr__(self) -> str:
        return f"MarkerTable(n_markers={self.n_markers}, n_individuals={self.n_individuals})"

    def summary(self) -> str:
        return f"MarkerTable with {self.n_markers} markers across {self.n_individuals} individuals"

    # ------------------------------------------------------------------ #
    # Core analysis methods (stubs that will grow)
    # ------------------------------------------------------------------ #

    def triage(
        self,
        *,
        popmap: Any,
        params: TriageParams | None = None,
        **kwargs: Any,
    ) -> TriageResult:
        """
        Run the hybrid strict + Bayesian triage.

        Supports both the ergonomic kwarg style and passing a TriageParams
        dataclass (more idiomatic when you want to version/share configs).
        """
        from .results import TriageResult
        from .params import TriageParams
        import tempfile
        from pathlib import Path
        import pandas as pd

        # Merge dataclass + kwargs (kwargs override)
        if params is None:
            p = TriageParams(**{k: v for k, v in kwargs.items() if hasattr(TriageParams, k)})
        else:
            p = params

        # Resolve to on-disk paths for the Rust streaming readers.
        # Output side is now pure in-memory Arrow IPC (no temp files on results).
        # Input serialization is still required when passing DataFrames because
        # the core readers (`MarkersTableStream` etc.) are path-based today.
        if self._df is not None:
            mpath = Path(tempfile.NamedTemporaryFile(suffix=".parquet", delete=False).name)
            from_narwhals(self._df, backend="pandas").to_parquet(mpath, index=False)
        else:
            mpath = self._path  # type: ignore[assignment]

        if is_dataframe_like(popmap):
            ppath = Path(tempfile.NamedTemporaryFile(suffix=".parquet", delete=False).name)
            from_narwhals(to_narwhals(popmap), backend="pandas").to_parquet(ppath, index=False)
        else:
            ppath = Path(popmap)

        outpath = Path(tempfile.NamedTemporaryFile(suffix="_triage.parquet", delete=False).name)

        import pyrsx as _pyrsx

        # Use the real Rust Arrow emitter for the result (pure in-memory IPC).
        # The input table is still written to a temp file for the Rust reader.
        if self._df is not None:
            # Preferred path: Rust produces the data directly as Arrow
            arrow_table = _pyrsx.triage_to_arrow(
                str(mpath),
                str(ppath),
                min_depth=p.min_depth,
                posterior_threshold=p.posterior_threshold,
                prior_probability=p.prior,
                linked_probability=p.linked_prob,
                group1=p.group1,
                group2=p.group2,
            )
            res_df = to_narwhals(arrow_table)   # thin adapter, no heavy Python work
        else:
            # Legacy path-based usage (still supported)
            _lowlevel_triage = _pyrsx.triage
            _lowlevel_triage(
                str(mpath),
                str(ppath),
                str(outpath),
                min_depth=p.min_depth,
                posterior_threshold=p.posterior_threshold,
                bayes_factor_threshold=p.bayes_factor_threshold,
                prior=p.prior,
                linked_prob=p.linked_prob,
                group1=p.group1,
                group2=p.group2,
            )
            res_df = to_narwhals(pd.read_parquet(outpath))  # still thin
        return TriageResult(
            _df=to_narwhals(res_df),
            params={
                "min_depth": p.min_depth,
                "posterior_threshold": p.posterior_threshold,
                "prior": p.prior,
                "linked_prob": p.linked_prob,
            },
            _input_backend=self._backend,
        )

    def pca(self, *, k: int = 2, **kwargs: Any) -> "PcaResult":
        """Compute streaming sample PCA (O(n_individuals²) memory)."""
        from .results import PcaResult
        import tempfile
        from pathlib import Path
        import pandas as pd

        import pyrsx as _pyrsx

        if self._df is not None:
            # Write input to temp for the Rust reader, then use the real in-memory
            # Arrow emitter (one-shot IPC for eigenvalues + loadings).
            mpath = Path(tempfile.NamedTemporaryFile(suffix=".tsv", delete=False).name)
            from_narwhals(self._df, backend="pandas").to_csv(mpath, sep="\t", index=False)

            arrow_res = _pyrsx.pca_to_arrow(
                str(mpath),
                min_depth=kwargs.get("min_depth", 1),
                n_components=k,
            )
            # loadings is the main table users plot (individuals in PC space)
            res_df = to_narwhals(arrow_res["loadings"])
            # We can also stash eigenvalues / provenance on the result later
            return PcaResult(
                _df=res_df,
                params={"k": k, "arrow": True, **kwargs},
                _input_backend=self._backend,
            )
        else:
            mpath = self._path  # type: ignore[assignment]

        outpath = Path(tempfile.NamedTemporaryFile(suffix="_pca.tsv", delete=False).name)

        _lowlevel_pca = getattr(_pyrsx, "pca", None)
        if _lowlevel_pca is None:
            raise NotImplementedError("Low-level pca not exposed yet in this build")

        _lowlevel_pca(str(mpath), str(outpath), k=k, **kwargs)

        res_df = pd.read_csv(outpath, sep="\t", comment="#")
        return PcaResult(
            _df=to_narwhals(res_df),
            params={"k": k, **kwargs},
            _input_backend=self._backend,
        )

    # Future methods: .depth_stats(), .distrib(), etc.

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
