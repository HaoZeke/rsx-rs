"""High-level MarkerTable class for the pyrsx Python API (Option B)."""

from __future__ import annotations

import io
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal

import narwhals as nw

from pyrsx._adapters import from_narwhals, is_dataframe_like, to_narwhals

if TYPE_CHECKING:
    from .params import TriageParams
    from .results import (
        DepthResult,
        DistribResult,
        FreqResult,
        PcaResult,
        SignifResult,
        TriageResult,
    )


def _table_to_ipc_bytes(table: Any) -> bytes:
    """Serialise a pyarrow.Table to Arrow IPC stream bytes (RAM only)."""
    import pyarrow as pa

    buf = io.BytesIO()
    with pa.ipc.new_stream(buf, table.schema) as writer:
        writer.write_table(table)
    return buf.getvalue()


def _core_tsv_skip_rows(path: str | Path) -> int:
    with open(path, encoding="utf-8") as f:
        first = f.readline()
    return 1 if first.startswith("#") else 0


def _read_core_tsv(path: str | Path) -> "nw.DataFrame":
    """Read a TSV produced by rsx core into a **narwhals DataFrame**.

    Implementation:
    - Uses `pyarrow.csv` and skips the optional leading metadata comment line
      when it is present, leaving the real TSV header for pyarrow.
    - This is a lightweight parse that does **not** require pandas or polars
      just for internal I/O of command outputs.

    The pyarrow Table is immediately passed to `to_narwhals(...)`. The returned
    object **is backend agnostic** — it is a narwhals DataFrame whose concrete
    backend is pyarrow by default. Callers (and end users) can convert on
    demand:

        result.to_pandas()
        result.to_polars()
        result.to_dataframe(backend="polars")
        # or pass the .df to siuba, plotnine, polars, etc.

    This is the **standard narwhals approach** used everywhere in the high-level
    API (see `_adapters.py`, `from_dataframe`, Arrow paths, and Result objects).

    Why this matters:
    - No "pandas fallback" (or any other backend) is ever an internal
      implementation detail for reading rsx artifacts.
    - Users stay in their preferred DataFrame library without the bindings
      pulling in extra deps.
    - Matches the design of the Arrow-based in-memory paths (Rust produces
      pyarrow Tables or batches → Python does `to_narwhals(...)`).

    See also:
    - The matching Rust helper `read_tsv_to_pyarrow_table` (in src/lib.rs)
      used by the `*_from_arrow` low-level functions.
    - `MarkerTable` path-backed methods and the various `TableResult` /
      `*Result` classes, which all end up with narwhals-backed data.
    """
    import pyarrow.csv as pa_csv

    read_options = pa_csv.ReadOptions(
        skip_rows=_core_tsv_skip_rows(path),
        use_threads=True,
    )
    parse_options = pa_csv.ParseOptions(delimiter="\t")
    table = pa_csv.read_csv(
        str(path),
        read_options=read_options,
        parse_options=parse_options,
    )
    return to_narwhals(table)


def _parse_marker_count(line: str) -> int | None:
    """Parse rsx marker-count metadata from a leading comment line."""
    if not line.startswith("#Number of markers"):
        return None
    try:
        return int(line.split(":", 1)[1].strip())
    except (IndexError, ValueError):
        return None


def _arrow_bytes_from(obj: Any) -> bytes:
    """Coerce a DataFrame-like or file path into Arrow IPC bytes.

    A DataFrame (any narwhals-compatible) goes through pyarrow → IPC.
    A str/Path is read as a TSV (`individual\tgroup` for popmaps,
    `id\tsequence\t...` for markers) and converted via pyarrow.csv.
    """
    import pyarrow as pa
    import pyarrow.csv as pa_csv

    if is_dataframe_like(obj):
        table = from_narwhals(to_narwhals(obj), backend="pyarrow")
        return _table_to_ipc_bytes(table)

    if isinstance(obj, (str, Path)):
        read_options = pa_csv.ReadOptions(use_threads=False)
        parse_options = pa_csv.ParseOptions(delimiter="\t")
        table = pa_csv.read_csv(
            str(obj),
            read_options=read_options,
            parse_options=parse_options,
        )
        return _table_to_ipc_bytes(table)

    raise TypeError(
        f"Expected a DataFrame-like or path, got {type(obj).__name__}"
    )


class MarkerTable:
    """
    High-level representation of a RAD-seq marker depth table.

    This is the central object in the pyrsx high-level API. It can be
    constructed from files (path-backed, for very large data) or any
    narwhals-compatible DataFrame (in-memory) and provides a fluent,
    Pythonic interface to the rsx analysis commands.

    Backend agnosticism (via narwhals):
    - When working with DataFrames (in-memory or results), everything goes
      through narwhals. The concrete backend (pandas / polars / pyarrow / ...)
      is chosen by the user or inferred ("auto").
    - Internal reading of rsx core TSV outputs (for path-backed results or
      the `*_from_arrow` helpers) is done with pyarrow + `to_narwhals(...)`.
      The exposed objects are always narwhals DataFrames — **backend
      agnostic**. See `_read_core_tsv` and the module-level docs in the
      Rust bindings for details. Pandas (or any other backend) is only a
      user-requested output conversion, never forced internally.

    All result objects (`TriageResult`, `TableResult`, etc.) delegate to
    their underlying narwhals DataFrame, so you get a great experience with
    siuba, plotnine, Polars expressions, etc.

    Examples
    --------
    >>> table = MarkerTable.from_path("markers.tsv")
    >>> triage = table.triage(popmap="popmap.tsv", min_depth=10)
    >>> df = triage.to_pandas()          # explicit backend
    >>> pl_df = triage.to_polars()
    >>> # or stay backend-agnostic
    >>> nw_df = triage.df
    >>> # works with plotnine / siuba on whatever the user has
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
        self._path_header: tuple[int, list[str]] | None = None

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

    def _read_path_header(self) -> tuple[int, list[str]]:
        if self._path_header is not None:
            return self._path_header
        if self._path is None:
            raise RuntimeError("MarkerTable has no path-backed source.")

        with open(self._path, encoding="utf-8") as f:
            first = f.readline()
            if not first:
                self._path_header = (0, [])
                return self._path_header

            marker_count = _parse_marker_count(first)
            if first.startswith("#"):
                header_line = f.readline()
            else:
                header_line = first

            columns = header_line.rstrip("\r\n").split("\t") if header_line else []
            if marker_count is None:
                marker_count = sum(1 for _ in f)

        self._path_header = (marker_count, columns)
        return self._path_header

    @property
    def n_markers(self) -> int:
        if self._df is not None:
            return len(self._df)
        count, _columns = self._read_path_header()
        return count

    @property
    def n_individuals(self) -> int:
        if self._df is not None:
            return len(self._df.columns) - 2
        _count, columns = self._read_path_header()
        return max(0, len(columns) - 2)

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

    def __arrow_c_stream__(self, requested_schema=None):
        """Arrow C stream protocol for direct zero-copy consumption of the marker table.

        If the table is in-memory (narwhals-backed), this can be very efficient.
        For path-backed, it will materialize via the best available backend.
        """
        try:
            import pyarrow as pa
            if self._df is not None:
                table = from_narwhals(self._df, backend="pyarrow")
                if hasattr(table, "__arrow_c_stream__"):
                    return table.__arrow_c_stream__(requested_schema)
                return pa.table(table).__arrow_c_stream__(requested_schema)
            # Path-backed tables require explicit materialization before Arrow export.
            raise NotImplementedError(
                "Arrow stream for path-backed MarkerTable requires loading first "
                "(use MarkerTable.from_dataframe after loading, or use low-level Arrow functions)"
            )
        except Exception as e:
            raise NotImplementedError(f"Arrow C stream unavailable: {e}") from e

    def summary(self) -> str:
        return f"MarkerTable with {self.n_markers} markers across {self.n_individuals} individuals"

    # ------------------------------------------------------------------ #
    # Core analysis methods
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
        import pyrsx as _pyrsx

        if params is None:
            p = TriageParams(**{k: v for k, v in kwargs.items() if hasattr(TriageParams, k)})
        else:
            p = params

        if self._df is not None:
            markers_bytes = _table_to_ipc_bytes(
                from_narwhals(self._df, backend="pyarrow")
            )
            popmap_bytes = _arrow_bytes_from(popmap)

            arrow_table = _pyrsx.triage_to_arrow_from_arrow(
                markers_bytes,
                popmap_bytes,
                min_depth=p.min_depth,
                posterior_threshold=p.posterior_threshold,
                prior_probability=p.prior,
                linked_probability=p.linked_prob,
                group1=p.group1,
                group2=p.group2,
            )
            return TriageResult(_df=to_narwhals(arrow_table), params=p)

        # Path-backed: delegate to the low-level triage CLI binding directly.
        import tempfile

        out = tempfile.NamedTemporaryFile(delete=False, suffix=".tsv")
        out.close()
        try:
            _pyrsx.triage(
                str(self._path),
                str(popmap),
                out.name,
                min_depth=p.min_depth,
                signif_threshold=0.05,
                posterior_threshold=p.posterior_threshold,
                bayes_factor_threshold=p.bayes_factor_threshold,
                prior_probability=p.prior,
                linked_probability=p.linked_prob,
                group1=p.group1,
                group2=p.group2,
            )
            res_df = _read_core_tsv(out.name)
        finally:
            Path(out.name).unlink(missing_ok=True)
        return TriageResult(_df=res_df, params=p)

    def pca(self, *, k: int = 2, min_depth: int = 1, **kwargs: Any) -> PcaResult:
        """Compute streaming sample PCA (O(n_individuals²) memory)."""
        from .results import PcaResult
        import pyrsx as _pyrsx

        if self._df is not None:
            markers_bytes = _table_to_ipc_bytes(
                from_narwhals(self._df, backend="pyarrow")
            )
            arrow_res = _pyrsx.pca_to_arrow_from_arrow(
                markers_bytes,
                min_depth=min_depth,
                n_components=k,
            )
            res_df = to_narwhals(arrow_res["loadings"])
            return PcaResult(
                _df=res_df,
                params={"k": k, "min_depth": min_depth, "arrow": True, **kwargs},
                _input_backend=self._backend,
            )

        # Path-backed: write to a hidden directory, read loadings back.
        import tempfile

        with tempfile.TemporaryDirectory() as outdir:
            _pyrsx.pca(
                str(self._path),
                outdir,
                min_depth=min_depth,
                n_components=k,
            )
            loadings_path = Path(outdir) / "loadings.tsv"
            res_df = _read_core_tsv(loadings_path)
        return PcaResult(
            _df=res_df,
            params={"k": k, "min_depth": min_depth, **kwargs},
            _input_backend=self._backend,
        )

    # ------------------------------------------------------------------ #
    # Tabular commands: freq / depth / distrib / signif
    # ------------------------------------------------------------------ #

    def freq(self, min_depth: int = 1, **kwargs: Any) -> "TableResult":
        """Compute marker frequency table."""
        from .results import TableResult
        import pyrsx as _pyrsx

        if self._df is not None:
            ipc_bytes = _table_to_ipc_bytes(
                from_narwhals(self._df, backend="pyarrow")
            )
            raw = _pyrsx.freq_from_arrow(ipc_bytes, min_depth=min_depth)
            res_df = to_narwhals(raw)
        else:
            import tempfile

            out = tempfile.NamedTemporaryFile(delete=False, suffix=".tsv")
            out.close()
            try:
                _pyrsx.freq(str(self._path), out.name, min_depth=min_depth)
                res_df = _read_core_tsv(out.name)
            finally:
                Path(out.name).unlink(missing_ok=True)

        return TableResult(
            _df=res_df, command="freq", params={"min_depth": min_depth, **kwargs}
        )

    def depth(
        self,
        popmap: Any,
        min_frequency: float = 0.75,
        **kwargs: Any,
    ) -> "TableResult":
        """Per-sample depth statistics (requires a popmap)."""
        from .results import TableResult
        import pyrsx as _pyrsx

        if self._df is not None:
            markers_bytes = _table_to_ipc_bytes(
                from_narwhals(self._df, backend="pyarrow")
            )
            popmap_bytes = _arrow_bytes_from(popmap)
            raw = _pyrsx.depth_from_arrow(
                markers_bytes,
                popmap_bytes,
                min_frequency=min_frequency,
            )
            res_df = to_narwhals(raw)
        else:
            import tempfile

            out = tempfile.NamedTemporaryFile(delete=False, suffix=".tsv")
            out.close()
            try:
                _pyrsx.depth(
                    str(self._path),
                    str(popmap),
                    out.name,
                    min_frequency=min_frequency,
                )
                res_df = _read_core_tsv(out.name)
            finally:
                Path(out.name).unlink(missing_ok=True)

        return TableResult(
            _df=res_df,
            command="depth",
            params={"min_frequency": min_frequency, **kwargs},
        )

    def distrib(
        self,
        popmap: Any,
        group1: str = "",
        group2: str = "",
        min_depth: int = 1,
        signif_threshold: float = 0.05,
        correction: str = "bonferroni",
        test: str = "chisq",
        **kwargs: Any,
    ) -> "TableResult":
        """Distribution of markers between two groups (requires a popmap)."""
        from .results import TableResult
        import pyrsx as _pyrsx

        if self._df is not None:
            markers_bytes = _table_to_ipc_bytes(
                from_narwhals(self._df, backend="pyarrow")
            )
            popmap_bytes = _arrow_bytes_from(popmap)
            raw = _pyrsx.distrib_from_arrow(
                markers_bytes,
                popmap_bytes,
                min_depth=min_depth,
                signif_threshold=signif_threshold,
                group1=group1,
                group2=group2,
                correction=correction,
                test=test,
            )
            res_df = to_narwhals(raw)
        else:
            import tempfile

            out = tempfile.NamedTemporaryFile(delete=False, suffix=".tsv")
            out.close()
            try:
                _pyrsx.distrib(
                    str(self._path),
                    str(popmap),
                    out.name,
                    min_depth=min_depth,
                    signif_threshold=signif_threshold,
                    group1=group1,
                    group2=group2,
                    correction=correction,
                    test=test,
                )
                res_df = _read_core_tsv(out.name)
            finally:
                Path(out.name).unlink(missing_ok=True)

        return TableResult(
            _df=res_df,
            command="distrib",
            params={
                "group1": group1,
                "group2": group2,
                "min_depth": min_depth,
                "signif_threshold": signif_threshold,
                "correction": correction,
                "test": test,
                **kwargs,
            },
        )

    def signif(
        self,
        popmap: Any,
        group1: str = "",
        group2: str = "",
        min_depth: int = 1,
        signif_threshold: float = 0.05,
        correction: str = "bonferroni",
        test: str = "chisq",
        output_fasta: bool = False,
        bayes: bool = False,
        **kwargs: Any,
    ) -> "TableResult":
        """Extract significantly associated markers (requires a popmap)."""
        from .results import TableResult
        import pyrsx as _pyrsx

        if self._df is not None:
            markers_bytes = _table_to_ipc_bytes(
                from_narwhals(self._df, backend="pyarrow")
            )
            popmap_bytes = _arrow_bytes_from(popmap)
            raw = _pyrsx.signif_from_arrow(
                markers_bytes,
                popmap_bytes,
                min_depth=min_depth,
                signif_threshold=signif_threshold,
                group1=group1,
                group2=group2,
                correction=correction,
                test=test,
                output_fasta=output_fasta,
                bayes=bayes,
            )
            res_df = to_narwhals(raw)
        else:
            import tempfile

            out = tempfile.NamedTemporaryFile(delete=False, suffix=".tsv")
            out.close()
            try:
                _pyrsx.signif(
                    str(self._path),
                    str(popmap),
                    out.name,
                    min_depth=min_depth,
                    signif_threshold=signif_threshold,
                    group1=group1,
                    group2=group2,
                    correction=correction,
                    test=test,
                    output_fasta=output_fasta,
                    bayes=bayes,
                )
                res_df = _read_core_tsv(out.name)
            finally:
                Path(out.name).unlink(missing_ok=True)

        return TableResult(
            _df=res_df,
            command="signif",
            params={
                "group1": group1,
                "group2": group2,
                "min_depth": min_depth,
                "signif_threshold": signif_threshold,
                "correction": correction,
                "test": test,
                "output_fasta": output_fasta,
                "bayes": bayes,
                **kwargs,
            },
        )

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
        raise RuntimeError(
            "MarkerTable was constructed from a path and has no in-memory "
            "DataFrame. Use `.to_path()` or load it with `from_dataframe`."
        )

    def to_pandas(self) -> Any:
        """Return the in-memory table as a pandas DataFrame."""
        return self.to_dataframe(backend="pandas")

    def to_polars(self) -> Any:
        """Return the in-memory table as a Polars DataFrame."""
        return self.to_dataframe(backend="polars")
