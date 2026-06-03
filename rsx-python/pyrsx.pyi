"""Type stubs for pyrsx (high-level + low-level API).

High-level API uses narwhals for backend-agnostic DataFrames.
Internal reading of rsx core TSV outputs (for path-backed results) is
done via pyarrow + to_narwhals (see _read_core_tsv); no pandas is
required internally. Results are always narwhals.DataFrame backed.
"""

# --------------------------------------------------------------------------- #
# High-level API (recommended)
# --------------------------------------------------------------------------- #

from dataclasses import dataclass
from typing import Any, Literal

import narwhals as nw  # type: ignore[import-untyped]

class MarkerTable:
    """High-level representation of a RAD-seq marker depth table.

    Data and results are narwhals DataFrames (backend agnostic).
    Use .to_pandas() / .to_polars() etc. for specific backends.
    Path-backed construction uses _read_core_tsv (narwhals via pyarrow)
    for outputs to avoid forcing any particular DF library.
    """

    @classmethod
    def from_path(
        cls,
        path: str | Any,  # PathLike
        *,
        backend: Literal["pandas", "polars", "pyarrow", "auto"] = "auto",
    ) -> MarkerTable: ...
    @classmethod
    def from_dataframe(
        cls,
        df: Any,
        *,
        backend: Literal["pandas", "polars", "pyarrow", "auto"] = "auto",
    ) -> MarkerTable: ...

    @property
    def n_markers(self) -> int: ...
    @property
    def n_individuals(self) -> int: ...

    @property
    def data(self) -> nw.DataFrame | None: ...
    @property
    def df(self) -> nw.DataFrame | None: ...

    def __len__(self) -> int: ...
    def __repr__(self) -> str: ...
    def summary(self) -> str: ...

    def to_dataframe(self, *, backend: Literal["pandas", "polars", "pyarrow", "auto"] = "auto") -> Any: ...
    def to_pandas(self) -> Any: ...
    def to_polars(self) -> Any: ...

    def __arrow_c_stream__(self, requested_schema: Any = None) -> Any: ...

    # Analysis methods (return narwhals-backed results)
    def triage(self, *, popmap: Any, params: "TriageParams | None" = None, **kwargs: Any) -> "TriageResult": ...
    def pca(self, *, k: int = 2, min_depth: int = 1, **kwargs: Any) -> "PcaResult": ...
    def freq(self, min_depth: int = 1, **kwargs: Any) -> "TableResult": ...
    def depth(self, popmap: Any, min_frequency: float = 0.75, **kwargs: Any) -> "TableResult": ...
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
    ) -> "TableResult": ...
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
    ) -> "TableResult": ...


class TriageResult:
    """Triage result (strict + Bayesian). Backed by narwhals DataFrame."""

    @property
    def df(self) -> nw.DataFrame: ...
    def to_dataframe(self, *, backend: str = "auto") -> Any: ...
    def to_pandas(self) -> Any: ...
    def to_polars(self) -> Any: ...
    def summary(self) -> str: ...
    def plot_evidence(self, **kwargs: Any) -> Any: ...
    def __arrow_c_stream__(self, requested_schema: Any = None) -> Any: ...
    def __getattr__(self, name: str) -> Any: ...


class PcaResult:
    """PCA result. Backed by narwhals DataFrame."""

    @property
    def df(self) -> nw.DataFrame: ...
    def to_dataframe(self, *, backend: str | None = None) -> Any: ...
    def to_pandas(self) -> Any: ...
    def to_polars(self) -> Any: ...
    def summary(self) -> str: ...
    def plot(self, **kwargs: Any) -> Any: ...
    def __arrow_c_stream__(self, requested_schema: Any = None) -> Any: ...
    def __getattr__(self, name: str) -> Any: ...


class TableResult:
    """Generic table result for freq/depth/distrib/signif etc. Narwhals backed."""

    command: str
    params: dict[str, Any]
    @property
    def df(self) -> nw.DataFrame: ...
    def to_dataframe(self, *, backend: str = "pandas") -> Any: ...
    def to_pandas(self) -> Any: ...
    def to_polars(self) -> Any: ...
    def __arrow_c_stream__(self, requested_schema: Any = None) -> Any: ...
    def __getattr__(self, name: str) -> Any: ...
    def __repr__(self) -> str: ...


@dataclass
class TriageParams:
    min_depth: int = 10
    posterior_threshold: float = 0.9
    prior: float = 0.01
    linked_prob: float = 0.9
    group1: str = "M"
    group2: str = "F"
    ...

# --------------------------------------------------------------------------- #
# Low-level direct bindings (for power users / compatibility)
# --------------------------------------------------------------------------- #

def process(
    input_dir: str,
    output_file: str,
    threads: int = 1,
    min_depth: int = 1,
    kmer_dedup: int | None = None,
) -> None:
    """Process demultiplexed FASTQ/FASTA reads into a marker depth table.

    Args:
        input_dir: Directory containing demultiplexed sequence files.
        output_file: Path to the output TSV marker depth table.
        threads: Number of threads for parallel processing.
        min_depth: Minimum depth in at least one individual to retain a marker.
        kmer_dedup: If set, group markers by canonical k-mer of this size.
    """
    ...

def distrib(
    table_path: str,
    popmap_path: str,
    output_file: str,
    min_depth: int = 1,
    signif_threshold: float = 0.05,
    group1: str = "",
    group2: str = "",
    correction: str = "bonferroni",
    test: str = "chisq",
) -> None:
    """Compute marker distribution between two groups.

    Args:
        table_path: Path to marker depth table from process.
        popmap_path: Path to population map (individual<tab>group).
        output_file: Path to output distribution table.
        min_depth: Minimum depth to consider a marker present.
        signif_threshold: P-value threshold for significance.
        group1: Name of first group (auto-detected if empty).
        group2: Name of second group.
        correction: Multiple testing correction: "bonferroni", "fdr", "none".
        test: Statistical test: "chisq", "fisher", "gtest".
    """
    ...

def signif(
    table_path: str,
    popmap_path: str,
    output_file: str,
    min_depth: int = 1,
    signif_threshold: float = 0.05,
    group1: str = "",
    group2: str = "",
    correction: str = "bonferroni",
    test: str = "chisq",
    output_fasta: bool = False,
    bayes: bool = False,
) -> None:
    """Extract markers significantly associated with a group.

    Args:
        table_path: Path to marker depth table.
        popmap_path: Path to population map.
        output_file: Path to output file.
        min_depth: Minimum depth to consider a marker present.
        signif_threshold: P-value threshold.
        group1: Name of first group.
        group2: Name of second group.
        correction: Correction method: "bonferroni", "fdr", "none".
        test: Test method: "chisq", "fisher", "gtest".
        output_fasta: If True, output FASTA instead of table.
        bayes: If True, add Bayes Factor and posterior columns.
    """
    ...

def triage(
    table_path: str,
    popmap_path: str,
    output_file: str,
    min_depth: int = 1,
    signif_threshold: float = 0.05,
    posterior_threshold: float = 0.9,
    bayes_factor_threshold: float = 10.0,
    prior_probability: float = 0.01,
    linked_probability: float = 0.9,
    group1: str = "",
    group2: str = "",
) -> None:
    """Rank strict and Bayesian marker candidates.

    Args:
        table_path: Path to marker depth table.
        popmap_path: Path to population map.
        output_file: Path to triage table.
        min_depth: Minimum depth to consider a marker present.
        signif_threshold: Family-wise p-value threshold for strict calls.
        posterior_threshold: Posterior P(sex-linked) threshold.
        bayes_factor_threshold: Bayes factor threshold.
        prior_probability: Prior probability that a marker is sex-linked.
        linked_probability: Expected marker prevalence in the linked sex.
        group1: Name of first group.
        group2: Name of second group.
    """
    ...

def freq(
    table_path: str,
    output_file: str,
    min_depth: int = 1,
) -> None:
    """Compute marker frequency distribution.

    Args:
        table_path: Path to marker depth table.
        output_file: Path to output frequency table.
        min_depth: Minimum depth to consider a marker present.
    """
    ...

def depth(
    table_path: str,
    popmap_path: str,
    output_file: str,
    min_frequency: float = 0.75,
) -> None:
    """Compute retained read statistics per individual.

    Args:
        table_path: Path to marker depth table.
        popmap_path: Path to population map.
        output_file: Path to output depth statistics.
        min_frequency: Minimum fraction of individuals for a marker to count.
    """
    ...

def merge(
    input_files: list[str],
    output_file: str,
    buffer_size: int = 2000000,
    output_parquet: bool = False,
) -> None:
    """Merge multiple marker depth tables by sequence identity.

    Uses bounded-memory external sort (configurable buffer).
    Handles 75M+ unique sequences without OOM.

    Args:
        input_files: List of paths to input marker depth tables.
        output_file: Path to output merged table.
        buffer_size: Number of entries to buffer before flushing to disk.
        output_parquet: If True, output Parquet instead of TSV.
    """
    ...

def pca(
    table_path: str,
    output_dir: str,
    min_depth: int = 1,
    n_components: int | None = None,
) -> None:
    """Streaming PCA of the depth matrix.

    Computes principal components via streaming Gram eigendecomposition.
    Memory: O(n_individuals^2), works on arbitrarily large tables.

    Args:
        table_path: Path to marker depth table.
        output_dir: Directory for eigenvalues.tsv, loadings.tsv, summary.txt.
        min_depth: Minimum depth to consider a marker present.
        n_components: Number of components to output (default: all).
    """
    ...

class PyrsxError(RuntimeError):
    """Base exception for all pyrsx errors (more specific than bare RuntimeError)."""
    ...
