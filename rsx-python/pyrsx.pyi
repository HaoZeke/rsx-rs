"""Type stubs for pyrsx."""

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
    """Streaming PCA of the depth matrix (Tucker mode-2 decomposition).

    Computes principal components via streaming Gram eigendecomposition.
    Memory: O(n_individuals^2), works on arbitrarily large tables.

    Args:
        table_path: Path to marker depth table.
        output_dir: Directory for eigenvalues.tsv, loadings.tsv, summary.txt.
        min_depth: Minimum depth to consider a marker present.
        n_components: Number of components to output (default: all).
    """
    ...
