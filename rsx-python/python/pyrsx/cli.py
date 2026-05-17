"""pyrsx CLI: command-line interface for rsx via click."""

import click

from . import depth as _depth
from . import distrib as _distrib
from . import freq as _freq
from . import merge as _merge
from . import pca as _pca
from . import process as _process
from . import signif as _signif
from . import triage as _triage


@click.group()
@click.version_option(package_name="pyrsx")
def main():
    """pyrsx: high-performance RAD-seq sex determination toolkit.

    Rust-powered analysis of RAD-seq data for sex-linked marker detection.
    All commands operate in bounded memory for arbitrarily large datasets.
    """


@main.command()
@click.option("-i", "--input-dir", required=True, help="Directory of FASTQ/FASTA files.")
@click.option("-o", "--output-file", required=True, help="Output marker depth table.")
@click.option("-T", "--threads", default=1, show_default=True, help="Number of threads.")
@click.option("-d", "--min-depth", default=1, show_default=True, help="Min depth to retain marker.")
@click.option("-k", "--kmer-dedup", default=None, type=int, help="K-mer size for deduplication.")
def process(input_dir, output_file, threads, min_depth, kmer_dedup):
    """Build a marker depth table from demultiplexed reads."""
    _process(input_dir, output_file, threads=threads, min_depth=min_depth, kmer_dedup=kmer_dedup)
    click.echo(f"Wrote {output_file}")


@main.command()
@click.option("-t", "--markers-table", required=True, help="Marker depth table from process.")
@click.option("-p", "--popmap", required=True, help="Population map (individual<tab>group).")
@click.option("-o", "--output-file", required=True, help="Output distribution table.")
@click.option("-d", "--min-depth", default=1, show_default=True)
@click.option("-S", "--signif-threshold", default=0.05, show_default=True)
@click.option("-G", "--groups", default="", help="Groups to compare (comma-separated).")
@click.option("--correction", default="bonferroni", show_default=True,
              type=click.Choice(["bonferroni", "fdr", "none"]))
@click.option("--test", "test_method", default="chisq", show_default=True,
              type=click.Choice(["chisq", "fisher", "gtest"]))
def distrib(markers_table, popmap, output_file, min_depth, signif_threshold,
            groups, correction, test_method):
    """Compute marker distribution between two groups."""
    g1, g2 = _parse_groups(groups)
    _distrib(markers_table, popmap, output_file, min_depth=min_depth,
             signif_threshold=signif_threshold, group1=g1, group2=g2,
             correction=correction, test=test_method)
    click.echo(f"Wrote {output_file}")


@main.command()
@click.option("-t", "--markers-table", required=True)
@click.option("-p", "--popmap", required=True)
@click.option("-o", "--output-file", required=True)
@click.option("-d", "--min-depth", default=1, show_default=True)
@click.option("-S", "--signif-threshold", default=0.05, show_default=True)
@click.option("-G", "--groups", default="")
@click.option("--correction", default="bonferroni", show_default=True,
              type=click.Choice(["bonferroni", "fdr", "none"]))
@click.option("--test", "test_method", default="chisq", show_default=True,
              type=click.Choice(["chisq", "fisher", "gtest"]))
@click.option("-a", "--output-fasta", is_flag=True, help="Output FASTA instead of table.")
@click.option("--bayes", is_flag=True, help="Add Bayes Factor + posterior columns.")
def signif(markers_table, popmap, output_file, min_depth, signif_threshold,
           groups, correction, test_method, output_fasta, bayes):
    """Extract markers significantly associated with a group."""
    g1, g2 = _parse_groups(groups)
    _signif(markers_table, popmap, output_file, min_depth=min_depth,
            signif_threshold=signif_threshold, group1=g1, group2=g2,
            correction=correction, test=test_method,
            output_fasta=output_fasta, bayes=bayes)
    click.echo(f"Wrote {output_file}")


@main.command()
@click.option("-t", "--markers-table", required=True)
@click.option("-p", "--popmap", required=True)
@click.option("-o", "--output-file", required=True)
@click.option("-d", "--min-depth", default=1, show_default=True)
@click.option("-S", "--signif-threshold", default=0.05, show_default=True)
@click.option("--posterior-threshold", default=0.9, show_default=True)
@click.option("--bayes-factor-threshold", default=10.0, show_default=True)
@click.option("--prior-probability", default=0.01, show_default=True)
@click.option("--linked-probability", default=0.9, show_default=True)
@click.option("-G", "--groups", default="")
def triage(markers_table, popmap, output_file, min_depth, signif_threshold,
           posterior_threshold, bayes_factor_threshold, prior_probability,
           linked_probability, groups):
    """Rank strict and Bayesian marker candidates for biological follow-up."""
    g1, g2 = _parse_groups(groups)
    _triage(markers_table, popmap, output_file, min_depth=min_depth,
            signif_threshold=signif_threshold,
            posterior_threshold=posterior_threshold,
            bayes_factor_threshold=bayes_factor_threshold,
            prior_probability=prior_probability,
            linked_probability=linked_probability,
            group1=g1, group2=g2)
    click.echo(f"Wrote {output_file}")


@main.command()
@click.option("-t", "--markers-table", required=True)
@click.option("-o", "--output-file", required=True)
@click.option("-d", "--min-depth", default=1, show_default=True)
def freq(markers_table, output_file, min_depth):
    """Compute marker frequency distribution."""
    _freq(markers_table, output_file, min_depth=min_depth)
    click.echo(f"Wrote {output_file}")


@main.command()
@click.option("-t", "--markers-table", required=True)
@click.option("-p", "--popmap", required=True)
@click.option("-o", "--output-file", required=True)
@click.option("-f", "--min-frequency", default=0.75, show_default=True)
def depth(markers_table, popmap, output_file, min_frequency):
    """Compute retained read statistics per individual."""
    _depth(markers_table, popmap, output_file, min_frequency=min_frequency)
    click.echo(f"Wrote {output_file}")


@main.command()
@click.option("-o", "--output-file", required=True, help="Output merged table.")
@click.option("-B", "--buffer-size", default=2000000, show_default=True,
              help="Entries to buffer before flushing to disk.")
@click.option("--output-parquet", is_flag=True, help="Output Parquet instead of TSV.")
@click.argument("input_files", nargs=-1, required=True)
def merge(output_file, buffer_size, output_parquet, input_files):
    """Merge multiple marker depth tables by sequence identity.

    Uses bounded-memory external sort. Handles 75M+ sequences.
    """
    _merge(list(input_files), output_file, buffer_size=buffer_size,
           output_parquet=output_parquet)
    click.echo(f"Wrote {output_file}")


@main.command()
@click.option("-t", "--markers-table", required=True)
@click.option("-o", "--output-dir", required=True, help="Directory for PCA results.")
@click.option("-d", "--min-depth", default=1, show_default=True)
@click.option("-r", "--components", default=None, type=int,
              help="Number of components (default: all).")
def pca(markers_table, output_dir, min_depth, components):
    """Streaming PCA of the depth matrix.

    Memory: O(n_individuals^2). Works on arbitrarily large tables.
    """
    _pca(markers_table, output_dir, min_depth=min_depth, n_components=components)
    click.echo(f"Results in {output_dir}/")


def _parse_groups(groups_str):
    """Parse comma-separated group names."""
    if not groups_str:
        return "", ""
    parts = groups_str.split(",")
    if len(parts) >= 2:
        return parts[0].strip(), parts[1].strip()
    return parts[0].strip(), ""
