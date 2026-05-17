"""Tests for pyrsx Python bindings and CLI."""

import os
import tempfile

import pytest
from click.testing import CliRunner

import pyrsx
from pyrsx.cli import main


@pytest.fixture
def test_data(tmp_path):
    """Create minimal test data."""
    # Markers table
    table = tmp_path / "markers.tsv"
    table.write_text(
        "#Number of markers : 4\n"
        "id\tsequence\tind1\tind2\tind3\tind4\tind5\n"
        "0\tATCGATCG\t10\t5\t8\t12\t7\n"
        "1\tGGGGAAAA\t15\t20\t10\t0\t0\n"
        "2\tCCCCTTTT\t0\t0\t0\t25\t30\n"
        "3\tAAAATTTT\t5\t0\t3\t8\t6\n"
    )

    # Popmap
    popmap = tmp_path / "popmap.tsv"
    popmap.write_text("ind1\tM\nind2\tM\nind3\tM\nind4\tF\nind5\tF\n")

    # FASTQ dir
    reads_dir = tmp_path / "reads"
    reads_dir.mkdir()
    for name in ["ind1", "ind2", "ind3", "ind4", "ind5"]:
        fq = reads_dir / f"{name}.fq"
        fq.write_text(
            f"@r1\nATCGATCGATCG\n+\nIIIIIIIIIIII\n"
            f"@r2\nGGGGAAAA\n+\nIIIIIIII\n"
        )

    return {
        "table": str(table),
        "popmap": str(popmap),
        "reads_dir": str(reads_dir),
        "tmp": str(tmp_path),
    }


def test_import():
    """pyrsx should be importable."""
    assert hasattr(pyrsx, "process")
    assert hasattr(pyrsx, "signif")
    assert hasattr(pyrsx, "pca")


def test_freq(test_data):
    """freq should produce frequency output."""
    out = os.path.join(test_data["tmp"], "freq.tsv")
    pyrsx.freq(test_data["table"], out, min_depth=1)
    assert os.path.exists(out)
    with open(out) as f:
        content = f.read()
    assert "Frequency\tCount" in content


def test_distrib(test_data):
    """distrib should produce distribution table."""
    out = os.path.join(test_data["tmp"], "distrib.tsv")
    pyrsx.distrib(test_data["table"], test_data["popmap"], out, min_depth=1)
    assert os.path.exists(out)
    with open(out) as f:
        content = f.read()
    assert "Markers\tP\tCorrectedP" in content


def test_signif_default(test_data):
    """signif with default settings."""
    out = os.path.join(test_data["tmp"], "signif.tsv")
    pyrsx.signif(
        test_data["table"], test_data["popmap"], out,
        min_depth=1, correction="none"
    )
    assert os.path.exists(out)


def test_signif_fisher_fdr(test_data):
    """signif with Fisher test + FDR correction."""
    out = os.path.join(test_data["tmp"], "signif_fisher.tsv")
    pyrsx.signif(
        test_data["table"], test_data["popmap"], out,
        min_depth=1, test="fisher", correction="fdr"
    )
    assert os.path.exists(out)


def test_signif_bayes(test_data):
    """signif with Bayesian output."""
    out = os.path.join(test_data["tmp"], "signif_bayes.tsv")
    pyrsx.signif(
        test_data["table"], test_data["popmap"], out,
        min_depth=1, correction="none", bayes=True
    )
    assert os.path.exists(out)
    with open(out) as f:
        header = f.readlines()[1]
    assert "Bayes_Factor" in header
    assert "Posterior_SexLinked" in header


def test_triage(test_data):
    """triage should produce marker-level candidate classes."""
    out = os.path.join(test_data["tmp"], "triage.tsv")
    pyrsx.triage(
        test_data["table"],
        test_data["popmap"],
        out,
        min_depth=1,
        group1="M",
        group2="F",
    )
    assert os.path.exists(out)
    with open(out) as f:
        content = f.read()
    assert "Candidate_Class" in content
    assert "Bayes_Factor" in content


def test_depth(test_data):
    """depth should produce per-individual stats."""
    out = os.path.join(test_data["tmp"], "depth.tsv")
    pyrsx.depth(test_data["table"], test_data["popmap"], out, min_frequency=0.5)
    assert os.path.exists(out)
    with open(out) as f:
        content = f.read()
    assert "Sample\tGroup" in content
    assert "ind1\tM" in content


def test_pca(test_data):
    """pca should produce eigenvalues and loadings."""
    out_dir = os.path.join(test_data["tmp"], "pca")
    pyrsx.pca(test_data["table"], out_dir, min_depth=1, n_components=3)
    assert os.path.exists(os.path.join(out_dir, "eigenvalues.tsv"))
    assert os.path.exists(os.path.join(out_dir, "loadings.tsv"))
    assert os.path.exists(os.path.join(out_dir, "summary.txt"))


def test_process(test_data):
    """process should create marker depth table from FASTQ."""
    out = os.path.join(test_data["tmp"], "processed.tsv")
    pyrsx.process(test_data["reads_dir"], out, threads=2, min_depth=1)
    assert os.path.exists(out)
    with open(out) as f:
        content = f.read()
    assert "#Number of markers" in content
    assert "id\tsequence" in content


def test_merge(test_data):
    """merge should combine two tables."""
    # Create second table
    table2 = os.path.join(test_data["tmp"], "table2.tsv")
    with open(table2, "w") as f:
        f.write("#Number of markers : 2\n")
        f.write("id\tsequence\tind6\tind7\n")
        f.write("0\tATCGATCG\t3\t9\n")
        f.write("1\tTTTTAAAA\t5\t7\n")

    out = os.path.join(test_data["tmp"], "merged.tsv")
    pyrsx.merge([test_data["table"], table2], out)
    assert os.path.exists(out)
    with open(out) as f:
        content = f.read()
    assert "#Number of markers" in content


def test_invalid_test_method(test_data):
    """Invalid test method should raise error."""
    out = os.path.join(test_data["tmp"], "bad.tsv")
    with pytest.raises(RuntimeError, match="Unknown test method"):
        pyrsx.signif(
            test_data["table"], test_data["popmap"], out,
            test="invalid_test"
        )


def test_invalid_correction(test_data):
    """Invalid correction should raise error."""
    out = os.path.join(test_data["tmp"], "bad.tsv")
    with pytest.raises(RuntimeError, match="Unknown correction"):
        pyrsx.signif(
            test_data["table"], test_data["popmap"], out,
            correction="invalid_corr"
        )


# === CLI tests ===

def test_cli_help():
    """CLI --help should show all commands."""
    runner = CliRunner()
    result = runner.invoke(main, ["--help"])
    assert result.exit_code == 0
    assert "process" in result.output
    assert "signif" in result.output
    assert "pca" in result.output
    assert "merge" in result.output


def test_cli_freq(test_data):
    """CLI freq command should work."""
    runner = CliRunner()
    out = os.path.join(test_data["tmp"], "cli_freq.tsv")
    result = runner.invoke(main, [
        "freq", "-t", test_data["table"], "-o", out, "-d", "1",
    ])
    assert result.exit_code == 0, result.output
    assert os.path.exists(out)


def test_cli_signif_bayes(test_data):
    """CLI signif with --bayes and --test fisher."""
    runner = CliRunner()
    out = os.path.join(test_data["tmp"], "cli_signif.tsv")
    result = runner.invoke(main, [
        "signif",
        "-t", test_data["table"],
        "-p", test_data["popmap"],
        "-o", out,
        "--correction", "none",
        "--test", "fisher",
        "--bayes",
    ])
    assert result.exit_code == 0, result.output
    assert os.path.exists(out)


def test_cli_triage(test_data):
    """CLI triage command."""
    runner = CliRunner()
    out = os.path.join(test_data["tmp"], "cli_triage.tsv")
    result = runner.invoke(main, [
        "triage",
        "-t", test_data["table"],
        "-p", test_data["popmap"],
        "-o", out,
        "-G", "M,F",
    ])
    assert result.exit_code == 0, result.output
    assert os.path.exists(out)


def test_cli_pca(test_data):
    """CLI pca command."""
    runner = CliRunner()
    out_dir = os.path.join(test_data["tmp"], "cli_pca")
    result = runner.invoke(main, [
        "pca", "-t", test_data["table"], "-o", out_dir, "-r", "3",
    ])
    assert result.exit_code == 0, result.output
    assert os.path.exists(os.path.join(out_dir, "eigenvalues.tsv"))
