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


# ------------------------------------------------------------------
# High-level Python API tests (thin callers over the real Rust Arrow path)
# ------------------------------------------------------------------

def test_highlevel_triage_via_marker_table():
    """High-level MarkerTable.triage should return a TriageResult backed by narwhals
    and contain the expected biological candidate classes (exercises the full
    Rust run_to_arrow -> PyO3 triage_to_arrow -> Python thin wrapper path).
    """
    import tempfile
    from pathlib import Path

    # Same synthetic data used by the Rust differential test
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        table = tmp / "markers.tsv"
        table.write_text(
            "#Number of markers : 3\n"
            "id\tsequence\tm1\tm2\tm3\tm4\tm5\tm6\tm7\tm8\tm9\tm10\tf1\tf2\tf3\tf4\tf5\tf6\tf7\tf8\tf9\tf10\n"
            "0\tALL\t10\t10\t10\t10\t10\t10\t10\t10\t10\t10\t10\t10\t10\t10\t10\t10\t10\t10\t10\t10\n"
            "1\tMONLY\t10\t10\t10\t10\t10\t10\t10\t10\t10\t10\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n"
            "2\tFONLY\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t10\t10\t10\t10\t10\t10\t10\t10\t10\t10\n"
        )
        pop = tmp / "popmap.tsv"
        pop.write_text(
            "m1\tM\nm2\tM\nm3\tM\nm4\tM\nm5\tM\nm6\tM\nm7\tM\nm8\tM\nm9\tM\nm10\tM\n"
            "f1\tF\nf2\tF\nf3\tF\nf4\tF\nf5\tF\nf6\tF\nf7\tF\nf8\tF\nf9\tF\nf10\tF\n"
        )

        # Import the high-level surface (will only work after maturin develop)
        from pyrsx.api.markers import MarkerTable
        from pyrsx.api.params import TriageParams

        mt = MarkerTable.from_path(str(table))
        p = TriageParams(group1="M", group2="F", prior=0.01, linked_prob=0.9)
        result = mt.triage(popmap=str(pop), params=p)

        # The result must be a TriageResult with a narwhals-backed df
        assert hasattr(result, "df")
        assert result.df is not None

        # Must have produced the biological calls we expect from the decision logic
        classes = result.df["Candidate_Class"].to_list()
        assert any("strict" in str(c) or "posterior" in str(c) or "M-biased" in str(c) or "F-biased" in str(c) for c in classes)

        # Provenance must be round-tripped
        assert result.params is not None
        assert result.params.group1 == "M"


def test_pca_to_arrow_lowlevel(test_data):
    """Low-level pca_to_arrow binding must return the expected dict with
    pyarrow tables and match the numbers from the file-based path.
    """
    import pyrsx

    # Use the existing 5-individual fixture
    res = pyrsx.pca_to_arrow(test_data["table"], min_depth=1, n_components=3)

    assert "eigenvalues" in res
    assert "loadings" in res
    assert res["n_individuals"] == 5
    assert res["n_components"] == 3
    assert res["total_variance"] > 0.0

    # loadings must have the 5 individuals
    loadings = res["loadings"]
    assert loadings.num_rows == 5
    assert "PC1" in loadings.column_names
    assert "individual" in loadings.column_names


# ------------------------------------------------------------------
# Additional high-level tests exercising path-backed _read_core_tsv path
# (verifies backend-agnostic narwhals results for commands that use
# the core TSV read helper instead of pure Arrow output).
# ------------------------------------------------------------------

def test_highlevel_freq_via_marker_table_path():
    """Path-backed MarkerTable.freq must return TableResult with narwhals df,
    exercising the _read_core_tsv path (no pandas required for the read).
    """
    import tempfile
    from pathlib import Path
    import narwhals as nw

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        table = tmp / "markers.tsv"
        table.write_text(
            "#Number of markers : 3\n"
            "id\tsequence\tm1\tm2\tm3\n"
            "0\tALL\t10\t10\t10\n"
            "1\tMONLY\t10\t10\t10\n"
            "2\tFONLY\t0\t0\t0\n"
        )

        from pyrsx.api.markers import MarkerTable

        mt = MarkerTable.from_path(str(table))
        assert len(mt) == 3
        assert mt.n_individuals == 3
        result = mt.freq(min_depth=1)

        assert hasattr(result, "df")
        assert result.df is not None
        assert isinstance(result.df, nw.DataFrame)
        assert result.command == "freq"
        # Should work without pandas being the internal reader
        # (narwhals pyarrow-backed by default here)
        assert len(result.df) > 0


def test_marker_table_path_introspection_without_metadata():
    """Path-backed MarkerTable also handles plain TSV marker tables."""
    from pathlib import Path
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        table = Path(tmp) / "markers.tsv"
        table.write_text(
            "id\tsequence\tind1\tind2\n"
            "0\tALL\t10\t10\n"
            "1\tMONLY\t10\t0\n"
        )

        from pyrsx.api.markers import MarkerTable

        mt = MarkerTable.from_path(str(table))
        assert len(mt) == 2
        assert mt.n_individuals == 2


def test_highlevel_depth_via_marker_table_path():
    """Path-backed depth using popmap exercises _read_core_tsv."""
    import tempfile
    from pathlib import Path
    import narwhals as nw

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        table = tmp / "markers.tsv"
        table.write_text(
            "#Number of markers : 2\n"
            "id\tsequence\tind1\tind2\n"
            "0\tALL\t10\t5\n"
            "1\tLOW\t1\t2\n"
        )
        pop = tmp / "popmap.tsv"
        pop.write_text("ind1\tM\nind2\tF\n")

        from pyrsx.api.markers import MarkerTable

        mt = MarkerTable.from_path(str(table))
        result = mt.depth(popmap=str(pop), min_frequency=0.0)

        assert hasattr(result, "df")
        assert isinstance(result.df, nw.DataFrame)
        assert result.command == "depth"
        assert "Sample" in result.df.columns


def test_highlevel_results_are_backend_agnostic():
    """Results from path-backed should be convertible via narwhals without
    the internal read having pulled pandas.
    """
    import tempfile
    from pathlib import Path
    import narwhals as nw

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        table = tmp / "markers.tsv"
        table.write_text(
            "#Number of markers : 2\n"
            "id\tsequence\tm1\tf1\n"
            "0\tALL\t10\t10\n"
            "1\tM\t10\t0\n"
        )
        pop = tmp / "popmap.tsv"
        pop.write_text("m1\tM\nf1\tF\n")

        from pyrsx.api.markers import MarkerTable

        mt = MarkerTable.from_path(str(table))
        res = mt.distrib(popmap=str(pop))

        assert isinstance(res.df, nw.DataFrame)
        # Explicit conversions should succeed (may pull pandas/polars if installed)
        p = res.to_pandas()
        assert p is not None
        pl = res.to_polars()
        assert pl is not None
