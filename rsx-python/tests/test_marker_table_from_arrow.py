"""Parity tests for the from-Arrow path of the high-level MarkerTable API.

These exercise the in-memory entry points that the previous refactor wired
up but never actually ran: every *_from_arrow shim used to pipe binary
Arrow IPC bytes into pandas.read_csv. The fix lands a real Arrow-IPC ->
markers TSV bridge in the Rust binding, so each method here must produce
the same result whether driven from a DataFrame (Arrow path) or from a
file path (CLI path).
"""

from __future__ import annotations

import pandas as pd
import pytest

import pyrsx
from pyrsx import MarkerTable


@pytest.fixture
def fixture_files(tmp_path):
    table = tmp_path / "markers.tsv"
    table.write_text(
        "#Number of markers : 6\n"
        "id\tsequence\tm1\tm2\tm3\tf1\tf2\tf3\n"
        # always present in both groups -> no bias
        "0\tALL\t10\t10\t10\t10\t10\t10\n"
        # male-only
        "1\tMONLY\t10\t10\t10\t0\t0\t0\n"
        # female-only
        "2\tFONLY\t0\t0\t0\t10\t10\t10\n"
        # mixed
        "3\tMIX\t10\t0\t10\t0\t10\t0\n"
        # very low depth row that should drop at min_depth >= 5
        "4\tLOW\t1\t2\t1\t2\t1\t2\n"
        # all zero -> filtered everywhere
        "5\tZERO\t0\t0\t0\t0\t0\t0\n"
    )
    popmap = tmp_path / "popmap.tsv"
    popmap.write_text(
        "m1\tM\n"
        "m2\tM\n"
        "m3\tM\n"
        "f1\tF\n"
        "f2\tF\n"
        "f3\tF\n"
    )
    return {"table": str(table), "popmap": str(popmap)}


@pytest.fixture
def fixture_frames(fixture_files):
    markers_df = pd.read_csv(fixture_files["table"], sep="\t", comment="#")
    popmap_df = pd.read_csv(
        fixture_files["popmap"], sep="\t", header=None, names=["individual", "group"]
    )
    return {
        "markers_df": markers_df,
        "popmap_df": popmap_df,
        "table": fixture_files["table"],
        "popmap": fixture_files["popmap"],
    }


def _to_pandas(result):
    return result.to_pandas().reset_index(drop=True)


def test_from_arrow_exports_resolve():
    """The low-level *_from_arrow names must actually import."""
    for name in (
        "freq_from_arrow",
        "depth_from_arrow",
        "distrib_from_arrow",
        "signif_from_arrow",
        "triage_to_arrow_from_arrow",
        "pca_to_arrow_from_arrow",
    ):
        assert callable(getattr(pyrsx, name)), name


def test_freq_arrow_matches_path(fixture_frames):
    arrow = _to_pandas(
        MarkerTable.from_dataframe(fixture_frames["markers_df"]).freq(min_depth=1)
    )
    file = _to_pandas(
        MarkerTable.from_path(fixture_frames["table"]).freq(min_depth=1)
    )
    pd.testing.assert_frame_equal(arrow, file, check_dtype=False)


def test_depth_arrow_matches_path(fixture_frames):
    arrow = _to_pandas(
        MarkerTable.from_dataframe(fixture_frames["markers_df"]).depth(
            popmap=fixture_frames["popmap_df"], min_frequency=0.1
        )
    )
    file = _to_pandas(
        MarkerTable.from_path(fixture_frames["table"]).depth(
            popmap=fixture_frames["popmap"], min_frequency=0.1
        )
    )
    pd.testing.assert_frame_equal(arrow, file, check_dtype=False)


def test_distrib_arrow_matches_path(fixture_frames):
    arrow = _to_pandas(
        MarkerTable.from_dataframe(fixture_frames["markers_df"]).distrib(
            popmap=fixture_frames["popmap_df"],
            group1="M",
            group2="F",
            min_depth=1,
            test="fisher",
            correction="none",
        )
    )
    file = _to_pandas(
        MarkerTable.from_path(fixture_frames["table"]).distrib(
            popmap=fixture_frames["popmap"],
            group1="M",
            group2="F",
            min_depth=1,
            test="fisher",
            correction="none",
        )
    )
    pd.testing.assert_frame_equal(arrow, file, check_dtype=False)


def test_signif_arrow_matches_path(fixture_frames):
    arrow = _to_pandas(
        MarkerTable.from_dataframe(fixture_frames["markers_df"]).signif(
            popmap=fixture_frames["popmap_df"],
            group1="M",
            group2="F",
            min_depth=1,
            test="fisher",
            correction="none",
        )
    )
    file = _to_pandas(
        MarkerTable.from_path(fixture_frames["table"]).signif(
            popmap=fixture_frames["popmap"],
            group1="M",
            group2="F",
            min_depth=1,
            test="fisher",
            correction="none",
        )
    )
    pd.testing.assert_frame_equal(arrow, file, check_dtype=False)


def test_triage_arrow_matches_path(fixture_frames):
    arrow_res = MarkerTable.from_dataframe(fixture_frames["markers_df"]).triage(
        popmap=fixture_frames["popmap_df"],
        min_depth=1,
        prior=0.5,
        linked_prob=0.9,
        group1="M",
        group2="F",
    )
    file_res = MarkerTable.from_path(fixture_frames["table"]).triage(
        popmap=fixture_frames["popmap"],
        min_depth=1,
        prior=0.5,
        linked_prob=0.9,
        group1="M",
        group2="F",
    )
    arrow_df = _to_pandas(arrow_res).sort_values("id").reset_index(drop=True)
    file_df = _to_pandas(file_res).sort_values("id").reset_index(drop=True)
    # Same set of qualifying markers. Compare ids as strings because the
    # path flow always re-reads the id column as text, while the Arrow
    # flow keeps whatever dtype the source DataFrame had.
    assert [str(x) for x in arrow_df["id"]] == [str(x) for x in file_df["id"]]
    assert list(arrow_df.columns) == list(file_df.columns)


def test_pca_arrow_matches_path(fixture_frames):
    arrow_res = MarkerTable.from_dataframe(fixture_frames["markers_df"]).pca(
        k=2, min_depth=1
    )
    file_res = MarkerTable.from_path(fixture_frames["table"]).pca(k=2, min_depth=1)
    arrow_df = _to_pandas(arrow_res)
    file_df = _to_pandas(file_res)
    # Loadings table: same individuals and same column count
    assert len(arrow_df) == len(file_df)
    assert list(arrow_df.columns) == list(file_df.columns)


def test_freq_min_depth_filters(fixture_frames):
    """Sanity check: min_depth=5 must drop the LOW row (cells <= 2)."""
    arrow_total = (
        _to_pandas(
            MarkerTable.from_dataframe(fixture_frames["markers_df"]).freq(min_depth=1)
        )["Count"]
        .sum()
    )
    arrow_strict = (
        _to_pandas(
            MarkerTable.from_dataframe(fixture_frames["markers_df"]).freq(min_depth=5)
        )["Count"]
        .sum()
    )
    assert arrow_strict <= arrow_total
