import tempfile
import unittest
from pathlib import Path

from benchmarks.merge_literature_results import merge_csv_files
from benchmarks.run_literature_bindings import count_data_rows, resolve_markers_table


class LiteratureBindingsFeatureTests(unittest.TestCase):
    def test_resolve_markers_table_prefers_main_workflow_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            dataset = Path(tmp)
            main = dataset / "markers_table.tsv"
            comparison = dataset / "comparison" / "rust" / "markers_table.tsv"
            comparison.parent.mkdir(parents=True)
            main.write_text("main\n")
            comparison.write_text("comparison\n")

            self.assertEqual(resolve_markers_table(dataset), main)

    def test_resolve_markers_table_accepts_comparison_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            dataset = Path(tmp)
            comparison = dataset / "comparison" / "rust" / "markers_table.tsv"
            comparison.parent.mkdir(parents=True)
            comparison.write_text("comparison\n")

            self.assertEqual(resolve_markers_table(dataset), comparison)

    def test_count_data_rows_skips_comments_and_header(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "table.tsv"
            path.write_text("#meta\nid\tsequence\ts1\n0\tAAAA\t1\n1\tCCCC\t2\n")

            self.assertEqual(count_data_rows(path), 2)

    def test_count_data_rows_counts_string_identifiers_after_header(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "depth.tsv"
            path.write_text("sample\tgroup\tdepth\nERR0001\tF\t17\nERR0002\tM\t23\n")

            self.assertEqual(count_data_rows(path), 2)

    def test_merge_csv_files_writes_single_header(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            a = root / "a.csv"
            b = root / "b.csv"
            out = root / "merged.csv"
            a.write_text("dataset,command\nalpha,process\n")
            b.write_text("dataset,command\nbeta,depth\n")

            merge_csv_files([a, b], out)

            self.assertEqual(out.read_text(), "dataset,command\nalpha,process\nbeta,depth\n")


if __name__ == "__main__":
    unittest.main()
