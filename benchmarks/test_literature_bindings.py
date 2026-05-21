import tempfile
import unittest
from pathlib import Path
import tomllib

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

    def test_merge_csv_files_keeps_latest_dataset_shard(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            old = root / "literature_speed_comparison_alpha_101.csv"
            new = root / "literature_speed_comparison_alpha_202.csv"
            beta = root / "literature_speed_comparison_beta_202.csv"
            out = root / "merged.csv"
            old.write_text("dataset,command,elapsed_seconds\nalpha,process,10\n")
            new.write_text("dataset,command,elapsed_seconds\nalpha,process,4\nalpha,depth,1\n")
            beta.write_text("dataset,command,elapsed_seconds\nbeta,process,7\n")

            old.touch()
            beta.touch()
            new.touch()
            merge_csv_files([old, new, beta], out)

            self.assertEqual(
                out.read_text(),
                "dataset,command,elapsed_seconds\n"
                "alpha,process,4\n"
                "alpha,depth,1\n"
                "beta,process,7\n",
            )

    def test_run_literature_bindings_task_covers_all_panel_datasets(self):
        pixi = tomllib.loads(Path("pixi.toml").read_text())
        command = pixi["feature"]["python"]["tasks"]["run-literature-bindings"]["cmd"]

        for dataset in [
            "danio_albolineatus",
            "notothenia_rossii",
            "plecoglossus_altivelis",
            "tinca_tinca",
        ]:
            self.assertIn(f"--dataset {dataset}", command)

    def test_paper_python_tasks_use_import_guard_before_building(self):
        pixi = tomllib.loads(Path("pixi.toml").read_text())
        tasks = pixi["feature"]["python"]["tasks"]

        self.assertIn("ensure-python", tasks)
        self.assertIn("import pyrsx", tasks["ensure-python"]["cmd"])
        self.assertEqual(tasks["run-literature-bindings"]["depends-on"], ["ensure-python"])
        self.assertEqual(tasks["analyze-literature-modes"]["depends-on"], ["ensure-python"])

    def test_binding_slurm_script_is_dataset_array_with_collatable_outputs(self):
        script = Path("benchmarks/slurm/literature_bindings_features.sbatch").read_text()

        self.assertIn("#SBATCH --array=0-3", script)
        self.assertIn("RSX_DATASETS:-danio_albolineatus notothenia_rossii plecoglossus_altivelis tinca_tinca", script)
        self.assertIn("SLURM_ARRAY_TASK_ID", script)
        self.assertIn("RSX_SLURM_RESULTS_DIR:-benchmarks/results/slurm", script)
        self.assertIn("literature_binding_results_${dataset}_${job_id}.csv", script)


if __name__ == "__main__":
    unittest.main()
