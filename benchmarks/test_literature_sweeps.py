import os
import tempfile
import unittest
from pathlib import Path

from benchmarks.collect_lowdepth_sweep import collect_shards
from benchmarks.plot_prior_sensitivity_heatmap import collect as collect_prior_sensitivity


TRIAGE_HEADER = "id\tPosterior_SexLinked\tBayes_Factor\tStrict_Call\n"


class LiteratureSweepCollectionTests(unittest.TestCase):
    def test_prior_sensitivity_prefers_newest_job_suffixed_cell(self):
        with tempfile.TemporaryDirectory() as tmp:
            slurm_dir = Path(tmp)
            canonical = slurm_dir / "triage_demo_species_pi0.01_psex0.9.tsv"
            suffixed = slurm_dir / "triage_demo_species_pi0.01_psex0.9_12345.tsv"
            canonical.write_text(TRIAGE_HEADER + "a\t0.95\t11\tfalse\n")
            suffixed.write_text(
                TRIAGE_HEADER
                + "a\t0.95\t11\tfalse\n"
                + "b\t0.91\t12\ttrue\n"
            )
            os.utime(canonical, (1, 1))
            os.utime(suffixed, (2, 2))

            summary = collect_prior_sensitivity(slurm_dir)

        self.assertEqual(len(summary), 1)
        row = summary.iloc[0]
        self.assertEqual(row["dataset"], "demo_species")
        self.assertEqual(row["n_posterior_gt_0_9"], 2)
        self.assertEqual(row["n_strict_call"], 1)

    def test_lowdepth_shards_are_combined_in_dataset_depth_mode_order(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "lowdepth_demo_species_d5.csv").write_text(
                "dataset,min_depth,mode,api_call\n"
                "demo_species,5,beta,call_b\n"
            )
            (root / "lowdepth_demo_species_d3.csv").write_text(
                "dataset,min_depth,mode,api_call\n"
                "demo_species,3,alpha,call_a\n"
            )

            combined = collect_shards(root)

        self.assertEqual(combined["min_depth"].tolist(), ["3", "5"])
        self.assertEqual(combined["mode"].tolist(), ["alpha", "beta"])

    def test_lowdepth_slurm_builds_python_extension_in_paper_environment(self):
        script = Path("benchmarks/slurm/literature_biology_low_depth.sbatch").read_text()

        self.assertIn('run -e paper build-python', script)
        self.assertIn('run -e paper python benchmarks/analyze_literature_modes.py', script)
        self.assertNotIn('|| true', script)


if __name__ == "__main__":
    unittest.main()
