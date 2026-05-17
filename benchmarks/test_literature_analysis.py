import tempfile
import unittest
import math
from pathlib import Path

from benchmarks.analyze_bayesian_evidence import (
    analyze_distrib,
    bayes_factor_2x2,
    posterior_sex_linked,
)
from benchmarks.plot_literature_benchmarks import (
    FALLBACK_RUHI_COLORS,
    candidate_recovery_rows,
    load_ruhi_colors,
    prepare_compute_phase_rows,
    summarize_dataset_rows,
)


class BayesianEvidenceTests(unittest.TestCase):
    def test_posterior_is_symmetric_for_either_enriched_group(self):
        group1 = posterior_sex_linked(10, 0, 10, 10, 0.01, 0.9)
        group2 = posterior_sex_linked(0, 10, 10, 10, 0.01, 0.9)

        self.assertGreater(group1, 0.9)
        self.assertGreater(group2, 0.9)
        self.assertAlmostEqual(group1, group2, places=12)

    def test_bayes_factor_elevates_extreme_distribution(self):
        self.assertGreater(bayes_factor_2x2(10, 0, 10, 10), 10.0)

    def test_analyze_distrib_weights_aggregate_marker_cells(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            popmap = root / "popmap.tsv"
            popmap.write_text("\n".join([f"m{i}\tmale" for i in range(10)] + [f"f{i}\tfemale" for i in range(10)]) + "\n")
            distrib = root / "distrib_1.tsv"
            distrib.write_text(
                "#source:rsx-distrib\n"
                "male\tfemale\tMarkers\tP\tCorrectedP\tSignif\tBias\n"
                "10\t0\t3\t0\t0\tTrue\t1\n"
                "5\t5\t7\t1\t1\tFalse\t0\n"
            )

            row = analyze_distrib("toy", 1, distrib, popmap, 0.01, 0.9)

            self.assertEqual(row["markers"], "10")
            self.assertEqual(row["markers_posterior_gt_0_9"], "3")
            self.assertEqual(row["top_cell_group1"], "10")
            self.assertEqual(row["top_cell_group2"], "0")


class LiteraturePlotTests(unittest.TestCase):
    def test_summarize_dataset_rows_tracks_process_and_analysis(self):
        rows = [
            {"dataset": "toy_species", "command": "metadata", "samples": "2", "total_spots": "10", "total_bases": "1000"},
            {"dataset": "toy_species", "command": "download", "elapsed_seconds": "3"},
            {"dataset": "toy_species", "command": "process", "elapsed_seconds": "2", "markers": "100"},
            {"dataset": "toy_species", "command": "freq", "elapsed_seconds": "0.5"},
            {"dataset": "toy_species", "command": "signif", "elapsed_seconds": "1.5", "significant_markers": "4"},
        ]

        summary = summarize_dataset_rows(rows)

        self.assertEqual(summary[0]["dataset"], "toy_species")
        self.assertEqual(summary[0]["markers"], "100")
        self.assertEqual(summary[0]["analysis_seconds"], "2.000")
        self.assertEqual(summary[0]["total_seconds"], "7.000")
        self.assertEqual(summary[0]["compute_seconds"], "4.000")
        self.assertEqual(summary[0]["markers_per_second"], "50.000")
        self.assertEqual(summary[0]["significant_markers"], "4")

    def test_compute_phase_rows_exclude_downloads(self):
        rows = [
            {"dataset": "toy_species", "dataset_label": "Toy Species", "command": "download", "elapsed_seconds": 30.0},
            {"dataset": "toy_species", "dataset_label": "Toy Species", "command": "process", "elapsed_seconds": 2.0},
            {"dataset": "toy_species", "dataset_label": "Toy Species", "command": "freq", "elapsed_seconds": 0.5},
            {"dataset": "toy_species", "dataset_label": "Toy Species", "command": "signif", "elapsed_seconds": 1.5},
        ]

        phases = prepare_compute_phase_rows(rows)

        self.assertEqual(
            phases,
            [
                {"dataset": "toy_species", "dataset_label": "Toy Species", "phase": "Process", "elapsed_seconds": 2.0},
                {"dataset": "toy_species", "dataset_label": "Toy Species", "phase": "Downstream analysis", "elapsed_seconds": 2.0},
            ],
        )

    def test_candidate_recovery_rows_compare_strict_and_posterior_counts(self):
        evidence = [
            {"dataset": "toy_species", "min_depth": "1", "markers_posterior_gt_0_9": "9", "markers_posterior_gt_0_5": "12"},
        ]
        benchmark = [
            {"dataset": "toy_species", "command": "signif", "min_depth": "1", "significant_markers": "0"},
        ]

        rows = candidate_recovery_rows(evidence, benchmark)

        self.assertEqual(rows[0]["metric"], "Strict Bonferroni FASTA")
        self.assertEqual(rows[0]["marker_count"], 0)
        self.assertEqual(rows[1]["metric"], "Posterior > 0.9")
        self.assertEqual(rows[1]["marker_count"], 9)
        self.assertAlmostEqual(rows[1]["log10_marker_count_plus_one"], math.log10(10.0))

    def test_ruhi_palette_falls_back_to_chemparseplot_colors(self):
        colors = load_ruhi_colors()

        self.assertEqual(colors["teal"], FALLBACK_RUHI_COLORS["teal"])
        self.assertEqual(colors["coral"], FALLBACK_RUHI_COLORS["coral"])


if __name__ == "__main__":
    unittest.main()
