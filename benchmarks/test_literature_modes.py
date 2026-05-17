import tempfile
import unittest
from pathlib import Path

from benchmarks.analyze_literature_modes import (
    summarize_depth,
    summarize_distrib,
    summarize_freq,
    summarize_pca,
    summarize_signif,
)


class LiteratureModeSummaryTests(unittest.TestCase):
    def test_summarize_freq_counts_singletons(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "freq.tsv"
            path.write_text(
                "#source:rsx-freq;min_depth:10\n"
                "Frequency\tCount\n"
                "1\t25\n"
                "2\t15\n"
                "3\t10\n"
            )
            summary = summarize_freq(path)
        self.assertEqual(summary["tested_markers"], "50")
        self.assertEqual(summary["output_rows"], "3")
        self.assertEqual(summary["singleton_fraction"], "0.5")

    def test_summarize_depth_reports_samples(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "depth.tsv"
            path.write_text(
                "Sample\tGroup\tReads\tMarkers\tRetained\tMin_depth\tMax_depth\tMedian_depth\tAverage_depth\n"
                "s1\tmale\t100\t10\t8\t0\t20\t2\t4\n"
                "s2\tfemale\t300\t12\t10\t0\t25\t3\t5\n"
            )
            summary = summarize_depth(path)
        self.assertEqual(summary["output_rows"], "2")
        self.assertIn("female:1", summary["summary"])
        self.assertIn("male:1", summary["summary"])

    def test_summarize_distrib_sums_significant_marker_cells(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "distrib.tsv"
            path.write_text(
                "#source:rsx-distrib;min_depth:10;signif_threshold:0.001;bonferroni:true;n_markers:1000\n"
                "male\tfemale\tMarkers\tP\tCorrectedP\tSignif\tBias\n"
                "10\t0\t4\t0.0001\t0.1\tTrue\t1\n"
                "9\t1\t6\t0.01\t1\tFalse\t0.8\n"
                "0\t10\t3\t0.0001\t0.1\tTrue\t-1\n"
            )
            summary = summarize_distrib(path)
        self.assertEqual(summary["tested_markers"], "1000")
        self.assertEqual(summary["output_rows"], "3")
        self.assertEqual(summary["significant_markers"], "7")

    def test_summarize_signif_counts_bayesian_columns(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "signif.tsv"
            path.write_text(
                "#source:rsx-signif;min_depth:10;signif_threshold:0.05;correction:none;test:chisq;n_markers:500\n"
                "id\tsequence\tBayes_Factor\tPosterior_SexLinked\n"
                "1\tAAAA\t12.5\t0.95\n"
                "2\tCCCC\t8.0\t0.75\n"
                "3\tGGGG\t2.0\t0.05\n"
            )
            summary = summarize_signif(path)
        self.assertEqual(summary["tested_markers"], "500")
        self.assertEqual(summary["output_rows"], "3")
        self.assertEqual(summary["posterior_gt_0_5"], "2")
        self.assertEqual(summary["posterior_gt_0_9"], "1")
        self.assertEqual(summary["bayes_factor_gt_10"], "1")

    def test_summarize_pca_computes_sex_loading_delta(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pca_dir = root / "pca"
            pca_dir.mkdir()
            (pca_dir / "eigenvalues.tsv").write_text(
                "component\teigenvalue\tvariance_fraction\tcumulative\n"
                "PC1\t8\t0.8\t0.8\n"
                "PC2\t1\t0.1\t0.9\n"
            )
            (pca_dir / "loadings.tsv").write_text(
                "individual\tPC1\tPC2\n"
                "m1\t1.0\t0.2\n"
                "m2\t0.8\t0.4\n"
                "f1\t0.2\t1.0\n"
                "f2\t0.0\t0.8\n"
            )
            popmap = root / "popmap.tsv"
            popmap.write_text("m1\tmale\nm2\tmale\nf1\tfemale\nf2\tfemale\n")
            summary = summarize_pca(pca_dir, popmap, "male", "female")
        self.assertEqual(summary["pc1_variance_fraction"], "0.8")
        self.assertEqual(summary["pc2_variance_fraction"], "0.1")
        self.assertEqual(summary["sex_loading_delta_pc1"], "0.8")
        self.assertEqual(summary["sex_loading_delta_pc2"], "0.6")


if __name__ == "__main__":
    unittest.main()
