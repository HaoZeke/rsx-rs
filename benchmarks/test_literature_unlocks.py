import tempfile
import unittest
from pathlib import Path

from benchmarks.plot_literature_benchmarks import prepare_bio_unlock_plot_rows
from benchmarks.summarize_literature_unlocks import infer_biological_signal, summarize_dataset


class LiteratureUnlockSummaryTests(unittest.TestCase):
    def test_infer_biological_signal_calls_confirmatory_xy_from_strict_male_bias(self):
        inference = infer_biological_signal(
            "demo_species",
            10,
            {
                "strict_candidates": "2",
                "posterior_gt_0_9": "2",
                "bayes_factor_only": "4",
                "singleton_fraction": "0.1",
            },
            [
                {"strict_call": True, "posterior_sex_linked": 0.99, "bias_direction": "male-biased"},
                {"strict_call": True, "posterior_sex_linked": 0.98, "bias_direction": "male-biased"},
                {"strict_call": False, "posterior_sex_linked": 0.95, "bias_direction": "female-biased"},
            ],
        )

        self.assertEqual(inference["evidence_class"], "confirmatory")
        self.assertEqual(inference["inferred_sex_system"], "XX/XY-supported")
        self.assertEqual(inference["strict_male_biased"], "2")
        self.assertIn("male-biased strict markers", inference["biological_inference"])

    def test_infer_biological_signal_reports_exploratory_zw_from_posterior_only_bias(self):
        inference = infer_biological_signal(
            "demo_species",
            10,
            {
                "strict_candidates": "0",
                "posterior_gt_0_9": "3",
                "bayes_factor_only": "1",
                "singleton_fraction": "0.2",
            },
            [
                {"strict_call": False, "posterior_sex_linked": 0.99, "bias_direction": "female-biased"},
                {"strict_call": False, "posterior_sex_linked": 0.96, "bias_direction": "female-biased"},
                {"strict_call": False, "posterior_sex_linked": 0.93, "bias_direction": "male-biased"},
            ],
        )

        self.assertEqual(inference["evidence_class"], "exploratory")
        self.assertEqual(inference["inferred_sex_system"], "ZZ/ZW-like")
        self.assertEqual(inference["posterior_female_biased"], "2")
        self.assertIn("source call remains strict-null", inference["biological_inference"])

    def test_infer_biological_signal_restrains_bayes_factor_only_rows(self):
        inference = infer_biological_signal(
            "demo_species",
            10,
            {
                "strict_candidates": "0",
                "posterior_gt_0_9": "0",
                "bayes_factor_only": "5",
                "singleton_fraction": "0.55",
            },
            [
                {"strict_call": False, "posterior_sex_linked": 0.8, "bias_direction": "female-biased"},
                {"strict_call": False, "posterior_sex_linked": 0.7, "bias_direction": "male-biased"},
            ],
        )

        self.assertEqual(inference["evidence_class"], "restrained_null")
        self.assertEqual(inference["inferred_sex_system"], "no high-posterior sex-system call")
        self.assertIn("not converted into a sex-system call", inference["biological_inference"])

    def test_prepare_bio_unlock_plot_rows_keeps_interpretable_classes(self):
        rows = prepare_bio_unlock_plot_rows(
            [
                {
                    "dataset": "demo_species",
                    "strict_and_posterior": "4",
                    "strict_only": "1",
                    "posterior_only": "2",
                    "bayes_factor_only": "3",
                }
            ]
        )

        by_class = {row["candidate_class"]: row["marker_count"] for row in rows}
        self.assertEqual(by_class["Strict + posterior"], 4)
        self.assertEqual(by_class["Strict only"], 1)
        self.assertEqual(by_class["Posterior only"], 2)
        self.assertEqual(by_class["Bayes-factor only"], 3)

    def test_summarize_dataset_counts_overlap_and_bayes_factor_only_rows(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            workdir = root / "workdir"
            mode_dir = workdir / "modes"
            dataset_dir = workdir / "demo_species"
            output_dir = mode_dir / "demo_species"
            dataset_dir.mkdir(parents=True)
            output_dir.mkdir(parents=True)
            (dataset_dir / "popmap.tsv").write_text(
                "m1\tmale\nm2\tmale\nm3\tmale\nf1\tfemale\nf2\tfemale\nf3\tfemale\n"
            )
            (output_dir / "signif_chisq_bonferroni_d10.tsv").write_text(
                "#source:rsx-signif;min_depth:10;signif_threshold:0.05;correction:bonferroni;test:chisq;n_markers:100\n"
                "id\tsequence\tm1\tm2\tm3\tf1\tf2\tf3\n"
                "1\tAAAA\t10\t10\t10\t0\t0\t0\n"
                "2\tCCCC\t10\t0\t10\t0\t0\t0\n"
            )
            (output_dir / "signif_bayes_chisq_none_d10.tsv").write_text(
                "#source:rsx-signif;min_depth:10;signif_threshold:0.05;correction:none;test:chisq;n_markers:100\n"
                "id\tsequence\tm1\tm2\tm3\tf1\tf2\tf3\tBayes_Factor\tPosterior_SexLinked\n"
                "1\tAAAA\t10\t10\t10\t0\t0\t0\t20.0\t0.95\n"
                "3\tGGGG\t10\t10\t10\t0\t0\t0\t15.0\t0.92\n"
                "4\tTTTT\t10\t10\t0\t0\t0\t0\t30.0\t0.80\n"
                "5\tACAC\t10\t0\t0\t0\t0\t0\t2.0\t0.40\n"
            )
            mode_effects = {
                "tested_markers": "100",
                "singleton_fraction": "0.25",
                "pc1_variance_fraction": "0.8",
                "sex_loading_delta_pc1": "0.1",
            }

            summary, candidates = summarize_dataset(
                workdir=workdir,
                mode_dir=mode_dir,
                dataset="demo_species",
                min_depth=10,
                group1="male",
                group2="female",
                mode_effects=mode_effects,
                top_candidates=10,
            )

        self.assertEqual(summary["strict_candidates"], "2")
        self.assertEqual(summary["posterior_gt_0_9"], "2")
        self.assertEqual(summary["strict_and_posterior"], "1")
        self.assertEqual(summary["strict_only"], "1")
        self.assertEqual(summary["posterior_only"], "1")
        self.assertEqual(summary["bayes_factor_gt_10"], "3")
        self.assertEqual(summary["bayes_factor_only"], "1")
        self.assertEqual(summary["unlock_class"], "strict recovery plus posterior expansion")
        self.assertEqual(len(candidates), 5)

    def test_summarize_dataset_reports_penetrance_and_candidate_classes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            workdir = root / "workdir"
            mode_dir = workdir / "modes"
            dataset_dir = workdir / "demo_species"
            output_dir = mode_dir / "demo_species"
            dataset_dir.mkdir(parents=True)
            output_dir.mkdir(parents=True)
            (dataset_dir / "popmap.tsv").write_text("m1\tmale\nm2\tmale\nf1\tfemale\nf2\tfemale\n")
            (output_dir / "signif_chisq_bonferroni_d10.tsv").write_text(
                "#source:rsx-signif;min_depth:10;signif_threshold:0.05;correction:bonferroni;test:chisq;n_markers:50\n"
                "id\tsequence\tm1\tm2\tf1\tf2\n"
                "1\tAAAA\t10\t10\t0\t0\n"
            )
            (output_dir / "signif_bayes_chisq_none_d10.tsv").write_text(
                "#source:rsx-signif;min_depth:10;signif_threshold:0.05;correction:none;test:chisq;n_markers:50\n"
                "id\tsequence\tm1\tm2\tf1\tf2\tBayes_Factor\tPosterior_SexLinked\n"
                "1\tAAAA\t10\t10\t0\t0\t25.0\t0.96\n"
                "2\tCCCC\t0\t0\t10\t10\t12.0\t0.93\n"
                "3\tGGGG\t10\t0\t0\t0\t30.0\t0.70\n"
            )

            _, candidates = summarize_dataset(
                workdir=workdir,
                mode_dir=mode_dir,
                dataset="demo_species",
                min_depth=10,
                group1="male",
                group2="female",
                mode_effects={},
                top_candidates=10,
            )

        by_id = {candidate["id"]: candidate for candidate in candidates}
        self.assertEqual(by_id["1"]["candidate_class"], "strict+posterior")
        self.assertEqual(by_id["1"]["male_penetrance"], "1")
        self.assertEqual(by_id["1"]["female_penetrance"], "0")
        self.assertEqual(by_id["1"]["bias_direction"], "male-biased")
        self.assertEqual(by_id["2"]["candidate_class"], "posterior_only")
        self.assertEqual(by_id["2"]["bias_direction"], "female-biased")
        self.assertEqual(by_id["3"]["candidate_class"], "bayes_factor_only")


if __name__ == "__main__":
    unittest.main()
