import math
import tempfile
import unittest
from pathlib import Path

from benchmarks.analyze_bayesian_evidence import (
    PrevalencePrior,
    posterior_sex_linked,
    posterior_sex_linked_with_model,
)
from benchmarks.analyze_posterior_families import analyze_profiles, read_profiles


class PosteriorFamilyTests(unittest.TestCase):
    def test_explicit_fixed_model_matches_compatibility_function(self):
        expected = posterior_sex_linked(10, 0, 10, 10, 0.01, 0.9)

        observed = posterior_sex_linked_with_model(
            10,
            0,
            10,
            10,
            linkage_prior=0.01,
            group1_linked_weight=0.5,
            linked=PrevalencePrior.fixed(0.9),
            null=PrevalencePrior.fixed(0.5),
        )

        self.assertAlmostEqual(observed, expected, places=12)

    def test_beta_model_uses_integrated_prevalence_evidence(self):
        observed = posterior_sex_linked_with_model(
            10,
            0,
            10,
            10,
            linkage_prior=0.01,
            group1_linked_weight=0.5,
            linked=PrevalencePrior.beta(9.0, 1.0),
            null=PrevalencePrior.beta(5.0, 5.0),
        )

        def log_beta(alpha, beta):
            return math.lgamma(alpha) + math.lgamma(beta) - math.lgamma(alpha + beta)

        linked_log_marginal = log_beta(29.0, 1.0) - log_beta(9.0, 1.0)
        null_log_marginal = log_beta(15.0, 15.0) - log_beta(5.0, 5.0)
        log_odds = linked_log_marginal - null_log_marginal + math.log(0.01 / 0.99)
        expected = 1.0 / (1.0 + math.exp(-log_odds))

        self.assertAlmostEqual(observed, expected, places=12)

    def test_toml_profiles_drive_weighted_real_table_summary(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            workdir = root / "workdir"
            dataset = workdir / "published_panel"
            dataset.mkdir(parents=True)
            (dataset / "popmap.tsv").write_text(
                "".join(f"m{i}\tmale\n" for i in range(10))
                + "".join(f"f{i}\tfemale\n" for i in range(10))
            )
            (dataset / "distrib_10.tsv").write_text(
                "male\tfemale\tMarkers\tP\tCorrectedP\tSignif\tBias\n"
                "10\t0\t3\t0\t0\tTrue\t1\n"
                "5\t5\t7\t1\t1\tFalse\t0\n"
            )
            profiles_path = root / "profiles.toml"
            profiles_path.write_text(
                "threshold = 0.9\n"
                "[[profile]]\n"
                "name = 'fixed-default'\n"
                "linkage_prior = 0.01\n"
                "group1_linked_weight = 0.5\n"
                "linked = { family = 'fixed', probability = 0.9 }\n"
                "null = { family = 'fixed', probability = 0.5 }\n"
                "[[profile]]\n"
                "name = 'beta-informed'\n"
                "linkage_prior = 0.01\n"
                "group1_linked_weight = 0.5\n"
                "linked = { family = 'beta', alpha = 9.0, beta = 1.0 }\n"
                "null = { family = 'beta', alpha = 5.0, beta = 5.0 }\n"
            )

            threshold, profiles = read_profiles(profiles_path)
            rows = analyze_profiles(workdir, ["published_panel"], [10], threshold, profiles)

        self.assertEqual([row["profile"] for row in rows], ["fixed-default", "beta-informed"])
        self.assertEqual(rows[0]["markers"], 10)
        self.assertEqual(rows[0]["markers_posterior_gt_threshold"], 3)
        self.assertEqual(rows[0]["threshold"], 0.9)
        self.assertEqual(rows[1]["linked_family"], "beta")
        self.assertEqual(rows[1]["source_table"], "published_panel/distrib_10.tsv")


if __name__ == "__main__":
    unittest.main()
