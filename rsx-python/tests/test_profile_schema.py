"""Contracts for the versioned Python run-profile schema."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
import pyrsx
from pydantic import ValidationError

from pyrsx.profile import RunProfile, parse_run_profile_toml


PROCESS_TOML = """
schema_version = 1
profile_name = "python-process-v1"
reproducibility_archive = "run-repro.zip"
write_hydrated_profile = "run.hydrated.toml"

[run]
command = "process"
input_dir = "reads"
output_file = "markers.tsv"
threads = 4
min_depth = 2
kmer_dedup = 31
"""


def test_process_profile_is_strict_and_command_discriminated() -> None:
    profile = parse_run_profile_toml(PROCESS_TOML)
    assert isinstance(profile, RunProfile)
    assert profile.run.command == "process"
    assert profile.run.threads == 4

    invalid = profile.model_dump()
    invalid["run"]["hidden_default"] = 7
    with pytest.raises(ValidationError, match="hidden_default"):
        RunProfile.model_validate(invalid)


def test_bayesian_commands_require_the_complete_directional_model() -> None:
    profile = RunProfile.model_validate(
        {
            "schema_version": 1,
            "profile_name": "python-distrib-v1",
            "run": {
                "command": "distrib",
                "markers_table": "markers.tsv",
                "popmap": "popmap.tsv",
                "output_file": "distrib.tsv",
                "min_depth": 1,
                "groups": ["M", "F"],
                "signif_threshold": 0.05,
                "disable_correction": False,
                "correction": "bonferroni",
                "test_method": "chisq",
                "output_bayes": True,
                "bayes_model": {
                    "linkage_prior": 0.01,
                    "linked_prevalence": 0.9,
                    "null_prevalence": 0.4,
                    "group1_linked_weight": 0.75,
                    "posterior": {
                        "linked": {
                            "family": "beta",
                            "alpha": 9.0,
                            "beta": 1.0,
                        },
                        "null": {
                            "family": "beta",
                            "alpha": 5.0,
                            "beta": 5.0,
                        },
                    },
                    "bayes_factor": {
                        "alternative_group1": {"alpha": 8.0, "beta": 2.0},
                        "alternative_group2": {"alpha": 2.0, "beta": 8.0},
                        "null": {"alpha": 10.0, "beta": 10.0},
                    },
                },
            },
        }
    )
    assert profile.run.bayes_model.null_prevalence == 0.4
    assert profile.run.bayes_model.group1_linked_weight == 0.75
    assert profile.run.bayes_model.posterior.linked.family == "beta"
    assert profile.run.bayes_model.posterior.linked.alpha == 9.0
    assert profile.run.bayes_model.bayes_factor.alternative_group1.alpha == 8.0

    invalid = profile.model_dump()
    invalid["run"]["bayes_model"]["bayes_factor"]["null"]["alpha"] = 0.0
    with pytest.raises(ValidationError, match="greater than 0"):
        RunProfile.model_validate(invalid)


def test_legacy_directional_model_hydrates_uniform_beta_priors() -> None:
    profile = RunProfile.model_validate(
        {
            "schema_version": 1,
            "profile_name": "python-legacy-v1",
            "run": {
                "command": "triage",
                "markers_table": "markers.tsv",
                "popmap": "popmap.tsv",
                "output_file": "triage.tsv",
                "min_depth": 1,
                "groups": ["M", "F"],
                "signif_threshold": 0.05,
                "posterior_threshold": 0.9,
                "bayes_factor_threshold": 10.0,
                "bayes_model": {
                    "linkage_prior": 0.01,
                    "linked_prevalence": 0.9,
                    "null_prevalence": 0.5,
                    "group1_linked_weight": 0.5,
                },
            },
        }
    )

    priors = profile.run.bayes_model.bayes_factor
    posterior = profile.run.bayes_model.posterior
    assert posterior.linked.family == "fixed"
    assert posterior.linked.probability == 0.9
    assert posterior.null.probability == 0.5
    assert priors.alternative_group1.alpha == 1.0
    assert priors.alternative_group2.beta == 1.0
    assert priors.null.alpha == 1.0
    assert "bayes_factor" in profile.model_dump()["run"]["bayes_model"]


def test_depth_profile_and_binding_expose_the_streaming_policy() -> None:
    profile = RunProfile.model_validate(
        {
            "schema_version": 1,
            "profile_name": "python-depth-v1",
            "run": {
                "command": "depth",
                "markers_table": "markers.tsv",
                "popmap": "popmap.tsv",
                "output_file": "depth.tsv",
                "min_frequency": 0.75,
                "streaming_mode": "memory",
                "streaming_threshold_bytes": 123,
            },
        }
    )
    assert profile.run.streaming_mode == "memory"
    assert profile.run.streaming_threshold_bytes == 123

    with pytest.raises(RuntimeError):
        pyrsx.depth(
            "missing.tsv",
            "missing-popmap.tsv",
            "unused.tsv",
            streaming=False,
            streaming_threshold_bytes=123,
        )


def test_checked_in_json_schema_matches_pydantic_model() -> None:
    schema_path = Path(__file__).parents[1] / "schema" / "run-profile-v1.schema.json"
    checked_in = json.loads(schema_path.read_text(encoding="utf-8"))
    generated = RunProfile.model_json_schema()

    assert checked_in == generated
    assert len(generated["$defs"]["RunCommand"]["oneOf"]) == 10


@pytest.mark.parametrize(
    "call",
    [
        lambda: pyrsx.distrib(
            "missing.tsv",
            "missing-popmap.tsv",
            "unused.tsv",
            bayes=True,
            prior_probability=0.02,
            linked_probability=0.85,
            null_prevalence=0.4,
            group1_linked_weight=0.7,
            bf_group1_alpha=8.0,
            bf_group1_beta=2.0,
            bf_group2_alpha=2.0,
            bf_group2_beta=8.0,
            bf_null_alpha=10.0,
            bf_null_beta=10.0,
            posterior_linked_family="beta",
            posterior_linked_alpha=9.0,
            posterior_linked_beta=1.0,
            posterior_null_family="beta",
            posterior_null_alpha=5.0,
            posterior_null_beta=5.0,
        ),
        lambda: pyrsx.signif(
            "missing.tsv",
            "missing-popmap.tsv",
            "unused.tsv",
            bayes=True,
            prior_probability=0.02,
            linked_probability=0.85,
            null_prevalence=0.4,
            group1_linked_weight=0.7,
            bf_group1_alpha=8.0,
            bf_group1_beta=2.0,
            bf_group2_alpha=2.0,
            bf_group2_beta=8.0,
            bf_null_alpha=10.0,
            bf_null_beta=10.0,
            posterior_linked_family="beta",
            posterior_linked_alpha=9.0,
            posterior_linked_beta=1.0,
            posterior_null_family="beta",
            posterior_null_alpha=5.0,
            posterior_null_beta=5.0,
        ),
        lambda: pyrsx.triage(
            "missing.tsv",
            "missing-popmap.tsv",
            "unused.tsv",
            prior_probability=0.02,
            linked_probability=0.85,
            null_prevalence=0.4,
            group1_linked_weight=0.7,
            bf_group1_alpha=8.0,
            bf_group1_beta=2.0,
            bf_group2_alpha=2.0,
            bf_group2_beta=8.0,
            bf_null_alpha=10.0,
            bf_null_beta=10.0,
            posterior_linked_family="beta",
            posterior_linked_alpha=9.0,
            posterior_linked_beta=1.0,
            posterior_null_family="beta",
            posterior_null_alpha=5.0,
            posterior_null_beta=5.0,
        ),
    ],
)
def test_low_level_bindings_accept_complete_bayesian_model(call) -> None:
    with pytest.raises(RuntimeError):
        call()


@pytest.mark.parametrize(
    "call",
    [
        lambda: pyrsx.triage_to_arrow(
            "missing.tsv",
            "missing-popmap.tsv",
            signif_threshold=0.123,
            bayes_factor_threshold=456.0,
        ),
        lambda: pyrsx.triage_to_arrow_from_arrow(
            b"not-arrow",
            b"not-arrow",
            signif_threshold=0.123,
            bayes_factor_threshold=456.0,
        ),
    ],
)
def test_arrow_bindings_accept_all_triage_decision_thresholds(call) -> None:
    with pytest.raises(RuntimeError):
        call()
