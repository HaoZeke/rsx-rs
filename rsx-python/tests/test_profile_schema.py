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
                },
            },
        }
    )
    assert profile.run.bayes_model.null_prevalence == 0.4
    assert profile.run.bayes_model.group1_linked_weight == 0.75


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
        ),
        lambda: pyrsx.triage(
            "missing.tsv",
            "missing-popmap.tsv",
            "unused.tsv",
            prior_probability=0.02,
            linked_probability=0.85,
            null_prevalence=0.4,
            group1_linked_weight=0.7,
        ),
    ],
)
def test_low_level_bindings_accept_complete_bayesian_model(call) -> None:
    with pytest.raises(RuntimeError):
        call()
