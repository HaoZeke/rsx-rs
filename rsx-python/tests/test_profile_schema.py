"""Contracts for the versioned Python run-profile schema."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
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


def test_checked_in_json_schema_matches_pydantic_model() -> None:
    schema_path = Path(__file__).parents[1] / "schema" / "run-profile-v1.schema.json"
    checked_in = json.loads(schema_path.read_text(encoding="utf-8"))
    generated = RunProfile.model_json_schema()

    assert checked_in == generated
    assert len(generated["$defs"]["RunCommand"]["oneOf"]) == 10
