"""Write the canonical Pydantic run-profile JSON Schema."""

from __future__ import annotations

import json
import importlib.util
import sys
from pathlib import Path


def load_run_profile(repository: Path):
    """Load the pure-Python schema without importing the native package module."""

    source = repository / "rsx-python" / "python" / "pyrsx" / "profile.py"
    spec = importlib.util.spec_from_file_location("pyrsx_profile_schema", source)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load profile schema module from {source}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module.RunProfile


def main() -> None:
    repository = Path(__file__).parents[1]
    output = repository / "rsx-python" / "schema"
    output.mkdir(parents=True, exist_ok=True)
    RunProfile = load_run_profile(repository)
    schema = RunProfile.model_json_schema()
    destination = output / "run-profile-v1.schema.json"
    destination.write_text(
        json.dumps(schema, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
