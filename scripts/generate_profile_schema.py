"""Write the canonical Pydantic run-profile JSON Schema."""

from __future__ import annotations

import json
from pathlib import Path

from pyrsx.profile import RunProfile


def main() -> None:
    output = Path(__file__).parents[1] / "rsx-python" / "schema"
    output.mkdir(parents=True, exist_ok=True)
    schema = RunProfile.model_json_schema()
    destination = output / "run-profile-v1.schema.json"
    destination.write_text(
        json.dumps(schema, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
