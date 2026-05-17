#!/usr/bin/env python3
"""Validate the literature-derived benchmark manifest."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


REQUIRED_COLUMNS = {
    "dataset",
    "source",
    "accession",
    "genome",
    "benchmark_role",
    "commands",
    "notes",
}

RADSEX_PANEL_DATASETS = {
    "cyprinus_carpio",
    "danio_aesculapii",
    "danio_albolineatus",
    "danio_choprae",
    "danio_kyathit",
    "gadus_morhua",
    "gymnocorymbus_ternetzi",
    "gymnotus_carapo",
    "hippocampus_abdominalis",
    "lepisosteus_oculatus",
    "notothenia_rossii",
    "plecoglossus_altivelis",
    "poecilia_sphenops",
    "sander_vitreus",
    "tinca_tinca",
}


def load_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise ValueError(f"{path} does not exist")
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = REQUIRED_COLUMNS.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} is missing columns: {', '.join(sorted(missing))}")
        return list(reader)


def validate(path: Path) -> list[str]:
    rows = load_rows(path)
    errors: list[str] = []
    by_dataset = {row["dataset"]: row for row in rows}

    if len(by_dataset) != len(rows):
        errors.append("dataset names must be unique")

    missing_panel = sorted(RADSEX_PANEL_DATASETS.difference(by_dataset))
    if missing_panel:
        errors.append(f"missing RADSex paper panel datasets: {', '.join(missing_panel)}")

    medaka = by_dataset.get("oryzias_latipes")
    if medaka is None:
        errors.append("missing RADSex medaka benchmark row: oryzias_latipes")
    elif medaka["accession"] != "PRJNA253959":
        errors.append("oryzias_latipes accession must be PRJNA253959")

    for dataset in RADSEX_PANEL_DATASETS:
        row = by_dataset.get(dataset)
        if row and row["accession"] != "PRJNA548074":
            errors.append(f"{dataset} accession must be PRJNA548074")

    for row in rows:
        commands = {command.strip() for command in row["commands"].split(",") if command.strip()}
        required = {"process", "distrib", "signif", "freq", "depth"}
        missing = sorted(required.difference(commands))
        if missing:
            errors.append(f"{row['dataset']} missing commands: {', '.join(missing)}")
        if not row["source"] or not row["benchmark_role"] or not row["notes"]:
            errors.append(f"{row['dataset']} has incomplete benchmark metadata")

    return errors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", default="benchmarks/literature_datasets.tsv", type=Path)
    args = parser.parse_args()

    errors = validate(args.manifest)
    if errors:
        for error in errors:
            print(error, file=sys.stderr)
        return 1

    print(f"{args.manifest} contains a complete literature benchmark manifest")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
