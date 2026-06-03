#!/usr/bin/env python3
"""Download and run literature-derived RADSex benchmark datasets."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import os
import shutil
import subprocess
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
HTTP_HEADERS = {"User-Agent": "rsx-rs-literature-benchmarks/0.1"}
RESULT_COLUMNS = [
    "dataset",
    "source",
    "accession",
    "command",
    "min_depth",
    "elapsed_seconds",
    "samples",
    "total_spots",
    "total_bases",
    "total_sra_bytes",
    "markers",
    "rows",
    "significant_markers",
    "output_bytes",
    "output_path",
]


@dataclass(frozen=True)
class Dataset:
    name: str
    source: str
    accession: str
    genome: str


@dataclass(frozen=True)
class Sample:
    name: str
    accession: str
    sex: str
    spots: int
    bases: int
    bytes: int


@dataclass(frozen=True)
class FastqRemote:
    url: str
    bytes: int
    md5: str


def read_manifest(path: Path) -> list[Dataset]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"dataset", "source", "accession", "genome"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} is missing required columns: {', '.join(sorted(missing))}")
        return [
            Dataset(
                name=row["dataset"],
                source=row["source"],
                accession=row["accession"],
                genome=row["genome"],
            )
            for row in reader
        ]


def ncbi_url(endpoint: str, params: dict[str, str]) -> str:
    return f"{NCBI_BASE}/{endpoint}.fcgi?{urllib.parse.urlencode(params)}"


def fetch_text(url: str, retries: int = 5, sleep_seconds: float = 2.0) -> str:
    request = urllib.request.Request(url, headers=HTTP_HEADERS)
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                return response.read().decode("utf-8")
        except Exception as error:  # noqa: BLE001 - network failures vary by urllib backend
            last_error = error
            if attempt == retries:
                break
            time.sleep(sleep_seconds)
    raise RuntimeError(f"failed to fetch {url}: {last_error}") from last_error


def query_sra_xml(dataset: Dataset, retries: int = 5) -> str:
    organism = dataset.name.replace("_", " ")
    term = f"{dataset.accession}[Bioproject] {organism}[Organism]"
    esearch = fetch_text(
        ncbi_url(
            "esearch",
            {
                "db": "sra",
                "term": term,
                "usehistory": "y",
                "retmax": "0",
                "tool": "rsx-rs-literature-benchmarks",
            },
        ),
        retries=retries,
    )
    root = ET.fromstring(esearch)
    count = int(root.findtext("Count") or "0")
    if count == 0:
        raise RuntimeError(f"NCBI SRA returned no runs for {dataset.name} ({dataset.accession})")
    query_key = root.findtext("QueryKey")
    webenv = root.findtext("WebEnv")
    if not query_key or not webenv:
        raise RuntimeError(f"NCBI SRA did not return query history for {dataset.name}")
    return fetch_text(
        ncbi_url(
            "efetch",
            {
                "db": "sra",
                "query_key": query_key,
                "WebEnv": webenv,
                "tool": "rsx-rs-literature-benchmarks",
            },
        ),
        retries=retries,
    )


def sample_sex(package: ET.Element) -> str:
    for attr in package.findall(".//SAMPLE_ATTRIBUTE"):
        tag = (attr.findtext("TAG") or "").strip().lower()
        if tag == "sex":
            return (attr.findtext("VALUE") or "unknown").strip()
    return "unknown"


def parse_int(value: str | None) -> int:
    if not value:
        return 0
    try:
        return int(value)
    except ValueError:
        return 0


def sample_name(package: ET.Element, run: ET.Element) -> str:
    experiment = package.find("EXPERIMENT")
    alias = experiment.attrib.get("alias") if experiment is not None else ""
    if not alias:
        library_name = package.findtext(".//LIBRARY_NAME") or ""
        alias = library_name.strip()
    if not alias:
        alias = run.attrib["accession"]
    name = alias.removeprefix("RAD_").replace("/", "_").replace(" ", "_")
    return name


def parse_sra_experiment_xml(xml_text: str) -> list[Sample]:
    root = ET.fromstring(xml_text)
    samples: list[Sample] = []
    seen: set[str] = set()
    for package in root.findall("EXPERIMENT_PACKAGE"):
        sex = sample_sex(package)
        for run in package.findall(".//RUN_SET/RUN"):
            accession = run.attrib.get("accession")
            if not accession:
                continue
            name = sample_name(package, run)
            if name in seen:
                name = f"{name}_{accession}"
            seen.add(name)
            samples.append(
                Sample(
                    name=name,
                    accession=accession,
                    sex=sex,
                    spots=parse_int(run.attrib.get("total_spots") or run.attrib.get("spots")),
                    bases=parse_int(run.attrib.get("total_bases") or run.attrib.get("bases")),
                    bytes=parse_int(run.attrib.get("size")),
                )
            )
    samples.sort(key=lambda sample: sample.name)
    if not samples:
        raise ValueError("SRA XML did not contain any RUN accessions")
    return samples


def write_dataset_metadata(dataset_dir: Path, samples: list[Sample]) -> None:
    dataset_dir.mkdir(parents=True, exist_ok=True)
    download_dir = dataset_dir / ".download"
    download_dir.mkdir(parents=True, exist_ok=True)
    with (dataset_dir / "popmap.tsv").open("w") as handle:
        for sample in samples:
            handle.write(f"{sample.name}\t{sample.sex}\n")
    with (dataset_dir / "samples.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample", "accession", "sex", "spots", "bases", "bytes"])
        for sample in samples:
            writer.writerow(
                [sample.name, sample.accession, sample.sex, sample.spots, sample.bases, sample.bytes]
            )
    for sample in samples:
        (download_dir / f"{sample.name}.accession").write_text(sample.accession)


def ena_report_url(accessions: Iterable[str]) -> str:
    return (
        "https://www.ebi.ac.uk/ena/portal/api/filereport?"
        + urllib.parse.urlencode(
            {
                "accession": ",".join(accessions),
                "result": "read_run",
                "fields": "run_accession,fastq_ftp,fastq_bytes,fastq_md5",
                "format": "tsv",
            }
        )
    )


def normalize_fastq_url(value: str) -> str:
    if value.startswith(("http://", "https://", "ftp://", "file://")):
        return value
    return f"https://{value}"


def parse_ena_fastq_report(report: str) -> dict[str, list[FastqRemote]]:
    rows: dict[str, list[FastqRemote]] = {}
    reader = csv.DictReader(report.splitlines(), delimiter="\t")
    for row in reader:
        accession = row.get("run_accession") or ""
        fastq_ftp = row.get("fastq_ftp") or ""
        if not accession or not fastq_ftp:
            continue
        urls = fastq_ftp.split(";")
        sizes = (row.get("fastq_bytes") or "").split(";")
        md5s = (row.get("fastq_md5") or "").split(";")
        rows[accession] = [
            FastqRemote(
                url=normalize_fastq_url(url),
                bytes=parse_int(sizes[index] if index < len(sizes) else ""),
                md5=md5s[index] if index < len(md5s) else "",
            )
            for index, url in enumerate(urls)
            if url
        ]
    return rows


def fetch_ena_fastq_report(
    samples: list[Sample],
    retries: int = 5,
    fetcher=fetch_text,
) -> dict[str, list[FastqRemote]]:
    fastqs: dict[str, list[FastqRemote]] = {}
    for sample in samples:
        report = fetcher(ena_report_url([sample.accession]), retries=retries)
        fastqs.update(parse_ena_fastq_report(report))
    return fastqs


def read_dataset_metadata(dataset_dir: Path) -> list[Sample]:
    samples_path = dataset_dir / "samples.tsv"
    if not samples_path.exists():
        return []
    with samples_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [
            Sample(
                name=row["sample"],
                accession=row["accession"],
                sex=row["sex"],
                spots=parse_int(row.get("spots")),
                bases=parse_int(row.get("bases")),
                bytes=parse_int(row.get("bytes")),
            )
            for row in reader
        ]


def total_file_size(paths: Iterable[Path]) -> int:
    return sum(path.stat().st_size for path in paths if path.exists())


def relpath(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(Path.cwd().resolve()))
    except ValueError:
        return str(path)


def result_base(dataset: Dataset, samples: list[Sample]) -> dict[str, str]:
    return {
        "dataset": dataset.name,
        "source": dataset.source,
        "accession": dataset.accession,
        "samples": str(len(samples)),
        "total_spots": str(sum(sample.spots for sample in samples)),
        "total_bases": str(sum(sample.bases for sample in samples)),
        "total_sra_bytes": str(sum(sample.bytes for sample in samples)),
    }


class ResultWriter:
    def __init__(self, path: Path) -> None:
        self.path = path
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.handle = self.path.open("a", newline="")
        self.writer = csv.DictWriter(self.handle, fieldnames=RESULT_COLUMNS)
        if self.path.stat().st_size == 0:
            self.writer.writeheader()
            self.handle.flush()

    def close(self) -> None:
        self.handle.close()

    def write(self, row: dict[str, str]) -> None:
        filled = {column: row.get(column, "") for column in RESULT_COLUMNS}
        self.writer.writerow(filled)
        self.handle.flush()


def prune_dataset_results(path: Path, dataset_names: set[str]) -> None:
    if not path.exists() or not dataset_names:
        return
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        rows = [row for row in reader if row.get("dataset") not in dataset_names]
        fieldnames = reader.fieldnames or RESULT_COLUMNS
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    tmp.replace(path)


def marker_count_from_table(path: Path) -> int:
    with path.open() as handle:
        first = handle.readline().strip()
    if first.startswith("#Number of markers"):
        return int(first.rsplit(":", maxsplit=1)[1].strip())
    return count_data_rows(path)


def count_data_rows(path: Path) -> int:
    rows = 0
    with open_text(path) as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            first = stripped.split("\t", maxsplit=1)[0]
            if first in {
                "id",
                "Frequency",
                "Sample",
                "male",
                "female",
                "M",
                "F",
                "Chromosome",
            }:
                continue
            rows += 1
    return rows


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def parse_source_markers(path: Path) -> str:
    with open_text(path) as handle:
        for line in handle:
            if not line.startswith("#source:"):
                continue
            for field in line.strip().split(";"):
                if field.startswith("n_markers:"):
                    return field.split(":", maxsplit=1)[1]
    return ""


def count_fasta_records(path: Path) -> int:
    with open_text(path) as handle:
        return sum(1 for line in handle if line.startswith(">"))


def summarize_output(command: str, path: Path) -> dict[str, str]:
    rows = count_data_rows(path) if path.exists() else 0
    summary = {
        "rows": str(rows),
        "output_bytes": str(path.stat().st_size if path.exists() else 0),
        "output_path": relpath(path),
    }
    if command == "process":
        summary["markers"] = str(marker_count_from_table(path))
    elif command in {"distrib", "signif", "map", "subset"}:
        markers = parse_source_markers(path)
        if markers:
            summary["markers"] = markers
    if command == "signif" and path.suffix in {".fa", ".fasta"}:
        summary["significant_markers"] = str(count_fasta_records(path))
    elif command == "signif":
        summary["significant_markers"] = str(rows)
    if command == "subset":
        summary["significant_markers"] = str(rows)
    return summary


def run_command(args: list[str], log_path: Path) -> float:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    start = time.perf_counter()
    with log_path.open("w") as log:
        subprocess.run(args, stdout=log, stderr=subprocess.STDOUT, check=True)
    return time.perf_counter() - start


def stream_url_to_handle(remote: FastqRemote, output) -> None:
    digest = hashlib.md5()
    total = 0
    request = urllib.request.Request(remote.url, headers=HTTP_HEADERS)
    with urllib.request.urlopen(request, timeout=1200) as response:
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
            output.write(chunk)
            total += len(chunk)
    if remote.bytes and total != remote.bytes:
        raise RuntimeError(f"{remote.url} expected {remote.bytes} bytes but downloaded {total}")
    if remote.md5 and digest.hexdigest() != remote.md5:
        raise RuntimeError(f"{remote.url} failed MD5 validation")


def download_fastq_urls(files: list[FastqRemote], output: Path) -> None:
    if not files:
        raise RuntimeError(f"no FASTQ URLs available for {output.name}")
    tmp = output.with_suffix(output.suffix + ".tmp")
    with tmp.open("wb") as handle:
        for remote in files:
            stream_url_to_handle(remote, handle)
    tmp.replace(output)


def download_samples(
    dataset_dir: Path,
    samples: list[Sample],
    fastq_dump: str,
    logs_dir: Path,
    force: bool,
    method: str,
) -> float:
    samples_dir = dataset_dir / "samples"
    samples_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)
    ena_files: dict[str, list[FastqRemote]] = {}
    if method == "ena":
        ena_files = fetch_ena_fastq_report(samples)
    start = time.perf_counter()
    for sample in samples:
        output = samples_dir / f"{sample.name}.fq.gz"
        if output.exists() and output.stat().st_size > 0 and not force:
            continue
        if method == "ena":
            download_fastq_urls(ena_files.get(sample.accession, []), output)
            continue
        tmp = output.with_suffix(output.suffix + ".tmp")
        log_path = logs_dir / f"download_{sample.name}.log"
        with tmp.open("wb") as out, log_path.open("w") as log:
            subprocess.run(
                [fastq_dump, "-Z", "--gzip", sample.accession],
                stdout=out,
                stderr=log,
                check=True,
            )
        tmp.replace(output)
    return time.perf_counter() - start


def command_output(dataset_dir: Path, command: str, min_depth: int | None = None) -> Path:
    if command == "process":
        return dataset_dir / "markers_table.tsv"
    if command == "depth":
        return dataset_dir / "depth.tsv"
    suffix = ".fa" if command == "signif" else ".tsv"
    if min_depth is None:
        return dataset_dir / f"{command}{suffix}"
    return dataset_dir / f"{command}_{min_depth}{suffix}"


def run_rsx_dataset(
    dataset: Dataset,
    dataset_dir: Path,
    samples: list[Sample],
    rsx: str,
    threads: int,
    min_depths: list[int],
    results: ResultWriter,
    force: bool,
) -> None:
    logs_dir = dataset_dir / "logs"
    base = result_base(dataset, samples)
    markers_table = command_output(dataset_dir, "process")
    if force or not markers_table.exists():
        elapsed = run_command(
            [
                rsx,
                "process",
                "-i",
                str(dataset_dir / "samples"),
                "-o",
                str(markers_table),
                "-T",
                str(threads),
                "-d",
                "1",
            ],
            logs_dir / "process.log",
        )
    else:
        elapsed = 0.0
    results.write(
        base
        | {
            "command": "process",
            "min_depth": "1",
            "elapsed_seconds": f"{elapsed:.6f}",
        }
        | summarize_output("process", markers_table)
    )

    depth_path = command_output(dataset_dir, "depth")
    if force or not depth_path.exists():
        elapsed = run_command(
            [
                rsx,
                "depth",
                "-t",
                str(markers_table),
                "-p",
                str(dataset_dir / "popmap.tsv"),
                "-o",
                str(depth_path),
            ],
            logs_dir / "depth.log",
        )
    else:
        elapsed = 0.0
    results.write(
        base
        | {
            "command": "depth",
            "elapsed_seconds": f"{elapsed:.6f}",
        }
        | summarize_output("depth", depth_path)
    )

    for min_depth in min_depths:
        for command in ("freq", "distrib", "signif"):
            output = command_output(dataset_dir, command, min_depth)
            args = [
                rsx,
                command,
                "-t",
                str(markers_table),
                "-o",
                str(output),
                "-d",
                str(min_depth),
            ]
            if command in {"distrib", "signif"}:
                args.extend(["-p", str(dataset_dir / "popmap.tsv"), "-G", "male,female"])
            if command == "signif":
                args.append("-a")
            if force or not output.exists():
                elapsed = run_command(args, logs_dir / f"{command}_{min_depth}.log")
            else:
                elapsed = 0.0
            results.write(
                base
                | {
                    "command": command,
                    "min_depth": str(min_depth),
                    "elapsed_seconds": f"{elapsed:.6f}",
                }
                | summarize_output(command, output)
            )


def summarize_results_csv(results: Path, summary: Path) -> None:
    by_dataset: dict[str, dict[str, str]] = {}
    elapsed: dict[str, float] = {}
    with results.open(newline="") as handle:
        for row in csv.DictReader(handle):
            dataset = row["dataset"]
            by_dataset.setdefault(dataset, row)
            try:
                elapsed[dataset] = elapsed.get(dataset, 0.0) + float(row.get("elapsed_seconds") or 0.0)
            except ValueError:
                pass
            if row["command"] == "process":
                by_dataset[dataset]["markers"] = row["markers"]
    lines = [
        "dataset\tsamples\tspots\tbases\tmarkers\telapsed_seconds",
    ]
    for dataset in sorted(by_dataset):
        row = by_dataset[dataset]
        lines.append(
            "\t".join(
                [
                    dataset,
                    row.get("samples", ""),
                    row.get("total_spots", ""),
                    row.get("total_bases", ""),
                    row.get("markers", ""),
                    f"{elapsed.get(dataset, 0.0):.3f}",
                ]
            )
        )
    summary.parent.mkdir(parents=True, exist_ok=True)
    summary.write_text("\n".join(lines) + "\n")


def parse_min_depths(value: str) -> list[int]:
    depths = [int(item) for item in value.split(",") if item.strip()]
    if not depths:
        raise argparse.ArgumentTypeError("at least one minimum depth is required")
    return depths


def resolve_binary(value: str) -> str:
    if os.sep in value:
        return value
    resolved = shutil.which(value)
    return resolved or value


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", default=Path("benchmarks/literature_datasets.tsv"), type=Path)
    parser.add_argument("--workdir", default=Path("benchmarks/literature-workdir"), type=Path)
    parser.add_argument("--results", default=Path("benchmarks/results/literature_benchmark_results.csv"), type=Path)
    parser.add_argument("--summary", default=Path("docs/figures/literature_benchmark_summary.tsv"), type=Path)
    parser.add_argument("--dataset", action="append", help="Dataset name to run; repeat for multiple datasets")
    parser.add_argument("--rsx", default=str(Path("target/release/rsx")))
    parser.add_argument("--fastq-dump", default="fastq-dump")
    parser.add_argument("--threads", default=4, type=int)
    parser.add_argument("--min-depths", default="1,2,5,10", type=parse_min_depths)
    parser.add_argument("--metadata-only", action="store_true")
    parser.add_argument("--download-only", action="store_true")
    parser.add_argument("--no-download", action="store_true")
    parser.add_argument("--download-method", default="ena", choices=["ena", "fastq-dump"])
    parser.add_argument("--append", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--retries", default=5, type=int)
    args = parser.parse_args()

    datasets = read_manifest(args.manifest)
    selected = set(args.dataset or [])
    if selected:
        datasets = [dataset for dataset in datasets if dataset.name in selected]
    missing = selected.difference({dataset.name for dataset in datasets})
    if missing:
        raise SystemExit(f"unknown dataset(s): {', '.join(sorted(missing))}")

    if not args.append:
        prune_dataset_results(args.results, {dataset.name for dataset in datasets})
    writer = ResultWriter(args.results)
    try:
        for dataset in datasets:
            dataset_dir = args.workdir / dataset.name
            samples = read_dataset_metadata(dataset_dir)
            if not samples or args.force:
                xml_text = query_sra_xml(dataset, retries=args.retries)
                samples = parse_sra_experiment_xml(xml_text)
                write_dataset_metadata(dataset_dir, samples)
            base = result_base(dataset, samples)
            writer.write(
                base
                | {
                    "command": "metadata",
                    "output_bytes": str((dataset_dir / "samples.tsv").stat().st_size),
                    "output_path": relpath(dataset_dir / "samples.tsv"),
                }
            )
            if args.metadata_only:
                continue

            if not args.no_download:
                elapsed = download_samples(
                    dataset_dir,
                    samples,
                    resolve_binary(args.fastq_dump),
                    dataset_dir / "logs",
                    args.force,
                    args.download_method,
                )
                writer.write(
                    base
                    | {
                        "command": "download",
                        "elapsed_seconds": f"{elapsed:.6f}",
                        "output_bytes": str(total_file_size((dataset_dir / "samples").glob("*.fq.gz"))),
                        "output_path": relpath(dataset_dir / "samples"),
                    }
                )
            if args.download_only:
                continue
            run_rsx_dataset(
                dataset,
                dataset_dir,
                samples,
                resolve_binary(args.rsx),
                args.threads,
                args.min_depths,
                writer,
                args.force,
            )
    finally:
        writer.close()
    summarize_results_csv(args.results, args.summary)
    print(f"Wrote {args.results}")
    print(f"Wrote {args.summary}")


if __name__ == "__main__":
    main()
