"""ASV benchmarks for rsx-rs binary performance."""

import shutil
import subprocess
import tempfile
from pathlib import Path

BENCH_DATA = Path(__file__).parent / "data"
RSX = shutil.which("rsx") or str(
    Path(__file__).parent.parent / "target" / "release" / "rsx"
)


def _run(args, **kwargs):
    """Run radsex with given arguments."""
    return subprocess.run(
        [RSX] + args,
        capture_output=True,
        check=True,
        **kwargs,
    )


class TimeFreqSmall:
    """Benchmark freq command on small dataset (1K markers, 10 individuals)."""

    timeout = 30
    repeat = 5
    number = 1
    warmup_time = 0

    def setup(self):
        self.tmpdir = tempfile.mkdtemp(prefix="asv_radsex_")
        self.table = str(BENCH_DATA / "small" / "markers.tsv")
        self.output = str(Path(self.tmpdir) / "freq.tsv")

    def teardown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def time_freq_small(self):
        _run(["freq", "-t", self.table, "-o", self.output, "-d", "5"])


class TimeFreqMedium:
    """Benchmark freq command on medium dataset (10K markers, 20 individuals)."""

    timeout = 30
    repeat = 5
    number = 1
    warmup_time = 0

    def setup(self):
        self.tmpdir = tempfile.mkdtemp(prefix="asv_radsex_")
        self.table = str(BENCH_DATA / "medium" / "markers.tsv")
        self.output = str(Path(self.tmpdir) / "freq.tsv")

    def teardown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def time_freq_medium(self):
        _run(["freq", "-t", self.table, "-o", self.output, "-d", "5"])


class TimeFreqLarge:
    """Benchmark freq command on large dataset (100K markers, 40 individuals)."""

    timeout = 60
    repeat = 5
    number = 1
    warmup_time = 0

    def setup(self):
        self.tmpdir = tempfile.mkdtemp(prefix="asv_radsex_")
        self.table = str(BENCH_DATA / "large" / "markers.tsv")
        self.output = str(Path(self.tmpdir) / "freq.tsv")

    def teardown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def time_freq_large(self):
        _run(["freq", "-t", self.table, "-o", self.output, "-d", "5"])


class TimeDistribLarge:
    """Benchmark distrib command on large dataset."""

    timeout = 60
    repeat = 5
    number = 1
    warmup_time = 0

    def setup(self):
        self.tmpdir = tempfile.mkdtemp(prefix="asv_radsex_")
        self.table = str(BENCH_DATA / "large" / "markers.tsv")
        self.popmap = str(BENCH_DATA / "large" / "popmap.tsv")
        self.output = str(Path(self.tmpdir) / "distrib.tsv")

    def teardown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def time_distrib_large(self):
        _run(["distrib", "-t", self.table, "-p", self.popmap,
              "-o", self.output, "-d", "5", "-G", "M,F"])


class TimeSignifLarge:
    """Benchmark signif command on large dataset."""

    timeout = 60
    repeat = 5
    number = 1
    warmup_time = 0

    def setup(self):
        self.tmpdir = tempfile.mkdtemp(prefix="asv_radsex_")
        self.table = str(BENCH_DATA / "large" / "markers.tsv")
        self.popmap = str(BENCH_DATA / "large" / "popmap.tsv")
        self.output = str(Path(self.tmpdir) / "signif.tsv")

    def teardown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def time_signif_large(self):
        _run(["signif", "-t", self.table, "-p", self.popmap,
              "-o", self.output, "-d", "5", "-G", "M,F"])


class TimeDepthLarge:
    """Benchmark depth command on large dataset."""

    timeout = 60
    repeat = 5
    number = 1
    warmup_time = 0

    def setup(self):
        self.tmpdir = tempfile.mkdtemp(prefix="asv_radsex_")
        self.table = str(BENCH_DATA / "large" / "markers.tsv")
        self.popmap = str(BENCH_DATA / "large" / "popmap.tsv")
        self.output = str(Path(self.tmpdir) / "depth.tsv")

    def teardown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def time_depth_large(self):
        _run(["depth", "-t", self.table, "-p", self.popmap,
              "-o", self.output])


class TimeProcessMedium:
    """Benchmark process command on medium dataset (20 individuals, FASTQ)."""

    timeout = 120
    repeat = 3
    number = 1
    warmup_time = 0

    def setup(self):
        self.tmpdir = tempfile.mkdtemp(prefix="asv_radsex_")
        self.reads = str(BENCH_DATA / "medium" / "reads")
        self.output = str(Path(self.tmpdir) / "markers.tsv")

    def teardown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def time_process_medium(self):
        _run(["process", "-i", self.reads, "-o", self.output, "-T", "4", "-d", "1"])


class TimeMapMedium:
    """Benchmark map command on medium dataset."""

    timeout = 120
    repeat = 3
    number = 1
    warmup_time = 0

    def setup(self):
        self.tmpdir = tempfile.mkdtemp(prefix="asv_radsex_")
        self.table = str(BENCH_DATA / "medium" / "markers.tsv")
        self.popmap = str(BENCH_DATA / "medium" / "popmap.tsv")
        self.genome = str(BENCH_DATA / "medium" / "genome.fa")
        self.output = str(Path(self.tmpdir) / "map.tsv")

    def teardown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def time_map_medium(self):
        _run(["map", "-t", self.table, "-p", self.popmap,
              "-g", self.genome, "-o", self.output,
              "-d", "5", "-G", "M,F", "-q", "0"])
