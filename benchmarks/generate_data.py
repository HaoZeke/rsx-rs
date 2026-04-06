#!/usr/bin/env python3
"""Generate synthetic RADSex benchmark data at multiple scales."""

import os
import random
import gzip
import argparse
from pathlib import Path

def generate_fastq_files(outdir, n_individuals, n_reads_per_ind, read_length=100):
    """Generate FASTQ files for n_individuals with n_reads_per_ind reads each."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Create a pool of ~1000 unique sequences to sample from
    # This simulates realistic RAD-seq where many sequences are shared
    random.seed(42)
    bases = "ATCG"
    n_unique = min(1000, n_reads_per_ind)
    seq_pool = ["".join(random.choices(bases, k=read_length)) for _ in range(n_unique)]

    # 10% of sequences are sex-biased
    male_only = seq_pool[:n_unique // 10]
    female_only = seq_pool[n_unique // 10 : n_unique // 5]
    common = seq_pool[n_unique // 5:]

    half = n_individuals // 2
    for i in range(1, n_individuals + 1):
        is_male = i <= half
        name = f"ind{i}"
        path = outdir / f"{name}.fq.gz"

        with gzip.open(path, "wt") as f:
            for r in range(n_reads_per_ind):
                # Pick a sequence: males get male_only, females get female_only
                if r < len(male_only) and is_male:
                    seq = male_only[r % len(male_only)]
                elif r < len(male_only) + len(female_only) and not is_male:
                    seq = female_only[r % len(female_only)]
                else:
                    seq = random.choice(common)

                # Add depth variation (1-5 copies)
                depth = random.randint(1, 5)
                for d in range(depth):
                    f.write(f"@{name}_r{r}_d{d}\n{seq}\n+\n{'I' * read_length}\n")

    # Create popmap
    popmap_path = outdir.parent / f"popmap_{n_individuals}.tsv"
    with open(popmap_path, "w") as f:
        for i in range(1, n_individuals + 1):
            sex = "M" if i <= half else "F"
            f.write(f"ind{i}\t{sex}\n")

    return str(outdir), str(popmap_path)


def generate_markers_table(outpath, n_markers, n_individuals):
    """Generate a markers depth table directly (for benchmarking downstream commands)."""
    random.seed(42)
    half = n_individuals // 2

    with open(outpath, "w") as f:
        f.write(f"#Number of markers : {n_markers}\n")
        header = "id\tsequence\t" + "\t".join(f"ind{i}" for i in range(1, n_individuals + 1))
        f.write(header + "\n")

        for m in range(n_markers):
            seq = "".join(random.choices("ATCG", k=100))
            depths = []
            for j in range(1, n_individuals + 1):
                is_male = j <= half
                if m < n_markers // 10:  # Male-biased
                    d = random.randint(5, 50) if is_male else random.randint(0, 2)
                elif m < n_markers // 5:  # Female-biased
                    d = random.randint(0, 2) if is_male else random.randint(5, 50)
                else:  # Common
                    d = random.randint(5, 50)
                depths.append(str(d))
            f.write(f"{m}\t{seq}\t" + "\t".join(depths) + "\n")


def generate_genome(outpath, n_contigs=40, contig_length=1_000_000):
    """Generate a synthetic genome FASTA."""
    random.seed(42)
    with open(outpath, "w") as f:
        for c in range(1, n_contigs + 1):
            f.write(f">chr{c}\n")
            seq = "".join(random.choices("ATCG", k=contig_length))
            # Write in 80-char lines
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate RADSex benchmark data")
    parser.add_argument("--outdir", default="benchmarks/data")
    parser.add_argument("--scales", nargs="+", default=["small", "medium", "large"])
    args = parser.parse_args()

    outdir = Path(args.outdir)

    # Scale definitions: (n_individuals, n_markers, n_reads_per_ind)
    scales = {
        "small":  (10,   1_000,   500),
        "medium": (20,  10_000,  2_000),
        "large":  (40, 100_000,  5_000),
    }

    for scale_name in args.scales:
        if scale_name not in scales:
            print(f"Unknown scale: {scale_name}")
            continue

        n_ind, n_markers, n_reads = scales[scale_name]
        scale_dir = outdir / scale_name
        scale_dir.mkdir(parents=True, exist_ok=True)

        print(f"=== {scale_name}: {n_ind} individuals, {n_markers} markers, {n_reads} reads/ind ===")

        # Generate markers table (for downstream commands)
        table_path = scale_dir / "markers.tsv"
        print(f"  Generating markers table: {table_path}")
        generate_markers_table(str(table_path), n_markers, n_ind)

        # Generate popmap
        popmap_path = scale_dir / "popmap.tsv"
        with open(popmap_path, "w") as f:
            for i in range(1, n_ind + 1):
                sex = "M" if i <= n_ind // 2 else "F"
                f.write(f"ind{i}\t{sex}\n")
        print(f"  Popmap: {popmap_path}")

        # Generate FASTQ files (only for small/medium to keep disk use reasonable)
        if scale_name in ("small", "medium"):
            reads_dir = scale_dir / "reads"
            print(f"  Generating FASTQ files: {reads_dir}")
            generate_fastq_files(str(reads_dir), n_ind, n_reads)

        # Generate genome (small, for map command)
        genome_path = scale_dir / "genome.fa"
        contig_len = 10_000 if scale_name == "small" else 100_000
        n_contigs = 10 if scale_name == "small" else 40
        print(f"  Generating genome: {genome_path} ({n_contigs} contigs x {contig_len}bp)")
        generate_genome(str(genome_path), n_contigs, contig_len)

    print("\nDone!")
