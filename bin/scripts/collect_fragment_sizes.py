#!/usr/bin/env python3
"""
collect_fragment_sizes.py - Collect cfDNA fragment size distribution from BAM file.

Extracts insert sizes from properly paired reads in a BAM file and generates
a fragment size distribution file for fetal fraction estimation.
Uses samtools via subprocess to avoid pysam dependency.
"""

import argparse
import json
import logging
import subprocess
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, Tuple

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)


def collect_insert_sizes(
    bam_path: str,
    min_size: int = 50,
    max_size: int = 600,
    max_reads: int = 10_000_000,
    region: str = None,
) -> Counter:
    """
    Collect insert sizes from a BAM file using samtools.
    Filters for properly paired, primary alignments with MAPQ >= 30.
    """
    logger.info("Collecting insert sizes from %s", bam_path)

    # samtools view flags:
    # -f 2: properly paired
    # -F 3852: exclude unmapped, mate unmapped, not primary, supplementary, duplicate, failed QC
    # -q 30: minimum mapping quality
    cmd = [
        "samtools", "view",
        "-f", "2",
        "-F", "3852",
        "-q", "30",
        bam_path,
    ]
    if region:
        cmd.append(region)

    logger.info("Running: %s", " ".join(cmd))

    size_counter = Counter()
    read_count = 0

    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)

    for line in proc.stdout:
        fields = line.split("\t")
        if len(fields) < 9:
            continue

        try:
            tlen = abs(int(fields[8]))
        except ValueError:
            continue

        if min_size <= tlen <= max_size:
            size_counter[tlen] += 1

        read_count += 1
        if read_count >= max_reads:
            logger.info("Reached max reads limit (%d)", max_reads)
            break

        if read_count % 1_000_000 == 0:
            logger.info("Processed %d reads...", read_count)

    proc.stdout.close()
    proc.wait()

    logger.info(
        "Collected %d fragment sizes from %d reads",
        sum(size_counter.values()),
        read_count,
    )
    return size_counter


def calculate_statistics(size_counter: Counter) -> Dict:
    """Calculate fragment size distribution statistics."""
    if not size_counter:
        return {}

    total = sum(size_counter.values())
    sizes = []
    for size, count in size_counter.items():
        sizes.extend([size] * count)
    sizes.sort()

    mean_size = sum(sizes) / len(sizes)
    median_size = sizes[len(sizes) // 2]

    # Mode
    mode_size = size_counter.most_common(1)[0][0]

    # Percentiles
    p10 = sizes[int(len(sizes) * 0.10)]
    p25 = sizes[int(len(sizes) * 0.25)]
    p75 = sizes[int(len(sizes) * 0.75)]
    p90 = sizes[int(len(sizes) * 0.90)]

    # Fragment categories for cfDNA
    short_frag = sum(c for s, c in size_counter.items() if s < 150)
    mono_nucleosomal = sum(c for s, c in size_counter.items() if 150 <= s <= 200)
    di_nucleosomal = sum(c for s, c in size_counter.items() if 300 <= s <= 400)
    long_frag = sum(c for s, c in size_counter.items() if s > 250)

    stats = {
        "total_fragments": total,
        "mean_size": round(mean_size, 1),
        "median_size": median_size,
        "mode_size": mode_size,
        "p10": p10,
        "p25": p25,
        "p75": p75,
        "p90": p90,
        "short_fragments_lt150": short_frag,
        "short_fragment_ratio": round(short_frag / total, 4) if total > 0 else 0,
        "mono_nucleosomal_150_200": mono_nucleosomal,
        "mono_nucleosomal_ratio": round(mono_nucleosomal / total, 4) if total > 0 else 0,
        "di_nucleosomal_300_400": di_nucleosomal,
        "long_fragments_gt250": long_frag,
        "long_fragment_ratio": round(long_frag / total, 4) if total > 0 else 0,
    }

    logger.info(
        "Fragment stats: mean=%.1f, median=%d, mode=%d, short_ratio=%.4f",
        mean_size,
        median_size,
        mode_size,
        stats["short_fragment_ratio"],
    )
    return stats


def count_chromosome_reads(bam_path: str) -> Tuple[int, int]:
    """Count reads on autosomes vs Y chromosome using samtools idxstats."""
    logger.info("Counting chromosome reads from %s", bam_path)

    result = subprocess.run(
        ["samtools", "idxstats", bam_path],
        capture_output=True,
        text=True,
    )

    autosome_reads = 0
    y_reads = 0

    for line in result.stdout.strip().split("\n"):
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        chrom = parts[0]
        mapped = int(parts[2])

        if chrom in ("chrY", "Y"):
            y_reads = mapped
        elif chrom.replace("chr", "").isdigit():
            chrom_num = int(chrom.replace("chr", ""))
            if 1 <= chrom_num <= 22:
                autosome_reads += mapped

    logger.info("Autosome reads: %d, Y reads: %d", autosome_reads, y_reads)
    return autosome_reads, y_reads


def main():
    parser = argparse.ArgumentParser(
        description="Collect cfDNA fragment size distribution from BAM file."
    )
    parser.add_argument("--bam", required=True, help="Input BAM file (sorted, indexed).")
    parser.add_argument("--sample-id", required=True, help="Sample identifier.")
    parser.add_argument("--output-sizes", required=True, help="Output fragment size distribution file.")
    parser.add_argument("--output-stats", required=True, help="Output statistics JSON file.")
    parser.add_argument("--min-size", type=int, default=50, help="Minimum fragment size.")
    parser.add_argument("--max-size", type=int, default=600, help="Maximum fragment size.")
    parser.add_argument("--max-reads", type=int, default=10_000_000, help="Maximum reads to process.")
    parser.add_argument("--region", default=None, help="Genomic region to restrict analysis.")
    parser.add_argument("--count-chroms", action="store_true", help="Also count chromosome reads for Y-chr FF.")
    args = parser.parse_args()

    # Collect insert sizes
    size_counter = collect_insert_sizes(
        args.bam,
        min_size=args.min_size,
        max_size=args.max_size,
        max_reads=args.max_reads,
        region=args.region,
    )

    # Write size distribution
    output_sizes_path = Path(args.output_sizes)
    output_sizes_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_sizes_path, "w") as fh:
        fh.write("#size\tcount\n")
        for size in sorted(size_counter.keys()):
            fh.write(f"{size}\t{size_counter[size]}\n")
    logger.info("Fragment size distribution written to %s", output_sizes_path)

    # Calculate statistics
    stats = calculate_statistics(size_counter)
    stats["sample_id"] = args.sample_id

    # Count chromosome reads if requested
    if args.count_chroms:
        autosome_reads, y_reads = count_chromosome_reads(args.bam)
        stats["autosome_reads"] = autosome_reads
        stats["y_chromosome_reads"] = y_reads

    # Write statistics
    output_stats_path = Path(args.output_stats)
    output_stats_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_stats_path, "w") as fh:
        json.dump(stats, fh, indent=2)
    logger.info("Fragment statistics written to %s", output_stats_path)


if __name__ == "__main__":
    main()
