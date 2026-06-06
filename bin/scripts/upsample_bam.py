#!/usr/bin/env python3
"""
upsample_bam.py — Upsample a BAM file by cloning all reads N times with renamed IDs.

Each copy gets a unique suffix appended to query_name so MarkDuplicates treats
them as independent fragments. The input BAM should already be post-markdup
(e.g. a .target.bam produced by the sgNIPT pipeline). The output BAM is
coordinate-sorted and indexed.

Usage:
    python3 upsample_bam.py \\
        --bam   input.target.bam \\
        --multiplier 6 \\
        --output upsampled.target.bam \\
        [--threads 4]

Multiplier:
    1 = no-op (output == input)
    6 = original + 5 copies  (6× total read count / depth)

Output:
    Coordinate-sorted, indexed BAM.
    Read names: original unchanged; copy k gets suffix  :dup:k  (k=1..N-1).
    Both reads of a pair get the same suffix → stays properly paired.
"""

import argparse
import os
import sys
import tempfile

import pysam


def upsample(bam_in: str, bam_out: str, multiplier: int, threads: int) -> None:
    if multiplier < 1:
        raise ValueError("--multiplier must be ≥ 1")

    with pysam.AlignmentFile(bam_in, "rb", threads=threads) as src:
        header = src.header.to_dict()
        total_reads = 0

        # Write to a temp unsorted BAM, then sort at the end
        tmp_unsorted = bam_out + ".tmp_unsorted.bam"
        try:
            with pysam.AlignmentFile(tmp_unsorted, "wb",
                                     header=header, threads=threads) as dst:

                for copy_idx in range(multiplier):
                    suffix = f":dup:{copy_idx}" if copy_idx > 0 else ""
                    label = suffix if suffix else "(original)"
                    print(f"[Pass {copy_idx + 1}/{multiplier}] Writing {label} …",
                          file=sys.stderr)

                    src.reset()
                    n = 0
                    for read in src:
                        if suffix:
                            read.query_name = read.query_name + suffix
                        dst.write(read)
                        n += 1
                        if n % 1_000_000 == 0:
                            print(f"  {label}: {n:,} reads …", file=sys.stderr)

                    print(f"  → {n:,} reads written", file=sys.stderr)
                    total_reads += n

            print(f"\nSorting → {bam_out} …", file=sys.stderr)
            pysam.sort("-o", bam_out, "-@", str(threads), tmp_unsorted)

            print("Indexing …", file=sys.stderr)
            pysam.index(bam_out, "-@", str(threads))

        finally:
            if os.path.exists(tmp_unsorted):
                os.remove(tmp_unsorted)

    print(f"\nDone. Total reads written: {total_reads:,}  "
          f"({multiplier}× of {total_reads // multiplier:,})", file=sys.stderr)
    print(f"Output: {bam_out}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description="Upsample a BAM file by cloning reads with unique names."
    )
    parser.add_argument("--bam",        required=True, help="Input BAM (coordinate-sorted)")
    parser.add_argument("--multiplier", required=True, type=int,
                        help="Total depth multiplier (e.g. 6 = original + 5 copies)")
    parser.add_argument("--output",     required=True, help="Output BAM path")
    parser.add_argument("--threads",    type=int, default=4,
                        help="SAM I/O threads (default 4)")
    args = parser.parse_args()

    if not os.path.exists(args.bam):
        sys.exit(f"ERROR: BAM not found: {args.bam}")

    print(f"Input BAM  : {args.bam}", file=sys.stderr)
    print(f"Multiplier : {args.multiplier}×  ({args.multiplier - 1} copies + original)",
          file=sys.stderr)
    print(f"Output BAM : {args.output}", file=sys.stderr)
    print(f"Threads    : {args.threads}", file=sys.stderr)

    upsample(args.bam, args.output, args.multiplier, args.threads)


if __name__ == "__main__":
    main()
