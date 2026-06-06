#!/usr/bin/env python3
"""
upsample_fastq.py — Upsample paired-end FASTQ by duplicating reads with renamed IDs.

Each copy gets a unique suffix appended to the read name so that downstream
MarkDuplicates tools treat them as independent fragments (different read names
→ different template IDs → not flagged as duplicates).

Usage:
    python3 upsample_fastq.py \\
        --r1 input_R1.fastq.gz \\
        --r2 input_R2.fastq.gz \\
        --multiplier 4 \\
        --output-r1 upsampled_R1.fastq.gz \\
        --output-r2 upsampled_R2.fastq.gz

Multiplier:
    1 = no-op (output == input)
    4 = original + 3 copies  (4× total read count)

Output:
    Gzip-compressed FASTQ files. Reads are written in copy order:
    all original reads first, then copy 1, copy 2, …
    Read names: original unchanged; copies get suffix  :dup:<N>  appended.
"""

import argparse
import gzip
import sys
import os
from typing import BinaryIO


def open_fastq(path: str, mode: str = "rb") -> BinaryIO:
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def copy_with_suffix(r1_in: str, r2_in: str,
                     r1_out: BinaryIO, r2_out: BinaryIO,
                     suffix: bytes, report_interval: int = 5_000_000) -> int:
    """Stream one pass over r1_in / r2_in, write to open output handles."""
    count = 0
    with open_fastq(r1_in) as f1, open_fastq(r2_in) as f2:
        while True:
            # Read 4 lines per read
            h1 = f1.readline()
            if not h1:
                break
            s1 = f1.readline()
            p1 = f1.readline()
            q1 = f1.readline()

            h2 = f2.readline()
            s2 = f2.readline()
            p2 = f2.readline()
            q2 = f2.readline()

            if not (s1 and p1 and q1 and h2 and s2 and p2 and q2):
                break

            if suffix:
                # Append suffix before any trailing whitespace / comment on the header
                # FASTQ header: @READ_NAME [optional comment]
                h1 = h1.rstrip(b"\n")
                h2 = h2.rstrip(b"\n")
                # Insert suffix after the read name (before first space if present)
                def inject(header: bytes, sfx: bytes) -> bytes:
                    assert header.startswith(b"@"), f"Bad header: {header!r}"
                    parts = header[1:].split(b" ", 1)
                    new_name = parts[0] + sfx
                    if len(parts) == 2:
                        return b"@" + new_name + b" " + parts[1] + b"\n"
                    return b"@" + new_name + b"\n"

                h1 = inject(h1, suffix)
                h2 = inject(h2, suffix)

            r1_out.write(h1)
            r1_out.write(s1)
            r1_out.write(p1)
            r1_out.write(q1)

            r2_out.write(h2)
            r2_out.write(s2)
            r2_out.write(p2)
            r2_out.write(q2)

            count += 1
            if count % report_interval == 0:
                label = suffix.decode() if suffix else "(original)"
                print(f"  [{label}] {count:,} reads written …", file=sys.stderr)

    return count


def main():
    parser = argparse.ArgumentParser(
        description="Upsample paired FASTQ by duplicating reads with unique names."
    )
    parser.add_argument("--r1",           required=True, help="Input R1 FASTQ(.gz)")
    parser.add_argument("--r2",           required=True, help="Input R2 FASTQ(.gz)")
    parser.add_argument("--multiplier",   required=True, type=int,
                        help="Total read multiplier (e.g. 4 = original + 3 copies)")
    parser.add_argument("--output-r1",    required=True, help="Output R1 FASTQ.gz")
    parser.add_argument("--output-r2",    required=True, help="Output R2 FASTQ.gz")
    parser.add_argument("--compress-level", type=int, default=1,
                        help="gzip compression level 1–9 (default 1 = fastest)")
    args = parser.parse_args()

    if args.multiplier < 1:
        parser.error("--multiplier must be ≥ 1")

    for f in [args.r1, args.r2]:
        if not os.path.exists(f):
            sys.exit(f"ERROR: file not found: {f}")

    print(f"Upsampling {os.path.basename(args.r1)} / {os.path.basename(args.r2)}", file=sys.stderr)
    print(f"  Multiplier : {args.multiplier}×  ({args.multiplier - 1} copies + original)", file=sys.stderr)
    print(f"  Output R1  : {args.output_r1}", file=sys.stderr)
    print(f"  Output R2  : {args.output_r2}", file=sys.stderr)

    total = 0
    with gzip.open(args.output_r1, "wb", compresslevel=args.compress_level) as out1, \
         gzip.open(args.output_r2, "wb", compresslevel=args.compress_level) as out2:

        # Pass 0: original reads (no suffix)
        print(f"\n[Pass 0/{ args.multiplier - 1 }] Writing original reads …", file=sys.stderr)
        n = copy_with_suffix(args.r1, args.r2, out1, out2, suffix=b"")
        total += n
        print(f"  → {n:,} read pairs written", file=sys.stderr)

        # Passes 1..M-1: copies with :dup:N suffix
        for copy_idx in range(1, args.multiplier):
            suffix = f":dup:{copy_idx}".encode()
            print(f"\n[Pass {copy_idx}/{args.multiplier - 1}] Writing copy {copy_idx} (suffix {suffix.decode()}) …",
                  file=sys.stderr)
            n = copy_with_suffix(args.r1, args.r2, out1, out2, suffix=suffix)
            total += n
            print(f"  → {n:,} read pairs written", file=sys.stderr)

    print(f"\nDone. Total read pairs written: {total:,}  ({args.multiplier}× of {total // args.multiplier:,})",
          file=sys.stderr)


if __name__ == "__main__":
    main()
