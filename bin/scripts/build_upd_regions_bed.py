#!/usr/bin/env python3
"""
build_upd_regions_bed.py — Build a BED file of analysis windows on the seven
UPD-relevant chromosomes (6, 7, 11, 14, 15, 16, 20).

Strategy
--------
The Twist UCL_SingleGeneNIPT panel covers CDS exons (~1 Mb total). For UPD
analysis we need many more informative maternal-heterozygous SNPs per
chromosome, so we extend the analysis window to the GENE locus
(transcription start to end) plus a configurable flanking region (default
50 kb on each side), and we MERGE the resulting per-gene intervals so that
adjacent / overlapping genes form one analysis window.

Inputs
------
--target-bed   The panel target BED (col 4 = "GENE_..." style).  We only use
               the gene name prefix — actual gene boundaries come from
               --gene-coords.
--gene-coords  TSV with columns: gene, chrom, start, end  (1-based or 0-based;
               we treat as 0-based half-open BED).  Easiest source:
                 mysql --user=genome --host=genome-mysql.soe.ucsc.edu -A -D hg38 \
                   -e "SELECT name2 AS gene, chrom, txStart, txEnd FROM ncbiRefSeq \
                       WHERE chrom IN ('chr6','chr7','chr11','chr14','chr15','chr16','chr20')" \
                   > hg38.upd_genes.tsv
               Or use refGene / GENCODE GFF parsed to {gene, chrom, min(txStart), max(txEnd)}.
--flank        Bases to extend each side of every gene (default 50000).
--chroms       Comma-separated chrom list (default chr6,chr7,chr11,chr14,chr15,chr16,chr20).
--output       Output BED path.

Output
------
A 4-column merged BED (chrom, start, end, gene-list) sorted by chrom, start.
Coordinates are clipped to ≥ 0 (no upper clip; downstream tools use --fai if needed).

Usage
-----
python3 bin/scripts/build_upd_regions_bed.py \
    --target-bed   data/target_bed/All_target_regions_*.bed \
    --gene-coords  data/upd/hg38.upd_genes.tsv \
    --flank        50000 \
    --output       data/upd/upd_analysis_regions.bed
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

DEFAULT_UPD_CHROMS = ["chr6", "chr7", "chr11", "chr14", "chr15", "chr16", "chr20"]


def parse_target_bed_genes(bed_path: str) -> set:
    """Return the set of gene-name prefixes found in column 4 of the panel BED."""
    genes: set = set()
    with open(bed_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 4:
                continue
            label = cols[3]
            gene = label.split("_", 1)[0].strip()
            if gene:
                genes.add(gene)
    return genes


def parse_gene_coords(tsv_path: str) -> dict:
    """
    Parse a 4-col TSV (gene, chrom, start, end). Aggregates per gene by taking
    the min(start) and max(end) so that all isoforms collapse into a single
    locus interval.
    """
    by_gene: dict[str, dict] = {}
    with open(tsv_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#") or line.lower().startswith("gene\t"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            gene, chrom, start_s, end_s = parts[0], parts[1], parts[2], parts[3]
            try:
                start, end = int(start_s), int(end_s)
            except ValueError:
                continue
            if not gene or not chrom or end <= start:
                continue
            entry = by_gene.get(gene)
            if entry is None:
                by_gene[gene] = {"chrom": chrom, "start": start, "end": end}
            else:
                if chrom != entry["chrom"]:
                    continue
                entry["start"] = min(entry["start"], start)
                entry["end"] = max(entry["end"], end)
    return by_gene


def merge_intervals(intervals: list[tuple[int, int, str]]) -> list[tuple[int, int, list[str]]]:
    """Merge overlapping/adjacent intervals; preserve gene names."""
    if not intervals:
        return []
    intervals = sorted(intervals, key=lambda x: (x[0], x[1]))
    merged: list[tuple[int, int, list[str]]] = []
    cur_s, cur_e, cur_names = intervals[0][0], intervals[0][1], [intervals[0][2]]
    for s, e, name in intervals[1:]:
        if s <= cur_e:
            cur_e = max(cur_e, e)
            if name not in cur_names:
                cur_names.append(name)
        else:
            merged.append((cur_s, cur_e, cur_names))
            cur_s, cur_e, cur_names = s, e, [name]
    merged.append((cur_s, cur_e, cur_names))
    return merged


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--target-bed", required=True)
    ap.add_argument("--gene-coords", required=True,
                    help="TSV: gene, chrom, start, end (BED-style 0-based)")
    ap.add_argument("--flank", type=int, default=50_000)
    ap.add_argument("--chroms", default=",".join(DEFAULT_UPD_CHROMS))
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    upd_chroms = {c.strip() for c in args.chroms.split(",") if c.strip()}
    panel_genes = parse_target_bed_genes(args.target_bed)
    gene_coords = parse_gene_coords(args.gene_coords)

    by_chrom: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    missing: list[str] = []
    for gene in panel_genes:
        info = gene_coords.get(gene)
        if info is None:
            missing.append(gene)
            continue
        if info["chrom"] not in upd_chroms:
            continue
        s = max(0, info["start"] - args.flank)
        e = info["end"] + args.flank
        by_chrom[info["chrom"]].append((s, e, gene))

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    total_bp = 0
    with out.open("w") as fh:
        for chrom in sorted(by_chrom.keys(),
                            key=lambda c: int(c.replace("chr", "")) if c.replace("chr", "").isdigit() else 99):
            merged = merge_intervals(by_chrom[chrom])
            for s, e, names in merged:
                fh.write(f"{chrom}\t{s}\t{e}\t{'|'.join(names)}\n")
                total_bp += e - s

    print(f"[build_upd_regions_bed] wrote {out}", file=sys.stderr)
    print(f"  chromosomes: {sorted(by_chrom.keys())}", file=sys.stderr)
    print(f"  total bp:    {total_bp:,}", file=sys.stderr)
    if missing:
        print(f"  WARN: {len(missing)} panel genes had no entry in --gene-coords "
              f"(first 10: {missing[:10]})", file=sys.stderr)


if __name__ == "__main__":
    main()
