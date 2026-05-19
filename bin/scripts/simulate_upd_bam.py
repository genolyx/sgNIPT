#!/usr/bin/env python3
"""
simulate_upd_bam.py — Synthetic BAM for UPD end-to-end BAM→mpileup→upd_detection testing.

Generates a BAM with reads at common SNP positions on a target chromosome.
Each read is a 150 bp fragment centered on the SNP; the base at the SNP
position encodes the desired allele.

VAF per site is chosen based on the specified UPD scenario:

  normal      : mat-het sites → VAF ≈ 0.50; mat-hom-ref → VAF ≈ 0.00
  matUPD_iso  : mat-het sites → VAF = 0.50 ± FF/2 (each site independently shifted)
  matUPD_het  : mat-het sites → VAF ≈ 0.50 (no shift; het-disomy = normal appearance)
  patUPD_iso  : mat-hom-ref sites → ~40 % show Mendelian error at VAF ≈ FF
  patUPD_het  : mat-hom-ref sites → ~20 % show Mendelian error at VAF ≈ FF/2

Usage (inside sgnipt docker):
  python3 simulate_upd_bam.py \\
      --vcf /data/upd/common_snps_maf05.upd.vcf.gz \\
      --chrom chr7 \\
      --scenario matUPD_iso \\
      --fetal-fraction 0.10 \\
      --depth 600 \\
      --n-snps 150 \\
      --output /tmp/matUPD7_FF10.bam
"""

import argparse
import gzip
import math
import os
import random
import sys

try:
    import pysam
except ImportError:
    sys.exit("[ERROR] pysam is required. Run inside the sgnipt Docker container.")


# hg38 chromosome lengths (needed for BAM header)
CHROM_LENGTHS = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559,
    "chr4": 190214555, "chr5": 181538259, "chr6": 170805979,
    "chr7": 159345973, "chr8": 145138636, "chr9": 138394717,
    "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189,
    "chr16": 90338345,  "chr17": 83257441,  "chr18": 80373285,
    "chr19": 58617616,  "chr20": 64444167,  "chr21": 46709983,
    "chr22": 50818468,  "chrX": 156040895,  "chrY": 57227415,
}


def parse_snps(vcf_path: str, chrom: str, n_snps: int,
               min_spacing: int = 200) -> list:
    """
    Load up to n_snps from VCF for the given chromosome.

    min_spacing: reject SNPs closer than this many bp to the previous accepted
    SNP.  This prevents read pile-ups from adjacent SNPs from diluting each
    other's VAF (a read of length 150 starting at pos-75 would overlap any
    SNP within ±75 bp, so spacing must exceed read_length = 150).
    Default 200 bp gives a comfortable margin.
    """
    open_fn = gzip.open if vcf_path.endswith(".gz") else open
    snps = []
    last_pos = -999999
    with open_fn(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip().split("\t")
            if len(cols) < 5:
                continue
            if cols[0] != chrom:
                continue
            ref, alt = cols[3], cols[4]
            if len(ref) != 1 or len(alt) != 1:
                continue
            pos = int(cols[1])
            if pos - last_pos < min_spacing:
                continue
            last_pos = pos
            # try to pull AF from INFO
            maf = 0.30
            for tag in cols[7].split(";"):
                if tag.upper().startswith("AF="):
                    try:
                        maf = float(tag.split("=")[1].split(",")[0])
                        if maf > 0.5:
                            maf = 1.0 - maf
                    except ValueError:
                        pass
                    break
            snps.append({"chrom": chrom, "pos1": pos, "ref": ref, "alt": alt, "maf": maf})
            if len(snps) >= n_snps:
                break
    return snps


def assign_maternal_gt(maf: float, rng: random.Random) -> str:
    """Random genotype drawn from Hardy-Weinberg."""
    r = rng.random()
    p_hom_ref = (1.0 - maf) ** 2
    p_het = 2.0 * maf * (1.0 - maf)
    if r < p_hom_ref:
        return "hom_ref"
    if r < p_hom_ref + p_het:
        return "het"
    return "hom_alt"


def target_vaf(mat_gt: str, scenario: str, ff: float, rng: random.Random) -> float:
    """
    Expected ALT VAF in cfDNA given maternal genotype and UPD scenario.

    cfDNA pool:  (1-FF) × maternal + FF × fetal allele frequencies.
    """
    if mat_gt == "hom_ref":
        if scenario == "normal":
            return 0.0
        if scenario == "patUPD_iso":
            # ~40 % of mat-hom-ref sites: father is het → fetal is hom_alt → Mendelian error
            return ff if rng.random() < 0.40 else 0.0
        if scenario == "patUPD_het":
            # ~20 % of mat-hom-ref sites: father is het → fetal is het → partial signal
            return (ff / 2) if rng.random() < 0.20 else 0.0
        return 0.0  # matUPD: fetal inherited maternal allele → also hom_ref

    if mat_gt == "het":
        if scenario == "normal":
            return 0.50
        if scenario == "matUPD_iso":
            # fetal inherits two copies of ONE maternal allele (iso)
            # 50% chance fetal is hom_ref → shift down; 50% hom_alt → shift up
            delta = ff / 2.0
            return (0.50 + delta) if rng.random() < 0.50 else (0.50 - delta)
        if scenario == "matUPD_het":
            return 0.50  # het-disomy: fetal inherits both maternal alleles → no shift
        if scenario == "patUPD_iso":
            # Fetal has no maternal allele; has two copies of one paternal allele
            # Paternal allele at a maternal-het site is unpredictable → treat as ~0.5 on average
            # Each site: 50% fetal is hom_ref, 50% hom_alt
            fetal_alt_frac = 1.0 if rng.random() < 0.5 else 0.0
            return (1.0 - ff) * 0.50 + ff * fetal_alt_frac
        if scenario == "patUPD_het":
            return 0.50  # paternal het-disomy behaves like normal at mat-het sites
        return 0.50

    # hom_alt
    if scenario == "normal":
        return 1.0
    if scenario == "matUPD_iso":
        return 1.0  # fetal inherited maternal alt allele → also hom_alt
    return 1.0


def make_read(
    read_id: int,
    chrom: str,
    read_start: int,
    snp_offset: int,
    allele_base: str,
    read_len: int,
    bam_header: pysam.AlignmentHeader,
    ref_fasta: "pysam.FastaFile",
) -> pysam.AlignedSegment:
    a = pysam.AlignedSegment(bam_header)
    a.query_name = f"updSim_{read_id:08d}"
    a.reference_id = bam_header.get_tid(chrom)
    a.reference_start = read_start
    a.mapping_quality = 60
    a.cigar = [(0, read_len)]  # nM

    # Use actual reference bases so only the SNP position differs from ref.
    chrom_len = CHROM_LENGTHS.get(chrom, 10**9)
    fetch_end = min(read_start + read_len, chrom_len)
    ref_seq = ref_fasta.fetch(chrom, read_start, fetch_end).upper()
    if len(ref_seq) < read_len:
        ref_seq = ref_seq + "N" * (read_len - len(ref_seq))

    seq = list(ref_seq)
    if 0 <= snp_offset < read_len:
        seq[snp_offset] = allele_base
    a.query_sequence = "".join(seq)
    a.query_qualities = pysam.qualitystring_to_array("I" * read_len)
    a.set_tag("RG", "sim")
    a.is_paired = False
    return a


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf",       required=True,  help="common_snps_maf05.upd.vcf.gz")
    ap.add_argument("--reference", required=True,  help="reference FASTA (must be indexed)")
    ap.add_argument("--chrom",  default="chr7", help="target chromosome")
    ap.add_argument("--scenario", default="matUPD_iso",
                    choices=["normal", "matUPD_iso", "matUPD_het",
                             "patUPD_iso", "patUPD_het"],
                    help="UPD scenario to simulate")
    ap.add_argument("--fetal-fraction", type=float, default=0.10,
                    help="fetal fraction (default: 0.10)")
    ap.add_argument("--depth",  type=int, default=600,
                    help="reads per SNP site (default: 600)")
    ap.add_argument("--n-snps", type=int, default=150,
                    help="number of SNPs to simulate (default: 150)")
    ap.add_argument("--min-snp-spacing", type=int, default=200,
                    help="minimum bp gap between successive SNPs (default: 200); "
                         "prevents read-overlap VAF dilution")
    ap.add_argument("--output", required=True, help="output BAM path")
    ap.add_argument("--seed",   type=int, default=42)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    ff  = args.fetal_fraction
    read_len = 150

    print(f"[INFO] Loading SNPs for {args.chrom} from {args.vcf}")
    snps = parse_snps(args.vcf, args.chrom, args.n_snps, args.min_snp_spacing)
    if not snps:
        sys.exit(f"[ERROR] No SNPs found for {args.chrom} in {args.vcf}")
    print(f"[INFO] Loaded {len(snps)} SNPs | scenario={args.scenario} | FF={ff:.2%} | depth={args.depth}")

    print(f"[INFO] Opening reference FASTA: {args.reference}")
    ref_fasta = pysam.FastaFile(args.reference)

    # Build BAM header with all hg38 chromosomes
    bam_header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": c, "LN": l} for c, l in CHROM_LENGTHS.items()],
        "RG": [{"ID": "sim", "SM": f"sim_{args.scenario}",
                "LB": "simLib", "PL": "ILLUMINA"}],
    })

    tmp_bam = args.output + ".tmp_unsorted.bam"
    total_reads = 0
    n_mat_het = n_mat_hom_ref = n_mat_hom_alt = 0

    with pysam.AlignmentFile(tmp_bam, "wb", header=bam_header) as out:
        for snp in snps:
            mat_gt = assign_maternal_gt(snp["maf"], rng)
            vaf    = target_vaf(mat_gt, args.scenario, ff, rng)

            pos0       = snp["pos1"] - 1          # 0-based
            snp_mid    = read_len // 2             # 75
            read_start = max(0, pos0 - snp_mid)
            snp_offset = pos0 - read_start

            # Sample observed allele counts from binomial distribution to
            # add realistic sequencing noise (avoids perfectly precise VAF).
            if 0.0 < vaf < 1.0:
                n_alt = rng.binomialvariate(args.depth, vaf)
            else:
                n_alt = round(args.depth * vaf)
            n_ref = args.depth - n_alt

            for _ in range(n_ref):
                out.write(make_read(total_reads, snp["chrom"], read_start,
                                    snp_offset, snp["ref"], read_len, bam_header,
                                    ref_fasta))
                total_reads += 1
            for _ in range(n_alt):
                out.write(make_read(total_reads, snp["chrom"], read_start,
                                    snp_offset, snp["alt"], read_len, bam_header,
                                    ref_fasta))
                total_reads += 1

            if mat_gt == "het":        n_mat_het     += 1
            elif mat_gt == "hom_ref":  n_mat_hom_ref += 1
            else:                      n_mat_hom_alt += 1

    print(f"[INFO] Genotype counts → het:{n_mat_het}  hom_ref:{n_mat_hom_ref}  hom_alt:{n_mat_hom_alt}")
    print(f"[INFO] Total reads written: {total_reads:,}")
    print(f"[INFO] Sorting → {args.output}")
    pysam.sort("-o", args.output, tmp_bam)
    pysam.index(args.output)
    os.remove(tmp_bam)
    print(f"[INFO] Done: {args.output}  ({os.path.getsize(args.output):,} bytes)")


if __name__ == "__main__":
    main()
