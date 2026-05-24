#!/usr/bin/env python3
"""
simulate_carrier_bam.py — Spike pathogenic variants into a WES BAM for carrier screening tests.

Unlike simulate_nipt_bam.py (cfDNA/fetal fraction), this script targets WES (germline) data:
  - Default VAF = 0.5  (heterozygous carrier)
  - No fetal-fraction SNP spiking
  - Checks maternal (source) genotype: skips positions already carrying the ALT allele

Usage:
    python3 simulate_carrier_bam.py \\
        --bam   NA12878.md.bam \\
        --output carrier_spiked.bam \\
        --variants carrier_variants.tsv \\
        --vaf 0.5 \\
        --seed 42

Note: Run inside a Docker image that has pysam (e.g., the sgnipt image).
"""

import argparse
import os
import random
import sys

try:
    import pysam
except ImportError:
    sys.exit("[ERROR] pysam is required. Run inside the sgnipt Docker container.")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def parse_variants_tsv(path: str) -> list:
    variants = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            chrom = parts[0]
            pos1  = int(parts[1])
            ref   = parts[2].upper()
            alt   = parts[3].upper()
            gene  = parts[4] if len(parts) > 4 else ""
            label = parts[5] if len(parts) > 5 else f"{chrom}:{pos1}_{ref}>{alt}"
            variants.append({
                "chrom": chrom, "pos0": pos1 - 1,
                "ref": ref, "alt": alt,
                "gene": gene, "label": label,
                "spike_base": alt,
            })
    return variants


def check_source_genotype(bam_in: pysam.AlignmentFile, chrom: str, pos0: int,
                           ref: str, alt: str) -> str:
    try:
        reads = [r for r in bam_in.fetch(chrom, pos0, pos0 + 1)
                 if not r.is_unmapped and not r.is_secondary and not r.is_supplementary]
    except ValueError:
        return "no_data"
    if not reads:
        return "no_data"

    ref_c = alt_c = other_c = 0
    for r in reads:
        pairs = {rp: qp for qp, rp in r.get_aligned_pairs(matches_only=True)}
        qp = pairs.get(pos0)
        if qp is None:
            continue
        base = r.query_sequence[qp].upper()
        if base == ref.upper():
            ref_c += 1
        elif base == alt.upper():
            alt_c += 1
        else:
            other_c += 1

    total = ref_c + alt_c + other_c
    if total == 0:
        return "no_data"
    ref_frac = ref_c / total
    if ref_frac >= 0.90:
        return "hom_ref"
    elif ref_frac <= 0.10:
        return "hom_alt"
    elif 0.30 <= ref_frac <= 0.70:
        return "het"
    return "other"


def spike_reads(bam_in: pysam.AlignmentFile, chrom: str, pos0: int,
                spike_base: str, target_vaf: float, rng: random.Random) -> list:
    try:
        existing = [r for r in bam_in.fetch(chrom, pos0, pos0 + 1)
                    if not r.is_unmapped and not r.is_secondary and not r.is_supplementary]
    except ValueError:
        print(f"  [WARN] Cannot fetch {chrom}:{pos0+1}", file=sys.stderr)
        return []

    if not existing:
        print(f"  [WARN] No reads at {chrom}:{pos0+1}", file=sys.stderr)
        return []

    n_existing = len(existing)
    # To reach target_vaf: n_add / (n_existing + n_add) = target_vaf
    # → n_add = n_existing * target_vaf / (1 - target_vaf)
    n_to_add = max(1, round(n_existing * target_vaf / max(1 - target_vaf, 1e-6)))

    new_reads = []
    for i in range(n_to_add):
        tmpl = rng.choice(existing)
        ref_to_query = {rp: qp for qp, rp in tmpl.get_aligned_pairs(matches_only=True)}
        qpos = ref_to_query.get(pos0)
        if qpos is None:
            continue

        new_r = pysam.AlignedSegment(bam_in.header)
        new_r.query_name           = f"{tmpl.query_name}_sim_{i}"
        new_r.flag                 = tmpl.flag & ~0x400
        new_r.reference_id         = tmpl.reference_id
        new_r.reference_start      = tmpl.reference_start
        new_r.mapping_quality      = tmpl.mapping_quality
        new_r.cigar                = tmpl.cigar
        new_r.next_reference_id    = tmpl.next_reference_id
        new_r.next_reference_start = tmpl.next_reference_start
        new_r.template_length      = tmpl.template_length

        seq = list(tmpl.query_sequence)
        seq[qpos] = spike_base
        new_r.query_sequence = "".join(seq)

        if tmpl.query_qualities is not None:
            new_r.query_qualities = tmpl.query_qualities
        new_r.set_tags(tmpl.get_tags())
        new_reads.append(new_r)

    return new_reads


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Spike pathogenic variants into a WES BAM for carrier screening tests",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument("--bam",       required=True, help="Input WES BAM (sorted + indexed)")
    ap.add_argument("--output",    required=True, help="Output BAM path")
    ap.add_argument("--variants",  required=True, help="TSV of variants to spike (chrom,pos,ref,alt,gene,label)")
    ap.add_argument("--vaf",       type=float, default=0.5,
                    help="Target variant allele fraction (default 0.5 = heterozygous carrier)")
    ap.add_argument("--seed",      type=int, default=42)
    args = ap.parse_args()

    if not (0 < args.vaf < 1.0):
        sys.exit(f"[ERROR] --vaf must be between 0 and 1 (got {args.vaf})")

    rng = random.Random(args.seed)
    variants = parse_variants_tsv(args.variants)
    print(f"[simulate_carrier_bam] Variants to spike: {len(variants)}  VAF={args.vaf:.1%}",
          file=sys.stderr)

    # ── Verify source genotype ──────────────────────────────────────────────
    print("[simulate_carrier_bam] Checking source genotype …", file=sys.stderr)
    bam_in = pysam.AlignmentFile(args.bam, "rb")
    skipped = []
    to_spike = []
    for v in variants:
        gt = check_source_genotype(bam_in, v["chrom"], v["pos0"], v["ref"], v["alt"])
        status = "OK" if gt == "hom_ref" else "WARN"
        print(f"  [check] {v['label']:50s} {v['chrom']}:{v['pos0']+1}  "
              f"maternal={gt}  {status}", file=sys.stderr)
        if gt == "hom_ref" or gt == "no_data":
            to_spike.append(v)
        else:
            skipped.append(v["label"])
    bam_in.close()

    if skipped:
        print(f"[simulate_carrier_bam] Skipped {len(skipped)} variants "
              f"(source already carries ALT): {skipped[:5]}{'...' if len(skipped)>5 else ''}",
              file=sys.stderr)

    # ── Spike reads ─────────────────────────────────────────────────────────
    bam_in   = pysam.AlignmentFile(args.bam, "rb")
    tmp_path = args.output + ".tmp_unsorted.bam"
    bam_out  = pysam.AlignmentFile(tmp_path, "wb", header=bam_in.header)

    all_new_reads = []
    for v in to_spike:
        new_reads = spike_reads(bam_in, v["chrom"], v["pos0"], v["spike_base"], args.vaf, rng)
        if new_reads:
            actual_vaf = len(new_reads) / (len(new_reads) + max(1, len(new_reads)))
            print(f"  [spike] {v['label']:50s} {v['chrom']}:{v['pos0']+1}  "
                  f"→{v['spike_base']}  added {len(new_reads)} reads (VAF≈{args.vaf:.1%})",
                  file=sys.stderr)
            all_new_reads.extend(new_reads)

    total_added = len(all_new_reads)
    for read in bam_in.fetch(until_eof=True):
        bam_out.write(read)
    for read in all_new_reads:
        bam_out.write(read)
    bam_in.close()
    bam_out.close()

    print(f"[simulate_carrier_bam] Writing {total_added} synthetic reads + all original reads …",
          file=sys.stderr)
    print(f"[simulate_carrier_bam] Sorting → {args.output}", file=sys.stderr)
    pysam.sort("-o", args.output, tmp_path)
    os.remove(tmp_path)
    print(f"[simulate_carrier_bam] Indexing …", file=sys.stderr)
    pysam.index(args.output)
    print(f"[simulate_carrier_bam] Done. Output: {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
