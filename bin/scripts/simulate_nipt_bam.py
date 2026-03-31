#!/usr/bin/env python3
"""
simulate_nipt_bam.py — Spike fetal reads into a maternal BAM to simulate cfDNA NIPT.

Strategy
--------
In real NIPT cfDNA, fetal DNA appears at ~FF/2 VAF at fetal-specific (paternal) loci,
because the cfDNA pool is FF% fetal and the fetal variant is heterozygous.

This script simulates that signal by cloning reads from an aligned maternal BAM and
changing the target base to the desired allele:

  1. Disease variants  (--variants TSV):
       Maternal is hom-ref (0/0) at pathogenic positions → spike ALT reads at FF/2 VAF.
       Simulates a fetal-specific (paternal) pathogenic allele.

  2. FF-informative SNPs  (--ff-vcf + --n-ff-snps):
       Uses homozygous-alt (GT=1/1) sites from an existing VCF → spike REF reads at FF/2.
       Simulates: maternal=1/1, fetal=0/1 → cfDNA sees REF at FF/2 as minor allele.
       fetal_fraction.py (SNP method) then detects these as informative sites.

Output
------
Modified BAM is sorted and indexed.  The original maternal reads are preserved;
synthetic reads are appended before sorting.  Each synthetic read gets a unique
name (original_name + _sim_<i>) to avoid marking-duplicate issues.

Usage
-----
python3 simulate_nipt_bam.py \\
    --bam     NA12878.dedup.bam \\
    --output  sim_FF10.bam \\
    --fetal-fraction 0.10 \\
    --variants data/test/test_variants.tsv \\
    --ff-vcf  NA12878.ff_variants.vcf \\
    --n-ff-snps 100 \\
    --seed 42

Note: Run inside the sgnipt Docker image (pysam is installed there).
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
    """
    Parse a TSV of variants to spike in.

    Expected columns (tab-separated, comment lines start with #):
      chrom  pos(1-based)  ref  alt  gene  label

    Returns list of dicts with keys: chrom, pos0(0-based), ref, alt, gene, label.
    """
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
            pos1  = int(parts[1])          # 1-based VCF/human convention
            ref   = parts[2].upper()
            alt   = parts[3].upper()
            gene  = parts[4] if len(parts) > 4 else ""
            label = parts[5] if len(parts) > 5 else f"{chrom}:{pos1}_{ref}>{alt}"
            variants.append({
                "chrom": chrom, "pos0": pos1 - 1,
                "ref": ref, "alt": alt,
                "gene": gene, "label": label,
                "spike_base": alt,          # what we inject
            })
    return variants


def parse_ff_snps(vcf_path: str, n_snps: int, min_depth: int = 30) -> list:
    """
    Extract up to n_snps homozygous-alt (GT=1/1) SNP sites from a VCF.

    At these sites the maternal genome is hom-alt; a fetal het contributes
    the REF allele at FF/2.  We spike REF reads to simulate that signal.

    Returns list of dicts: chrom, pos0, ref, alt, spike_base(=ref), label.
    """
    snps = []
    seen = set()
    with open(vcf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 10:
                continue
            chrom, pos1, ref, alt = parts[0], int(parts[1]), parts[3], parts[4]

            # SNPs only
            if len(ref) != 1 or len(alt) != 1:
                continue
            key = (chrom, pos1)
            if key in seen:
                continue
            seen.add(key)

            fmt_keys = parts[8].split(":")
            fmt_vals = parts[9].split(":")
            fmt = dict(zip(fmt_keys, fmt_vals))

            gt = fmt.get("GT", "./.").replace("|", "/")
            if gt not in ("1/1",):
                continue

            ad = fmt.get("AD", "0,0").split(",")
            try:
                depth = int(ad[0]) + int(ad[1])
            except (ValueError, IndexError):
                continue
            if depth < min_depth:
                continue

            snps.append({
                "chrom": chrom, "pos0": pos1 - 1,
                "ref": ref, "alt": alt,
                "gene": "", "label": f"ff_snp_{chrom}:{pos1}",
                "spike_base": ref,          # inject REF to simulate fetal contribution
            })
            if len(snps) >= n_snps:
                break

    return snps


def spike_reads(bam_in: pysam.AlignmentFile,
                chrom: str,
                pos0: int,
                spike_base: str,
                target_vaf: float,
                rng: random.Random) -> list:
    """
    Clone reads overlapping pos0 and change the target base to spike_base.

    Number of reads to add = round(existing_depth * vaf / (1 - vaf)).
    Returns list of new AlignedSegment objects.
    """
    try:
        existing = [
            r for r in bam_in.fetch(chrom, pos0, pos0 + 1)
            if not r.is_unmapped
            and not r.is_secondary
            and not r.is_supplementary
        ]
    except ValueError:
        print(f"  [WARN] Cannot fetch {chrom}:{pos0+1} — contig not in BAM index", file=sys.stderr)
        return []

    if not existing:
        print(f"  [WARN] No reads at {chrom}:{pos0+1}", file=sys.stderr)
        return []

    n_existing = len(existing)
    n_to_add   = max(1, round(n_existing * target_vaf / max(1 - target_vaf, 1e-6)))

    new_reads = []
    for i in range(n_to_add):
        tmpl = rng.choice(existing)

        # Find which query position corresponds to this reference position
        ref_to_query = {rp: qp for qp, rp in tmpl.get_aligned_pairs(matches_only=True)}
        qpos = ref_to_query.get(pos0)
        if qpos is None:
            continue                        # read doesn't cover exact pos (soft-clip etc.)

        new_r = pysam.AlignedSegment(bam_in.header)
        new_r.query_name          = f"{tmpl.query_name}_sim_{i}"
        new_r.flag                = tmpl.flag & ~0x400  # clear duplicate flag
        new_r.reference_id        = tmpl.reference_id
        new_r.reference_start     = tmpl.reference_start
        new_r.mapping_quality     = tmpl.mapping_quality
        new_r.cigar               = tmpl.cigar
        new_r.next_reference_id   = tmpl.next_reference_id
        new_r.next_reference_start = tmpl.next_reference_start
        new_r.template_length     = tmpl.template_length

        # Modify the base
        seq       = list(tmpl.query_sequence)
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
        description="Simulate NIPT cfDNA BAM with spiked fetal alleles",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument("--bam",            required=True, help="Input BAM (sorted + indexed)")
    ap.add_argument("--output",         required=True, help="Output BAM path")
    ap.add_argument("--fetal-fraction", required=True, type=float,
                    help="Fetal fraction to simulate (0.0–0.5), e.g. 0.10 for 10%%")
    ap.add_argument("--variants",       default=None,
                    help="TSV of disease variants to spike (chrom,pos,ref,alt,gene,label)")
    ap.add_argument("--ff-vcf",         default=None,
                    help="VCF file with hom-alt sites for FF-informative SNP spiking")
    ap.add_argument("--n-ff-snps",      type=int, default=100,
                    help="Number of hom-alt sites to use for FF simulation (default: 100)")
    ap.add_argument("--ff-min-depth",   type=int, default=30,
                    help="Minimum read depth at FF-informative sites (default: 30)")
    ap.add_argument("--seed",           type=int, default=42)
    args = ap.parse_args()

    ff  = args.fetal_fraction
    if not (0 < ff < 0.5):
        sys.exit(f"[ERROR] --fetal-fraction must be between 0 and 0.5 (got {ff})")

    target_vaf = ff / 2.0           # fetal het variant → appears at FF/2 in cfDNA
    rng        = random.Random(args.seed)

    print(f"[simulate_nipt_bam] Fetal fraction: {ff:.1%}  →  spike VAF: {target_vaf:.1%}",
          file=sys.stderr)

    # ---- load variant lists ----
    disease_variants = []
    if args.variants:
        disease_variants = parse_variants_tsv(args.variants)
        print(f"[simulate_nipt_bam] Disease variants to spike: {len(disease_variants)}",
              file=sys.stderr)
    else:
        print("[simulate_nipt_bam] No --variants file provided; skipping disease variant spiking",
              file=sys.stderr)

    ff_snps = []
    if args.ff_vcf:
        ff_snps = parse_ff_snps(args.ff_vcf, args.n_ff_snps, args.ff_min_depth)
        print(f"[simulate_nipt_bam] FF-informative SNP sites loaded: {len(ff_snps)}",
              file=sys.stderr)
    else:
        print("[simulate_nipt_bam] No --ff-vcf provided; FF estimation will not be improved",
              file=sys.stderr)

    if not disease_variants and not ff_snps:
        sys.exit("[ERROR] Nothing to spike. Provide --variants and/or --ff-vcf.")

    # ---- spike reads ----
    bam_in    = pysam.AlignmentFile(args.bam, "rb")
    tmp_path  = args.output + ".tmp_unsorted.bam"

    total_added = 0
    extra_reads = []

    for item in disease_variants + ff_snps:
        chrom      = item["chrom"]
        pos0       = item["pos0"]
        spike_base = item["spike_base"]
        label      = item["label"]

        added = spike_reads(bam_in, chrom, pos0, spike_base, target_vaf, rng)
        extra_reads.extend(added)
        total_added += len(added)
        print(f"  [spike] {label}  {chrom}:{pos0+1}  →{spike_base}  "
              f"added {len(added)} reads (VAF≈{target_vaf:.1%})", file=sys.stderr)

    # ---- write output ----
    print(f"[simulate_nipt_bam] Writing {total_added} synthetic reads + all original reads …",
          file=sys.stderr)

    with pysam.AlignmentFile(tmp_path, "wb", header=bam_in.header) as out:
        bam_in.reset()
        for read in bam_in:
            out.write(read)
        for read in extra_reads:
            out.write(read)

    bam_in.close()

    # ---- sort + index ----
    print(f"[simulate_nipt_bam] Sorting → {args.output}", file=sys.stderr)
    pysam.sort("-o", args.output, tmp_path)
    os.remove(tmp_path)

    print(f"[simulate_nipt_bam] Indexing …", file=sys.stderr)
    pysam.index(args.output)

    print(f"[simulate_nipt_bam] Done. Output: {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
