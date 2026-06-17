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
      chrom  pos(1-based)  ref  alt  gene  label  [vaf]

    Column 7 (vaf) is optional. When provided, it overrides the global --fetal-fraction
    target VAF for that variant (direct VAF, 0.0–1.0, e.g. 0.5 for heterozygous,
    1.0 for homozygous).

    Returns list of dicts with keys: chrom, pos0(0-based), ref, alt, gene, label, vaf.
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
            if chrom == "chrom":  # skip header row
                continue
            pos1  = int(parts[1])          # 1-based VCF/human convention
            ref   = parts[2].upper()
            alt   = parts[3].upper()
            gene  = parts[4] if len(parts) > 4 else ""
            label = parts[5] if len(parts) > 5 else f"{chrom}:{pos1}_{ref}>{alt}"
            per_vaf = None
            if len(parts) > 6 and parts[6].strip():
                try:
                    per_vaf = float(parts[6].strip())
                except ValueError:
                    pass
            variants.append({
                "chrom": chrom, "pos0": pos1 - 1,
                "ref": ref, "alt": alt,
                "gene": gene, "label": label,
                "spike_base": alt,          # what we inject
                "per_vaf": per_vaf,         # optional per-variant VAF override
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


def check_maternal_genotype(bam_in: pysam.AlignmentFile,
                           chrom: str,
                           pos0: int,
                           ref: str,
                           alt: str) -> str:
    """
    Check the maternal (source BAM) genotype at a position.

    Returns 'hom_ref' if ≥90% of reads match ref, 'het' if 30-70% alt,
    'hom_alt' if ≥90% match alt, or 'other'.
    """
    try:
        reads = [
            r for r in bam_in.fetch(chrom, pos0, pos0 + 1)
            if not r.is_unmapped and not r.is_secondary and not r.is_supplementary
        ]
    except ValueError:
        return "no_data"
    if not reads:
        return "no_data"

    ref_count = alt_count = other_count = 0
    for r in reads:
        pairs = {rp: qp for qp, rp in r.get_aligned_pairs(matches_only=True)}
        qp = pairs.get(pos0)
        if qp is None:
            continue
        base = r.query_sequence[qp].upper()
        if base == ref.upper():
            ref_count += 1
        elif base == alt.upper():
            alt_count += 1
        else:
            other_count += 1

    total = ref_count + alt_count + other_count
    if total == 0:
        return "no_data"
    ref_frac = ref_count / total
    if ref_frac >= 0.90:
        return "hom_ref"
    elif ref_frac <= 0.10:
        return "hom_alt"
    elif 0.30 <= ref_frac <= 0.70:
        return "het"
    return "other"


def spike_reads(bam_in: pysam.AlignmentFile,
                chrom: str,
                pos0: int,
                spike_base: str,
                target_vaf: float,
                rng: random.Random) -> list:
    """
    Clone reads overlapping pos0 and change the target base to spike_base.

    For target_vaf < 0.95 (additive mode):
        Number of reads to add = round(existing_depth * vaf / (1 - vaf)).
    For target_vaf >= 0.95 (replace mode):
        Replace all REF-carrying reads with ALT reads in-place (no new reads added).
        Returns an empty list; the BAM is modified via the returned modified_reads dict.

    Returns list of new AlignedSegment objects (additive mode only).
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

    def _make_alt_read(tmpl, i):
        ref_to_query = {rp: qp for qp, rp in tmpl.get_aligned_pairs(matches_only=True)}
        qpos = ref_to_query.get(pos0)
        if qpos is None:
            return None
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
        seq       = list(tmpl.query_sequence)
        seq[qpos] = spike_base
        new_r.query_sequence = "".join(seq)
        if tmpl.query_qualities is not None:
            new_r.query_qualities = tmpl.query_qualities
        new_r.set_tags(tmpl.get_tags())
        return new_r

    if target_vaf >= 0.95:
        # Replace mode: add enough ALT copies so REF is overwhelmed (homozygous simulation)
        # Strategy: for each non-ALT read, add one ALT clone; keep everything
        new_reads = []
        for i, tmpl in enumerate(existing):
            ref_to_query = {rp: qp for qp, rp in tmpl.get_aligned_pairs(matches_only=True)}
            qpos = ref_to_query.get(pos0)
            if qpos is None:
                continue
            base = tmpl.query_sequence[qpos].upper()
            if base != spike_base:
                r = _make_alt_read(tmpl, i)
                if r:
                    new_reads.append(r)
        return new_reads

    # Additive mode
    n_to_add = max(1, round(n_existing * target_vaf / max(1 - target_vaf, 1e-6)))
    new_reads = []
    for i in range(n_to_add):
        tmpl = rng.choice(existing)
        r = _make_alt_read(tmpl, i)
        if r:
            new_reads.append(r)
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
    ap.add_argument("--fetal-fraction", required=False, type=float, default=0.0,
                    help="Fetal fraction to simulate (0.0–0.5). Use 0 (default) for "
                         "direct-VAF mode where each variant's VAF is read from TSV column 7.")
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
    # Allow ff=0 as a sentinel meaning "use per-variant VAF from TSV column 7"
    direct_vaf_mode = (ff == 0.0)
    if not direct_vaf_mode and not (0 < ff < 0.5):
        sys.exit(f"[ERROR] --fetal-fraction must be between 0 and 0.5, or 0 for direct-VAF mode (got {ff})")

    target_vaf = ff / 2.0 if not direct_vaf_mode else None
    rng        = random.Random(args.seed)

    if direct_vaf_mode:
        print("[simulate_nipt_bam] Direct-VAF mode: per-variant VAF from TSV column 7",
              file=sys.stderr)
    else:
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

    # Pre-flight: verify genotype at disease variant positions
    if disease_variants:
        print("[simulate_nipt_bam] Verifying genotype at disease positions …",
              file=sys.stderr)
        skipped = []
        for v in disease_variants:
            gt = check_maternal_genotype(bam_in, v["chrom"], v["pos0"], v["ref"], v["alt"])
            per_vaf = v.get("per_vaf")
            if direct_vaf_mode and per_vaf is not None:
                # In direct-VAF mode, skip if no data or already at/above target
                if gt == "no_data":
                    skipped.append(v)
                    print(f"  [SKIP]  {v['label']} — no reads at position", file=sys.stderr)
                elif (per_vaf >= 0.95 and gt == "hom_alt") or \
                     (per_vaf < 0.95 and gt == "het"):
                    # Already at target — no extra spiking needed
                    v["_skip_spike"] = True
                    print(f"  [check] {v['label']}  {v['chrom']}:{v['pos0']+1}  "
                          f"ref={v['ref']} alt={v['alt']}  current_gt={gt}  "
                          f"target_vaf={per_vaf:.0%}  SKIP (already at target)",
                          file=sys.stderr)
                else:
                    print(f"  [check] {v['label']}  {v['chrom']}:{v['pos0']+1}  "
                          f"ref={v['ref']} alt={v['alt']}  current_gt={gt}  "
                          f"target_vaf={per_vaf:.0%}  OK",
                          file=sys.stderr)
            else:
                status = "OK" if gt == "hom_ref" else f"WARN ({gt})"
                print(f"  [check] {v['label']}  {v['chrom']}:{v['pos0']+1}  "
                      f"ref={v['ref']} alt={v['alt']}  maternal={gt}  {status}",
                      file=sys.stderr)
                if gt != "hom_ref":
                    skipped.append(v)
                    print(f"  [SKIP]  {v['label']} — maternal is NOT hom-ref, spike would be meaningless",
                          file=sys.stderr)
        for v in skipped:
            disease_variants.remove(v)
        if skipped:
            print(f"[simulate_nipt_bam] Skipped {len(skipped)} variants",
                  file=sys.stderr)

    total_added = 0
    extra_reads = []

    for item in disease_variants + ff_snps:
        chrom      = item["chrom"]
        pos0       = item["pos0"]
        spike_base = item["spike_base"]
        label      = item["label"]

        # Skip positions already at their target VAF
        if item.get("_skip_spike"):
            print(f"  [skip]  {label}  {chrom}:{pos0+1}  already at target VAF",
                  file=sys.stderr)
            continue

        # Resolve effective VAF: per-variant override > global target_vaf
        per_vaf = item.get("per_vaf")
        if per_vaf is not None and direct_vaf_mode:
            effective_vaf = per_vaf
        elif target_vaf is not None:
            effective_vaf = target_vaf
        else:
            print(f"  [SKIP]  {label} — no VAF defined (set --fetal-fraction or add vaf column)",
                  file=sys.stderr)
            continue

        added = spike_reads(bam_in, chrom, pos0, spike_base, effective_vaf, rng)
        extra_reads.extend(added)
        total_added += len(added)
        print(f"  [spike] {label}  {chrom}:{pos0+1}  →{spike_base}  "
              f"added {len(added)} reads (VAF≈{effective_vaf:.1%})", file=sys.stderr)

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
