#!/usr/bin/env python3
"""
simulate_dark_gene_bam.py — Simulate dark gene dosage/allele changes in cfDNA BAM.

Supported diseases and scenarios
─────────────────────────────────
alpha_thal  (HBA1/HBA2  chr16:172,000-178,000)
    --scenario carrier      fetal 3 copies  (-α/αα)
    --scenario affected     fetal 2 copies  (--/αα or -α/-α)
    --scenario hbh          fetal 1 copy    (--/-α, HbH disease)
    --scenario hydrops      fetal 0 copies  (--/--, Hb Bart's)

sma  (SMN1/SMN2  chr5:70,920,000-70,970,000, c.840 C→T at 70,951,946)
    --scenario carrier      fetal 1×SMN1 + 2×SMN2
    --scenario affected     fetal 0×SMN1 + 4×SMN2

gaucher  (GBA1/GBAP1  chr1:155,210,000-155,245,000)
    --scenario carrier      fetal 1 copy
    --scenario affected     fetal 0 copies

Mechanism
─────────
alpha_thal / gaucher: Remove reads randomly from the dark gene window.
  n_remove = round(region_reads × FF × (1 − fetal_cn / normal_cn))

sma: Change existing C reads to T at chr5:70,951,946.
  fetal_smn1_frac = smn1_copies / (smn1_copies + smn2_copies)
  Assuming maternal smn1_frac ≈ 0.5 (normal 2×SMN1 + 2×SMN2).
  Δ = (maternal_frac − expected_observed_frac) = FF × (maternal_frac − fetal_frac)
  n_to_change = round(total_c840_reads × Δ / current_C_frac)

Usage:
    python3 simulate_dark_gene_bam.py \\
        --bam           input.bam \\
        --disease       alpha_thal \\
        --scenario      carrier \\
        --fetal-fraction 0.10 \\
        --output        output.bam \\
        [--threads 4] [--seed 42]
"""

import argparse
import os
import random
import sys
from typing import Set

import pysam

# ── Dark gene region definitions (GRCh38) ─────────────────────────────────────

DARK_GENE_REGIONS = {
    "alpha_thal": {
        "chrom": "chr16",
        "start": 172000,   # 0-based
        "end":   178000,
        "normal_cn": 4,
        "label": "HBA1/HBA2",
    },
    "sma": {
        "chrom": "chr5",
        "start": 70920000,
        "end":   70970000,
        "normal_cn": 4,   # 2×SMN1 + 2×SMN2
        "label": "SMN1/SMN2",
        "c840_pos0": 70951945,  # 0-based
        "smn1_base": "C",
        "smn2_base": "T",
    },
    "gaucher": {
        "chrom": "chr1",
        "start": 155210000,
        "end":   155245000,
        "normal_cn": 2,
        "label": "GBA1/GBAP1",
    },
}

SCENARIOS = {
    "alpha_thal": {
        "carrier":  3,
        "affected": 2,
        "hbh":      1,
        "hydrops":  0,
    },
    "sma": {
        # (smn1_copies, smn2_copies)
        "carrier":  (1, 2),
        "affected": (0, 4),
    },
    "gaucher": {
        "carrier":  1,
        "affected": 0,
    },
}


# ── Helpers ───────────────────────────────────────────────────────────────────

def get_aligned_base(read: pysam.AlignedSegment, ref_pos0: int) -> tuple[int, str]:
    """Return (query_index, base) for a given 0-based ref position, or (-1, '')."""
    try:
        pairs = read.get_aligned_pairs(matches_only=True)
    except Exception:
        return -1, ""
    for qpos, rpos in pairs:
        if rpos == ref_pos0:
            if qpos is None:
                return -1, ""
            seq = read.query_sequence
            if seq and 0 <= qpos < len(seq):
                return qpos, seq[qpos].upper()
            return -1, ""
    return -1, ""


def subsample_reads_in_region(
    bam_in: pysam.AlignmentFile,
    chrom: str,
    start: int,
    end: int,
    n_remove: int,
    rng: random.Random,
) -> Set[str]:
    """Return a set of (query_name, is_read1) pairs to exclude."""
    candidates = []
    for read in bam_in.fetch(chrom, start, end):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        candidates.append((read.query_name, read.is_read1))
    if n_remove >= len(candidates):
        print(f"  WARNING: n_remove={n_remove} >= available={len(candidates)}; removing all",
              file=sys.stderr)
        return set(map(tuple, candidates))
    to_remove = rng.sample(candidates, n_remove)
    return set(map(tuple, to_remove))


# ── Simulation modes ──────────────────────────────────────────────────────────

def simulate_depth_reduction(
    bam_in_path: str,
    bam_out_path: str,
    region: dict,
    fetal_cn: int,
    fetal_fraction: float,
    threads: int,
    rng: random.Random,
) -> None:
    """Remove reads from dark gene window to simulate fetal copy loss."""

    chrom  = region["chrom"]
    start  = region["start"]
    end    = region["end"]
    norm_cn = region["normal_cn"]
    label  = region["label"]

    with pysam.AlignmentFile(bam_in_path, "rb", threads=threads) as bam:
        region_reads = sum(
            1 for r in bam.fetch(chrom, start, end)
            if not (r.is_unmapped or r.is_secondary or r.is_supplementary)
        )

    print(f"  Region reads in {chrom}:{start}-{end} ({label}): {region_reads:,}",
          file=sys.stderr)

    # How many reads to remove
    # observed_ratio = 1.0×(1-FF) + (fetal_cn/norm_cn)×FF
    # depth_factor   = observed_ratio = 1 - FF×(1 - fetal_cn/norm_cn)
    depth_factor = 1.0 - fetal_fraction * (1.0 - fetal_cn / norm_cn)
    n_remove = round(region_reads * (1.0 - depth_factor))
    print(f"  Fetal CN {fetal_cn}/{norm_cn}, FF={fetal_fraction:.2f} "
          f"→ depth_factor={depth_factor:.4f}, removing {n_remove} reads",
          file=sys.stderr)

    with pysam.AlignmentFile(bam_in_path, "rb", threads=threads) as bam:
        to_remove = subsample_reads_in_region(bam, chrom, start, end, n_remove, rng)

    print(f"  Selected {len(to_remove)} reads to remove", file=sys.stderr)

    tmp_unsorted = bam_out_path + ".tmp_unsorted.bam"
    try:
        with pysam.AlignmentFile(bam_in_path, "rb", threads=threads) as src, \
             pysam.AlignmentFile(tmp_unsorted, "wb",
                                 header=src.header.to_dict(), threads=threads) as dst:
            kept = removed = 0
            for read in src:
                key = (read.query_name, read.is_read1)
                if key in to_remove:
                    removed += 1
                    to_remove.discard(key)  # remove only once per read
                else:
                    dst.write(read)
                    kept += 1

        print(f"  Kept {kept:,} / Removed {removed:,} reads", file=sys.stderr)
        print(f"  Sorting …", file=sys.stderr)
        pysam.sort("-o", bam_out_path, "-@", str(threads), tmp_unsorted)
        print(f"  Indexing …", file=sys.stderr)
        pysam.index(bam_out_path, "-@", str(threads))
    finally:
        if os.path.exists(tmp_unsorted):
            os.remove(tmp_unsorted)


def simulate_smn(
    bam_in_path: str,
    bam_out_path: str,
    smn1_copies: int,
    smn2_copies: int,
    fetal_fraction: float,
    threads: int,
    rng: random.Random,
) -> None:
    """Modify C→T at SMN1 c.840 position to simulate SMN1 copy loss in fetus."""

    region  = DARK_GENE_REGIONS["sma"]
    chrom   = region["chrom"]
    pos0    = region["c840_pos0"]  # 0-based

    total_smn = smn1_copies + smn2_copies
    fetal_smn1_frac = smn1_copies / total_smn if total_smn > 0 else 0.0
    maternal_smn1_frac = 0.5  # assumed normal maternal

    # expected observed fraction
    # obs = maternal*(1-FF) + fetal*FF
    expected_obs = maternal_smn1_frac * (1 - fetal_fraction) + fetal_smn1_frac * fetal_fraction
    delta = maternal_smn1_frac - expected_obs  # fraction of C reads to change to T

    print(f"  SMN1 copies: {smn1_copies}, SMN2 copies: {smn2_copies}", file=sys.stderr)
    print(f"  Fetal SMN1 frac: {fetal_smn1_frac:.3f}, FF={fetal_fraction:.2f}",
          file=sys.stderr)
    print(f"  Expected observed frac: {expected_obs:.3f} (Δ={delta:.4f})",
          file=sys.stderr)

    # Collect reads overlapping c.840
    with pysam.AlignmentFile(bam_in_path, "rb", threads=threads) as bam:
        c_reads = []   # reads showing C (SMN1-like)
        for read in bam.fetch(chrom, pos0, pos0 + 1):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            qpos, base = get_aligned_base(read, pos0)
            if base == "C" and qpos >= 0:
                c_reads.append((read.query_name, read.is_read1, qpos))

    total_c = len(c_reads)
    n_to_change = round(total_c * delta / max(maternal_smn1_frac, 1e-6))
    print(f"  C reads at c.840: {total_c}, changing {n_to_change} C→T", file=sys.stderr)

    if n_to_change <= 0:
        print("  Nothing to change; copying BAM as-is", file=sys.stderr)
        pysam.sort("-o", bam_out_path, "-@", str(threads), bam_in_path)
        pysam.index(bam_out_path)
        return

    selected = rng.sample(c_reads, min(n_to_change, total_c))
    change_keys = {(name, is_r1): qpos for name, is_r1, qpos in selected}

    tmp_unsorted = bam_out_path + ".tmp_unsorted.bam"
    try:
        modified = 0
        with pysam.AlignmentFile(bam_in_path, "rb", threads=threads) as src, \
             pysam.AlignmentFile(tmp_unsorted, "wb",
                                 header=src.header.to_dict(), threads=threads) as dst:
            for read in src:
                key = (read.query_name, read.is_read1)
                if key in change_keys:
                    qpos = change_keys[key]
                    seq = bytearray(read.query_sequence.encode())
                    seq[qpos] = ord("T")
                    read.query_sequence = seq.decode()
                    modified += 1
                dst.write(read)

        print(f"  Modified {modified} reads (C→T at c.840)", file=sys.stderr)
        pysam.sort("-o", bam_out_path, "-@", str(threads), tmp_unsorted)
        pysam.index(bam_out_path)
    finally:
        if os.path.exists(tmp_unsorted):
            os.remove(tmp_unsorted)


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Simulate dark gene dosage/allele changes in cfDNA BAM."
    )
    parser.add_argument("--bam",             required=True, help="Input BAM")
    parser.add_argument("--disease",         required=True,
                        choices=["alpha_thal", "sma", "gaucher"])
    parser.add_argument("--scenario",        required=True,
                        help="carrier / affected / hbh / hydrops  (see --help)")
    parser.add_argument("--fetal-fraction",  required=True, type=float)
    parser.add_argument("--output",          required=True, help="Output BAM path")
    parser.add_argument("--threads",         type=int, default=4)
    parser.add_argument("--seed",            type=int, default=42)
    args = parser.parse_args()

    if not os.path.exists(args.bam):
        sys.exit(f"ERROR: BAM not found: {args.bam}")

    rng = random.Random(args.seed)

    print(f"Input  : {args.bam}", file=sys.stderr)
    print(f"Disease: {args.disease}  Scenario: {args.scenario}", file=sys.stderr)
    print(f"FF     : {args.fetal_fraction:.3f}", file=sys.stderr)
    print(f"Output : {args.output}", file=sys.stderr)

    disease   = args.disease
    scenario  = args.scenario
    ff        = args.fetal_fraction

    valid_scenarios = list(SCENARIOS[disease].keys())
    if scenario not in valid_scenarios:
        sys.exit(f"ERROR: '{scenario}' not valid for {disease}. Choose: {valid_scenarios}")

    if disease in ("alpha_thal", "gaucher"):
        fetal_cn = SCENARIOS[disease][scenario]
        region   = DARK_GENE_REGIONS[disease]
        print(f"\nSimulating {region['label']} — fetal CN={fetal_cn} "
              f"(normal={region['normal_cn']})", file=sys.stderr)
        simulate_depth_reduction(
            args.bam, args.output, region, fetal_cn, ff, args.threads, rng
        )

    elif disease == "sma":
        smn1, smn2 = SCENARIOS["sma"][scenario]
        print(f"\nSimulating SMN1/SMN2 — fetal: {smn1}×SMN1 + {smn2}×SMN2",
              file=sys.stderr)
        simulate_smn(
            args.bam, args.output, smn1, smn2, ff, args.threads, rng
        )

    print("\nDone.", file=sys.stderr)


if __name__ == "__main__":
    main()
