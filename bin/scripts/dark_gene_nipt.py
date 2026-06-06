#!/usr/bin/env python3
"""
dark_gene_nipt.py — Dark Gene CNV Analysis for Single Gene NIPT

Estimates copy-number dosage for pseudogene/paralog-heavy regions
(SMN1/2, HBA1/2, GBA1/GBAP1) from cfDNA BAM data.

Key principle — cfDNA mixing model:
  observed_ratio = maternal_ratio * (1 - FF) + fetal_ratio * FF

  Where ratio = region_mean_depth / autosome_background_mean_depth

For a diploid autosome background:
  normal region ratio ≈ 1.0 (if region has same copy number as background)

Adapted from gx-exome dark gene analysis (fallback.nf, coverage.nf)
with cfDNA-specific fetal fraction correction.

Outputs:
  <sample_id>.dark_gene_report.json   — structured results
  <sample_id>.dark_gene_depth.tsv     — per-region depth table
"""
from __future__ import annotations

import argparse
import json
import logging
import math
import os
import subprocess
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)

# ─────────────────────────────────────────────────────────────────────────────
# Region definitions (GRCh38 / hg38)
# ─────────────────────────────────────────────────────────────────────────────
DARK_GENE_REGIONS = {
    "SMN1_SMN2": {
        "chrom": "chr5",
        "start": 70920000,
        "end": 70970000,
        "name": "SMN1_SMN2_Region",
        "gene": "SMN1/SMN2",
        "disease": "Spinal Muscular Atrophy (SMA)",
        "normal_copy_number": 4,   # 2×SMN1 + 2×SMN2 (diploid combined)
        "interpretation": "smn",
    },
    "GBA1_GBAP1": {
        "chrom": "chr1",
        "start": 155210000,
        "end": 155245000,
        "name": "GBA1_GBAP1_Region",
        "gene": "GBA1/GBAP1",
        "disease": "Gaucher Disease",
        "normal_copy_number": 2,
        "interpretation": "depth_ratio",
    },
    "HBA1_HBA2": {
        "chrom": "chr16",
        "start": 172000,
        "end": 178000,
        "name": "HBA1_HBA2_Region",
        "gene": "HBA1/HBA2",
        "disease": "Alpha Thalassemia",
        "normal_copy_number": 4,   # 2×HBA1 + 2×HBA2 (diploid alpha globin locus)
        "interpretation": "hba",
    },
}

# SMN1-specific discriminating position: c.840 (g.27134 in SMN1) — GRCh38 1-based
# SMN1: chr5:70951946 (c.840 C→T change not present in SMN2, i.e. SMN1=C, SMN2=T)
# This is the position that distinguishes SMN1 from SMN2 reads
SMN1_C840_POS = 70951946   # 1-based GRCh38 position
SMN2_C840_POS = 70076526   # SMN2 corresponding position (mirror)

# Reference allele at SMN1 c.840 position (SMN1-specific = C; SMN2 = T)
SMN1_SPECIFIC_BASE = "C"
SMN2_SPECIFIC_BASE = "T"


# ─────────────────────────────────────────────────────────────────────────────
# Data classes
# ─────────────────────────────────────────────────────────────────────────────
@dataclass
class RegionDepth:
    name: str
    chrom: str
    start: int
    end: int
    total_bases: int
    covered_bases: int
    mean_depth: float
    median_depth: float


@dataclass
class DarkGeneResult:
    region: str
    gene: str
    disease: str
    observed_ratio: float          # region_depth / background_depth
    background_depth: float
    region_depth: float
    fetal_fraction: float
    ff_confidence: str
    # FF-aware interpretation
    expected_ratio_normal: float   # expected if fetus is normal
    expected_ratio_1del: float     # expected if fetus has 1-copy deletion
    expected_ratio_2del: float     # expected if fetus has 2-copy deletion
    expected_ratio_affected: float # expected if fetus fully affected
    risk_classification: str       # NORMAL / LOW_RISK / MODERATE_RISK / HIGH_RISK
    interpretation_notes: str
    # SMN-specific fields (only for SMN1/2)
    smn1_c840_depth: int = 0
    smn1_c840_c_count: int = 0
    smn1_c840_t_count: int = 0
    smn1_frac: float = 0.0         # SMN1 fraction at c.840 (C / (C+T))
    smn1_expected_normal: float = 0.0
    smn1_flag: str = ""


# ─────────────────────────────────────────────────────────────────────────────
# samtools helpers
# ─────────────────────────────────────────────────────────────────────────────
def run_cmd(cmd: List[str], check: bool = True) -> subprocess.CompletedProcess:
    logger.debug("CMD: %s", " ".join(cmd))
    return subprocess.run(cmd, capture_output=True, text=True, check=check)


def samtools_bedcov(bed_path: str, bam_path: str) -> Dict[str, float]:
    """
    Run samtools bedcov and return {region_name: mean_depth}.
    bedcov outputs sum of read bases per region; divide by length for mean.
    """
    result = run_cmd(["samtools", "bedcov", bed_path, bam_path])
    depths = {}
    for line in result.stdout.splitlines():
        parts = line.strip().split("\t")
        if len(parts) < 4:
            continue
        chrom, start, end = parts[0], int(parts[1]), int(parts[2])
        name = parts[3] if len(parts) > 3 else f"{chrom}:{start}-{end}"
        read_sum = int(parts[-1])
        length = end - start
        mean = read_sum / length if length > 0 else 0.0
        depths[name] = mean
    return depths


def samtools_idxstats(bam_path: str) -> Dict[str, Tuple[int, int]]:
    """Return {chrom: (length, mapped_reads)} from samtools idxstats."""
    result = run_cmd(["samtools", "idxstats", bam_path])
    stats = {}
    for line in result.stdout.splitlines():
        parts = line.strip().split("\t")
        if len(parts) < 4:
            continue
        chrom, length, mapped, unmapped = parts[0], int(parts[1]), int(parts[2]), int(parts[3])
        if length > 0:
            stats[chrom] = (length, mapped)
    return stats


def compute_autosome_background_depth(
    bam_path: str,
    exclude_chroms: Optional[List[str]] = None,
    target_bed: Optional[str] = None,
) -> float:
    """
    Compute mean depth across autosome background regions.

    Two modes:
    1. target_bed provided (panel/capture BAM):
       Use samtools bedcov on target BED regions (excluding dark gene chroms).
       This gives per-base mean depth that is directly comparable to dark gene
       bedcov output.  Required for capture-based panel BAMs where the read
       distribution is non-uniform genome-wide.

    2. target_bed not provided (WGS or low-pass cfDNA):
       Use samtools idxstats (reads/chrom_length × read_length).
       Only valid when reads are approximately uniformly distributed.
    """
    skip_chroms = {"chrX", "chrY", "X", "Y", "chrM", "MT", "chrMT", "*"}
    if exclude_chroms:
        skip_chroms.update(exclude_chroms)

    # ── Panel mode: bedcov on target BED, skip dark gene chromosomes ──────────
    if target_bed:
        total_sum  = 0
        total_bp   = 0
        with open(target_bed) as fh:
            for raw in fh:
                raw = raw.strip()
                if not raw or raw.startswith("#"):
                    continue
                parts = raw.split("\t")
                if len(parts) < 3:
                    continue
                chrom, start, end = parts[0], int(parts[1]), int(parts[2])
                if chrom in skip_chroms:
                    continue
                if "_" in chrom or chrom.startswith("chrUn") or chrom.startswith("HLA"):
                    continue
                total_bp += end - start

        if total_bp == 0:
            logger.warning("target_bed has no usable autosome regions after filtering; "
                           "falling back to idxstats background.")
        else:
            # Write a temporary filtered BED and run bedcov once
            import tempfile
            with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as tmp:
                tmp_path = tmp.name
                with open(target_bed) as fh:
                    for raw in fh:
                        raw = raw.strip()
                        if not raw or raw.startswith("#"):
                            continue
                        parts = raw.split("\t")
                        if len(parts) < 3:
                            continue
                        chrom, start, end = parts[0], int(parts[1]), int(parts[2])
                        if chrom in skip_chroms:
                            continue
                        if "_" in chrom or chrom.startswith("chrUn") or chrom.startswith("HLA"):
                            continue
                        tmp.write(f"{chrom}\t{start}\t{end}\n")

            result = run_cmd(["samtools", "bedcov", tmp_path, bam_path], check=False)
            os.unlink(tmp_path)

            if result.returncode == 0:
                total_cov = 0
                total_len = 0
                for line in result.stdout.splitlines():
                    cols = line.strip().split("\t")
                    if len(cols) < 4:
                        continue
                    length = int(cols[2]) - int(cols[1])
                    cov_sum = int(cols[-1])
                    total_cov += cov_sum
                    total_len += length
                if total_len > 0:
                    bg = total_cov / total_len
                    logger.info("Background depth from target BED bedcov: %.2f×", bg)
                    return bg

            logger.warning("bedcov on target BED failed; falling back to idxstats.")

    # ── WGS / fallback mode: idxstats ─────────────────────────────────────────
    # Returns reads/bp which is then multiplied by an assumed read length (~150)
    # during output formatting.  The ratio region_depth/bg_depth is only valid
    # here when both are in the same depth units — so we multiply bg by the
    # effective read length derived from the first chromosome with meaningful coverage.
    stats = samtools_idxstats(bam_path)

    total_reads  = 0
    total_length = 0
    total_cov    = 0   # sum of (mapped * read_len) as proxy for total base coverage
    READ_LEN_PROXY = 150

    for chrom, (length, mapped) in stats.items():
        if chrom in skip_chroms:
            continue
        if "_" in chrom or chrom.startswith("chrUn") or chrom.startswith("HLA"):
            continue
        total_reads  += mapped
        total_length += length
        total_cov    += mapped * READ_LEN_PROXY

    if total_length == 0:
        return 0.0

    bg = total_cov / total_length   # approximate mean depth (×) across genome
    logger.info("Background depth from idxstats (WGS proxy): %.4f×", bg)
    return bg


def mpileup_at_position(
    bam_path: str,
    chrom: str,
    pos_1based: int,
    ref_fasta: Optional[str] = None,
    min_mapq: int = 0,
    min_bq: int = 0,
) -> Dict[str, int]:
    """
    Run samtools mpileup at a single 1-based position.
    Returns {A: n, C: n, G: n, T: n, total: n}.
    """
    region = f"{chrom}:{pos_1based}-{pos_1based}"
    cmd = ["samtools", "mpileup", "-r", region,
           "-q", str(min_mapq), "-Q", str(min_bq),
           "--no-BAQ", "--max-depth", "1000000"]
    if ref_fasta:
        cmd += ["-f", ref_fasta]
    cmd.append(bam_path)

    result = run_cmd(cmd, check=False)
    counts = {"A": 0, "C": 0, "G": 0, "T": 0, "total": 0}

    for line in result.stdout.splitlines():
        parts = line.strip().split("\t")
        if len(parts) < 5:
            continue
        pileup_str = parts[4].upper()
        for base in pileup_str:
            if base in ("A", "C", "G", "T"):
                counts[base] = counts.get(base, 0) + 1
        counts["total"] = sum(counts[b] for b in ("A", "C", "G", "T"))
    return counts


# ─────────────────────────────────────────────────────────────────────────────
# FF-aware interpretation
# ─────────────────────────────────────────────────────────────────────────────
def compute_expected_ratios(
    normal_cn: int,
    deletion_cn_1: int,
    deletion_cn_2: int,
    affected_cn: int,
    fetal_fraction: float,
    maternal_cn: int,
) -> Tuple[float, float, float, float]:
    """
    Compute expected depth ratios for different fetal CNV scenarios.

    Model:
      observed_ratio = maternal_ratio * (1 - FF) + fetal_ratio * FF

    Where ratio = region_cn / diploid_background_cn.
    For diploid background (CN=2), ratio = cn / 2.
    maternal_ratio = maternal_cn / 2
    """
    def cn_to_ratio(cn: int) -> float:
        return cn / 2.0  # relative to diploid background

    mat_ratio = cn_to_ratio(maternal_cn)
    ff = fetal_fraction

    r_normal   = mat_ratio * (1 - ff) + cn_to_ratio(normal_cn)   * ff
    r_1del     = mat_ratio * (1 - ff) + cn_to_ratio(deletion_cn_1) * ff
    r_2del     = mat_ratio * (1 - ff) + cn_to_ratio(deletion_cn_2) * ff
    r_affected = mat_ratio * (1 - ff) + cn_to_ratio(affected_cn)  * ff

    return r_normal, r_1del, r_2del, r_affected


def classify_hba_risk(
    observed: float,
    expected_normal: float,
    expected_1del: float,
    expected_2del: float,
    expected_affected: float,
    fetal_fraction: float,
) -> Tuple[str, str]:
    """
    Classify alpha thalassemia risk from HBA depth ratio.
    Returns (risk_level, notes).

    Scenarios (assuming maternal normal = 4 copies):
      normal (4-copy fetal):   ratio ~ expected_normal    → NORMAL
      -α/αα  (3-copy fetal):   ratio ~ expected_1del      → LOW_RISK (silent carrier)
      --/αα or -α/-α (2-copy): ratio ~ expected_2del      → MODERATE_RISK (trait)
      --/-α  (1-copy):         ratio ~ expected_affected  → HIGH_RISK (HbH risk)
      0 copies:                very low ratio             → CRITICAL
    """
    if fetal_fraction < 0.03:
        return "INDETERMINATE", (
            f"Fetal fraction too low ({fetal_fraction:.1%}) for reliable dark gene CNV call. "
            "Minimum recommended FF ≥ 3%."
        )

    # Tolerance: ±1.5 sigma based on FF-dependent signal
    # Signal per copy = FF / 2 (for diploid background)
    copy_unit = fetal_fraction / 2.0
    tol = copy_unit * 0.8  # 80% of one copy-unit tolerance

    if observed >= expected_normal - tol:
        risk = "NORMAL"
        notes = (
            f"HBA ratio {observed:.3f} consistent with fetal 4-copy (normal). "
            f"Expected normal: {expected_normal:.3f}."
        )
    elif abs(observed - expected_1del) < tol:
        risk = "LOW_RISK"
        notes = (
            f"HBA ratio {observed:.3f} ~ expected 3-copy (-α/αα; e.g. -α3.7/-α4.2). "
            f"Expected 1-deletion: {expected_1del:.3f}. "
            "Possible silent carrier. Confirmatory testing recommended."
        )
    elif abs(observed - expected_2del) < tol or (
        expected_2del - tol <= observed < expected_1del - tol
    ):
        risk = "MODERATE_RISK"
        notes = (
            f"HBA ratio {observed:.3f} ~ expected 2-copy (--/αα or -α/-α). "
            f"Expected 2-deletion: {expected_2del:.3f}. "
            "Alpha thalassemia trait (at risk for HbH offspring). "
            "Confirmatory testing strongly recommended."
        )
    elif observed < expected_2del - tol:
        risk = "HIGH_RISK"
        notes = (
            f"HBA ratio {observed:.3f} below 2-copy threshold. "
            f"Expected HbH-risk (--/-α): {expected_affected:.3f}. "
            "High-risk for alpha thalassemia major / HbH disease. "
            "Urgent confirmatory testing required."
        )
    else:
        risk = "INDETERMINATE"
        notes = (
            f"HBA ratio {observed:.3f} falls between thresholds. "
            "Cannot reliably distinguish copy number with current FF. "
            f"Expected ranges: normal={expected_normal:.3f}, "
            f"3-copy={expected_1del:.3f}, 2-copy={expected_2del:.3f}."
        )
    return risk, notes


def classify_smn_risk(
    smn1_frac: float,
    smn1_depth: int,
    fetal_fraction: float,
    depth_ratio: float,
    expected_normal: float,
) -> Tuple[str, str]:
    """
    Classify SMA risk from SMN1/SMN2 c.840 pileup in cfDNA.

    In cfDNA (mixed maternal+fetal):
      measured_smn1_frac = maternal_smn1_frac * (1 - FF) + fetal_smn1_frac * FF

    Normal germline: SMN1_frac ~ 0.5 (2 SMN1 out of 4 total SMN copies)
    Carrier (1 SMN1 copy): SMN1_frac ~ 0.25
    Affected (0 SMN1 copy): SMN1_frac ~ 0.0

    In cfDNA, if fetal is affected (0 copies) and mother is carrier (1 copy):
      measured = 0.25 * (1 - FF) + 0 * FF = 0.25 * (1 - FF)
    """
    if smn1_depth < 20:
        return "INSUFFICIENT_DATA", (
            f"Insufficient reads at SMN1 c.840 position (depth={smn1_depth}). "
            "Minimum 20 reads required for SMN analysis."
        )

    if fetal_fraction < 0.03:
        return "INDETERMINATE", (
            f"Fetal fraction too low ({fetal_fraction:.1%}). "
            "SMN1 analysis requires FF ≥ 3%."
        )

    ff = fetal_fraction
    # Expected SMN1 fraction under different scenarios
    # Assuming mother is normal (SMN1_frac = 0.5):
    exp_normal        = 0.5 * (1 - ff) + 0.50 * ff  # fetal normal (2 copies)
    exp_carrier_fetal = 0.5 * (1 - ff) + 0.25 * ff  # fetal carrier (1 copy)
    exp_affected_fetal= 0.5 * (1 - ff) + 0.00 * ff  # fetal affected (0 copies)
    # Assuming mother is carrier:
    exp_carrier_mat_normal_fetal = 0.25 * (1 - ff) + 0.50 * ff
    exp_carrier_mat_carrier_fetal= 0.25 * (1 - ff) + 0.25 * ff
    exp_carrier_mat_affected_fetal=0.25 * (1 - ff) + 0.00 * ff

    # Tolerance: half copy-unit (1 copy changes frac by ~0.25 for germline)
    tol = (0.25 * ff) * 0.8

    if smn1_frac >= exp_normal - tol:
        risk = "NORMAL"
        notes = (
            f"SMN1 fraction at c.840: {smn1_frac:.3f} (depth={smn1_depth}). "
            f"Consistent with normal fetal SMN1 copies. "
            f"Expected normal: {exp_normal:.3f}."
        )
    elif abs(smn1_frac - exp_carrier_fetal) < tol:
        risk = "MODERATE_RISK"
        notes = (
            f"SMN1 fraction {smn1_frac:.3f} ~ carrier pattern "
            f"(expected {exp_carrier_fetal:.3f} if fetal has 1 SMN1 copy). "
            "Possible fetal SMA carrier. Confirmatory testing recommended."
        )
    elif smn1_frac <= exp_affected_fetal + tol:
        risk = "HIGH_RISK"
        notes = (
            f"SMN1 fraction {smn1_frac:.3f} very low, near expected 0-copy pattern "
            f"(expected {exp_affected_fetal:.3f} if fetal affected). "
            "High-risk for SMA. Urgent confirmatory testing required."
        )
    else:
        risk = "INDETERMINATE"
        notes = (
            f"SMN1 fraction {smn1_frac:.3f} falls in ambiguous range. "
            f"Expected: normal={exp_normal:.3f}, carrier={exp_carrier_fetal:.3f}, "
            f"affected={exp_affected_fetal:.3f}. "
            "Consider higher-depth re-sequencing or orthogonal method."
        )
    return risk, notes


def classify_gba_risk(
    observed: float,
    expected_normal: float,
    fetal_fraction: float,
) -> Tuple[str, str]:
    """
    GBA1 depth ratio classification.
    GBA1/GBAP1 homology makes depth interpretation unreliable for point variants,
    but gross deletion can be detected by depth.
    """
    if fetal_fraction < 0.03:
        return "INDETERMINATE", (
            "Fetal fraction too low for GBA1 depth analysis."
        )

    copy_unit = fetal_fraction / 2.0
    tol = copy_unit * 1.2

    if observed >= expected_normal - tol:
        risk = "NORMAL"
        notes = (
            f"GBA1/GBAP1 depth ratio {observed:.3f} consistent with normal copy number. "
            "NOTE: Depth analysis cannot detect point mutations or small indels; "
            "these require specific variant calling."
        )
    elif observed >= expected_normal - 2 * tol:
        risk = "LOW_RISK"
        notes = (
            f"GBA1/GBAP1 ratio {observed:.3f} slightly below normal. "
            "Possible GBA1 dosage reduction. "
            "NOTE: GBA1/GBAP1 pseudogene homology reduces specificity; "
            "confirmatory testing with long-read or gene-conversion assay recommended."
        )
    else:
        risk = "MODERATE_RISK"
        notes = (
            f"GBA1/GBAP1 ratio {observed:.3f} significantly below normal. "
            "Possible fetal GBA1 deletion or gene conversion. "
            "Confirmatory testing strongly recommended."
        )
    return risk, notes


# ─────────────────────────────────────────────────────────────────────────────
# Main analysis
# ─────────────────────────────────────────────────────────────────────────────
def analyse_dark_genes(
    bam_path: str,
    dark_gene_bed: str,
    ff_json_path: str,
    sample_id: str,
    ref_fasta: Optional[str] = None,
    min_background_depth: float = 5.0,
    target_bed: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Main entry point.
    Returns full analysis dict (for JSON output).
    """
    # ── Load fetal fraction ──────────────────────────────────────────────────
    with open(ff_json_path) as fh:
        ff_data = json.load(fh)

    primary_ff = ff_data.get("primary_fetal_fraction")
    ff_confidence = ff_data.get("ff_confidence", "unknown")

    if primary_ff is None:
        logger.warning("Fetal fraction not available; using default 0.05 for thresholds")
        primary_ff = 0.05
        ff_confidence = "unavailable"
    fetal_fraction = float(primary_ff)

    logger.info("Fetal fraction: %.4f (confidence: %s)", fetal_fraction, ff_confidence)

    # ── Compute background depth ─────────────────────────────────────────────
    logger.info("Computing autosome background depth...")
    # Exclude dark gene chromosomes (chr1/chr5/chr16) to avoid circular bias.
    # target_bed (panel BED) is strongly preferred for capture-based panels —
    # bedcov on target regions gives depth in the same units as dark gene bedcov.
    bg_depth = compute_autosome_background_depth(
        bam_path,
        exclude_chroms=["chr5", "chr1", "chr16"],
        target_bed=target_bed,
    )
    if bg_depth < 1e-9:
        logger.error("Background depth is zero. BAM may be empty or wrong chromosome names.")
        bg_depth = 1e-9

    # ── bedcov per region ────────────────────────────────────────────────────
    logger.info("Running samtools bedcov on dark gene regions...")
    region_depths = samtools_bedcov(dark_gene_bed, bam_path)
    logger.info("Region depths: %s", region_depths)

    # ── SMN1/SMN2 c.840 pileup ───────────────────────────────────────────────
    logger.info("Pileup at SMN1 c.840 discriminating position (chr5:%d)...", SMN1_C840_POS)
    c840_counts = mpileup_at_position(
        bam_path, "chr5", SMN1_C840_POS,
        ref_fasta=ref_fasta, min_mapq=0, min_bq=0,
    )
    logger.info("SMN1 c.840 counts: %s", c840_counts)

    smn1_c_count = c840_counts.get("C", 0)
    smn1_t_count = c840_counts.get("T", 0)
    smn1_total   = smn1_c_count + smn1_t_count
    smn1_frac    = smn1_c_count / smn1_total if smn1_total > 0 else 0.0

    # ── Per-region analysis ──────────────────────────────────────────────────
    results: List[DarkGeneResult] = []

    for region_key, region_meta in DARK_GENE_REGIONS.items():
        region_name = region_meta["name"]
        region_depth = region_depths.get(region_name)

        if region_depth is None:
            logger.warning("No depth found for region %s (name not in bedcov output). "
                           "Available keys: %s", region_name, list(region_depths.keys()))
            # Try partial match
            for key in region_depths:
                if region_key.split("_")[0] in key:
                    region_depth = region_depths[key]
                    break
            if region_depth is None:
                region_depth = 0.0

        normal_cn  = region_meta["normal_copy_number"]
        # Assume maternal is normal diploid unless specified otherwise
        maternal_cn = normal_cn

        # Expected ratios (background is diploid = 2 copies)
        r_normal, r_1del, r_2del, r_affected = compute_expected_ratios(
            normal_cn=normal_cn,
            deletion_cn_1=normal_cn - 1,
            deletion_cn_2=normal_cn - 2,
            affected_cn=max(0, normal_cn - normal_cn),  # fully affected = 0
            fetal_fraction=fetal_fraction,
            maternal_cn=maternal_cn,
        )

        observed_ratio = (region_depth / bg_depth) if bg_depth > 0 else 0.0

        interp = region_meta["interpretation"]

        if interp == "hba":
            risk, notes = classify_hba_risk(
                observed_ratio, r_normal,
                r_normal - (fetal_fraction / 2),    # -1 copy
                r_normal - fetal_fraction,            # -2 copies
                r_normal - (3 * fetal_fraction / 2), # -3 copies (HbH)
                fetal_fraction,
            )
            smn_fields = {}
        elif interp == "smn":
            # Primary: c.840 pileup
            risk, notes = classify_smn_risk(
                smn1_frac, smn1_total, fetal_fraction,
                observed_ratio, r_normal,
            )
            smn_fields = {
                "smn1_c840_depth": smn1_total,
                "smn1_c840_c_count": smn1_c_count,
                "smn1_c840_t_count": smn1_t_count,
                "smn1_frac": round(smn1_frac, 4),
                "smn1_expected_normal": round(
                    0.5 * (1 - fetal_fraction) + 0.5 * fetal_fraction, 4
                ),
                "smn1_flag": (
                    "LOW_SMN1_FRAC" if smn1_frac < 0.5 * (1 - fetal_fraction) - 0.05
                    else "NORMAL_SMN1_FRAC"
                ),
            }
        else:  # depth_ratio / GBA1
            risk, notes = classify_gba_risk(
                observed_ratio, r_normal, fetal_fraction,
            )
            smn_fields = {}

        res = DarkGeneResult(
            region=region_key,
            gene=region_meta["gene"],
            disease=region_meta["disease"],
            observed_ratio=round(observed_ratio, 4),
            background_depth=round(bg_depth, 2),
            region_depth=round(region_depth, 2),
            fetal_fraction=round(fetal_fraction, 4),
            ff_confidence=ff_confidence,
            expected_ratio_normal=round(r_normal, 4),
            expected_ratio_1del=round(r_normal - fetal_fraction / 2, 4),
            expected_ratio_2del=round(r_normal - fetal_fraction, 4),
            expected_ratio_affected=round(r_normal - 3 * fetal_fraction / 2, 4),
            risk_classification=risk,
            interpretation_notes=notes,
            **smn_fields,
        )
        results.append(res)

    # ── Build output dict ────────────────────────────────────────────────────
    has_high_risk   = any(r.risk_classification == "HIGH_RISK"     for r in results)
    has_moderate    = any(r.risk_classification == "MODERATE_RISK" for r in results)
    has_indeterminate = any(r.risk_classification in ("INDETERMINATE", "INSUFFICIENT_DATA") for r in results)

    if has_high_risk:
        overall_status = "HIGH_RISK"
    elif has_moderate:
        overall_status = "MODERATE_RISK"
    elif has_indeterminate:
        overall_status = "INDETERMINATE"
    else:
        overall_status = "NORMAL"

    output = {
        "sample_id": sample_id,
        "module": "dark_gene_nipt",
        "fetal_fraction": round(fetal_fraction, 4),
        "ff_confidence": ff_confidence,
        "overall_status": overall_status,
        "caveat": (
            "Dark gene CNV analysis is based on read-depth ratios from cfDNA. "
            "Results are SCREENING-ONLY and require confirmatory diagnostic testing "
            "(MLPA, long-read sequencing, or gene-specific assay) for clinical decisions. "
            "SMN/HBA/GBA1 pseudogene homology reduces specificity of short-read depth analysis. "
            f"Analysis uses fetal fraction = {fetal_fraction:.1%} "
            f"(confidence: {ff_confidence})."
        ),
        "regions": [asdict(r) for r in results],
        "summary": {
            "total_regions": len(results),
            "normal": sum(1 for r in results if r.risk_classification == "NORMAL"),
            "low_risk": sum(1 for r in results if r.risk_classification == "LOW_RISK"),
            "moderate_risk": sum(1 for r in results if r.risk_classification == "MODERATE_RISK"),
            "high_risk": sum(1 for r in results if r.risk_classification == "HIGH_RISK"),
            "indeterminate": sum(1 for r in results if r.risk_classification in (
                "INDETERMINATE", "INSUFFICIENT_DATA",
            )),
        },
    }
    return output


def write_depth_tsv(results: List[DarkGeneResult], out_path: str) -> None:
    """Write per-region depth summary TSV."""
    header = [
        "region", "gene", "disease",
        "observed_ratio", "expected_ratio_normal",
        "expected_ratio_1del", "expected_ratio_2del",
        "fetal_fraction", "ff_confidence",
        "risk_classification", "smn1_frac", "smn1_c840_depth",
    ]
    with open(out_path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in results:
            row = [
                r.region, r.gene, r.disease,
                f"{r.observed_ratio:.4f}", f"{r.expected_ratio_normal:.4f}",
                f"{r.expected_ratio_1del:.4f}", f"{r.expected_ratio_2del:.4f}",
                f"{r.fetal_fraction:.4f}", r.ff_confidence,
                r.risk_classification,
                f"{r.smn1_frac:.4f}" if hasattr(r, "smn1_frac") else "NA",
                str(r.smn1_c840_depth) if hasattr(r, "smn1_c840_depth") else "NA",
            ]
            fh.write("\t".join(row) + "\n")


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="Dark Gene CNV Analysis for Single Gene NIPT (cfDNA)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--bam",           required=True, help="Input BAM file (dedup, pre-target-extract recommended)")
    parser.add_argument("--dark-gene-bed", required=True, help="BED file with dark gene windows (SMN/HBA/GBA)")
    parser.add_argument("--ff-json",       required=True, help="Fetal fraction JSON from ESTIMATE_FETAL_FRACTION")
    parser.add_argument("--sample-id",     required=True, help="Sample identifier")
    parser.add_argument("--ref-fasta",     default=None,  help="Reference FASTA (optional, improves mpileup quality)")
    parser.add_argument("--target-bed",    default=None,
                        help="Panel target BED (recommended for capture panels). "
                             "Used for background depth via bedcov — more accurate than idxstats "
                             "when reads are non-uniformly distributed across the genome.")
    parser.add_argument("--output-json",   required=True, help="Output JSON report path")
    parser.add_argument("--output-tsv",    required=True, help="Output depth TSV path")
    parser.add_argument("--min-ff",        type=float, default=0.03,
                        help="Minimum FF for reliable analysis (default: 0.03)")
    args = parser.parse_args()

    logger.info("=== Dark Gene NIPT Analysis: %s ===", args.sample_id)

    result = analyse_dark_genes(
        bam_path=args.bam,
        dark_gene_bed=args.dark_gene_bed,
        ff_json_path=args.ff_json,
        sample_id=args.sample_id,
        ref_fasta=args.ref_fasta,
        min_background_depth=args.min_ff,
        target_bed=args.target_bed,
    )

    Path(args.output_json).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output_json, "w") as fh:
        json.dump(result, fh, indent=2)
    logger.info("Report written: %s", args.output_json)

    regions = [DarkGeneResult(**r) for r in result["regions"]]
    write_depth_tsv(regions, args.output_tsv)
    logger.info("Depth TSV written: %s", args.output_tsv)

    # ── Print summary ──────────────────────────────────────────────────────
    print(f"\n{'='*65}")
    print(f"DARK GENE NIPT — {args.sample_id}")
    print(f"{'='*65}")
    print(f"Fetal Fraction  : {result['fetal_fraction']:.1%} ({result['ff_confidence']})")
    print(f"Overall Status  : {result['overall_status']}")
    for r in result["regions"]:
        smn_note = ""
        if r.get("smn1_frac", 0) > 0:
            smn_note = f" | SMN1_frac={r['smn1_frac']:.3f}"
        print(
            f"  {r['gene']:20s} | ratio={r['observed_ratio']:.3f}"
            f" (exp_normal={r['expected_ratio_normal']:.3f}) "
            f"| {r['risk_classification']}{smn_note}"
        )
    print(f"\nCAVEAT: {result['caveat'][:100]}...")
    print(f"{'='*65}\n")


if __name__ == "__main__":
    main()
