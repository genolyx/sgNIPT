#!/usr/bin/env python3
"""
fetal_fraction.py - Fetal Fraction Estimation Module for Single Gene NIPT Pipeline

Estimates fetal fraction from cfDNA using SNP-based methods.
Uses informative SNPs where the mother is homozygous and the fetus carries
a different allele (detectable as minor allele in cfDNA).

Methods implemented:
1. Informative SNP method (primary): FF = 2 * mean(minor_AF) at informative loci
2. Y-chromosome method (supplementary, male fetus only)
3. Fragment size distribution method (supplementary)
"""

import argparse
import csv
import json
import logging
import math
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Configuration defaults
# ---------------------------------------------------------------------------
DEFAULT_CONFIG: Dict[str, Any] = {
    # Minimum read depth at a SNP site to be considered
    "min_depth": 50,
    # Minimum minor allele frequency to consider as informative (noise filter)
    "min_minor_af": 0.005,
    # Maximum minor allele frequency for informative SNP (should be < 50%)
    "max_minor_af": 0.25,
    # Minimum number of informative SNPs required for reliable estimation
    "min_informative_snps": 10,
    # Minimum acceptable fetal fraction for clinical use
    "min_fetal_fraction": 0.02,
    # Maximum expected fetal fraction (sanity check)
    "max_fetal_fraction": 0.50,
    # Number of bootstrap iterations for confidence interval
    "bootstrap_iterations": 1000,
    # Outlier removal: number of standard deviations
    "outlier_sd_threshold": 3.0,
    # Minimum base quality for allele counting
    "min_base_quality": 20,
    # Minimum mapping quality for read counting
    "min_mapping_quality": 30,
}


def load_config(config_file: Optional[str]) -> Dict[str, Any]:
    """Load configuration from JSON file, falling back to defaults."""
    config = DEFAULT_CONFIG.copy()
    if config_file and Path(config_file).exists():
        logger.info("Loading fetal fraction config from %s", config_file)
        with open(config_file, "r") as fh:
            custom = json.load(fh)
        config.update(custom)
    return config


def parse_vcf_for_snps(
    vcf_path: str,
    config: Dict[str, Any],
) -> List[Dict[str, Any]]:
    """
    Parse VCF file to extract SNP allele frequencies at target loci.
    Expects a VCF with allele depth (AD) or allele frequency (AF) information.
    """
    logger.info("Parsing VCF for informative SNPs: %s", vcf_path)
    snps = []

    with open(vcf_path, "r") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 10:
                continue

            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]
            qual = fields[5]
            filter_field = fields[6]
            info = fields[7]
            format_field = fields[8]
            sample_data = fields[9]

            # Skip non-SNP variants
            if len(ref) != 1 or any(len(a) != 1 for a in alt.split(",")):
                continue

            # Parse FORMAT and sample fields
            fmt_keys = format_field.split(":")
            fmt_values = sample_data.split(":")
            fmt_dict = dict(zip(fmt_keys, fmt_values))

            # Extract allele depths
            ad = None
            dp = None
            af = None

            if "AD" in fmt_dict:
                try:
                    ad = [int(x) for x in fmt_dict["AD"].split(",")]
                    dp = sum(ad)
                except (ValueError, IndexError):
                    pass

            if "DP" in fmt_dict:
                try:
                    dp = int(fmt_dict["DP"])
                except ValueError:
                    pass

            if "AF" in fmt_dict:
                try:
                    af = float(fmt_dict["AF"].split(",")[0])
                except (ValueError, IndexError):
                    pass

            if dp is None or dp < config["min_depth"]:
                continue

            # Calculate minor allele frequency
            if ad and len(ad) >= 2:
                ref_count = ad[0]
                alt_count = sum(ad[1:])
                minor_af = alt_count / dp if dp > 0 else 0.0
            elif af is not None:
                minor_af = af
                ref_count = int(dp * (1 - af))
                alt_count = dp - ref_count
            else:
                continue

            snps.append({
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "depth": dp,
                "ref_count": ref_count,
                "alt_count": alt_count,
                "minor_af": round(minor_af, 6),
            })

    logger.info("Parsed %d SNP sites from VCF", len(snps))
    return snps


def parse_pileup_for_snps(
    pileup_path: str,
    target_snps_path: str,
    config: Dict[str, Any],
) -> List[Dict[str, Any]]:
    """
    Parse mpileup output at target SNP positions.
    target_snps_path: TSV file with columns: chrom, pos, ref, alt
    """
    logger.info("Parsing pileup for informative SNPs: %s", pileup_path)

    # Load target SNP positions
    target_positions = set()
    target_info = {}
    with open(target_snps_path, "r") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            chrom = row.get("chrom", row.get("chr", row.get("#chrom", "")))
            pos = int(row.get("pos", row.get("position", 0)))
            key = (chrom, pos)
            target_positions.add(key)
            target_info[key] = {
                "ref": row.get("ref", ""),
                "alt": row.get("alt", ""),
            }

    logger.info("Loaded %d target SNP positions", len(target_positions))

    snps = []
    with open(pileup_path, "r") as fh:
        for line in fh:
            fields = line.strip().split("\t")
            if len(fields) < 6:
                continue

            chrom = fields[0]
            pos = int(fields[1])
            key = (chrom, pos)

            if key not in target_positions:
                continue

            ref_base = fields[2].upper()
            depth = int(fields[3])
            bases = fields[4]
            quals = fields[5]

            if depth < config["min_depth"]:
                continue

            # Count alleles from pileup bases string
            allele_counts = count_pileup_alleles(bases, quals, ref_base, config["min_base_quality"])
            total = sum(allele_counts.values())

            if total < config["min_depth"]:
                continue

            info = target_info.get(key, {})
            alt_base = info.get("alt", "").upper()

            ref_count = allele_counts.get(ref_base, 0)
            alt_count = allele_counts.get(alt_base, 0) if alt_base else 0

            # If alt not specified, use the most common non-ref allele
            if not alt_base:
                non_ref = {k: v for k, v in allele_counts.items() if k != ref_base}
                if non_ref:
                    alt_base = max(non_ref, key=non_ref.get)
                    alt_count = non_ref[alt_base]

            minor_af = alt_count / total if total > 0 else 0.0

            snps.append({
                "chrom": chrom,
                "pos": pos,
                "ref": ref_base,
                "alt": alt_base,
                "depth": total,
                "ref_count": ref_count,
                "alt_count": alt_count,
                "minor_af": round(minor_af, 6),
            })

    logger.info("Parsed %d SNP sites from pileup", len(snps))
    return snps


def count_pileup_alleles(
    bases: str, quals: str, ref: str, min_bq: int
) -> Dict[str, int]:
    """Count alleles from a pileup bases string, filtering by base quality."""
    counts: Dict[str, int] = defaultdict(int)
    i = 0
    qi = 0  # quality index

    while i < len(bases):
        c = bases[i]

        if c == "^":
            # Start of read segment, skip mapping quality char
            i += 2
            continue
        elif c == "$":
            # End of read segment
            i += 1
            continue
        elif c in "+-":
            # Indel: skip the indel bases
            i += 1
            num_str = ""
            while i < len(bases) and bases[i].isdigit():
                num_str += bases[i]
                i += 1
            if num_str:
                i += int(num_str)
            continue
        elif c == "*":
            # Deletion placeholder
            i += 1
            qi += 1
            continue
        elif c in ".,ACGTNacgtn":
            # Actual base
            if qi < len(quals):
                bq = ord(quals[qi]) - 33
            else:
                bq = min_bq  # assume passing if quality string is short

            if bq >= min_bq:
                if c in ".,":
                    counts[ref] += 1
                else:
                    counts[c.upper()] += 1

            qi += 1
            i += 1
        else:
            i += 1

    return dict(counts)


def identify_informative_snps(
    snps: List[Dict[str, Any]],
    config: Dict[str, Any],
) -> List[Dict[str, Any]]:
    """
    Identify informative SNPs for fetal fraction estimation.

    Informative SNPs are those where:
    - The mother is likely homozygous for the reference allele
    - The fetus carries a paternal allele (visible as minor allele in cfDNA)
    - Minor AF is in the expected range for fetal contribution
    """
    informative = []

    for snp in snps:
        maf = snp["minor_af"]

        # Filter: minor AF should be in the range expected for fetal alleles
        # For a heterozygous fetal allele in cfDNA: expected AF ≈ FF/2
        # Typical FF range: 2-30%, so expected minor AF: 1-15%
        if maf < config["min_minor_af"]:
            continue
        if maf > config["max_minor_af"]:
            continue

        # The major allele fraction should be high (mother is homozygous)
        # Major AF should be > 0.75 (i.e., 1 - max_minor_af)
        major_af = 1.0 - maf
        if major_af < (1.0 - config["max_minor_af"]):
            continue

        snp["is_informative"] = True
        informative.append(snp)

    logger.info(
        "Identified %d informative SNPs out of %d total",
        len(informative),
        len(snps),
    )
    return informative


def estimate_fetal_fraction_snp(
    informative_snps: List[Dict[str, Any]],
    config: Dict[str, Any],
) -> Dict[str, Any]:
    """
    Estimate fetal fraction using the informative SNP method.

    For each informative SNP where mother is AA and fetus is AB:
    FF_i = 2 * (alt_count / total_depth)

    The overall FF is the median of all FF_i values (robust to outliers).
    """
    if len(informative_snps) < config["min_informative_snps"]:
        logger.warning(
            "Insufficient informative SNPs (%d < %d). Cannot reliably estimate FF.",
            len(informative_snps),
            config["min_informative_snps"],
        )
        return {
            "method": "snp_informative",
            "fetal_fraction": None,
            "confidence_interval": None,
            "num_informative_snps": len(informative_snps),
            "status": "INSUFFICIENT_DATA",
            "message": f"Only {len(informative_snps)} informative SNPs found "
                       f"(minimum {config['min_informative_snps']} required).",
        }

    # Calculate per-SNP fetal fraction estimates
    ff_estimates = []
    for snp in informative_snps:
        ff_i = 2.0 * snp["minor_af"]
        snp["ff_estimate"] = round(ff_i, 6)
        ff_estimates.append(ff_i)

    # Remove outliers using IQR method
    ff_estimates_clean = remove_outliers_iqr(ff_estimates)

    if len(ff_estimates_clean) < config["min_informative_snps"]:
        ff_estimates_clean = ff_estimates  # Fall back to all estimates

    # Calculate statistics
    ff_mean = sum(ff_estimates_clean) / len(ff_estimates_clean)
    ff_sorted = sorted(ff_estimates_clean)
    ff_median = ff_sorted[len(ff_sorted) // 2]

    # Standard deviation
    variance = sum((x - ff_mean) ** 2 for x in ff_estimates_clean) / len(ff_estimates_clean)
    ff_sd = math.sqrt(variance)

    # Standard error
    ff_se = ff_sd / math.sqrt(len(ff_estimates_clean))

    # 95% confidence interval (using normal approximation)
    ci_lower = max(0.0, ff_median - 1.96 * ff_se)
    ci_upper = min(1.0, ff_median + 1.96 * ff_se)

    # Bootstrap confidence interval
    bootstrap_ci = bootstrap_confidence_interval(
        ff_estimates_clean, config["bootstrap_iterations"]
    )

    # Determine status
    if ff_median < config["min_fetal_fraction"]:
        status = "LOW_FF"
        message = (
            f"Estimated fetal fraction ({ff_median:.4f}) is below minimum "
            f"threshold ({config['min_fetal_fraction']:.4f}). Results may be unreliable."
        )
    elif ff_median > config["max_fetal_fraction"]:
        status = "ABNORMAL_HIGH_FF"
        message = (
            f"Estimated fetal fraction ({ff_median:.4f}) is abnormally high "
            f"(>{config['max_fetal_fraction']:.4f}). Check for sample issues."
        )
    else:
        status = "PASS"
        message = "Fetal fraction estimation successful."

    result = {
        "method": "snp_informative",
        "fetal_fraction": round(ff_median, 4),
        "fetal_fraction_mean": round(ff_mean, 4),
        "fetal_fraction_sd": round(ff_sd, 4),
        "fetal_fraction_se": round(ff_se, 4),
        "confidence_interval_95": {
            "lower": round(ci_lower, 4),
            "upper": round(ci_upper, 4),
        },
        "bootstrap_ci_95": {
            "lower": round(bootstrap_ci[0], 4),
            "upper": round(bootstrap_ci[1], 4),
        },
        "num_informative_snps": len(informative_snps),
        "num_snps_after_outlier_removal": len(ff_estimates_clean),
        "per_snp_estimates": ff_estimates,
        "status": status,
        "message": message,
    }

    logger.info(
        "Fetal Fraction: %.4f (median), %.4f (mean), SD=%.4f, 95%% CI=[%.4f, %.4f]",
        ff_median,
        ff_mean,
        ff_sd,
        ci_lower,
        ci_upper,
    )
    logger.info("Status: %s - %s", status, message)
    return result


def remove_outliers_iqr(values: List[float]) -> List[float]:
    """Remove outliers using the IQR method."""
    if len(values) < 4:
        return values

    sorted_vals = sorted(values)
    q1_idx = len(sorted_vals) // 4
    q3_idx = 3 * len(sorted_vals) // 4
    q1 = sorted_vals[q1_idx]
    q3 = sorted_vals[q3_idx]
    iqr = q3 - q1

    lower_bound = q1 - 1.5 * iqr
    upper_bound = q3 + 1.5 * iqr

    cleaned = [v for v in values if lower_bound <= v <= upper_bound]
    removed = len(values) - len(cleaned)
    if removed > 0:
        logger.info("Removed %d outliers using IQR method", removed)
    return cleaned


def bootstrap_confidence_interval(
    values: List[float], n_iterations: int = 1000, ci: float = 0.95
) -> Tuple[float, float]:
    """Calculate bootstrap confidence interval for the median."""
    import random

    if len(values) < 2:
        if values:
            return (values[0], values[0])
        return (0.0, 0.0)

    random.seed(42)  # Reproducibility
    medians = []
    n = len(values)

    for _ in range(n_iterations):
        sample = [random.choice(values) for _ in range(n)]
        sample.sort()
        medians.append(sample[len(sample) // 2])

    medians.sort()
    alpha = (1 - ci) / 2
    lower_idx = int(alpha * len(medians))
    upper_idx = int((1 - alpha) * len(medians)) - 1

    return (medians[lower_idx], medians[upper_idx])


def estimate_ff_y_chromosome(
    y_reads: int,
    autosome_reads: int,
    total_y_length: int = 57_227_415,
    total_autosome_length: int = 2_875_001_522,
) -> Dict[str, Any]:
    """
    Estimate fetal fraction using Y-chromosome read proportion.
    Only applicable for male fetuses.

    FF ≈ 2 * (Y_reads / autosome_reads) * (autosome_length / Y_length)
    """
    if autosome_reads == 0 or y_reads == 0:
        return {
            "method": "y_chromosome",
            "fetal_fraction": None,
            "status": "NOT_APPLICABLE",
            "message": "No Y-chromosome reads detected (likely female fetus).",
        }

    # Normalize by chromosome length
    y_density = y_reads / total_y_length
    auto_density = autosome_reads / total_autosome_length

    # FF = 2 * normalized Y proportion
    ff = 2.0 * y_density / auto_density

    if ff > 0.01:
        predicted_sex = "male"
        status = "PASS"
    else:
        predicted_sex = "female"
        status = "LIKELY_FEMALE"

    return {
        "method": "y_chromosome",
        "fetal_fraction": round(ff, 4),
        "y_reads": y_reads,
        "autosome_reads": autosome_reads,
        "predicted_fetal_sex": predicted_sex,
        "status": status,
    }


def estimate_ff_fragment_size(
    fragment_sizes_file: str,
    short_threshold: int = 150,
    long_threshold: int = 250,
) -> Dict[str, Any]:
    """
    Estimate fetal fraction using cfDNA fragment size distribution.
    Fetal cfDNA fragments are typically shorter (~143bp) than maternal (~166bp).

    FF ≈ proportion of short fragments (< threshold)
    This is a rough estimate and should be used as supplementary evidence.
    """
    logger.info("Estimating FF from fragment size distribution: %s", fragment_sizes_file)

    sizes = []
    with open(fragment_sizes_file, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                try:
                    size = int(parts[0])
                    count = int(parts[1])
                    sizes.extend([size] * count)
                except ValueError:
                    continue
            else:
                try:
                    sizes.append(int(parts[0]))
                except ValueError:
                    continue

    if not sizes:
        return {
            "method": "fragment_size",
            "fetal_fraction": None,
            "status": "NO_DATA",
        }

    total = len(sizes)
    short_fragments = sum(1 for s in sizes if s < short_threshold)
    long_fragments = sum(1 for s in sizes if s >= long_threshold)

    # Simple ratio-based estimate
    short_ratio = short_fragments / total
    mean_size = sum(sizes) / total
    median_size = sorted(sizes)[total // 2]

    # Empirical relationship: FF ≈ f(short_ratio)
    # This is a simplified model; real implementations use trained models
    ff_estimate = short_ratio * 0.5  # Rough approximation

    return {
        "method": "fragment_size",
        "fetal_fraction": round(ff_estimate, 4),
        "mean_fragment_size": round(mean_size, 1),
        "median_fragment_size": median_size,
        "short_fragment_ratio": round(short_ratio, 4),
        "total_fragments": total,
        "status": "SUPPLEMENTARY",
        "message": "Fragment size-based FF is a rough estimate. Use SNP-based method as primary.",
    }


def main():
    parser = argparse.ArgumentParser(
        description="Estimate fetal fraction from cfDNA using SNP-based methods."
    )
    parser.add_argument("--sample-id", required=True, help="Sample identifier.")
    parser.add_argument(
        "--vcf",
        default=None,
        help="VCF file with allele depths at target SNP positions.",
    )
    parser.add_argument(
        "--pileup",
        default=None,
        help="Pileup file at target SNP positions (alternative to VCF).",
    )
    parser.add_argument(
        "--target-snps",
        default=None,
        help="TSV file with target SNP positions (chrom, pos, ref, alt). Used with --pileup.",
    )
    parser.add_argument(
        "--y-reads",
        type=int,
        default=None,
        help="Number of reads mapped to Y chromosome.",
    )
    parser.add_argument(
        "--autosome-reads",
        type=int,
        default=None,
        help="Number of reads mapped to autosomes.",
    )
    parser.add_argument(
        "--fragment-sizes",
        default=None,
        help="File with fragment size distribution (size\\tcount per line).",
    )
    parser.add_argument(
        "--config",
        default=None,
        help="JSON configuration file for fetal fraction estimation.",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output JSON file path.",
    )
    args = parser.parse_args()

    config = load_config(args.config)

    results = {
        "sample_id": args.sample_id,
        "methods": {},
        "primary_fetal_fraction": None,
        "primary_method": None,
        "status": None,
    }

    # Method 1: SNP-based (primary)
    snps = []
    if args.vcf:
        snps = parse_vcf_for_snps(args.vcf, config)
    elif args.pileup and args.target_snps:
        snps = parse_pileup_for_snps(args.pileup, args.target_snps, config)

    if snps:
        informative = identify_informative_snps(snps, config)
        snp_result = estimate_fetal_fraction_snp(informative, config)
        results["methods"]["snp_informative"] = snp_result

        # Set as primary if successful
        if snp_result["fetal_fraction"] is not None:
            results["primary_fetal_fraction"] = snp_result["fetal_fraction"]
            results["primary_method"] = "snp_informative"
            results["status"] = snp_result["status"]

    # Method 2: Y-chromosome (supplementary)
    if args.y_reads is not None and args.autosome_reads is not None:
        y_result = estimate_ff_y_chromosome(args.y_reads, args.autosome_reads)
        results["methods"]["y_chromosome"] = y_result

        # Use as primary if SNP method failed
        if results["primary_fetal_fraction"] is None and y_result["fetal_fraction"] is not None:
            results["primary_fetal_fraction"] = y_result["fetal_fraction"]
            results["primary_method"] = "y_chromosome"
            results["status"] = y_result["status"]

    # Method 3: Fragment size (supplementary)
    if args.fragment_sizes:
        frag_result = estimate_ff_fragment_size(args.fragment_sizes)
        results["methods"]["fragment_size"] = frag_result

    # Final status
    if results["primary_fetal_fraction"] is None:
        results["status"] = "FAILED"
        logger.error("Failed to estimate fetal fraction for sample %s", args.sample_id)
    else:
        logger.info(
            "Primary fetal fraction for %s: %.4f (method: %s)",
            args.sample_id,
            results["primary_fetal_fraction"],
            results["primary_method"],
        )

    # Write output
    # Remove per-SNP estimates from output to keep file size manageable
    output_results = json.loads(json.dumps(results))
    if "snp_informative" in output_results.get("methods", {}):
        snp_data = output_results["methods"]["snp_informative"]
        if "per_snp_estimates" in snp_data:
            snp_data["per_snp_estimates_count"] = len(snp_data["per_snp_estimates"])
            # Keep only summary stats, not all individual estimates
            del snp_data["per_snp_estimates"]

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as fh:
        json.dump(output_results, fh, indent=2)
    logger.info("Fetal fraction results written to %s", output_path)

    # Also write detailed per-SNP results for debugging
    detail_path = output_path.with_suffix(".detail.json")
    with open(detail_path, "w") as fh:
        json.dump(results, fh, indent=2)
    logger.info("Detailed results written to %s", detail_path)


if __name__ == "__main__":
    main()
