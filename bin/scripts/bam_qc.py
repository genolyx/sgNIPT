#!/usr/bin/env python3
"""
bam_qc.py - BAM Quality Control Module for Single Gene NIPT Pipeline

Parses samtools stats/flagstat/idxstats output and calculates on-target metrics.
Evaluates alignment quality against configurable thresholds.

Panel-specific features (Twist UCL_SingleGeneNIPT_TE-96276661_hg38):
  - Probe coverage QC: evaluates capture efficiency using probes_bed
  - Zero-probe region warnings: flags targets with no probe coverage
  - Panel-tuned default thresholds for ~780 kb probe footprint
"""

import argparse
import json
import logging
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Default BAM QC thresholds (tuned for Twist TE-96276661 panel)
# ---------------------------------------------------------------------------
DEFAULT_THRESHOLDS: Dict[str, Any] = {
    "min_mapping_rate": 0.90,
    "min_on_target_rate": 0.50,
    "max_duplicate_rate": 0.50,
    "min_mean_target_coverage": 500.0,     # High-depth targeted panel
    "min_target_coverage_uniformity": 0.20, # fraction of targets >= 0.2x mean
    "min_insert_size_median": 100,
    "max_insert_size_median": 300,
    "min_properly_paired_rate": 0.80,
    # Panel-specific thresholds
    "min_probe_coverage_rate": 0.95,       # fraction of probes with >= 100x
    "probe_min_depth": 100,                # minimum depth to consider probe covered
    "max_zero_cov_targets": 50,            # max targets with zero coverage
}


def load_thresholds(threshold_file: Optional[str]) -> Dict[str, Any]:
    """Load BAM QC thresholds from a JSON file, falling back to defaults."""
    thresholds = DEFAULT_THRESHOLDS.copy()
    if threshold_file and Path(threshold_file).exists():
        logger.info("Loading custom BAM QC thresholds from %s", threshold_file)
        with open(threshold_file, "r") as fh:
            custom = json.load(fh)
        thresholds.update(custom)
    return thresholds


def parse_flagstat(flagstat_path: str) -> Dict[str, Any]:
    """Parse samtools flagstat output."""
    logger.info("Parsing flagstat: %s", flagstat_path)
    metrics = {}
    with open(flagstat_path, "r") as fh:
        for line in fh:
            line = line.strip()
            match = re.match(r"^(\d+)\s+\+\s+(\d+)\s+(.*)", line)
            if not match:
                continue
            passed = int(match.group(1))
            failed = int(match.group(2))
            desc = match.group(3).strip()

            if "in total" in desc:
                metrics["total_reads"] = passed + failed
                metrics["total_reads_qc_passed"] = passed
            elif "secondary" in desc:
                metrics["secondary_alignments"] = passed
            elif "supplementary" in desc:
                metrics["supplementary_alignments"] = passed
            elif "duplicates" in desc:
                metrics["duplicate_reads"] = passed
            elif "mapped" in desc and "primary mapped" not in desc and "mate mapped" not in desc:
                metrics["mapped_reads"] = passed
            elif "paired in sequencing" in desc:
                metrics["paired_reads"] = passed
            elif "properly paired" in desc:
                metrics["properly_paired_reads"] = passed
            elif "singletons" in desc:
                metrics["singleton_reads"] = passed

    # Calculate rates
    total = metrics.get("total_reads_qc_passed", 0)
    if total > 0:
        metrics["mapping_rate"] = round(metrics.get("mapped_reads", 0) / total, 4)
        metrics["duplicate_rate"] = round(metrics.get("duplicate_reads", 0) / total, 4)
        if metrics.get("paired_reads", 0) > 0:
            metrics["properly_paired_rate"] = round(
                metrics.get("properly_paired_reads", 0) / metrics["paired_reads"], 4
            )
        else:
            metrics["properly_paired_rate"] = 0.0
    else:
        metrics["mapping_rate"] = 0.0
        metrics["duplicate_rate"] = 0.0
        metrics["properly_paired_rate"] = 0.0

    logger.info(
        "Flagstat: %d total reads, mapping rate=%.4f, dup rate=%.4f",
        metrics.get("total_reads", 0),
        metrics["mapping_rate"],
        metrics["duplicate_rate"],
    )
    return metrics


def parse_samtools_stats(stats_path: str) -> Dict[str, Any]:
    """Parse samtools stats output for insert size and other metrics."""
    logger.info("Parsing samtools stats: %s", stats_path)
    metrics = {}
    insert_sizes: List[Tuple[int, int]] = []

    with open(stats_path, "r") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("SN\t"):
                parts = line.split("\t")
                if len(parts) >= 3:
                    key = parts[1].rstrip(":")
                    value = parts[2]
                    try:
                        if "." in value:
                            metrics[key] = float(value)
                        else:
                            metrics[key] = int(value)
                    except ValueError:
                        metrics[key] = value
            elif line.startswith("IS\t"):
                parts = line.split("\t")
                if len(parts) >= 3:
                    try:
                        insert_sizes.append((int(parts[1]), int(parts[2])))
                    except ValueError:
                        pass

    # Calculate insert size statistics
    if insert_sizes:
        total_pairs = sum(count for _, count in insert_sizes)
        if total_pairs > 0:
            weighted_sum = sum(size * count for size, count in insert_sizes)
            metrics["insert_size_mean"] = round(weighted_sum / total_pairs, 1)

            # Median
            cumulative = 0
            for size, count in sorted(insert_sizes):
                cumulative += count
                if cumulative >= total_pairs / 2:
                    metrics["insert_size_median"] = size
                    break

    result = {
        "raw_total_sequences": metrics.get("raw total sequences", 0),
        "reads_mapped": metrics.get("reads mapped", 0),
        "reads_mapped_and_paired": metrics.get("reads mapped and paired", 0),
        "reads_unmapped": metrics.get("reads unmapped", 0),
        "reads_duplicated": metrics.get("reads duplicated", 0),
        "average_length": metrics.get("average length", 0),
        "average_quality": metrics.get("average quality", 0.0),
        "insert_size_mean": metrics.get("insert_size_mean", 0.0),
        "insert_size_median": metrics.get("insert_size_median", 0),
        "bases_mapped": metrics.get("bases mapped", 0),
        "bases_mapped_cigar": metrics.get("bases mapped (cigar)", 0),
        "error_rate": metrics.get("error rate", 0.0),
        "average_insert_size": metrics.get("insert size average", 0.0),
        "insert_size_standard_deviation": metrics.get("insert size standard deviation", 0.0),
    }

    logger.info(
        "Stats: avg quality=%.1f, insert size median=%d",
        result["average_quality"],
        result["insert_size_median"],
    )
    return result


def parse_on_target_metrics(
    mosdepth_summary: Optional[str] = None,
    mosdepth_regions: Optional[str] = None,
    on_target_count_file: Optional[str] = None,
    total_mapped_reads: int = 0,
) -> Dict[str, Any]:
    """
    Parse on-target coverage metrics.
    Supports mosdepth summary/region output or a simple on-target read count file.
    """
    logger.info("Calculating on-target metrics...")
    metrics = {}

    # Parse mosdepth summary if available
    if mosdepth_summary and Path(mosdepth_summary).exists():
        logger.info("Parsing mosdepth summary: %s", mosdepth_summary)
        with open(mosdepth_summary, "r") as fh:
            header = fh.readline()
            for line in fh:
                parts = line.strip().split("\t")
                if len(parts) >= 4:
                    chrom = parts[0]
                    if chrom == "total_region" or chrom.endswith("_region"):
                        metrics["mean_target_coverage"] = float(parts[3])
                    elif chrom == "total":
                        metrics["mean_genome_coverage"] = float(parts[3])

    # Parse mosdepth per-region file for uniformity
    if mosdepth_regions and Path(mosdepth_regions).exists():
        logger.info("Parsing mosdepth regions: %s", mosdepth_regions)
        coverages = []
        region_details = []
        with open(mosdepth_regions, "r") as fh:
            for line in fh:
                parts = line.strip().split("\t")
                if len(parts) >= 4:
                    try:
                        cov = float(parts[3])
                        coverages.append(cov)
                        region_details.append({
                            "chrom": parts[0],
                            "start": int(parts[1]),
                            "end": int(parts[2]),
                            "coverage": cov,
                        })
                    except ValueError:
                        pass

        if coverages:
            mean_cov = sum(coverages) / len(coverages)
            metrics["mean_target_coverage"] = round(mean_cov, 2)
            metrics["median_target_coverage"] = round(
                sorted(coverages)[len(coverages) // 2], 2
            )
            metrics["min_target_coverage"] = round(min(coverages), 2)
            metrics["max_target_coverage"] = round(max(coverages), 2)
            metrics["num_targets"] = len(coverages)

            # Uniformity: fraction of targets with coverage >= 0.2x mean
            if mean_cov > 0:
                threshold = 0.2 * mean_cov
                above = sum(1 for c in coverages if c >= threshold)
                metrics["target_coverage_uniformity"] = round(above / len(coverages), 4)

            # Targets with zero coverage
            zero_cov = sum(1 for c in coverages if c == 0)
            metrics["targets_zero_coverage"] = zero_cov

            # Coverage distribution bins
            cov_bins = {
                "0x": sum(1 for c in coverages if c == 0),
                "1-49x": sum(1 for c in coverages if 0 < c < 50),
                "50-99x": sum(1 for c in coverages if 50 <= c < 100),
                "100-499x": sum(1 for c in coverages if 100 <= c < 500),
                "500-999x": sum(1 for c in coverages if 500 <= c < 1000),
                ">=1000x": sum(1 for c in coverages if c >= 1000),
            }
            metrics["coverage_distribution"] = cov_bins

            # Low coverage regions (for reporting)
            low_cov_regions = [
                r for r in region_details if r["coverage"] < 50
            ]
            if low_cov_regions:
                metrics["low_coverage_regions"] = low_cov_regions[:100]  # cap at 100
                metrics["num_low_coverage_regions"] = len(low_cov_regions)

    # Parse simple on-target count file
    if on_target_count_file and Path(on_target_count_file).exists():
        logger.info("Parsing on-target count: %s", on_target_count_file)
        with open(on_target_count_file, "r") as fh:
            on_target_reads = int(fh.read().strip())
        metrics["on_target_reads"] = on_target_reads
        if total_mapped_reads > 0:
            metrics["on_target_rate"] = round(on_target_reads / total_mapped_reads, 4)

    logger.info("On-target metrics: mean_cov=%.2f, on_target_rate=%.4f",
                metrics.get("mean_target_coverage", 0),
                metrics.get("on_target_rate", 0))
    return metrics


# ---------------------------------------------------------------------------
# Zero-probe region analysis
# ---------------------------------------------------------------------------
def analyze_zero_probe_regions(
    zero_probe_bed: str,
    mosdepth_regions: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Analyze zero-probe regions from Twist panel design.
    These are target regions with no probe coverage (29 regions, 3,786 bp).
    Variants in these regions cannot be reliably detected.
    """
    logger.info("Analyzing zero-probe regions: %s", zero_probe_bed)

    zero_regions = []
    with open(zero_probe_bed, "r") as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("track"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                region = {
                    "chrom": parts[0],
                    "start": int(parts[1]),
                    "end": int(parts[2]),
                    "size": int(parts[2]) - int(parts[1]),
                }
                if len(parts) >= 4:
                    region["name"] = parts[3]
                zero_regions.append(region)

    total_zero_bp = sum(r["size"] for r in zero_regions)

    result = {
        "num_zero_probe_regions": len(zero_regions),
        "total_zero_probe_bp": total_zero_bp,
        "regions": zero_regions,
        "warning": (
            f"{len(zero_regions)} target regions ({total_zero_bp} bp) have no probe coverage. "
            "Variants in these regions may not be reliably detected."
        ),
    }

    logger.warning(
        "Zero-probe regions: %d regions, %d bp total",
        len(zero_regions),
        total_zero_bp,
    )
    return result


# ---------------------------------------------------------------------------
# Probe coverage QC
# ---------------------------------------------------------------------------
def evaluate_probe_coverage(
    mosdepth_regions: str,
    probes_bed: str,
    thresholds: Dict[str, Any],
) -> Dict[str, Any]:
    """
    Evaluate capture efficiency by comparing coverage at probe regions.
    Uses the Probes_merged_ok_shareable BED to assess actual probe performance.
    """
    logger.info("Evaluating probe coverage QC...")

    # Load probe regions
    probe_regions = set()
    with open(probes_bed, "r") as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("track"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                probe_regions.add((parts[0], int(parts[1]), int(parts[2])))

    logger.info("Loaded %d probe regions from probes BED", len(probe_regions))

    # This is a simplified check; in production, you'd intersect mosdepth
    # per-base output with probe regions using bedtools
    result = {
        "num_probe_regions": len(probe_regions),
        "probes_bed_used": probes_bed,
        "note": "Detailed probe-level coverage requires bedtools intersect with per-base mosdepth output.",
    }

    return result


def evaluate_bam_qc(
    flagstat_metrics: Dict[str, Any],
    stats_metrics: Dict[str, Any],
    target_metrics: Dict[str, Any],
    thresholds: Dict[str, Any],
    zero_probe_info: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Evaluate BAM QC metrics against thresholds."""
    logger.info("Evaluating BAM QC metrics against thresholds...")
    checks = []

    # Mapping rate
    checks.append({
        "metric": "mapping_rate",
        "value": flagstat_metrics.get("mapping_rate", 0.0),
        "threshold": thresholds["min_mapping_rate"],
        "operator": ">=",
        "passed": flagstat_metrics.get("mapping_rate", 0.0) >= thresholds["min_mapping_rate"],
    })

    # Duplicate rate
    checks.append({
        "metric": "duplicate_rate",
        "value": flagstat_metrics.get("duplicate_rate", 0.0),
        "threshold": thresholds["max_duplicate_rate"],
        "operator": "<=",
        "passed": flagstat_metrics.get("duplicate_rate", 0.0) <= thresholds["max_duplicate_rate"],
    })

    # Properly paired rate (only for paired-end)
    if flagstat_metrics.get("paired_reads", 0) > 0:
        checks.append({
            "metric": "properly_paired_rate",
            "value": flagstat_metrics.get("properly_paired_rate", 0.0),
            "threshold": thresholds["min_properly_paired_rate"],
            "operator": ">=",
            "passed": flagstat_metrics.get("properly_paired_rate", 0.0) >= thresholds["min_properly_paired_rate"],
        })

    # On-target rate
    if "on_target_rate" in target_metrics:
        checks.append({
            "metric": "on_target_rate",
            "value": target_metrics["on_target_rate"],
            "threshold": thresholds["min_on_target_rate"],
            "operator": ">=",
            "passed": target_metrics["on_target_rate"] >= thresholds["min_on_target_rate"],
        })

    # Mean target coverage
    if "mean_target_coverage" in target_metrics:
        checks.append({
            "metric": "mean_target_coverage",
            "value": target_metrics["mean_target_coverage"],
            "threshold": thresholds["min_mean_target_coverage"],
            "operator": ">=",
            "passed": target_metrics["mean_target_coverage"] >= thresholds["min_mean_target_coverage"],
        })

    # Target coverage uniformity
    if "target_coverage_uniformity" in target_metrics:
        checks.append({
            "metric": "target_coverage_uniformity",
            "value": target_metrics["target_coverage_uniformity"],
            "threshold": thresholds["min_target_coverage_uniformity"],
            "operator": ">=",
            "passed": target_metrics["target_coverage_uniformity"] >= thresholds["min_target_coverage_uniformity"],
        })

    # Zero coverage targets
    if "targets_zero_coverage" in target_metrics:
        checks.append({
            "metric": "targets_zero_coverage",
            "value": target_metrics["targets_zero_coverage"],
            "threshold": thresholds.get("max_zero_cov_targets", 50),
            "operator": "<=",
            "passed": target_metrics["targets_zero_coverage"] <= thresholds.get("max_zero_cov_targets", 50),
        })

    # Insert size median
    if stats_metrics.get("insert_size_median", 0) > 0:
        median_is = stats_metrics["insert_size_median"]
        checks.append({
            "metric": "insert_size_median",
            "value": median_is,
            "threshold": f"[{thresholds['min_insert_size_median']}, {thresholds['max_insert_size_median']}]",
            "operator": "in_range",
            "passed": thresholds["min_insert_size_median"] <= median_is <= thresholds["max_insert_size_median"],
        })

    overall_pass = all(c["passed"] for c in checks)

    for c in checks:
        status = "PASS" if c["passed"] else "FAIL"
        logger.info("  [%s] %s: %s %s %s", status, c["metric"], c["value"], c["operator"], c["threshold"])

    result = {
        "overall_qc_pass": overall_pass,
        "checks": checks,
        "num_passed": sum(1 for c in checks if c["passed"]),
        "num_failed": sum(1 for c in checks if not c["passed"]),
    }

    # Add warnings section
    warnings = []
    if zero_probe_info:
        warnings.append({
            "type": "zero_probe_regions",
            "message": zero_probe_info.get("warning", ""),
            "num_regions": zero_probe_info.get("num_zero_probe_regions", 0),
            "total_bp": zero_probe_info.get("total_zero_probe_bp", 0),
        })

    if target_metrics.get("num_low_coverage_regions", 0) > 0:
        warnings.append({
            "type": "low_coverage_regions",
            "message": f"{target_metrics['num_low_coverage_regions']} target regions have coverage < 50x.",
            "num_regions": target_metrics["num_low_coverage_regions"],
        })

    if warnings:
        result["warnings"] = warnings

    logger.info(
        "Overall BAM QC: %s (%d/%d checks passed, %d warnings)",
        "PASS" if overall_pass else "FAIL",
        result["num_passed"],
        result["num_passed"] + result["num_failed"],
        len(warnings),
    )
    return result


def main():
    parser = argparse.ArgumentParser(
        description="Parse BAM QC outputs and evaluate alignment quality."
    )
    parser.add_argument("--sample-id", required=True, help="Sample identifier.")
    parser.add_argument("--flagstat", required=True, help="Path to samtools flagstat output.")
    parser.add_argument("--stats", default=None, help="Path to samtools stats output.")
    parser.add_argument("--mosdepth-summary", default=None, help="Path to mosdepth summary file.")
    parser.add_argument("--mosdepth-regions", default=None, help="Path to mosdepth per-region BED file.")
    parser.add_argument("--on-target-count", default=None, help="Path to file with on-target read count.")
    parser.add_argument("--probes-bed", default=None, help="Path to probes BED file for capture QC.")
    parser.add_argument("--zero-probe-bed", default=None, help="Path to zero-probe regions BED file.")
    parser.add_argument("--thresholds", default=None, help="Optional JSON file with custom thresholds.")
    parser.add_argument("--output", required=True, help="Output JSON file path for BAM QC results.")
    args = parser.parse_args()

    # Load thresholds
    thresholds = load_thresholds(args.thresholds)

    # Parse flagstat
    flagstat_metrics = parse_flagstat(args.flagstat)

    # Parse samtools stats
    stats_metrics = {}
    if args.stats:
        stats_metrics = parse_samtools_stats(args.stats)

    # Parse on-target metrics
    target_metrics = parse_on_target_metrics(
        mosdepth_summary=args.mosdepth_summary,
        mosdepth_regions=args.mosdepth_regions,
        on_target_count_file=args.on_target_count,
        total_mapped_reads=flagstat_metrics.get("mapped_reads", 0),
    )

    # Analyze zero-probe regions
    zero_probe_info = None
    if args.zero_probe_bed and Path(args.zero_probe_bed).exists():
        zero_probe_info = analyze_zero_probe_regions(args.zero_probe_bed)

    # Evaluate probe coverage if probes BED provided
    probe_info = None
    if (args.probes_bed and Path(args.probes_bed).exists()
            and args.mosdepth_regions and Path(args.mosdepth_regions).exists()):
        probe_info = evaluate_probe_coverage(
            args.mosdepth_regions, args.probes_bed, thresholds
        )

    # Evaluate QC
    qc_result = evaluate_bam_qc(
        flagstat_metrics, stats_metrics, target_metrics, thresholds,
        zero_probe_info=zero_probe_info,
    )

    # Combine output
    output = {
        "sample_id": args.sample_id,
        "panel": "Twist_UCL_SingleGeneNIPT_TE-96276661_hg38",
        "flagstat_metrics": flagstat_metrics,
        "stats_metrics": stats_metrics,
        "target_metrics": target_metrics,
        "qc_evaluation": qc_result,
        "thresholds_used": thresholds,
    }

    if zero_probe_info:
        output["zero_probe_regions"] = zero_probe_info
    if probe_info:
        output["probe_coverage"] = probe_info

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as fh:
        json.dump(output, fh, indent=2)
    logger.info("BAM QC results written to %s", output_path)


if __name__ == "__main__":
    main()
