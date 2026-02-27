#!/usr/bin/env python3
"""
fastq_qc.py - FASTQ Quality Control Module for Single Gene NIPT Pipeline

Parses fastp JSON output and generates a structured QC report.
Evaluates sequencing quality metrics against configurable thresholds.
"""

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Any, Dict, Optional

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Default QC thresholds (can be overridden via --thresholds JSON file)
# ---------------------------------------------------------------------------
DEFAULT_THRESHOLDS: Dict[str, Any] = {
    "min_total_reads": 1_000_000,
    "min_q30_rate": 0.80,
    "min_q20_rate": 0.90,
    "max_duplication_rate": 0.50,
    "min_mean_quality_after_filter": 25.0,
    "max_adapter_trimmed_rate": 0.30,
    "min_reads_after_filter": 500_000,
}


def load_thresholds(threshold_file: Optional[str]) -> Dict[str, Any]:
    """Load QC thresholds from a JSON file, falling back to defaults."""
    thresholds = DEFAULT_THRESHOLDS.copy()
    if threshold_file and Path(threshold_file).exists():
        logger.info("Loading custom thresholds from %s", threshold_file)
        with open(threshold_file, "r") as fh:
            custom = json.load(fh)
        thresholds.update(custom)
    return thresholds


def parse_fastp_json(fastp_json_path: str) -> Dict[str, Any]:
    """Parse fastp JSON report and extract key QC metrics."""
    logger.info("Parsing fastp JSON: %s", fastp_json_path)
    with open(fastp_json_path, "r") as fh:
        data = json.load(fh)

    summary = data.get("summary", {})
    before = summary.get("before_filtering", {})
    after = summary.get("after_filtering", {})
    filtering = data.get("filtering_result", {})
    duplication = data.get("duplication", {})
    adapter = data.get("adapter_cutting", {})

    total_reads_before = before.get("total_reads", 0)
    total_reads_after = after.get("total_reads", 0)
    total_bases_before = before.get("total_bases", 0)
    total_bases_after = after.get("total_bases", 0)

    # Q20/Q30 rates after filtering
    q20_rate = after.get("q20_rate", 0.0)
    q30_rate = after.get("q30_rate", 0.0)

    # Handle percentage strings (e.g., "97.5%") or float values
    if isinstance(q20_rate, str):
        q20_rate = float(q20_rate.replace("%", "")) / 100.0
    if isinstance(q30_rate, str):
        q30_rate = float(q30_rate.replace("%", "")) / 100.0

    # Mean quality scores
    quality_before_r1 = before.get("read1_mean_length", 0)
    quality_after = after.get("quality_curves", {})

    metrics = {
        "sample_id": Path(fastp_json_path).stem.replace(".fastp", "").replace("_fastp", ""),
        "total_reads_before_filter": total_reads_before,
        "total_reads_after_filter": total_reads_after,
        "total_bases_before_filter": total_bases_before,
        "total_bases_after_filter": total_bases_after,
        "reads_passed_filter": filtering.get("passed_filter_reads", total_reads_after),
        "reads_low_quality": filtering.get("low_quality_reads", 0),
        "reads_too_many_N": filtering.get("too_many_N_reads", 0),
        "reads_too_short": filtering.get("too_short_reads", 0),
        "q20_rate_after_filter": round(q20_rate, 4),
        "q30_rate_after_filter": round(q30_rate, 4),
        "gc_content_before": round(before.get("gc_content", 0.0), 4),
        "gc_content_after": round(after.get("gc_content", 0.0), 4),
        "duplication_rate": round(duplication.get("rate", 0.0), 4),
        "read1_mean_length_before": before.get("read1_mean_length", 0),
        "read1_mean_length_after": after.get("read1_mean_length", 0),
        "adapter_trimmed_reads": adapter.get("adapter_trimmed_reads", 0),
        "adapter_trimmed_bases": adapter.get("adapter_trimmed_bases", 0),
    }

    # Calculate adapter trimmed rate
    if total_reads_before > 0:
        metrics["adapter_trimmed_rate"] = round(
            metrics["adapter_trimmed_reads"] / total_reads_before, 4
        )
    else:
        metrics["adapter_trimmed_rate"] = 0.0

    # Read2 metrics if paired-end
    if "read2_mean_length" in before:
        metrics["read2_mean_length_before"] = before.get("read2_mean_length", 0)
        metrics["read2_mean_length_after"] = after.get("read2_mean_length", 0)
        metrics["sequencing_mode"] = "paired-end"
    else:
        metrics["sequencing_mode"] = "single-end"

    # Insert size (if available, paired-end only)
    insert_size = data.get("insert_size", {})
    if insert_size:
        metrics["insert_size_peak"] = insert_size.get("peak", 0)
        metrics["insert_size_mean"] = insert_size.get("mean", 0.0)
        metrics["insert_size_std"] = insert_size.get("std", 0.0)

    logger.info(
        "Parsed metrics: %d reads before → %d reads after filtering",
        total_reads_before,
        total_reads_after,
    )
    return metrics


def evaluate_qc(metrics: Dict[str, Any], thresholds: Dict[str, Any]) -> Dict[str, Any]:
    """Evaluate QC metrics against thresholds and return pass/fail status."""
    logger.info("Evaluating QC metrics against thresholds...")
    checks = []

    # Check total reads
    check_total = {
        "metric": "total_reads_before_filter",
        "value": metrics["total_reads_before_filter"],
        "threshold": thresholds["min_total_reads"],
        "operator": ">=",
        "passed": metrics["total_reads_before_filter"] >= thresholds["min_total_reads"],
    }
    checks.append(check_total)

    # Check reads after filter
    check_after = {
        "metric": "total_reads_after_filter",
        "value": metrics["total_reads_after_filter"],
        "threshold": thresholds["min_reads_after_filter"],
        "operator": ">=",
        "passed": metrics["total_reads_after_filter"] >= thresholds["min_reads_after_filter"],
    }
    checks.append(check_after)

    # Check Q30 rate
    check_q30 = {
        "metric": "q30_rate_after_filter",
        "value": metrics["q30_rate_after_filter"],
        "threshold": thresholds["min_q30_rate"],
        "operator": ">=",
        "passed": metrics["q30_rate_after_filter"] >= thresholds["min_q30_rate"],
    }
    checks.append(check_q30)

    # Check Q20 rate
    check_q20 = {
        "metric": "q20_rate_after_filter",
        "value": metrics["q20_rate_after_filter"],
        "threshold": thresholds["min_q20_rate"],
        "operator": ">=",
        "passed": metrics["q20_rate_after_filter"] >= thresholds["min_q20_rate"],
    }
    checks.append(check_q20)

    # Check duplication rate
    check_dup = {
        "metric": "duplication_rate",
        "value": metrics["duplication_rate"],
        "threshold": thresholds["max_duplication_rate"],
        "operator": "<=",
        "passed": metrics["duplication_rate"] <= thresholds["max_duplication_rate"],
    }
    checks.append(check_dup)

    # Check adapter trimmed rate
    check_adapter = {
        "metric": "adapter_trimmed_rate",
        "value": metrics["adapter_trimmed_rate"],
        "threshold": thresholds["max_adapter_trimmed_rate"],
        "operator": "<=",
        "passed": metrics["adapter_trimmed_rate"] <= thresholds["max_adapter_trimmed_rate"],
    }
    checks.append(check_adapter)

    overall_pass = all(c["passed"] for c in checks)

    for c in checks:
        status = "PASS" if c["passed"] else "FAIL"
        logger.info(
            "  [%s] %s: %.4f %s %.4f",
            status,
            c["metric"],
            float(c["value"]),
            c["operator"],
            float(c["threshold"]),
        )

    result = {
        "overall_qc_pass": overall_pass,
        "checks": checks,
        "num_passed": sum(1 for c in checks if c["passed"]),
        "num_failed": sum(1 for c in checks if not c["passed"]),
    }

    logger.info(
        "Overall QC: %s (%d/%d checks passed)",
        "PASS" if overall_pass else "FAIL",
        result["num_passed"],
        result["num_passed"] + result["num_failed"],
    )
    return result


def main():
    parser = argparse.ArgumentParser(
        description="Parse fastp JSON output and evaluate FASTQ QC metrics."
    )
    parser.add_argument(
        "--fastp-json",
        required=True,
        help="Path to fastp JSON report file.",
    )
    parser.add_argument(
        "--thresholds",
        default=None,
        help="Optional JSON file with custom QC thresholds.",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output JSON file path for QC results.",
    )
    parser.add_argument(
        "--sample-id",
        default=None,
        help="Override sample ID (default: inferred from filename).",
    )
    args = parser.parse_args()

    # Load thresholds
    thresholds = load_thresholds(args.thresholds)

    # Parse fastp JSON
    metrics = parse_fastp_json(args.fastp_json)
    if args.sample_id:
        metrics["sample_id"] = args.sample_id

    # Evaluate QC
    qc_result = evaluate_qc(metrics, thresholds)

    # Combine output
    output = {
        "sample_id": metrics["sample_id"],
        "metrics": metrics,
        "qc_evaluation": qc_result,
        "thresholds_used": thresholds,
    }

    # Write output
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as fh:
        json.dump(output, fh, indent=2)
    logger.info("QC results written to %s", output_path)

    # Exit with non-zero if QC failed
    if not qc_result["overall_qc_pass"]:
        logger.warning("FASTQ QC FAILED for sample %s", metrics["sample_id"])
        sys.exit(0)  # Don't fail the pipeline, just flag it


if __name__ == "__main__":
    main()
