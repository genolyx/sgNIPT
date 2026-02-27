#!/usr/bin/env python3
"""
generate_result_json.py - Daemon-Compatible Result JSON Generator

Aggregates all per-sample pipeline outputs into a single Order-level result JSON
that is compatible with the nipt-daemon monitoring pattern.

The daemon monitors for {ORDER_ID}.json in the output directory.
This JSON contains:
  - Order metadata
  - Per-sample results (QC, fetal fraction, variants)
  - Overall order status
"""

import argparse
import json
import logging
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)


def load_json_safe(path: Path) -> Optional[Dict[str, Any]]:
    """Safely load a JSON file, returning None on failure."""
    if not path.exists():
        logger.warning("File not found: %s", path)
        return None
    try:
        with open(path, "r") as fh:
            return json.load(fh)
    except (json.JSONDecodeError, IOError) as e:
        logger.error("Failed to load %s: %s", path, e)
        return None


def find_sample_outputs(output_dir: Path, sample_id: str) -> Dict[str, Optional[Path]]:
    """Locate all output files for a given sample within the output directory."""
    # Search patterns for each output type
    patterns = {
        "fastq_qc": [
            output_dir / sample_id / f"{sample_id}_fastq_qc.json",
            output_dir / "fastq_qc" / f"{sample_id}_fastq_qc.json",
            output_dir / f"{sample_id}_fastq_qc.json",
        ],
        "bam_qc": [
            output_dir / sample_id / f"{sample_id}_bam_qc.json",
            output_dir / "bam_qc" / f"{sample_id}_bam_qc.json",
            output_dir / f"{sample_id}_bam_qc.json",
        ],
        "fetal_fraction": [
            output_dir / sample_id / f"{sample_id}_fetal_fraction.json",
            output_dir / "fetal_fraction" / f"{sample_id}_fetal_fraction.json",
            output_dir / f"{sample_id}_fetal_fraction.json",
        ],
        "variant_report": [
            output_dir / sample_id / f"{sample_id}_variant_report.json",
            output_dir / "variant_analysis" / f"{sample_id}_variant_report.json",
            output_dir / f"{sample_id}_variant_report.json",
        ],
        "final_report": [
            output_dir / sample_id / f"{sample_id}_final_report.json",
            output_dir / "reports" / f"{sample_id}_final_report.json",
            output_dir / f"{sample_id}_final_report.json",
        ],
    }

    found = {}
    for key, candidates in patterns.items():
        found[key] = None
        for candidate in candidates:
            if candidate.exists():
                found[key] = candidate
                break
        # Fallback: glob search
        if found[key] is None:
            matches = list(output_dir.rglob(f"{sample_id}*{key}*.json"))
            if not matches:
                matches = list(output_dir.rglob(f"*{sample_id}*{key.replace('_', '*')}*.json"))
            if matches:
                found[key] = matches[0]

    return found


def build_sample_result(sample_id: str, output_dir: Path) -> Dict[str, Any]:
    """Build a result dictionary for a single sample."""
    logger.info("Processing sample: %s", sample_id)

    files = find_sample_outputs(output_dir, sample_id)

    result = {
        "sample_id": sample_id,
        "status": "PASS",
        "status_flags": [],
    }

    # FASTQ QC
    if files["fastq_qc"]:
        data = load_json_safe(files["fastq_qc"])
        if data:
            result["fastq_qc"] = {
                "total_reads": data.get("metrics", {}).get("total_reads_before_filter"),
                "reads_after_filter": data.get("metrics", {}).get("total_reads_after_filter"),
                "q30_rate": data.get("metrics", {}).get("q30_rate_after_filter"),
                "gc_content": data.get("metrics", {}).get("gc_content_after_filter"),
                "overall_pass": data.get("qc_evaluation", {}).get("overall_qc_pass", True),
            }
            if not result["fastq_qc"]["overall_pass"]:
                result["status_flags"].append("FASTQ_QC_FAIL")
    else:
        result["fastq_qc"] = None

    # BAM QC
    if files["bam_qc"]:
        data = load_json_safe(files["bam_qc"])
        if data:
            result["bam_qc"] = {
                "mapped_reads": data.get("metrics", {}).get("mapped_reads"),
                "mapping_rate": data.get("metrics", {}).get("mapping_rate"),
                "mean_coverage": data.get("metrics", {}).get("mean_coverage"),
                "target_coverage": data.get("metrics", {}).get("target_mean_coverage"),
                "on_target_rate": data.get("metrics", {}).get("on_target_rate"),
                "dup_rate": data.get("metrics", {}).get("duplicate_rate"),
                "overall_pass": data.get("qc_evaluation", {}).get("overall_qc_pass", True),
            }
            if not result["bam_qc"]["overall_pass"]:
                result["status_flags"].append("BAM_QC_FAIL")
    else:
        result["bam_qc"] = None

    # Fetal Fraction
    if files["fetal_fraction"]:
        data = load_json_safe(files["fetal_fraction"])
        if data:
            result["fetal_fraction"] = {
                "primary_ff": data.get("primary_fetal_fraction"),
                "primary_method": data.get("primary_method"),
                "status": data.get("status"),
                "num_informative_snps": data.get("methods", {}).get(
                    "snp_informative", {}
                ).get("num_informative_snps"),
                "confidence_interval": data.get("methods", {}).get(
                    "snp_informative", {}
                ).get("confidence_interval_95"),
            }
            ff_status = data.get("status", "")
            if ff_status == "FAILED":
                result["status_flags"].append("FF_ESTIMATION_FAILED")
            elif ff_status == "LOW_FF":
                result["status_flags"].append("LOW_FETAL_FRACTION")
    else:
        result["fetal_fraction"] = None

    # Variant Analysis
    if files["variant_report"]:
        data = load_json_safe(files["variant_report"])
        if data:
            variants = data.get("variants", [])
            pathogenic = [
                v for v in variants
                if v.get("classification") in ("PATHOGENIC", "LIKELY_PATHOGENIC")
            ]
            result["variant_analysis"] = {
                "total_target_loci": data.get("summary", {}).get("total_target_loci"),
                "covered_loci": data.get("summary", {}).get("covered_loci"),
                "total_variants": data.get("summary", {}).get("total_variants_detected"),
                "fetal_variants": data.get("summary", {}).get("fetal_specific_variants"),
                "pathogenic_variants": len(pathogenic),
                "pathogenic_details": pathogenic[:10],  # Top 10 for summary
            }
            if pathogenic:
                result["status_flags"].append("PATHOGENIC_VARIANT_DETECTED")
    else:
        result["variant_analysis"] = None

    # Final report path
    if files["final_report"]:
        result["final_report_path"] = str(files["final_report"])

    # Determine overall status
    if "FF_ESTIMATION_FAILED" in result["status_flags"]:
        result["status"] = "NO_CALL"
    elif "FASTQ_QC_FAIL" in result["status_flags"] or "BAM_QC_FAIL" in result["status_flags"]:
        result["status"] = "WARNING"
    elif "PATHOGENIC_VARIANT_DETECTED" in result["status_flags"]:
        result["status"] = "POSITIVE"
    elif "LOW_FETAL_FRACTION" in result["status_flags"]:
        result["status"] = "WARNING"
    else:
        result["status"] = "NEGATIVE"

    return result


def build_order_result(
    order_id: str,
    output_dir: Path,
    sample_names: List[str],
    pipeline_version: str,
) -> Dict[str, Any]:
    """Build the complete Order-level result JSON."""
    logger.info("Building order result for: %s", order_id)

    sample_results = []
    for sample_id in sample_names:
        sample_result = build_sample_result(sample_id, output_dir)
        sample_results.append(sample_result)

    # Determine order-level status
    statuses = [s["status"] for s in sample_results]
    if "NO_CALL" in statuses:
        order_status = "NO_CALL"
    elif "POSITIVE" in statuses:
        order_status = "POSITIVE"
    elif "WARNING" in statuses:
        order_status = "WARNING"
    else:
        order_status = "NEGATIVE"

    order_result = {
        "order_id": order_id,
        "pipeline": "single_gene_nipt",
        "pipeline_version": pipeline_version,
        "analysis_type": "single_gene_disorder_screening",
        "generated_at": datetime.now().isoformat(),
        "status": order_status,
        "sample_count": len(sample_results),
        "samples": sample_results,
    }

    return order_result


def main():
    parser = argparse.ArgumentParser(
        description="Generate daemon-compatible result JSON for Single Gene NIPT"
    )
    parser.add_argument("--order_id", required=True, help="Order ID")
    parser.add_argument("--output_dir", required=True, help="Pipeline output directory")
    parser.add_argument("--sample_names", required=True, help="Comma-separated sample names")
    parser.add_argument("--pipeline_version", default="1.0.0", help="Pipeline version")
    parser.add_argument("--output", required=True, help="Output JSON file path")

    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    sample_names = [s.strip() for s in args.sample_names.split(",") if s.strip()]

    if not output_dir.exists():
        logger.error("Output directory does not exist: %s", output_dir)
        sys.exit(1)

    if not sample_names:
        logger.error("No sample names provided")
        sys.exit(1)

    # Build result
    result = build_order_result(
        order_id=args.order_id,
        output_dir=output_dir,
        sample_names=sample_names,
        pipeline_version=args.pipeline_version,
    )

    # Write JSON
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w", encoding="utf-8") as fh:
        json.dump(result, fh, indent=2, ensure_ascii=False, default=str)

    logger.info("Result JSON written to: %s", output_path)
    logger.info("Order status: %s", result["status"])
    logger.info("Samples processed: %d", len(result["samples"]))


if __name__ == "__main__":
    main()
