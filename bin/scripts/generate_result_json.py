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
from typing import Any, Dict, List, Optional, Tuple

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
    """Locate all output files for a given sample within the output directory.

    Prefers carrier-screening-style flat layout under outdir:
      qc/, variant/, summary/, fetal_fraction/
    Older sgNIPT nested layout (outdir/sample_id/...) is still accepted.
    """
    patterns = {
        "fastq_qc": [
            output_dir / "qc" / f"{sample_id}.fastq_qc.json",
            output_dir / sample_id / "qc" / f"{sample_id}.fastq_qc.json",
            output_dir / sample_id / f"{sample_id}_fastq_qc.json",
            output_dir / "fastq_qc" / f"{sample_id}_fastq_qc.json",
            output_dir / f"{sample_id}_fastq_qc.json",
            output_dir / f"{sample_id}.fastq_qc.json",
        ],
        "bam_qc": [
            output_dir / "qc" / f"{sample_id}.bam_qc.json",
            output_dir / sample_id / "qc" / f"{sample_id}.bam_qc.json",
            output_dir / sample_id / f"{sample_id}_bam_qc.json",
            output_dir / "bam_qc" / f"{sample_id}_bam_qc.json",
            output_dir / f"{sample_id}_bam_qc.json",
            output_dir / f"{sample_id}.bam_qc.json",
        ],
        "fetal_fraction": [
            output_dir / "fetal_fraction" / f"{sample_id}.fetal_fraction.json",
            output_dir / sample_id / "fetal_fraction" / f"{sample_id}.fetal_fraction.json",
            output_dir / sample_id / f"{sample_id}_fetal_fraction.json",
            output_dir / "fetal_fraction" / f"{sample_id}_fetal_fraction.json",
            output_dir / f"{sample_id}_fetal_fraction.json",
            output_dir / f"{sample_id}.fetal_fraction.json",
        ],
        "variant_report": [
            output_dir / "variant" / f"{sample_id}.variant_report.json",
            output_dir / sample_id / "variants" / f"{sample_id}.variant_report.json",
            output_dir / sample_id / f"{sample_id}_variant_report.json",
            output_dir / "variant_analysis" / f"{sample_id}_variant_report.json",
            output_dir / f"{sample_id}_variant_report.json",
            output_dir / f"{sample_id}.variant_report.json",
        ],
        "final_report": [
            output_dir / "summary" / f"{sample_id}_final_report.json",
            output_dir / sample_id / "report" / f"{sample_id}_final_report.json",
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


def _parse_bam_qc_blob(data: Dict[str, Any]) -> Dict[str, Any]:
    """Map pipeline bam_qc.json (flagstat/stats/target) to portal summary fields."""
    legacy = data.get("metrics") if isinstance(data.get("metrics"), dict) else {}
    flagstat = data.get("flagstat_metrics") if isinstance(data.get("flagstat_metrics"), dict) else {}
    stats = data.get("stats_metrics") if isinstance(data.get("stats_metrics"), dict) else {}
    target = data.get("target_metrics") if isinstance(data.get("target_metrics"), dict) else {}

    def pick(*keys: str):
        for key in keys:
            for src in (legacy, flagstat, stats, target):
                if key in src and src[key] is not None:
                    return src[key]
        return None

    dup = pick("duplicate_rate", "dup_rate")
    return {
        "total_reads": pick("total_reads", "raw_total_sequences"),
        "mapped_reads": pick("mapped_reads", "reads_mapped"),
        "mapping_rate": pick("mapping_rate"),
        "mean_coverage": pick("mean_coverage", "mean_target_coverage", "mean_genome_coverage"),
        "target_coverage": pick("target_mean_coverage", "mean_target_coverage"),
        "on_target_rate": pick("on_target_rate"),
        "dup_rate": dup,
        "mean_mapping_quality": pick("mean_mapping_quality", "average_quality"),
        "overall_pass": (data.get("qc_evaluation") or {}).get("overall_qc_pass", True),
    }


_FETAL_ORIGINS = frozenset({"fetal_specific", "fetal", "FETAL", "FETAL_SPECIFIC"})
_PATHO_CLNSIGS = frozenset({
    "Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic",
    "PATHOGENIC", "LIKELY_PATHOGENIC",
})


def _compact_vep(vep: Optional[Dict[str, Any]]) -> Optional[Dict[str, Any]]:
    """Return a slim VEP dict with only the fields relevant to the daemon."""
    if not vep:
        return None
    clinvar = vep.get("clinvar") or {}
    return {
        "hgvsc":       vep.get("hgvsc", ""),
        "hgvsp":       vep.get("hgvsp", ""),
        "consequence": vep.get("consequence", ""),
        "impact":      vep.get("impact", ""),
        "symbol":      vep.get("symbol", ""),
        "mane":        vep.get("mane", ""),
        "gnomad_af":   vep.get("gnomad_af", ""),
        "existing":    vep.get("existing", ""),
        "sift":        vep.get("sift", ""),
        "polyphen":    vep.get("polyphen", ""),
        "clinvar": {
            "clnsig":     clinvar.get("clnsig", ""),
            "clndn":      clinvar.get("clndn", ""),
            "clnrevstat": clinvar.get("clnrevstat", ""),
            "clnhgvs":    clinvar.get("clnhgvs", ""),
        },
    }


def _is_pathogenic_finding(finding: Dict[str, Any]) -> bool:
    """Return True when a finding should be treated as pathogenic / likely-pathogenic."""
    # pathogenic_variant field (from variant_analysis.py)
    if finding.get("pathogenic_variant"):
        return True
    # classification field (legacy)
    cls = (finding.get("classification") or "").upper()
    if cls in ("PATHOGENIC", "LIKELY_PATHOGENIC"):
        return True
    # VEP ClinVar signal
    vep = finding.get("vep") or {}
    clnsig = (vep.get("clinvar") or {}).get("clnsig", "")
    if clnsig in _PATHO_CLNSIGS:
        return True
    return False


def _extract_vep_from_final_report(
    final_data: Dict[str, Any],
) -> Tuple[bool, Dict[str, Dict[str, Any]]]:
    """
    Extract VEP annotation map {CHROM:POS:REF:ALT → vep_dict} from final_report.
    Returns (vep_annotated: bool, vep_map: dict).
    """
    vep_annotated = bool(
        (final_data.get("report_metadata") or {}).get("vep_annotated")
    )
    vep_map: Dict[str, Dict[str, Any]] = {}
    if not vep_annotated:
        return False, vep_map

    va = final_data.get("variant_analysis") or {}
    for finding in va.get("clinical_findings") or []:
        vep = finding.get("vep")
        if not vep:
            continue
        key = (
            f"{finding.get('chrom')}:{finding.get('pos')}"
            f":{finding.get('ref')}:{finding.get('alt')}"
        )
        vep_map[key] = vep

    logger.info(
        "VEP annotation loaded from final_report: %d entries (vep_annotated=%s)",
        len(vep_map),
        vep_annotated,
    )
    return vep_annotated, vep_map


def _build_finding_key(finding: Dict[str, Any]) -> str:
    return (
        f"{finding.get('chrom')}:{finding.get('pos')}"
        f":{finding.get('ref')}:{finding.get('alt')}"
    )


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
            result["bam_qc"] = _parse_bam_qc_blob(data)
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

    # ── Final report: load once, used for VEP + UPD ─────────────────────────
    vep_annotated = False
    vep_map: Dict[str, Dict[str, Any]] = {}
    final_data: Optional[Dict[str, Any]] = None

    if files["final_report"]:
        result["final_report_path"] = str(files["final_report"])
        final_data = load_json_safe(files["final_report"])
        if final_data:
            vep_annotated, vep_map = _extract_vep_from_final_report(final_data)
            upd = final_data.get("upd_analysis")
            if upd:
                result["upd_analysis"] = upd
                summary = (upd.get("summary_call") or "").strip()
                if summary and summary not in ("normal", "INSUFFICIENT_INFORMATION", ""):
                    result["status_flags"].append(f"UPD_SUSPECTED:{summary}")
                elif summary == "INSUFFICIENT_INFORMATION":
                    result["status_flags"].append("UPD_INSUFFICIENT_INFO")

    # ── Variant Analysis ─────────────────────────────────────────────────────
    if files["variant_report"]:
        data = load_json_safe(files["variant_report"])
        if data:
            summary_block = data.get("summary", {})
            # clinical_findings is the primary list (variant_analysis.py output)
            findings: List[Dict[str, Any]] = data.get("clinical_findings", [])

            # Overlay VEP annotation from final_report when available
            if vep_map:
                for f in findings:
                    key = _build_finding_key(f)
                    if key in vep_map and "vep" not in f:
                        f["vep"] = vep_map[key]

            # Deduplicate findings by CHROM:POS:REF:ALT (same locus may appear
            # in multiple target regions)
            seen_keys: set = set()
            deduped: List[Dict[str, Any]] = []
            for f in findings:
                k = _build_finding_key(f)
                if k not in seen_keys:
                    seen_keys.add(k)
                    deduped.append(f)
            findings = deduped

            # Fetal-specific variants (paternal de novo / fetal-origin)
            fetal_findings = [
                f for f in findings if f.get("origin") in _FETAL_ORIGINS
            ]

            # Pathogenic variants (any origin)
            pathogenic = [f for f in findings if _is_pathogenic_finding(f)]

            # Build compact fetal detail (with VEP)
            fetal_detail = []
            for f in fetal_findings[:20]:  # cap at 20
                entry = {
                    "chrom":             f.get("chrom"),
                    "pos":               f.get("pos"),
                    "ref":               f.get("ref"),
                    "alt":               f.get("alt"),
                    "vaf":               f.get("vaf"),
                    "depth":             f.get("depth"),
                    "gene":              f.get("gene", ""),
                    "disease":           f.get("disease", ""),
                    "pathogenic_variant":f.get("pathogenic_variant", ""),
                    "fetal_genotype":    f.get("fetal_genotype", ""),
                    "maternal_genotype": f.get("maternal_genotype", ""),
                    "confidence":        f.get("confidence", ""),
                    "details":           f.get("details", ""),
                }
                vep_compact = _compact_vep(f.get("vep"))
                if vep_compact:
                    entry["vep"] = vep_compact
                fetal_detail.append(entry)

            # Build compact pathogenic detail (with VEP)
            pathogenic_detail = []
            for f in pathogenic[:10]:  # cap at 10
                entry = {
                    "chrom":             f.get("chrom"),
                    "pos":               f.get("pos"),
                    "ref":               f.get("ref"),
                    "alt":               f.get("alt"),
                    "vaf":               f.get("vaf"),
                    "depth":             f.get("depth"),
                    "origin":            f.get("origin", ""),
                    "gene":              f.get("gene", ""),
                    "disease":           f.get("disease", ""),
                    "pathogenic_variant":f.get("pathogenic_variant", ""),
                    "fetal_genotype":    f.get("fetal_genotype", ""),
                    "maternal_genotype": f.get("maternal_genotype", ""),
                    "confidence":        f.get("confidence", ""),
                }
                vep_compact = _compact_vep(f.get("vep"))
                if vep_compact:
                    entry["vep"] = vep_compact
                pathogenic_detail.append(entry)

            result["variant_analysis"] = {
                "vep_annotated":      vep_annotated,
                "total_target_loci":  summary_block.get("total_target_loci"),
                "covered_loci":       summary_block.get("covered_loci"),
                "total_variants":     summary_block.get("total_variants_analyzed"),
                "fetal_variants":     summary_block.get("fetal_specific_variants"),
                "pathogenic_variants":len(pathogenic),
                "pathogenic_details": pathogenic_detail,
                "fetal_specific_detail": fetal_detail,
            }

            if pathogenic:
                result["status_flags"].append("PATHOGENIC_VARIANT_DETECTED")
    else:
        result["variant_analysis"] = None

    # Determine overall status
    if "FF_ESTIMATION_FAILED" in result["status_flags"]:
        result["status"] = "NO_CALL"
    elif "FASTQ_QC_FAIL" in result["status_flags"] or "BAM_QC_FAIL" in result["status_flags"]:
        result["status"] = "WARNING"
    elif "PATHOGENIC_VARIANT_DETECTED" in result["status_flags"]:
        result["status"] = "POSITIVE"
    elif any(f.startswith("UPD_SUSPECTED:") for f in result["status_flags"]):
        result["status"] = "WARNING"
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
