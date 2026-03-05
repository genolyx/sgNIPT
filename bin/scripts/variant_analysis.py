#!/usr/bin/env python3
"""
variant_analysis.py - Variant Analysis Module for Single Gene NIPT Pipeline

Analyzes variants at target loci to determine fetal genotype status.
Classifies variants as maternal, fetal-specific, or shared based on
allele frequencies and fetal fraction estimates.

Key concepts:
- Maternal homozygous (AA): expected AF ~ 0% (ref) or ~ 100% (alt)
- Maternal heterozygous (AB): expected AF ~ 50%
- Fetal-specific (paternal) allele: expected AF ~ FF/2
- Fetal heterozygous from maternal carrier: expected AF ~ 50% + FF/2 or 50% - FF/2

Panel-specific features (Twist UCL_SingleGeneNIPT_TE-96276661_hg38):
  - Zero-probe region flagging: variants in 29 zero-probe regions are flagged
  - Gene list validation: cross-checks detected genes against singleGene_NIPT.xlsx
  - 182 target genes across 5,138 regions
"""

import argparse
import csv
import json
import logging
import math
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Default configuration
# ---------------------------------------------------------------------------
DEFAULT_CONFIG: Dict[str, Any] = {
    # Minimum read depth for variant analysis
    "min_depth": 100,
    # Minimum variant allele frequency to call a variant present
    "min_vaf": 0.005,
    # Maternal homozygous threshold: AF > this -> likely maternal hom alt
    "maternal_hom_alt_threshold": 0.85,
    # Maternal homozygous ref threshold: AF < this -> likely maternal hom ref
    "maternal_hom_ref_threshold": 0.15,
    # Maternal heterozygous range
    "maternal_het_lower": 0.35,
    "maternal_het_upper": 0.65,
    # Fetal-specific allele detection: AF should be around FF/2
    # We use a window of FF/2 +/- tolerance
    "fetal_af_tolerance_factor": 2.0,  # multiplier for expected FF/2
    # Minimum base quality
    "min_base_quality": 20,
    # Minimum mapping quality
    "min_mapping_quality": 30,
    # Binomial test p-value threshold for fetal allele detection
    "binomial_pvalue_threshold": 0.01,
}


def load_config(config_file: Optional[str]) -> Dict[str, Any]:
    """Load configuration from JSON file."""
    config = DEFAULT_CONFIG.copy()
    if config_file and Path(config_file).exists():
        logger.info("Loading variant analysis config from %s", config_file)
        with open(config_file, "r") as fh:
            custom = json.load(fh)
        config.update(custom)
    return config


# ---------------------------------------------------------------------------
# Gene list loading and validation
# ---------------------------------------------------------------------------
def load_gene_list(gene_list_path: str) -> Dict[str, Dict[str, Any]]:
    """
    Load the target gene list (singleGene_NIPT.xlsx exported as TSV/CSV,
    or a simple text file with one gene per line).

    Returns dict: { gene_symbol: { "disease": ..., "inheritance": ..., ... } }
    """
    logger.info("Loading gene list: %s", gene_list_path)
    genes = {}

    path = Path(gene_list_path)
    suffix = path.suffix.lower()

    if suffix in (".tsv", ".csv", ".txt"):
        delimiter = "\t" if suffix == ".tsv" else ","
        with open(gene_list_path, "r") as fh:
            # Try to detect if there's a header
            first_line = fh.readline().strip()
            fh.seek(0)

            if any(kw in first_line.lower() for kw in ["gene", "symbol", "name", "disease"]):
                reader = csv.DictReader(fh, delimiter=delimiter)
                for row in reader:
                    gene = (
                        row.get("Gene", "")
                        or row.get("gene", "")
                        or row.get("Gene Symbol", "")
                        or row.get("symbol", "")
                    ).strip()
                    if gene:
                        genes[gene] = {
                            "disease": row.get("Disease", row.get("disease", "")),
                            "inheritance": row.get("Inheritance", row.get("inheritance", "")),
                            "omim": row.get("OMIM", row.get("omim", "")),
                        }
            else:
                # Simple one-gene-per-line format
                for line in fh:
                    gene = line.strip()
                    if gene and not gene.startswith("#"):
                        genes[gene] = {}
    elif suffix in (".json",):
        with open(gene_list_path, "r") as fh:
            data = json.load(fh)
        if isinstance(data, list):
            for item in data:
                if isinstance(item, str):
                    genes[item] = {}
                elif isinstance(item, dict):
                    gene = item.get("gene", item.get("Gene", ""))
                    if gene:
                        genes[gene] = item
        elif isinstance(data, dict):
            genes = data

    logger.info("Loaded %d genes from gene list", len(genes))
    return genes


def validate_gene_coverage(
    targets: List[Dict[str, Any]],
    gene_list: Dict[str, Dict[str, Any]],
) -> Dict[str, Any]:
    """
    Validate that all genes in the gene list are covered by target regions.
    Returns coverage report.
    """
    logger.info("Validating gene coverage against gene list...")

    # Extract genes from target BED
    target_genes = set()
    gene_to_regions = {}
    for t in targets:
        gene = t.get("gene", "")
        name = t.get("name", "")
        # Try to extract gene from name field if gene field is empty
        if not gene and name:
            # Common BED name formats: "GENE_exon1", "GENE:NM_xxx", etc.
            gene = name.split("_")[0].split(":")[0].split("|")[0]
        if gene:
            target_genes.add(gene)
            if gene not in gene_to_regions:
                gene_to_regions[gene] = []
            gene_to_regions[gene].append(t)

    # Check coverage
    covered = set()
    missing = set()
    for gene in gene_list:
        if gene in target_genes:
            covered.add(gene)
        else:
            missing.add(gene)

    # Extra genes in panel not in gene list
    extra = target_genes - set(gene_list.keys())

    result = {
        "total_genes_in_list": len(gene_list),
        "genes_covered_by_panel": len(covered),
        "genes_missing_from_panel": len(missing),
        "extra_genes_in_panel": len(extra),
        "coverage_rate": round(len(covered) / len(gene_list), 4) if gene_list else 0.0,
        "missing_genes": sorted(missing),
        "extra_genes": sorted(extra)[:50],  # cap output
    }

    if missing:
        logger.warning(
            "Gene coverage: %d/%d genes covered. Missing: %s",
            len(covered),
            len(gene_list),
            ", ".join(sorted(missing)[:10]),
        )
    else:
        logger.info("Gene coverage: All %d genes covered by panel", len(gene_list))

    return result


# ---------------------------------------------------------------------------
# Zero-probe region handling
# ---------------------------------------------------------------------------
def load_zero_probe_regions(zero_probe_bed: str) -> List[Dict[str, Any]]:
    """Load zero-probe regions from BED file."""
    regions = []
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
                }
                if len(parts) >= 4:
                    region["name"] = parts[3]
                regions.append(region)
    logger.info("Loaded %d zero-probe regions", len(regions))
    return regions


def check_variant_in_zero_probe(
    variant: Dict[str, Any],
    zero_probe_regions: List[Dict[str, Any]],
) -> bool:
    """Check if a variant falls within a zero-probe region."""
    chrom = variant["chrom"]
    pos = variant["pos"]
    for region in zero_probe_regions:
        if (region["chrom"] == chrom
                and region["start"] <= pos <= region["end"]):
            return True
    return False


# ---------------------------------------------------------------------------
# Target BED and VCF parsing
# ---------------------------------------------------------------------------
def load_target_bed(bed_path: str) -> List[Dict[str, Any]]:
    """
    Load target BED file with disease/gene annotations.
    Expected format: chrom, start, end, name (gene/disease info), [additional columns]
    """
    logger.info("Loading target BED: %s", bed_path)
    targets = []

    with open(bed_path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                continue

            target = {
                "chrom": fields[0],
                "start": int(fields[1]),
                "end": int(fields[2]),
                "name": fields[3] if len(fields) > 3 else f"{fields[0]}:{fields[1]}-{fields[2]}",
            }

            # Optional additional columns
            if len(fields) > 4:
                target["score"] = fields[4]
            if len(fields) > 5:
                target["strand"] = fields[5]
            if len(fields) > 6:
                target["gene"] = fields[6]
            if len(fields) > 7:
                target["disease"] = fields[7]
            if len(fields) > 8:
                target["pathogenic_variant"] = fields[8]

            targets.append(target)

    logger.info("Loaded %d target regions from BED file", len(targets))
    return targets


def parse_vcf_variants(vcf_path: str, config: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Parse VCF file and extract variant information with allele depths."""
    logger.info("Parsing VCF for variant analysis: %s", vcf_path)
    variants = []

    with open(vcf_path, "r") as fh:
        header_fields = []
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header_fields = line.strip().split("\t")
                continue

            fields = line.strip().split("\t")
            if len(fields) < 10:
                continue

            chrom = fields[0]
            pos = int(fields[1])
            var_id = fields[2]
            ref = fields[3]
            alt_alleles = fields[4].split(",")
            qual = fields[5]
            filt = fields[6]
            info = fields[7]
            fmt = fields[8]
            sample = fields[9]

            # Parse FORMAT and sample
            fmt_keys = fmt.split(":")
            fmt_values = sample.split(":")
            fmt_dict = dict(zip(fmt_keys, fmt_values))

            # Extract genotype
            gt = fmt_dict.get("GT", "./.")

            # Extract allele depths
            ad = None
            dp = None
            if "AD" in fmt_dict:
                try:
                    ad = [int(x) for x in fmt_dict["AD"].split(",")]
                except ValueError:
                    ad = None

            if "DP" in fmt_dict:
                try:
                    dp = int(fmt_dict["DP"])
                except ValueError:
                    dp = None

            if dp is None and ad:
                dp = sum(ad)

            if dp is None or dp < config["min_depth"]:
                continue

            # Parse INFO field
            info_dict = {}
            for item in info.split(";"):
                if "=" in item:
                    k, v = item.split("=", 1)
                    info_dict[k] = v
                else:
                    info_dict[item] = True

            for i, alt in enumerate(alt_alleles):
                if alt == ".":
                    continue

                ref_count = ad[0] if ad and len(ad) > 0 else 0
                alt_count = ad[i + 1] if ad and len(ad) > i + 1 else 0
                vaf = alt_count / dp if dp > 0 else 0.0

                if vaf < config["min_vaf"]:
                    continue

                variant = {
                    "chrom": chrom,
                    "pos": pos,
                    "id": var_id,
                    "ref": ref,
                    "alt": alt,
                    "qual": qual,
                    "filter": filt,
                    "genotype": gt,
                    "depth": dp,
                    "ref_count": ref_count,
                    "alt_count": alt_count,
                    "vaf": round(vaf, 6),
                    "info": info_dict,
                }
                variants.append(variant)

    logger.info("Parsed %d variants from VCF", len(variants))
    return variants


# ---------------------------------------------------------------------------
# Variant classification
# ---------------------------------------------------------------------------
def classify_variant_origin(
    variant: Dict[str, Any],
    fetal_fraction: float,
    config: Dict[str, Any],
) -> Dict[str, Any]:
    """
    Classify a variant's origin based on allele frequency and fetal fraction.

    Classification logic:
    1. VAF ~ 50%: Maternal heterozygous (could also be shared with fetus)
    2. VAF ~ 100%: Maternal homozygous alt
    3. VAF ~ FF/2 (low): Fetal-specific (paternal origin, fetus heterozygous)
    4. VAF ~ 50% + FF/2: Maternal carrier + fetal inherited
    5. VAF ~ 50% - FF/2: Maternal carrier, fetus did NOT inherit
    """
    vaf = variant["vaf"]
    depth = variant["depth"]
    alt_count = variant["alt_count"]

    expected_fetal_af = fetal_fraction / 2.0
    tolerance = expected_fetal_af * config["fetal_af_tolerance_factor"]

    classification = {
        "origin": "unknown",
        "confidence": "low",
        "fetal_genotype_prediction": "unknown",
        "maternal_genotype_prediction": "unknown",
        "details": "",
    }

    # Case 1: Very high VAF -> maternal homozygous alt
    if vaf >= config["maternal_hom_alt_threshold"]:
        classification["origin"] = "maternal"
        classification["maternal_genotype_prediction"] = "hom_alt"
        classification["confidence"] = "high"
        classification["fetal_genotype_prediction"] = "at_least_het"
        classification["details"] = (
            f"VAF={vaf:.4f} indicates maternal homozygous alt. "
            "Fetus inherits at least one alt allele from mother."
        )

    # Case 2: VAF in heterozygous range -> maternal heterozygous
    elif config["maternal_het_lower"] <= vaf <= config["maternal_het_upper"]:
        classification["origin"] = "maternal_het"
        classification["maternal_genotype_prediction"] = "het"

        deviation = vaf - 0.5
        if abs(deviation) < tolerance:
            classification["confidence"] = "medium"
            classification["fetal_genotype_prediction"] = "uncertain"
            classification["details"] = (
                f"VAF={vaf:.4f} consistent with maternal heterozygous. "
                "Cannot determine if fetus inherited the variant allele."
            )
        elif deviation > 0 and abs(deviation - expected_fetal_af) < tolerance:
            classification["confidence"] = "medium"
            classification["fetal_genotype_prediction"] = "likely_het"
            classification["details"] = (
                f"VAF={vaf:.4f} slightly above 50%, consistent with maternal het + "
                f"fetal inheritance (expected shift ~ {expected_fetal_af:.4f})."
            )
        elif deviation < 0 and abs(deviation + expected_fetal_af) < tolerance:
            classification["confidence"] = "medium"
            classification["fetal_genotype_prediction"] = "likely_ref"
            classification["details"] = (
                f"VAF={vaf:.4f} slightly below 50%, consistent with maternal het + "
                "fetal NOT inheriting the variant allele."
            )
        else:
            classification["confidence"] = "low"
            classification["fetal_genotype_prediction"] = "uncertain"

    # Case 3: Low VAF around FF/2 -> fetal-specific (paternal allele)
    elif abs(vaf - expected_fetal_af) < tolerance and expected_fetal_af > 0:
        p_value = binomial_test(alt_count, depth, 0.001)

        if p_value < config["binomial_pvalue_threshold"]:
            classification["origin"] = "fetal_specific"
            classification["maternal_genotype_prediction"] = "hom_ref"
            classification["fetal_genotype_prediction"] = "het"
            classification["confidence"] = "high"
            classification["details"] = (
                f"VAF={vaf:.4f} ~ FF/2={expected_fetal_af:.4f}. "
                f"Binomial p-value={p_value:.2e}. "
                "Likely fetal-specific variant (paternal origin)."
            )
        else:
            classification["origin"] = "possible_fetal"
            classification["confidence"] = "low"
            classification["fetal_genotype_prediction"] = "uncertain"
            classification["details"] = (
                f"VAF={vaf:.4f} near expected fetal AF but binomial test "
                f"not significant (p={p_value:.2e})."
            )

    # Case 4: Very low VAF -> possible noise or very low fetal signal
    elif vaf < config["min_vaf"] * 2:
        classification["origin"] = "noise"
        classification["confidence"] = "low"
        classification["details"] = f"VAF={vaf:.4f} too low, likely sequencing noise."

    # Case 5: Intermediate VAF
    else:
        classification["origin"] = "ambiguous"
        classification["confidence"] = "low"
        classification["details"] = (
            f"VAF={vaf:.4f} does not clearly fit maternal or fetal-specific patterns."
        )

    # Add statistical support
    classification["binomial_pvalue"] = binomial_test(alt_count, depth, expected_fetal_af)

    return classification


def binomial_test(successes: int, trials: int, expected_prob: float) -> float:
    """
    Simple binomial test (one-sided, greater).
    Returns approximate p-value using normal approximation for large n.
    """
    if trials == 0 or expected_prob <= 0:
        return 1.0

    observed_prob = successes / trials
    se = math.sqrt(expected_prob * (1 - expected_prob) / trials)

    if se == 0:
        return 1.0

    z = (observed_prob - expected_prob) / se
    p_value = 0.5 * math.erfc(z / math.sqrt(2))
    return p_value


def annotate_with_targets(
    variants: List[Dict[str, Any]],
    targets: List[Dict[str, Any]],
) -> List[Dict[str, Any]]:
    """Annotate variants with target region information."""
    logger.info("Annotating variants with target regions...")

    # Build interval lookup
    target_lookup = {}
    for t in targets:
        chrom = t["chrom"]
        if chrom not in target_lookup:
            target_lookup[chrom] = []
        target_lookup[chrom].append(t)

    annotated = 0
    for var in variants:
        chrom = var["chrom"]
        pos = var["pos"]
        var["in_target"] = False
        var["target_info"] = None

        if chrom in target_lookup:
            for t in target_lookup[chrom]:
                if t["start"] <= pos <= t["end"]:
                    var["in_target"] = True
                    var["target_info"] = {
                        "name": t.get("name", ""),
                        "gene": t.get("gene", ""),
                        "disease": t.get("disease", ""),
                        "pathogenic_variant": t.get("pathogenic_variant", ""),
                    }
                    annotated += 1
                    break

    logger.info("Annotated %d variants within target regions", annotated)
    return variants


def generate_variant_report(
    variants: List[Dict[str, Any]],
    fetal_fraction: float,
    sample_id: str,
    gene_coverage: Optional[Dict[str, Any]] = None,
    zero_probe_stats: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Generate a structured variant analysis report."""
    logger.info("Generating variant analysis report...")

    # Summary statistics
    total_variants = len(variants)
    in_target = [v for v in variants if v.get("in_target", False)]
    fetal_specific = [v for v in variants if v.get("classification", {}).get("origin") == "fetal_specific"]
    maternal_het = [v for v in variants if v.get("classification", {}).get("origin") == "maternal_het"]
    maternal_hom = [v for v in variants if v.get("classification", {}).get("origin") == "maternal"]
    in_zero_probe = [v for v in variants if v.get("in_zero_probe_region", False)]

    # Clinically relevant findings
    clinical_findings = []
    for var in variants:
        if not var.get("in_target"):
            continue

        cls = var.get("classification", {})
        target = var.get("target_info", {})

        finding = {
            "chrom": var["chrom"],
            "pos": var["pos"],
            "ref": var["ref"],
            "alt": var["alt"],
            "vaf": var["vaf"],
            "depth": var["depth"],
            "origin": cls.get("origin", "unknown"),
            "fetal_genotype": cls.get("fetal_genotype_prediction", "unknown"),
            "maternal_genotype": cls.get("maternal_genotype_prediction", "unknown"),
            "confidence": cls.get("confidence", "low"),
            "gene": target.get("gene", ""),
            "disease": target.get("disease", ""),
            "target_name": target.get("name", ""),
            "pathogenic_variant": target.get("pathogenic_variant", ""),
            "details": cls.get("details", ""),
            "in_zero_probe_region": var.get("in_zero_probe_region", False),
        }

        # Add warning for zero-probe region variants
        if var.get("in_zero_probe_region", False):
            finding["warning"] = (
                "This variant falls in a zero-probe region. "
                "Detection reliability is reduced."
            )

        clinical_findings.append(finding)

    # Sort by confidence (high first) then by chromosome and position
    confidence_order = {"high": 0, "medium": 1, "low": 2}
    clinical_findings.sort(
        key=lambda x: (
            confidence_order.get(x["confidence"], 3),
            x["chrom"],
            x["pos"],
        )
    )

    report = {
        "sample_id": sample_id,
        "panel": "Twist_UCL_SingleGeneNIPT_TE-96276661_hg38",
        "fetal_fraction_used": fetal_fraction,
        "summary": {
            "total_variants_analyzed": total_variants,
            "variants_in_target": len(in_target),
            "fetal_specific_variants": len(fetal_specific),
            "maternal_heterozygous_variants": len(maternal_het),
            "maternal_homozygous_variants": len(maternal_hom),
            "variants_in_zero_probe_regions": len(in_zero_probe),
            "clinical_findings_count": len(clinical_findings),
        },
        "clinical_findings": clinical_findings,
        "all_target_variants": [
            {
                "chrom": v["chrom"],
                "pos": v["pos"],
                "ref": v["ref"],
                "alt": v["alt"],
                "vaf": v["vaf"],
                "depth": v["depth"],
                "origin": v.get("classification", {}).get("origin", "unknown"),
                "fetal_genotype": v.get("classification", {}).get("fetal_genotype_prediction", "unknown"),
                "confidence": v.get("classification", {}).get("confidence", "low"),
                "target_name": v.get("target_info", {}).get("name", "") if v.get("target_info") else "",
                "in_zero_probe_region": v.get("in_zero_probe_region", False),
            }
            for v in in_target
        ],
    }

    # Add gene coverage validation if available
    if gene_coverage:
        report["gene_coverage_validation"] = gene_coverage

    # Add zero-probe region stats
    if zero_probe_stats:
        report["zero_probe_region_stats"] = zero_probe_stats

    # Warnings section
    warnings = []
    if in_zero_probe:
        warnings.append({
            "type": "zero_probe_variants",
            "message": f"{len(in_zero_probe)} variant(s) detected in zero-probe regions. "
                       "These may have reduced detection reliability.",
            "count": len(in_zero_probe),
        })
    if gene_coverage and gene_coverage.get("genes_missing_from_panel", 0) > 0:
        warnings.append({
            "type": "missing_genes",
            "message": f"{gene_coverage['genes_missing_from_panel']} gene(s) from the target list "
                       "are not covered by the panel.",
            "missing": gene_coverage.get("missing_genes", []),
        })
    if warnings:
        report["warnings"] = warnings

    logger.info(
        "Report: %d total variants, %d in target, %d fetal-specific, %d clinical findings",
        total_variants,
        len(in_target),
        len(fetal_specific),
        len(clinical_findings),
    )
    return report


def main():
    parser = argparse.ArgumentParser(
        description="Analyze variants at target loci for fetal genotype determination."
    )
    parser.add_argument("--sample-id", required=True, help="Sample identifier.")
    parser.add_argument("--vcf", required=True, help="VCF file with variants at target loci.")
    parser.add_argument("--target-bed", required=True, help="Target BED file with disease annotations.")
    parser.add_argument(
        "--fetal-fraction",
        type=float,
        required=True,
        help="Estimated fetal fraction (0-1).",
    )
    parser.add_argument(
        "--gene-list",
        default=None,
        help="Gene list file (TSV/CSV/JSON) for coverage validation.",
    )
    parser.add_argument(
        "--zero-probe-bed",
        default=None,
        help="BED file with zero-probe regions for variant flagging.",
    )
    parser.add_argument("--config", default=None, help="JSON configuration file.")
    parser.add_argument("--output", required=True, help="Output JSON report file.")
    args = parser.parse_args()

    config = load_config(args.config)

    # Load targets
    targets = load_target_bed(args.target_bed)

    # Parse variants
    variants = parse_vcf_variants(args.vcf, config)

    # Annotate with target regions
    variants = annotate_with_targets(variants, targets)

    # Load and check zero-probe regions
    zero_probe_regions = []
    zero_probe_stats = None
    if args.zero_probe_bed and Path(args.zero_probe_bed).exists():
        zero_probe_regions = load_zero_probe_regions(args.zero_probe_bed)
        zero_probe_stats = {
            "num_regions": len(zero_probe_regions),
            "total_bp": sum(r["end"] - r["start"] for r in zero_probe_regions),
        }
        # Flag variants in zero-probe regions
        zp_count = 0
        for var in variants:
            var["in_zero_probe_region"] = check_variant_in_zero_probe(var, zero_probe_regions)
            if var["in_zero_probe_region"]:
                zp_count += 1
        if zp_count > 0:
            logger.warning(
                "%d variants fall within zero-probe regions", zp_count
            )

    # Gene list validation
    gene_coverage = None
    if args.gene_list and Path(args.gene_list).exists():
        gene_list = load_gene_list(args.gene_list)
        gene_coverage = validate_gene_coverage(targets, gene_list)

    # Classify each variant
    logger.info("Classifying variant origins (FF=%.4f)...", args.fetal_fraction)
    for var in variants:
        var["classification"] = classify_variant_origin(var, args.fetal_fraction, config)

    # Generate report
    report = generate_variant_report(
        variants,
        args.fetal_fraction,
        args.sample_id,
        gene_coverage=gene_coverage,
        zero_probe_stats=zero_probe_stats,
    )

    # Write output
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as fh:
        json.dump(report, fh, indent=2)
    logger.info("Variant analysis report written to %s", output_path)

    # Summary to stdout
    summary = report["summary"]
    print(f"\n{'='*60}")
    print(f"VARIANT ANALYSIS SUMMARY - {args.sample_id}")
    print(f"{'='*60}")
    print(f"Panel: Twist UCL_SingleGeneNIPT TE-96276661 (hg38)")
    print(f"Fetal Fraction: {args.fetal_fraction:.4f}")
    print(f"Total variants analyzed: {summary['total_variants_analyzed']}")
    print(f"Variants in target regions: {summary['variants_in_target']}")
    print(f"Fetal-specific variants: {summary['fetal_specific_variants']}")
    print(f"Variants in zero-probe regions: {summary['variants_in_zero_probe_regions']}")
    print(f"Clinical findings: {summary['clinical_findings_count']}")

    if gene_coverage:
        print(f"\nGene Coverage: {gene_coverage['genes_covered_by_panel']}/{gene_coverage['total_genes_in_list']} "
              f"({gene_coverage['coverage_rate']*100:.1f}%)")
        if gene_coverage["missing_genes"]:
            print(f"Missing genes: {', '.join(gene_coverage['missing_genes'][:10])}")

    print(f"{'='*60}\n")

    if report["clinical_findings"]:
        print("CLINICAL FINDINGS:")
        for i, finding in enumerate(report["clinical_findings"], 1):
            zp_flag = " [ZERO-PROBE]" if finding.get("in_zero_probe_region") else ""
            print(
                f"  {i}. {finding['chrom']}:{finding['pos']} {finding['ref']}>{finding['alt']} "
                f"| VAF={finding['vaf']:.4f} | Origin={finding['origin']} "
                f"| Fetal GT={finding['fetal_genotype']} "
                f"| Gene={finding['gene']} | Disease={finding['disease']} "
                f"| Confidence={finding['confidence']}{zp_flag}"
            )
        print()


if __name__ == "__main__":
    main()
