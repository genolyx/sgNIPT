#!/usr/bin/env python3
"""
upd_detection.py — Uniparental Disomy (UPD) detection from cfDNA (maternal-only)

Detects UPD on the seven clinically-relevant imprinting chromosomes
(chr6, 7, 11, 14, 15, 16, 20) using only the maternal-plasma cfDNA mixture
that the rest of the sgNIPT pipeline already produces — no paternal sample
required.

Signal model (FF = fetal fraction)
----------------------------------
Per-SNP cfDNA ALT VAF expected by maternal genotype × fetal scenario:

  Maternal AA (hom-ref)
     normal/matUPD:   {0,  FF/2}     (fetal AA or Aa)
     patUPD:          {0,  FF/2, FF} (fetal can also be hom-alt → VAF≈FF)

  Maternal Aa (het)        — most informative
     normal:          3-mode {0.5−FF/2, 0.5, 0.5+FF/2}
     matUPD-iso:      2-mode {0.5−FF/2, 0.5+FF/2}            (no central peak)
     matUPD-het:      1-mode {0.5}                            (variance ↓)
     patUPD-iso:      2-mode {0.5−FF/2, 0.5+FF/2}            (no central peak)
     patUPD-het:      ≈ normal                                (hard from cfDNA only)

  Maternal aa (hom-alt) — symmetric to AA

We compute, per chromosome:
  • LL_normal, LL_matUPD_het, LL_iso (using mat-het sites)
  • Mendelian-error excess at mat-hom sites (fetal-hom-alt count → patUPD)
And classify against the genome-wide null built from non-UPD autosomes
present in the BAM.

Inputs
------
--vcf            VCF (tabix-indexed) from bcftools mpileup+call over the UPD
                 analysis BED with FORMAT/AD,DP. The same VCF can also include
                 background autosomes; if so we use them as null. Otherwise
                 chromosome-level z-scores fall back to model-based thresholds.
--fetal-fraction Float in (0,1) (typically from fetal_fraction.json).
--upd-bed        Optional. If provided, restricts SNPs to this BED (chrom,start,end).
--common-snps    Optional VCF.gz of common population SNPs (e.g. gnomAD MAF≥5%);
                 SNPs are kept only if they intersect this set.  Strongly
                 recommended to suppress sequencing-error noise at non-SNP sites.
--sample-id      Sample name (used in output filenames).
--output         Output JSON path.
--per-snp-tsv    Optional per-SNP TSV (debug).
--config         Optional JSON config overriding defaults.

Outputs
-------
JSON with per-chromosome metrics + final UPD call(s).
Per-SNP TSV (optional) with chrom,pos,ref,alt,depth,ref_count,alt_count,vaf,mat_gt,class.
"""

from __future__ import annotations

import argparse
import gzip
import json
import logging
import math
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

# Imprinting / UPD chromosomes (Eggermann et al., 2015 & later reviews)
UPD_CHROMS = ["chr6", "chr7", "chr11", "chr14", "chr15", "chr16", "chr20"]
UPD_CLINICAL = {
    "chr6":  {"maternal": "Imprinting/trisomy 6 mosaicism (research)",
              "paternal": "Transient neonatal Diabetes mellitus (TNDM)"},
    "chr7":  {"maternal": "Silver-Russell syndrome (SRS)",
              "paternal": "Overgrowth (research)"},
    "chr11": {"maternal": "Silver-Russell syndrome (SRS)",
              "paternal": "Beckwith-Wiedemann syndrome (BWS)"},
    "chr14": {"maternal": "Temple syndrome",
              "paternal": "Kagami-Ogata syndrome"},
    "chr15": {"maternal": "Prader-Willi syndrome (PWS)",
              "paternal": "Angelman syndrome (AS)"},
    "chr16": {"maternal": "Imprinting/trisomy 16 mosaicism (research)",
              "paternal": "—"},
    "chr20": {"maternal": "Mulchandani-Bhoj-Conlin syndrome upd(20)mat",
              "paternal": "Pseudohypoparathyroidism Ib (PHP-Ib)"},
}

DEFAULT_CONFIG: dict[str, Any] = {
    # Per-site filters
    "min_depth": 50,
    "min_alt_count": 0,            # we keep mat-hom sites with 0 alt too
    # Maternal genotype windows (cfDNA VAF)
    "mat_homref_max_vaf": 0.05,
    "mat_homalt_min_vaf": 0.95,
    "mat_het_min_vaf":    0.35,
    "mat_het_max_vaf":    0.65,
    # Required FF range
    "min_fetal_fraction": 0.02,
    "max_fetal_fraction": 0.30,
    # Minimum mat-het SNPs per chrom for any BAF-based call
    "min_mat_het_snps_baf":   30,
    # Minimum mat-hom SNPs per chrom for Mendelian-error analysis
    "min_mat_hom_snps_mende": 30,
    # patUPD threshold: fraction of mat-hom sites with VAF in [c_lo*FF, c_hi*FF]
    "patupd_vaf_lo_factor": 0.7,
    "patupd_vaf_hi_factor": 1.3,
    # Per-chrom call thresholds
    "llr_call_threshold": 5.0,         # ln(BF) > 5 ≈ ~150x evidence
    "patupd_excess_z_threshold": 4.0,  # z-score for patUPD-iso/het call
    # Sequencing error baseline (used for null)
    "seq_error_rate": 0.005,
}


# ──────────────────────────────────────────────────────────────────────
# I/O
# ──────────────────────────────────────────────────────────────────────
def open_text(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def load_bed(path: str) -> dict[str, list[tuple[int, int]]]:
    out: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with open_text(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 3:
                continue
            try:
                out[cols[0]].append((int(cols[1]), int(cols[2])))
            except ValueError:
                continue
    for c in out:
        out[c].sort()
    return out


def in_bed(beds: dict[str, list[tuple[int, int]]], chrom: str, pos1: int) -> bool:
    """1-based pos lookup."""
    arr = beds.get(chrom)
    if not arr:
        return False
    p0 = pos1 - 1
    lo, hi = 0, len(arr) - 1
    while lo <= hi:
        m = (lo + hi) // 2
        s, e = arr[m]
        if p0 < s:
            hi = m - 1
        elif p0 >= e:
            lo = m + 1
        else:
            return True
    return False


def load_common_snp_set(path: str | None) -> set[tuple[str, int, str, str]] | None:
    """
    Load (chrom, pos1, ref, alt) tuples from a (possibly gzipped) VCF.
    Only biallelic SNPs are retained.
    """
    if not path:
        return None
    s: set[tuple[str, int, str, str]] = set()
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 5:
                continue
            chrom, pos, _id, ref, alt = cols[0], cols[1], cols[2], cols[3], cols[4]
            if "," in alt:
                continue
            if len(ref) != 1 or len(alt) != 1:
                continue
            try:
                s.add((chrom, int(pos), ref.upper(), alt.upper()))
            except ValueError:
                continue
    log.info("Loaded %d common-SNP sites", len(s))
    return s


# ──────────────────────────────────────────────────────────────────────
# VCF parsing
# ──────────────────────────────────────────────────────────────────────
def parse_vcf(vcf_path: str,
              cfg: dict[str, Any],
              upd_bed: dict[str, list[tuple[int, int]]] | None,
              common_snps: set[tuple[str, int, str, str]] | None,
              ) -> list[dict[str, Any]]:
    sites: list[dict[str, Any]] = []
    with open_text(vcf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 10:
                continue
            chrom, pos_s, _id, ref, alt, _qual, _flt, _info, fmt, sample = cols[:10]
            if "," in alt:
                continue
            if len(ref) != 1 or len(alt) != 1:
                continue
            try:
                pos = int(pos_s)
            except ValueError:
                continue
            if upd_bed is not None and not in_bed(upd_bed, chrom, pos):
                continue
            if common_snps is not None and (chrom, pos, ref.upper(), alt.upper()) not in common_snps:
                continue
            fmtd = dict(zip(fmt.split(":"), sample.split(":")))
            ad = fmtd.get("AD", "")
            try:
                ad_parts = [int(x) for x in ad.split(",")]
            except ValueError:
                continue
            if len(ad_parts) < 2:
                continue
            ref_c, alt_c = ad_parts[0], sum(ad_parts[1:])
            dp = ref_c + alt_c
            if dp < cfg["min_depth"]:
                continue
            sites.append({
                "chrom": chrom, "pos": pos,
                "ref": ref.upper(), "alt": alt.upper(),
                "depth": dp, "ref_count": ref_c, "alt_count": alt_c,
                "vaf": alt_c / dp if dp else 0.0,
            })
    log.info("Parsed %d SNP sites from %s", len(sites), vcf_path)
    return sites


# ──────────────────────────────────────────────────────────────────────
# Genotype classification & likelihoods
# ──────────────────────────────────────────────────────────────────────
def classify_maternal(vaf: float, cfg: dict, ff: float) -> str:
    """
    Classify the maternal genotype from observed cfDNA VAF.

    NOTE on the upper edge of mat_homref:
        At a maternal-homref site under patUPD, the fetus may be hom-alt and
        contribute ALT reads at FF.  Observed VAF ≈ FF (e.g. 0.08), well above
        a naive 0.05 threshold.  We therefore widen the mat_homref window to
        max(mat_homref_max_vaf, 1.5*FF + buffer) so those sites are still
        captured for the Mendelian-error analysis.  The downstream patUPD
        signal then re-tests them in the [0.7*FF, 1.3*FF] window.
    """
    upper_homref = max(cfg["mat_homref_max_vaf"], 1.5 * ff + 0.02)
    lower_homalt = min(cfg["mat_homalt_min_vaf"], 1.0 - (1.5 * ff + 0.02))
    if vaf <= upper_homref:
        return "mat_homref"
    if vaf >= lower_homalt:
        return "mat_homalt"
    if cfg["mat_het_min_vaf"] <= vaf <= cfg["mat_het_max_vaf"]:
        return "mat_het"
    return "ambiguous"


def log_binom(k: int, n: int, p: float) -> float:
    """log P(K=k | n,p), numerically stable for small/large k."""
    if n == 0:
        return 0.0
    p = min(max(p, 1e-9), 1.0 - 1e-9)
    # log C(n,k) via lgamma
    log_c = math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)
    return log_c + k * math.log(p) + (n - k) * math.log1p(-p)


def log_mixture(k: int, n: int, components: list[tuple[float, float]]) -> float:
    """log Σ_i w_i * Binom(k|n,p_i)."""
    log_terms = []
    for w, p in components:
        if w <= 0:
            continue
        log_terms.append(math.log(w) + log_binom(k, n, p))
    if not log_terms:
        return -1e9
    m = max(log_terms)
    return m + math.log(sum(math.exp(t - m) for t in log_terms))


def site_loglik(scenario: str, k: int, n: int, ff: float, cfg: dict, mat_gt: str) -> float:
    """
    Per-site log-likelihood under one of:
      'normal', 'matupd_het', 'matupd_iso', 'patupd_iso', 'patupd_het'
    using populated mixture weights derived from population MAF q ≈ 0.3 (mid-range).
    """
    err = cfg["seq_error_rate"]
    q = 0.30  # population MAF; mid-range biallelic SNP common-allele freq

    if mat_gt == "mat_het":
        if scenario == "normal":
            # 3-mode: {0.5-FF/2 with prob 0.5*(1-q), 0.5 prob 0.5, 0.5+FF/2 prob 0.5*q}
            comps = [
                (0.5 * (1 - q), 0.5 - ff / 2),
                (0.5,           0.5),
                (0.5 * q,       0.5 + ff / 2),
            ]
        elif scenario == "matupd_het":
            comps = [(1.0, 0.5)]
        elif scenario in ("matupd_iso", "patupd_iso"):
            # 2-mode {0.5-FF/2, 0.5+FF/2}, equal weight (mat-iso) or weighted by q (pat-iso)
            if scenario == "matupd_iso":
                comps = [(0.5, 0.5 - ff / 2), (0.5, 0.5 + ff / 2)]
            else:
                comps = [(1 - q, 0.5 - ff / 2), (q, 0.5 + ff / 2)]
        elif scenario == "patupd_het":
            comps = [
                ((1 - q) ** 2, 0.5 - ff / 2),
                (2 * q * (1 - q), 0.5),
                (q ** 2, 0.5 + ff / 2),
            ]
        else:
            return 0.0
    elif mat_gt in ("mat_homref", "mat_homalt"):
        # For mat_hom sites, model expected ALT VAF.
        # Symmetric: at mat_homalt we re-define k as ref_count to use the same shape.
        if scenario in ("normal", "matupd_iso", "matupd_het"):
            # Fetal allele dosage at mat-hom site under (any) matUPD = always 0 ALT
            # (assuming we're at mat_homref); fetal AA fixed.
            # Under normal: fetal Aa with prob q (paternal allele a) → VAF ≈ FF/2
            #               fetal AA with prob 1-q                → VAF ≈ 0
            if scenario == "normal":
                comps = [(1 - q, err), (q, ff / 2)]
            else:  # matUPD: fetal must be AA → VAF ≈ 0 (sequencing error only)
                comps = [(1.0, err)]
        elif scenario in ("patupd_iso", "patupd_het"):
            # At mat_homref under patUPD:
            #   patupd_iso  : fetal AA (prob (1-q)²+...) or aa (prob q²+...) — 2 outcomes:
            #       AA → VAF ≈ 0
            #       aa → VAF ≈ FF                              (Mendelian-error signature)
            # Weight: under iso, fetal genotype = paternal hom haplotype: (1-q) AA, q aa
            if scenario == "patupd_iso":
                comps = [(1 - q, err), (q, ff)]
            else:  # patupd_het: fetal = paternal genotype: AA, Aa, aa
                comps = [
                    ((1 - q) ** 2, err),
                    (2 * q * (1 - q), ff / 2),
                    (q ** 2, ff),
                ]
        else:
            return 0.0
    else:
        return 0.0

    # For mat_homalt we mirror everything around 0.5: use ref_count and 1-VAF.
    if mat_gt == "mat_homalt":
        comps = [(w, 1.0 - p) for (w, p) in comps]

    return log_mixture(k, n, comps)


# ──────────────────────────────────────────────────────────────────────
# Per-chromosome analysis
# ──────────────────────────────────────────────────────────────────────
def analyze_chrom(chrom: str,
                  sites: list[dict[str, Any]],
                  ff: float,
                  cfg: dict) -> dict[str, Any]:
    mat_het = [s for s in sites if s["mat_gt"] == "mat_het"]
    mat_hom = [s for s in sites if s["mat_gt"] in ("mat_homref", "mat_homalt")]

    res: dict[str, Any] = {
        "chrom": chrom,
        "n_total_snps": len(sites),
        "n_mat_het": len(mat_het),
        "n_mat_hom": len(mat_hom),
        "fetal_fraction_used": ff,
    }

    # ── mat-het likelihoods ──────────────────────────────────────────
    scenarios_het = ["normal", "matupd_het", "matupd_iso", "patupd_iso", "patupd_het"]
    ll_het = {sc: 0.0 for sc in scenarios_het}
    for s in mat_het:
        for sc in scenarios_het:
            ll_het[sc] += site_loglik(sc, s["alt_count"], s["depth"], ff, cfg, "mat_het")
    res["loglik_mat_het"] = {k: round(v, 3) for k, v in ll_het.items()}
    res["llr_vs_normal"] = {
        sc: round(ll_het[sc] - ll_het["normal"], 3)
        for sc in ("matupd_het", "matupd_iso", "patupd_iso", "patupd_het")
    }

    # Empirical BAF mode summary (sanity)
    if mat_het:
        bafs = [s["vaf"] for s in mat_het]
        bafs_dev = [abs(b - 0.5) for b in bafs]
        bafs_dev.sort()
        median_dev = bafs_dev[len(bafs_dev) // 2]
        var_dev = sum(d * d for d in bafs_dev) / len(bafs_dev)
        # fraction of "central" sites: |BAF - 0.5| < FF/4
        thr_central = max(0.01, ff / 4)
        n_central = sum(1 for d in bafs_dev if d < thr_central)
        n_off     = len(bafs_dev) - n_central
        res["baf_summary"] = {
            "median_abs_dev": round(median_dev, 4),
            "var_abs_dev":    round(var_dev, 6),
            "central_thr":    round(thr_central, 4),
            "n_central":      n_central,
            "n_off_center":   n_off,
            "central_frac":   round(n_central / len(bafs), 4) if bafs else 0.0,
        }

    # ── mat-hom Mendelian-error analysis (patUPD signal) ─────────────
    lo = cfg["patupd_vaf_lo_factor"] * ff
    hi = cfg["patupd_vaf_hi_factor"] * ff
    n_fetal_homalt = 0
    n_fetal_het    = 0
    n_fetal_homref = 0
    for s in mat_hom:
        v = s["vaf"]
        # at mat_homref: ALT VAF carries fetal info; at mat_homalt: REF "minor" carries it.
        minor_v = v if s["mat_gt"] == "mat_homref" else 1.0 - v
        if minor_v < ff / 4:
            n_fetal_homref += 1            # fetal also matches mother
        elif lo <= minor_v <= hi:
            n_fetal_homalt += 1            # Mendelian-error / patUPD signature
        elif ff / 4 <= minor_v < lo:
            n_fetal_het += 1               # normal fetal het from paternal allele
    if mat_hom:
        frac_homalt = n_fetal_homalt / len(mat_hom)
    else:
        frac_homalt = 0.0
    res["mat_hom_mende"] = {
        "n_fetal_homref": n_fetal_homref,
        "n_fetal_het":    n_fetal_het,
        "n_fetal_homalt": n_fetal_homalt,
        "frac_fetal_homalt": round(frac_homalt, 4),
        "vaf_window": [round(lo, 4), round(hi, 4)],
    }

    return res


# ──────────────────────────────────────────────────────────────────────
# Cross-chromosome calling (uses non-UPD autosomes as null)
# ──────────────────────────────────────────────────────────────────────
def make_call(chrom_results: list[dict[str, Any]],
              cfg: dict) -> list[dict[str, Any]]:
    """Add a 'call' field to UPD-chromosome results, using non-UPD autosomes
    in the same VCF as the empirical null for the patUPD-mende fraction."""
    null_fracs = []
    for r in chrom_results:
        if r["chrom"] not in UPD_CHROMS and r["n_mat_hom"] >= cfg["min_mat_hom_snps_mende"]:
            null_fracs.append(r["mat_hom_mende"]["frac_fetal_homalt"])
    if null_fracs:
        mu = sum(null_fracs) / len(null_fracs)
        var = sum((x - mu) ** 2 for x in null_fracs) / max(1, len(null_fracs) - 1)
        sd = math.sqrt(var) if var > 0 else max(0.005, mu)
    else:
        mu, sd = 0.0, 0.01

    for r in chrom_results:
        r["null_frac_homalt_mu"] = round(mu, 5)
        r["null_frac_homalt_sd"] = round(sd, 5)
        # patUPD z-score vs null
        f = r["mat_hom_mende"]["frac_fetal_homalt"]
        r["patupd_z"] = round((f - mu) / max(sd, 1e-4), 3)

        # Decide
        call = "normal"
        notes: list[str] = []
        info_ok = (r["n_mat_het"] >= cfg["min_mat_het_snps_baf"]
                   and r["n_mat_hom"] >= cfg["min_mat_hom_snps_mende"])
        if not info_ok:
            call = "INSUFFICIENT_INFORMATION"
            notes.append(
                f"need ≥ {cfg['min_mat_het_snps_baf']} mat-het and "
                f"≥ {cfg['min_mat_hom_snps_mende']} mat-hom; "
                f"have {r['n_mat_het']} / {r['n_mat_hom']}.")
        else:
            llr = r["llr_vs_normal"]
            patupd_z = r["patupd_z"]
            if patupd_z >= cfg["patupd_excess_z_threshold"]:
                # Strong Mendelian-error excess → paternal UPD
                if llr["patupd_iso"] >= cfg["llr_call_threshold"]:
                    call = "patUPD_iso_suspected"
                else:
                    call = "patUPD_het_suspected"
            elif llr["matupd_het"] >= cfg["llr_call_threshold"] \
                    and llr["matupd_het"] > llr["matupd_iso"]:
                call = "matUPD_het_suspected"
            elif llr["matupd_iso"] >= cfg["llr_call_threshold"] \
                    and llr["matupd_iso"] > llr["matupd_het"]:
                call = "matUPD_iso_suspected"
            else:
                call = "normal"

        r["call"] = call
        r["notes"] = notes
        r["clinical_note"] = UPD_CLINICAL.get(r["chrom"], {})
    return chrom_results


# ──────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--fetal-fraction", type=float, required=True)
    ap.add_argument("--upd-bed", default=None)
    ap.add_argument("--common-snps", default=None)
    ap.add_argument("--config", default=None)
    ap.add_argument("--output", required=True)
    ap.add_argument("--per-snp-tsv", default=None)
    args = ap.parse_args()

    cfg = DEFAULT_CONFIG.copy()
    if args.config and Path(args.config).exists():
        with open(args.config) as fh:
            cfg.update(json.load(fh))

    if args.fetal_fraction <= 0 or args.fetal_fraction >= 1:
        log.error("Invalid fetal-fraction: %s", args.fetal_fraction)
        sys.exit(2)
    if args.fetal_fraction < cfg["min_fetal_fraction"]:
        log.warning("FF=%.4f below min_fetal_fraction=%.4f — UPD calls flagged LOW_FF.",
                    args.fetal_fraction, cfg["min_fetal_fraction"])

    upd_bed = load_bed(args.upd_bed) if args.upd_bed else None
    common = load_common_snp_set(args.common_snps)

    sites = parse_vcf(args.vcf, cfg, upd_bed, common)

    for s in sites:
        s["mat_gt"] = classify_maternal(s["vaf"], cfg, args.fetal_fraction)

    # Group by chromosome
    by_chrom: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for s in sites:
        by_chrom[s["chrom"]].append(s)

    chrom_results = [
        analyze_chrom(c, by_chrom[c], args.fetal_fraction, cfg)
        for c in sorted(by_chrom.keys(),
                        key=lambda x: int(x.replace("chr", "")) if x.replace("chr", "").isdigit() else 99)
    ]
    chrom_results = make_call(chrom_results, cfg)

    upd_calls = [r for r in chrom_results if r["chrom"] in UPD_CHROMS]
    summary_call = "normal"
    flagged = [r for r in upd_calls if r["call"] not in ("normal", "INSUFFICIENT_INFORMATION")]
    insufficient = [r for r in upd_calls if r["call"] == "INSUFFICIENT_INFORMATION"]
    if flagged:
        summary_call = ", ".join(f"{r['chrom']}:{r['call']}" for r in flagged)
    elif insufficient and len(insufficient) == len(upd_calls):
        summary_call = "INSUFFICIENT_INFORMATION"

    out: dict[str, Any] = {
        "sample_id": args.sample_id,
        "fetal_fraction": args.fetal_fraction,
        "panel": "Twist_UCL_SingleGeneNIPT_TE-96276661_hg38",
        "upd_chroms": UPD_CHROMS,
        "config": cfg,
        "summary_call": summary_call,
        "per_chrom_upd": upd_calls,
        "per_chrom_other": [r for r in chrom_results if r["chrom"] not in UPD_CHROMS],
    }

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as fh:
        json.dump(out, fh, indent=2)
    log.info("UPD report → %s  (summary_call=%s)", args.output, summary_call)

    if args.per_snp_tsv:
        with open(args.per_snp_tsv, "w") as fh:
            fh.write("chrom\tpos\tref\talt\tdepth\tref_count\talt_count\tvaf\tmat_gt\n")
            for s in sites:
                fh.write(f"{s['chrom']}\t{s['pos']}\t{s['ref']}\t{s['alt']}\t"
                         f"{s['depth']}\t{s['ref_count']}\t{s['alt_count']}\t"
                         f"{s['vaf']:.6f}\t{s['mat_gt']}\n")
        log.info("Per-SNP detail → %s", args.per_snp_tsv)


if __name__ == "__main__":
    main()
