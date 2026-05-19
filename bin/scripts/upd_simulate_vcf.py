#!/usr/bin/env python3
"""
upd_simulate_vcf.py — Synthesise a cfDNA-style VCF for UPD detection unit tests.

Generates a VCF with FORMAT/AD,DP that mimics what bcftools mpileup+call would
produce on a maternal-cfDNA mixture, for one or more UPD scenarios on one or
more chromosomes.

Useful for end-to-end testing of upd_detection.py without needing a real BAM:

  python3 bin/scripts/upd_simulate_vcf.py \
      --output /tmp/sim.vcf \
      --fetal-fraction 0.08 \
      --depth 800 \
      --n-snps-per-chrom 400 \
      --scenario chr15:matupd_iso \
      --scenario chr7:matupd_het \
      --scenario chr11:patupd_iso \
      --background-chroms chr1,chr2,chr3,chr4,chr5

  python3 bin/scripts/upd_detection.py \
      --sample-id sim_upd \
      --vcf /tmp/sim.vcf \
      --fetal-fraction 0.08 \
      --output /tmp/sim.upd.json \
      --per-snp-tsv /tmp/sim.upd.tsv

Scenarios
---------
  normal        biparental (default for chroms not in --scenario)
  matupd_iso    maternal isodisomy
  matupd_het    maternal heterodisomy
  patupd_iso    paternal isodisomy
  patupd_het    paternal heterodisomy
"""

from __future__ import annotations

import argparse
import random
import sys
from pathlib import Path

SCENARIOS = ("normal", "matupd_iso", "matupd_het", "patupd_iso", "patupd_het")


def sample_mat_genotype(rng: random.Random, q: float) -> str:
    """Pick maternal genotype from HWE with MAF q."""
    r = rng.random()
    if r < (1 - q) ** 2:
        return "AA"
    if r < (1 - q) ** 2 + 2 * q * (1 - q):
        return "Aa"
    return "aa"


def fetal_dosage_from_scenario(scenario: str, mat_gt: str,
                               q: float, rng: random.Random) -> tuple[float, float]:
    """
    Return (fetal_alt_dosage_in_[0,1], expected ALT VAF central value).

    fetal_alt_dosage = mean ALT dosage of the fetal contribution at this site
                       (0=hom-ref, 0.5=het, 1=hom-alt for the fetal genotype).
    """
    if scenario == "normal":
        # Fetal: one allele from mom (random), one from population (q)
        mat_alleles = {"AA": (0, 0), "Aa": (0, 1), "aa": (1, 1)}[mat_gt]
        f_mat = rng.choice(mat_alleles)
        f_pat = 1 if rng.random() < q else 0
        return (f_mat + f_pat) / 2, None  # type: ignore
    if scenario == "matupd_iso":
        # Fetal = homozygous for one of mom's alleles
        mat_alleles = {"AA": (0, 0), "Aa": (0, 1), "aa": (1, 1)}[mat_gt]
        f_a = rng.choice(mat_alleles)
        return float(f_a), None  # type: ignore
    if scenario == "matupd_het":
        mat_alleles = {"AA": (0, 0), "Aa": (0, 1), "aa": (1, 1)}[mat_gt]
        # Fetal = same as mom (both homologs)
        return (mat_alleles[0] + mat_alleles[1]) / 2, None  # type: ignore
    if scenario == "patupd_iso":
        f_a = 1 if rng.random() < q else 0
        return float(f_a), None  # type: ignore
    if scenario == "patupd_het":
        a1 = 1 if rng.random() < q else 0
        a2 = 1 if rng.random() < q else 0
        return (a1 + a2) / 2, None  # type: ignore
    raise ValueError(scenario)


def cfdna_alt_count(mat_gt: str, fetal_alt_dose: float, ff: float, depth: int,
                    err: float, rng: random.Random) -> int:
    """Sample cfDNA ALT read count given maternal GT, fetal dosage, FF, depth."""
    # Maternal alt fraction
    mat_alt_frac = {"AA": 0.0, "Aa": 0.5, "aa": 1.0}[mat_gt]
    # cfDNA alt fraction
    p = (1 - ff) * mat_alt_frac + ff * fetal_alt_dose
    # Add small symmetric error noise: small chance of flipping
    p = p * (1 - err) + (1 - p) * err
    p = min(max(p, 1e-6), 1 - 1e-6)
    # Sample alt count
    return sum(1 for _ in range(depth) if rng.random() < p)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", required=True, help="Output VCF path")
    ap.add_argument("--sample-id", default="sim_cfDNA")
    ap.add_argument("--fetal-fraction", type=float, default=0.08)
    ap.add_argument("--depth", type=int, default=600,
                    help="Mean depth per site (Poisson-jittered)")
    ap.add_argument("--maf", type=float, default=0.30,
                    help="Per-site population minor-allele freq (HWE)")
    ap.add_argument("--n-snps-per-chrom", type=int, default=400)
    ap.add_argument("--seq-error", type=float, default=0.005)
    ap.add_argument("--scenario", action="append", default=[],
                    help="<chrom>:<scenario>, repeatable; e.g. chr15:matupd_iso")
    ap.add_argument("--background-chroms", default="chr1,chr2,chr3,chr4,chr5",
                    help="Comma-separated chroms simulated as 'normal' (null)")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    rng = random.Random(args.seed)

    chrom_scenarios: dict[str, str] = {}
    for spec in args.scenario:
        if ":" not in spec:
            sys.exit(f"[ERROR] --scenario must be CHROM:NAME, got {spec!r}")
        c, name = spec.split(":", 1)
        if name not in SCENARIOS:
            sys.exit(f"[ERROR] unknown scenario {name!r} (choose from {SCENARIOS})")
        chrom_scenarios[c] = name
    for c in args.background_chroms.split(","):
        c = c.strip()
        if c and c not in chrom_scenarios:
            chrom_scenarios[c] = "normal"

    # Synthesise positions: 1 SNP every (chrom_length / n) bp; chrom_length=10 Mb default.
    chrom_lengths = {c: 10_000_000 for c in chrom_scenarios}

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write('##INFO=<ID=SCN,Number=1,Type=String,Description="Simulated scenario">\n')
        fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype (maternal)">\n')
        fh.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">\n')
        fh.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n')
        fh.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{args.sample_id}\n")

        for chrom, scenario in chrom_scenarios.items():
            n = args.n_snps_per_chrom
            length = chrom_lengths[chrom]
            step = max(1, length // (n + 1))
            pos = step
            for _ in range(n):
                # Sample
                mat_gt = sample_mat_genotype(rng, args.maf)
                fetal_dose, _ = fetal_dosage_from_scenario(scenario, mat_gt, args.maf, rng)
                # Depth ~ Poisson-ish (uniform jitter)
                d = max(20, int(args.depth * (0.7 + 0.6 * rng.random())))
                alt_c = cfdna_alt_count(mat_gt, fetal_dose, args.fetal_fraction,
                                        d, args.seq_error, rng)
                ref_c = d - alt_c
                # GT field: maternal (caller would emit dominant)
                gt_str = {"AA": "0/0", "Aa": "0/1", "aa": "1/1"}[mat_gt]
                fh.write(
                    f"{chrom}\t{pos}\t.\tA\tG\t.\tPASS\tSCN={scenario}\t"
                    f"GT:AD:DP\t{gt_str}:{ref_c},{alt_c}:{d}\n"
                )
                pos += step

    print(f"[upd_simulate_vcf] wrote {out}", file=sys.stderr)
    print("  scenarios:", chrom_scenarios, file=sys.stderr)


if __name__ == "__main__":
    main()
