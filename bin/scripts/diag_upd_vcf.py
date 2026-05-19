#!/usr/bin/env python3
"""Diagnose VAF distribution in a bcftools-called VCF for UPD analysis."""
import sys, subprocess, json

def scan_vcf(vcf_path, ff=0.10):
    mat_hom_ref, mat_het, mat_hom_alt, other = [], [], [], []
    homref_max = max(0.05, 1.5 * ff + 0.02)

    result = subprocess.run(
        ["bcftools", "view", "-H", vcf_path],
        capture_output=True, text=True
    )
    for line in result.stdout.splitlines():
        cols = line.split("\t")
        if len(cols) < 10:
            continue
        fmt = dict(zip(cols[8].split(":"), cols[9].split(":")))
        ad = [int(x) for x in fmt.get("AD", "0,0").split(",")]
        dp = sum(ad)
        if dp < 50:
            continue
        vaf = ad[1] / dp if len(ad) > 1 else 0.0
        if vaf <= homref_max:
            mat_hom_ref.append(vaf)
        elif vaf >= (1.0 - 0.05):
            mat_hom_alt.append(vaf)
        elif 0.35 <= vaf <= 0.65:
            mat_het.append(vaf)
        else:
            other.append(vaf)

    print(f"  VCF: {vcf_path}")
    print(f"  mat_hom_ref (VAF<=0.{homref_max:.2f}): n={len(mat_hom_ref):3d}  "
          f"median={sorted(mat_hom_ref)[len(mat_hom_ref)//2]:.3f}" if mat_hom_ref else
          f"  mat_hom_ref: n=  0")
    print(f"  mat_het     (0.35-0.65):  n={len(mat_het):3d}  "
          f"median={sorted(mat_het)[len(mat_het)//2]:.3f}" if mat_het else
          f"  mat_het    : n=  0")
    print(f"  mat_hom_alt (VAF>=0.95):  n={len(mat_hom_alt):3d}")
    print(f"  ambiguous:                n={len(other):3d}  "
          f"examples={sorted(other)[:5]}")
    return len(mat_hom_ref), len(mat_het)

if __name__ == "__main__":
    vcf = sys.argv[1]
    ff  = float(sys.argv[2]) if len(sys.argv) > 2 else 0.10
    scan_vcf(vcf, ff)
