#!/usr/bin/env bash
# fetch_upd_common_snps.sh — Build the common-SNP VCF that upd_detection.py
# uses to suppress sequencing-error noise.
#
# Streams gnomAD genomes per-chromosome via remote tabix, restricts to the
# UPD analysis BED, keeps only biallelic SNPs with INFO/AF in [0.05, 0.95],
# and concatenates into one tabix-indexed VCF.gz.
#
# Output:  data/upd/common_snps_maf05.upd.vcf.gz (+ .tbi)
#
# Usage:
#     bin/scripts/fetch_upd_common_snps.sh \
#         [data/upd/upd_analysis_regions.bed] \
#         [data/upd/common_snps_maf05.upd.vcf.gz]
#
# Requires: bcftools (≥1.10) with htslib remote-https tabix support.

set -euo pipefail

BED="${1:-data/upd/upd_analysis_regions.bed}"
OUT="${2:-data/upd/common_snps_maf05.upd.vcf.gz}"

if [[ ! -f "$BED" ]]; then
    echo "[ERROR] BED not found: $BED" >&2
    exit 1
fi

GNOMAD_BASE="https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes"
CHROMS=("chr6" "chr7" "chr11" "chr14" "chr15" "chr16" "chr20")

mkdir -p "$(dirname "$OUT")"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

echo "[fetch_upd_common_snps] BED: $BED  →  OUT: $OUT"

for chr in "${CHROMS[@]}"; do
    if ! awk -v c="$chr" '$1==c {found=1; exit} END{exit !found}' "$BED"; then
        echo "  [skip] $chr (no intervals in BED)"
        continue
    fi
    awk -v c="$chr" '$1==c {print $1"\t"$2"\t"$3}' "$BED" > "$TMP/${chr}.bed"
    URL="${GNOMAD_BASE}/gnomad.genomes.v3.1.2.sites.${chr}.vcf.bgz"
    echo "  [fetch] $chr  ←  $URL"
    bcftools view \
        -m2 -M2 -v snps \
        -i 'INFO/AF>=0.05 && INFO/AF<=0.95' \
        -R "$TMP/${chr}.bed" \
        "$URL" \
        -Oz -o "$TMP/${chr}.vcf.gz"
    bcftools index -t "$TMP/${chr}.vcf.gz"
    n=$(bcftools view -H "$TMP/${chr}.vcf.gz" | wc -l)
    echo "         kept ${n} common SNPs"
done

echo "[fetch_upd_common_snps] concatenating ..."
ls "$TMP"/chr*.vcf.gz >/dev/null 2>&1 || { echo "[ERROR] no per-chrom slices created"; exit 1; }
bcftools concat -Oz -o "$OUT" \
    $(ls "$TMP"/chr*.vcf.gz | sort -V)
bcftools index -t "$OUT"

TOTAL=$(bcftools view -H "$OUT" | wc -l)
SIZE=$(du -h "$OUT" | cut -f1)
echo "[fetch_upd_common_snps] wrote $OUT  (${TOTAL} sites, ${SIZE})"
