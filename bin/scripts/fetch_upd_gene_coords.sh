#!/usr/bin/env bash
# fetch_upd_gene_coords.sh — Build the gene-coordinate TSV that
# build_upd_regions_bed.py expects.
#
# Pulls UCSC refFlat (hg38) from hgdownload.soe.ucsc.edu, keeps only the
# 7 UPD-relevant chromosomes, and aggregates per-gene to (gene, chrom,
# min(txStart), max(txEnd)) so that all isoforms collapse into a single
# locus.
#
# Output format (tab-separated, BED-style 0-based coordinates):
#     gene<TAB>chrom<TAB>start<TAB>end
#
# Usage:
#     bin/scripts/fetch_upd_gene_coords.sh [output_path]
# Default output: data/upd/hg38.upd_genes.tsv
#
# Requires: curl (or wget), gunzip, awk, sort.

set -euo pipefail

OUT="${1:-data/upd/hg38.upd_genes.tsv}"
URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz"
CHROMS_RE='^(chr6|chr7|chr11|chr14|chr15|chr16|chr20)$'

mkdir -p "$(dirname "$OUT")"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

echo "[fetch_upd_gene_coords] downloading: $URL"
if command -v curl >/dev/null 2>&1; then
    curl -fSL --retry 3 -o "$TMP/refFlat.txt.gz" "$URL"
elif command -v wget >/dev/null 2>&1; then
    wget -q -O "$TMP/refFlat.txt.gz" "$URL"
else
    echo "[ERROR] need curl or wget" >&2
    exit 1
fi

gunzip -c "$TMP/refFlat.txt.gz" > "$TMP/refFlat.txt"
echo "[fetch_upd_gene_coords] refFlat lines: $(wc -l < "$TMP/refFlat.txt")"

# refFlat columns:
#   1=geneName 2=name 3=chrom 4=strand 5=txStart 6=txEnd
#   7=cdsStart 8=cdsEnd 9=exonCount 10=exonStarts 11=exonEnds
#
# Aggregate per (geneName, chrom): min(txStart), max(txEnd)
# Drops alt-haplotype contigs by chromosome regex.
awk -v OFS='\t' -v re="$CHROMS_RE" '
    $3 ~ re {
        key = $1 "\t" $3
        s = $5 + 0
        e = $6 + 0
        if (!(key in mn) || s < mn[key]) mn[key] = s
        if (!(key in mx) || e > mx[key]) mx[key] = e
    }
    END {
        for (k in mn) {
            split(k, a, "\t")
            print a[1], a[2], mn[k], mx[k]
        }
    }
' "$TMP/refFlat.txt" \
    | sort -k2,2V -k3,3n -k1,1 \
    | (echo -e "gene\tchrom\tstart\tend"; cat) \
    > "$OUT"

N=$(awk 'NR>1' "$OUT" | wc -l)
echo "[fetch_upd_gene_coords] wrote $OUT  (genes: $N across UPD chroms)"
