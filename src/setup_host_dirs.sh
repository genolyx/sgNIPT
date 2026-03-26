#!/bin/bash
###############################################################################
#  Single Gene NIPT - Host Directory Setup Script
###############################################################################
#  Creates the required directory structure on the host machine.
#
#  Usage:
#    ./setup_host_dirs.sh /data/SgNIPT
###############################################################################

set -euo pipefail

ROOT_DIR="${1:?'Usage: ./setup_host_dirs.sh <ROOT_DIR>'}"

echo "Creating SgNIPT directory structure under: ${ROOT_DIR}"

# Main directories
mkdir -p "${ROOT_DIR}/fastq"
mkdir -p "${ROOT_DIR}/data/reference"
mkdir -p "${ROOT_DIR}/data/bwa_index"
mkdir -p "${ROOT_DIR}/data/target_bed"
mkdir -p "${ROOT_DIR}/data/annotation"
mkdir -p "${ROOT_DIR}/config"
mkdir -p "${ROOT_DIR}/analysis"
mkdir -p "${ROOT_DIR}/output"
mkdir -p "${ROOT_DIR}/log"

echo ""
echo "Directory structure created:"
echo ""
tree "${ROOT_DIR}" 2>/dev/null || find "${ROOT_DIR}" -type d | sort | sed 's|[^/]*/|  |g'
echo ""
echo "Next steps:"
echo "  1. Place hg38 reference FASTA in:  ${ROOT_DIR}/data/reference/"
echo "  2. Place BWA-MEM2 index files in:  ${ROOT_DIR}/data/bwa_index/"
echo "  3. Place target BED file in:       ${ROOT_DIR}/data/target_bed/"
echo "  4. Place FASTQ files in:           ${ROOT_DIR}/fastq/{YYMM}/{SAMPLE}/"
echo "     (e.g. fastq/2603/NA12878_Twist_Exome/*.fastq.gz)"
echo "  5. (Optional) Place ClinVar VCF:   ${ROOT_DIR}/data/annotation/clinvar.vcf.gz"
echo "  6. (Optional) Custom config:       ${ROOT_DIR}/config/params.yaml"
echo "  7. (Optional) QC thresholds:       ${ROOT_DIR}/config/qc_thresholds.json"
echo ""
