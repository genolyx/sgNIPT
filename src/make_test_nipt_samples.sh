#!/usr/bin/env bash
# =============================================================================
# make_test_nipt_samples.sh — Create artificial NIPT samples for pipeline testing
#
# Overview
# --------
# Takes the target.bam from a completed sgNIPT run (the "maternal" BAM) and
# creates simulated cfDNA BAMs with spiked disease variants and FF-informative
# SNPs for each requested fetal fraction (5%, 10%, 15% by default).
#
# The simulated BAMs are used directly (no FASTQ conversion!) via the pipeline's
# BAM input mode (--input_bam):
#
#   FASTQ mode (default):   FASTQ → FASTP → BWA-MEM2 → MARK_DUPLICATES → EXTRACT_TARGET → …
#   BAM mode (--input_bam): [skips above]  →  BAM QC → FF estimation → Variant calling → …
#
# Why target.bam instead of dedup.bam?
#   After EXTRACT_TARGET_READS, the pipeline's ch_bam IS the target.bam.
#   Spiking reads into the target.bam means they bypass duplicate marking,
#   which would otherwise remove cloned reads with identical 5' coordinates.
#
# Usage
# -----
#   src/make_test_nipt_samples.sh [OPTIONS]
#
# Options
#   --source-order   ORDER_ID of completed run to steal target.bam from
#                    (default: Test_sgNIPT)
#   --work-id        YYMM work ID (default: 2603)
#   --ff             Space-separated FF percentages (default: "5 10 15")
#   --seed           Random seed for pysam spiking (default: 42)
#   --docker-image   Docker image name (default: sgnipt)
#   --dry-run        Print commands without executing
#   -h / --help      Show this message
#
# Output
# ------
#   data/test/sim_bam/FF<N>_sim_<SOURCE>.target.bam   — simulated BAM per FF
#   data/test/sim_bam/FF<N>_sim_<SOURCE>.bam_samplesheet.csv  — Nextflow BAM CSV
#
# Run commands printed at the end of execution.
# =============================================================================
set -euo pipefail

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
SOURCE_ORDER="Test_sgNIPT"
WORK_ID="2603"
FF_VALUES="5 10 15"
DOCKER_IMAGE="sgnipt"
DRY_RUN=false
SEED=42

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case $1 in
        --source-order) SOURCE_ORDER="$2"; shift 2 ;;
        --work-id)      WORK_ID="$2";      shift 2 ;;
        --ff)           FF_VALUES="$2";    shift 2 ;;
        --seed)         SEED="$2";         shift 2 ;;
        --dry-run)      DRY_RUN=true;      shift   ;;
        --docker-image) DOCKER_IMAGE="$2"; shift 2 ;;
        -h|--help)
            sed -n '2,40p' "$0"
            exit 0 ;;
        *)
            echo "[ERROR] Unknown argument: $1" >&2; exit 1 ;;
    esac
done

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
# Use the target.bam (already extracted to panel regions) from a previous run.
# This BAM enters the pipeline at Step 5 (BAM QC) via --input_bam mode,
# bypassing FASTP, alignment, dedup, and extract_target steps.
SOURCE_BAM="${REPO_ROOT}/output/${WORK_ID}/${SOURCE_ORDER}/${SOURCE_ORDER}/alignment/${SOURCE_ORDER}.target.bam"
SOURCE_BAI="${SOURCE_BAM}.bai"
FF_VCF="${REPO_ROOT}/output/${WORK_ID}/${SOURCE_ORDER}/${SOURCE_ORDER}/fetal_fraction/${SOURCE_ORDER}.ff_variants.vcf"
VARIANTS_TSV="${REPO_ROOT}/data/test/test_variants.tsv"
SIM_SCRIPT="${REPO_ROOT}/bin/scripts/simulate_nipt_bam.py"
SIM_BAM_DIR="${REPO_ROOT}/data/test/sim_bam"

echo "============================================================"
echo " sgNIPT Artificial Sample Generator"
echo "============================================================"
echo " Source order    : ${SOURCE_ORDER}"
echo " Source BAM      : ${SOURCE_BAM}"
echo " FF values       : ${FF_VALUES}"
echo " Variants TSV    : ${VARIANTS_TSV}"
echo " FF SNPs VCF     : ${FF_VCF}"
echo " Docker image    : ${DOCKER_IMAGE}"
echo " Dry run         : ${DRY_RUN}"
echo "============================================================"

if [[ ! -f "${SOURCE_BAM}" ]]; then
    echo "[ERROR] Source target.bam not found: ${SOURCE_BAM}"
    echo "        Run a full sgNIPT pipeline first on '${SOURCE_ORDER}'."
    exit 1
fi

if [[ ! -f "${SOURCE_BAI}" ]]; then
    echo "[WARN] BAM index not found: ${SOURCE_BAI}"
    echo "       Indexing source BAM …"
    if ! ${DRY_RUN}; then
        docker run --rm --user "$(id -u):$(id -g)" \
            -v "${SOURCE_BAM}:/target.bam" \
            -v "$(dirname "${SOURCE_BAM}"):/bamdir" \
            "${DOCKER_IMAGE}" \
            samtools index /bamdir/"$(basename "${SOURCE_BAM}")"
    fi
fi

FF_VCF_MOUNT=""
FF_VCF_ARG=""
if [[ -f "${FF_VCF}" ]]; then
    FF_VCF_MOUNT="-v ${FF_VCF}:/ff_variants.vcf:ro"
    FF_VCF_ARG="--ff-vcf /ff_variants.vcf --n-ff-snps 100"
else
    echo "[WARN] FF VCF not found: ${FF_VCF}"
    echo "       FF-informative SNP spiking will be skipped."
    echo "       The SNP-based fetal fraction method will still report INSUFFICIENT_DATA."
fi

run_or_echo() {
    if ${DRY_RUN}; then
        echo "[DRY-RUN] $*"
    else
        "$@"
    fi
}

run_or_echo mkdir -p "${SIM_BAM_DIR}"

# ---------------------------------------------------------------------------
# Create one simulated BAM per FF value
# ---------------------------------------------------------------------------
GENERATED_SAMPLES=()

for FF_PCT in ${FF_VALUES}; do
    FF_FLOAT=$(python3 -c "print(${FF_PCT}/100)")
    SAMPLE_NAME="FF${FF_PCT}_sim_${SOURCE_ORDER}"
    SIM_BAM="${SIM_BAM_DIR}/${SAMPLE_NAME}.target.bam"
    BAM_CSV="${SIM_BAM_DIR}/${SAMPLE_NAME}.bam_samplesheet.csv"

    echo ""
    echo "------------------------------------------------------------"
    echo " Generating: ${SAMPLE_NAME}  (FF=${FF_PCT}%,  target VAF=$(python3 -c "print(f'{${FF_PCT}/200:.1%}')") per allele)"
    echo "------------------------------------------------------------"

    # ── Step 1: spike reads into target.bam (runs inside Docker) ──────────
    echo "[Step 1/2] Spiking reads with simulate_nipt_bam.py …"

    DOCKER_CMD=(
        docker run --rm
            --user "$(id -u):$(id -g)"
            -v "${SOURCE_BAM}:/input.bam:ro"
            -v "${SOURCE_BAI}:/input.bam.bai:ro"
            -v "${VARIANTS_TSV}:/variants.tsv:ro"
            -v "${SIM_SCRIPT}:/simulate_nipt_bam.py:ro"
            -v "${SIM_BAM_DIR}:/out"
    )
    if [[ -n "${FF_VCF_MOUNT}" ]]; then
        DOCKER_CMD+=(${FF_VCF_MOUNT})
    fi
    DOCKER_CMD+=(
        "${DOCKER_IMAGE}"
        python3 /simulate_nipt_bam.py
            --bam     /input.bam
            --output  "/out/${SAMPLE_NAME}.target.bam"
            --fetal-fraction "${FF_FLOAT}"
            --variants /variants.tsv
            ${FF_VCF_ARG}
            --seed "${SEED}"
    )

    if ${DRY_RUN}; then
        echo "[DRY-RUN] ${DOCKER_CMD[*]}"
    else
        "${DOCKER_CMD[@]}"
    fi

    # ── Step 2: create BAM samplesheet CSV (with container-internal paths) ────
    # HOST_DATA_DIR is mounted as /Work/SgNIPT/data inside the container.
    # sim_bam/ lives under data/test/sim_bam/ → accessible at the same relative path.
    echo "[Step 2/2] Writing BAM samplesheet: ${BAM_CSV}"
    CONTAINER_DATA_DIR="/Work/SgNIPT/data"
    # Compute container-internal BAM path by replacing the repo data/ prefix
    CONTAINER_SIM_BAM="${CONTAINER_DATA_DIR}/test/sim_bam/${SAMPLE_NAME}.target.bam"
    CONTAINER_SIM_BAI="${CONTAINER_SIM_BAM}.bai"

    if ${DRY_RUN}; then
        echo "[DRY-RUN] Write ${BAM_CSV}:"
        echo "          sample_id,bam,bai"
        echo "          ${SAMPLE_NAME},${CONTAINER_SIM_BAM},${CONTAINER_SIM_BAI}"
    else
        {
            echo "sample_id,bam,bai"
            echo "${SAMPLE_NAME},${CONTAINER_SIM_BAM},${CONTAINER_SIM_BAI}"
        } > "${BAM_CSV}"
    fi

    echo "[Done] ${SAMPLE_NAME}"
    GENERATED_SAMPLES+=("${SAMPLE_NAME}")
done

# ---------------------------------------------------------------------------
# Print run commands
# ---------------------------------------------------------------------------
echo ""
echo "============================================================"
echo " Simulated samples ready. Run commands (BAM mode):"
echo "============================================================"
echo ""
echo " The --input_bam flag tells the pipeline to skip:"
echo "   FASTP → BWA_MEM2_ALIGN → MARK_DUPLICATES → EXTRACT_TARGET_READS"
echo " and start directly at:"
echo "   BAM_QC → FF Estimation → Variant Calling → VEP → Report"
echo ""

for SAMPLE_NAME in "${GENERATED_SAMPLES[@]}"; do
    BAM_CSV="${SIM_BAM_DIR}/${SAMPLE_NAME}.bam_samplesheet.csv"
    # The BAM samplesheet will be mounted into the container at /Work/SgNIPT/analysis/bam_samplesheet.csv
    echo "  # ${SAMPLE_NAME}"
    echo "  src/run_sgnipt.sh \\"
    echo "      --order-id   ${SAMPLE_NAME} \\"
    echo "      --work-id    ${WORK_ID} \\"
    echo "      --input-bam  ${BAM_CSV} \\"
    echo "      --fresh"
    echo ""
done

echo "============================================================"
echo " After each run, verify:"
echo "   • fetal_fraction.json → primary_fetal_fraction ≈ expected FF"
echo "   • variant_report.json → disease variants at VAF ≈ FF/2"
echo "   • report PDF → fetal variants correctly classified"
echo "============================================================"
