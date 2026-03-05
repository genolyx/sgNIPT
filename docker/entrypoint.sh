#!/bin/bash
###############################################################################
#  Single Gene NIPT Pipeline - Docker Entrypoint
###############################################################################
#  This script is the ENTRYPOINT for the Docker container.
#  It receives Order-level arguments, generates a samplesheet,
#  configures paths, and launches the Nextflow pipeline.
#
#  Compatible with nipt-daemon monitoring pattern:
#    - Writes {ORDER_ID}.json + {ORDER_ID}.output.tar on success
#    - Writes {ORDER_ID}_progress.txt for real-time progress tracking
#    - daemon polls for JSON+TAR presence to detect completion
#
#  Panel: Twist UCL_SingleGeneNIPT_TE-96276661_hg38
#    - Auto-detects BED files in data directory:
#      * target_bed/All_target_regions*.bed  → Variant calling & QC
#      * target_bed/ff_snps*.bed             → Fetal Fraction SNP list
#      * target_bed/Probes_merged*.bed       → Probe coverage QC
#      * target_bed/Target_regions_with_zero_probes*.bed → Zero-probe flagging
#      * target_bed/gene_list*.tsv           → Gene list validation
#
#  Usage (inside docker run):
#    --order_id        ORDER_ID        (required)
#    --sample_name     SAMPLE_NAME     (required, comma-separated for multi-sample)
#    --fastq_r1        FASTQ_R1_FILES  (required, comma-separated)
#    --fastq_r2        FASTQ_R2_FILES  (optional, comma-separated)
#    --target_bed      TARGET_BED_NAME (optional, default: auto-detect)
#    --variant_caller  CALLER          (optional, default: bcftools)
#    --extra_args      "..."           (optional, extra Nextflow arguments)
###############################################################################

set -euo pipefail

# ── Constants ────────────────────────────────────────────────────────────────
PIPELINE_DIR="/Work/SgNIPT/pipeline"
FASTQ_DIR="/Work/SgNIPT/fastq"
DATA_DIR="/Work/SgNIPT/data"
CONFIG_DIR="/Work/SgNIPT/config"
ANALYSIS_DIR="/Work/SgNIPT/analysis"
OUTPUT_DIR="/Work/SgNIPT/output"
LOG_DIR="/Work/SgNIPT/log"

TIMESTAMP=$(date '+%Y%m%d_%H%M%S')
PIPELINE_VERSION="1.1.0"

# ── Logging functions ────────────────────────────────────────────────────────
log_info()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INFO]  $*" | tee -a "${LOG_DIR}/entrypoint.log"; }
log_warn()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [WARN]  $*" | tee -a "${LOG_DIR}/entrypoint.log"; }
log_error() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ERROR] $*" | tee -a "${LOG_DIR}/entrypoint.log"; }

# ── Progress tracking (daemon-compatible) ────────────────────────────────────
update_progress() {
    local step="$1"
    local message="$2"
    local progress_file="${ORDER_OUTPUT_DIR}/${ORDER_ID}_progress.txt"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [${step}] ${message}" >> "${progress_file}"
}

# ── Parse arguments ──────────────────────────────────────────────────────────
ORDER_ID=""
SAMPLE_NAMES=""
FASTQ_R1_FILES=""
FASTQ_R2_FILES=""
TARGET_BED_NAME=""
VARIANT_CALLER="bcftools"
FF_VARIANT_CALLER="bcftools"
EXTRA_NF_ARGS=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --order_id)          ORDER_ID="$2";          shift 2 ;;
        --sample_name)       SAMPLE_NAMES="$2";      shift 2 ;;
        --fastq_r1)          FASTQ_R1_FILES="$2";    shift 2 ;;
        --fastq_r2)          FASTQ_R2_FILES="$2";    shift 2 ;;
        --target_bed)        TARGET_BED_NAME="$2";   shift 2 ;;
        --variant_caller)    VARIANT_CALLER="$2";    shift 2 ;;
        --ff_variant_caller) FF_VARIANT_CALLER="$2"; shift 2 ;;
        --extra_args)        EXTRA_NF_ARGS="$2";     shift 2 ;;
        *)
            log_error "Unknown argument: $1"
            exit 1
            ;;
    esac
done

# ── Validate required arguments ──────────────────────────────────────────────
if [[ -z "$ORDER_ID" ]]; then
    log_error "--order_id is required"
    exit 1
fi
if [[ -z "$SAMPLE_NAMES" ]]; then
    log_error "--sample_name is required"
    exit 1
fi
if [[ -z "$FASTQ_R1_FILES" ]]; then
    log_error "--fastq_r1 is required"
    exit 1
fi

# ── Prepare output directories early (for progress tracking) ─────────────────
ORDER_OUTPUT_DIR="${OUTPUT_DIR}/${ORDER_ID}"
ORDER_LOG_DIR="${LOG_DIR}/${ORDER_ID}"
WORK_DIR="${ANALYSIS_DIR}/${ORDER_ID}_work"

mkdir -p "${ORDER_OUTPUT_DIR}" "${ORDER_LOG_DIR}" "${WORK_DIR}"

log_info "============================================================"
log_info "  Single Gene NIPT Pipeline v${PIPELINE_VERSION} - Starting"
log_info "  Panel: Twist UCL_SingleGeneNIPT_TE-96276661 (hg38)"
log_info "============================================================"
log_info "Order ID       : ${ORDER_ID}"
log_info "Sample Names   : ${SAMPLE_NAMES}"
log_info "Variant Caller : ${VARIANT_CALLER}"
log_info "FF Caller      : ${FF_VARIANT_CALLER}"
log_info "Timestamp      : ${TIMESTAMP}"
log_info "============================================================"

update_progress "INIT" "Pipeline starting for order ${ORDER_ID}"

# ── Detect reference files ───────────────────────────────────────────────────
log_info "Detecting reference files in ${DATA_DIR} ..."
update_progress "INIT" "Detecting reference files"

# Reference FASTA
REF_FASTA=$(find "${DATA_DIR}/reference" -name "*.fa" -o -name "*.fasta" 2>/dev/null | head -1)
if [[ -z "$REF_FASTA" ]]; then
    log_error "No reference FASTA found in ${DATA_DIR}/reference/"
    update_progress "ERROR" "No reference FASTA found"
    exit 1
fi
log_info "Reference FASTA: ${REF_FASTA}"

# Check .fai index
if [[ ! -f "${REF_FASTA}.fai" ]]; then
    log_warn "Reference .fai index not found. Creating with samtools..."
    samtools faidx "${REF_FASTA}"
fi

# BWA-MEM2 index directory
BWA_INDEX_DIR="${DATA_DIR}/bwa_index"
if [[ ! -d "$BWA_INDEX_DIR" ]]; then
    log_error "BWA-MEM2 index directory not found: ${BWA_INDEX_DIR}"
    update_progress "ERROR" "BWA-MEM2 index not found"
    exit 1
fi
log_info "BWA Index Dir  : ${BWA_INDEX_DIR}"

# ── Detect BED files (Twist panel-specific) ──────────────────────────────────
log_info "Detecting Twist panel BED files in ${DATA_DIR}/target_bed/ ..."
update_progress "INIT" "Detecting Twist panel BED files"

BED_DIR="${DATA_DIR}/target_bed"

# 1. Target BED (All_target_regions) — primary target for variant calling & QC
if [[ -n "$TARGET_BED_NAME" ]]; then
    TARGET_BED="${BED_DIR}/${TARGET_BED_NAME}"
else
    TARGET_BED=$(find "${BED_DIR}" -name "All_target_regions*.bed" 2>/dev/null | head -1)
    if [[ -z "$TARGET_BED" ]]; then
        TARGET_BED=$(find "${BED_DIR}" -name "*.bed" -not -name "Probes*" -not -name "Target_bases*" -not -name "Target_regions_with_zero*" -not -name "ff_snps*" 2>/dev/null | head -1)
    fi
fi
if [[ -z "$TARGET_BED" ]] || [[ ! -f "$TARGET_BED" ]]; then
    log_error "Target BED file not found in ${BED_DIR}/"
    update_progress "ERROR" "Target BED file not found"
    exit 1
fi
log_info "Target BED     : ${TARGET_BED}"

# 2. FF SNPs BED — common SNP positions for Fetal Fraction estimation
FF_SNPS_BED=$(find "${BED_DIR}" -name "ff_snps*.bed" 2>/dev/null | head -1)
FF_SNPS_ARG=""
if [[ -n "$FF_SNPS_BED" ]] && [[ -f "$FF_SNPS_BED" ]]; then
    FF_SNPS_ARG="--ff_snps_bed ${FF_SNPS_BED}"
    log_info "FF SNPs BED    : ${FF_SNPS_BED}"
else
    log_warn "FF SNPs BED not found. Using target_bed for FF estimation."
    log_warn "  → For better FF accuracy, create ff_snps_<name>.bed with common SNPs (MAF>5%)"
fi

# 3. Probes BED — probe coverage QC
PROBES_BED=$(find "${BED_DIR}" -name "Probes_merged*.bed" 2>/dev/null | head -1)
PROBES_ARG=""
if [[ -n "$PROBES_BED" ]] && [[ -f "$PROBES_BED" ]]; then
    PROBES_ARG="--probes_bed ${PROBES_BED}"
    log_info "Probes BED     : ${PROBES_BED}"
else
    log_warn "Probes BED not found. Probe coverage QC will be skipped."
fi

# 4. Zero-probe BED — flagging uncovered target regions
ZERO_PROBE_BED=$(find "${BED_DIR}" -name "Target_regions_with_zero_probes*.bed" 2>/dev/null | head -1)
ZERO_PROBE_ARG=""
if [[ -n "$ZERO_PROBE_BED" ]] && [[ -f "$ZERO_PROBE_BED" ]]; then
    ZERO_PROBE_ARG="--zero_probe_bed ${ZERO_PROBE_BED}"
    log_info "Zero-probe BED : ${ZERO_PROBE_BED}"
    ZP_COUNT=$(grep -c -v "^#\|^track\|^browser" "${ZERO_PROBE_BED}" 2>/dev/null || echo "0")
    log_info "  → ${ZP_COUNT} zero-probe regions will be flagged in QC/variant reports"
else
    log_warn "Zero-probe BED not found. Zero-probe region flagging will be skipped."
fi

# 5. Gene list — for gene coverage validation
GENE_LIST=$(find "${BED_DIR}" -name "gene_list*" -o -name "singleGene*" 2>/dev/null | head -1)
if [[ -z "$GENE_LIST" ]]; then
    GENE_LIST=$(find "${DATA_DIR}" -name "gene_list*" -o -name "singleGene*" 2>/dev/null | head -1)
fi
GENE_LIST_ARG=""
if [[ -n "$GENE_LIST" ]] && [[ -f "$GENE_LIST" ]]; then
    GENE_LIST_ARG="--gene_list ${GENE_LIST}"
    log_info "Gene list      : ${GENE_LIST}"
else
    log_warn "Gene list not found. Gene coverage validation will be skipped."
fi

update_progress "INIT" "BED files detected successfully"

# ── Generate samplesheet ────────────────────────────────────────────────────
SAMPLESHEET="${ANALYSIS_DIR}/${ORDER_ID}_samplesheet.csv"
log_info "Generating samplesheet: ${SAMPLESHEET}"
update_progress "INIT" "Generating samplesheet"

echo "sample_id,fastq_1,fastq_2" > "${SAMPLESHEET}"

IFS=',' read -ra NAMES <<< "${SAMPLE_NAMES}"
IFS=',' read -ra R1_FILES <<< "${FASTQ_R1_FILES}"

if [[ -n "$FASTQ_R2_FILES" ]]; then
    IFS=',' read -ra R2_FILES <<< "${FASTQ_R2_FILES}"
else
    R2_FILES=()
fi

for i in "${!NAMES[@]}"; do
    SAMPLE="${NAMES[$i]}"
    R1="${FASTQ_DIR}/${R1_FILES[$i]}"
    R2=""

    if [[ ${#R2_FILES[@]} -gt $i ]] && [[ -n "${R2_FILES[$i]}" ]]; then
        R2="${FASTQ_DIR}/${R2_FILES[$i]}"
    fi

    # Validate FASTQ files exist
    if [[ ! -f "$R1" ]]; then
        log_error "FASTQ R1 not found: ${R1}"
        update_progress "ERROR" "FASTQ R1 not found: ${R1}"
        exit 1
    fi
    if [[ -n "$R2" ]] && [[ ! -f "$R2" ]]; then
        log_error "FASTQ R2 not found: ${R2}"
        update_progress "ERROR" "FASTQ R2 not found: ${R2}"
        exit 1
    fi

    echo "${SAMPLE},${R1},${R2}" >> "${SAMPLESHEET}"
    log_info "  Sample: ${SAMPLE}  R1: ${R1}  R2: ${R2:-N/A}"
done

log_info "Samplesheet generated with ${#NAMES[@]} sample(s)"
update_progress "INIT" "Samplesheet generated with ${#NAMES[@]} sample(s)"

# ── Check for user-provided config overrides ─────────────────────────────────
NF_PARAMS_FILE=""
if [[ -f "${CONFIG_DIR}/params.yaml" ]]; then
    NF_PARAMS_FILE="-params-file ${CONFIG_DIR}/params.yaml"
    log_info "Using custom params file: ${CONFIG_DIR}/params.yaml"
fi

QC_THRESHOLDS_ARG=""
if [[ -f "${CONFIG_DIR}/qc_thresholds.json" ]]; then
    QC_THRESHOLDS_ARG="--fastq_qc_thresholds ${CONFIG_DIR}/qc_thresholds.json --bam_qc_thresholds ${CONFIG_DIR}/qc_thresholds.json"
    log_info "Using custom QC thresholds: ${CONFIG_DIR}/qc_thresholds.json"
fi

ANNOTATION_ARGS=""
if [[ -f "${DATA_DIR}/annotation/clinvar.vcf.gz" ]]; then
    ANNOTATION_ARGS="--run_annotation true --annotation_tool snpsift --clinvar_vcf ${DATA_DIR}/annotation/clinvar.vcf.gz"
    log_info "ClinVar annotation enabled"
fi

# ── Run Nextflow pipeline ────────────────────────────────────────────────────
log_info "============================================================"
log_info "  Launching Nextflow Pipeline"
log_info "============================================================"

update_progress "PIPELINE" "Launching Nextflow pipeline"

NF_LOG="${ORDER_LOG_DIR}/nextflow_${TIMESTAMP}.log"

cd "${PIPELINE_DIR}"

# Build the Nextflow command
NF_CMD="nextflow run main.nf \
    -profile docker_internal \
    ${NF_PARAMS_FILE} \
    --input ${SAMPLESHEET} \
    --target_bed ${TARGET_BED} \
    ${FF_SNPS_ARG} \
    ${PROBES_ARG} \
    ${ZERO_PROBE_ARG} \
    ${GENE_LIST_ARG} \
    --reference_fasta ${REF_FASTA} \
    --bwa_index_dir ${BWA_INDEX_DIR} \
    --outdir ${ORDER_OUTPUT_DIR} \
    --variant_caller ${VARIANT_CALLER} \
    --ff_variant_caller ${FF_VARIANT_CALLER} \
    ${QC_THRESHOLDS_ARG} \
    ${ANNOTATION_ARGS} \
    ${EXTRA_NF_ARGS} \
    -work-dir ${WORK_DIR} \
    -with-trace ${ORDER_LOG_DIR}/trace.txt \
    -with-timeline ${ORDER_LOG_DIR}/timeline.html \
    -with-report ${ORDER_LOG_DIR}/report.html \
    -with-dag ${ORDER_LOG_DIR}/dag.html \
    -ansi-log false"

log_info "Command: ${NF_CMD}"
log_info "Log file: ${NF_LOG}"

# Execute
eval ${NF_CMD} 2>&1 | tee "${NF_LOG}"
NF_EXIT_CODE=${PIPESTATUS[0]}

# ── Post-execution ───────────────────────────────────────────────────────────
if [[ $NF_EXIT_CODE -eq 0 ]]; then
    log_info "============================================================"
    log_info "  Nextflow Pipeline completed SUCCESSFULLY"
    log_info "============================================================"
    update_progress "PIPELINE" "Nextflow pipeline completed successfully"

    # ── Generate daemon-compatible result JSON ───────────────────────────
    # This is the file that nipt-daemon monitors: {ORDER_ID}.json
    update_progress "RESULTS" "Generating result JSON"

    RESULT_JSON="${ORDER_OUTPUT_DIR}/${ORDER_ID}.json"
    python3 "${PIPELINE_DIR}/scripts/generate_result_json.py" \
        --order_id "${ORDER_ID}" \
        --output_dir "${ORDER_OUTPUT_DIR}" \
        --sample_names "${SAMPLE_NAMES}" \
        --pipeline_version "${PIPELINE_VERSION}" \
        --output "${RESULT_JSON}" \
        2>&1 | tee -a "${ORDER_LOG_DIR}/result_json.log"

    if [[ ! -f "${RESULT_JSON}" ]]; then
        log_error "Failed to generate result JSON"
        update_progress "ERROR" "Failed to generate result JSON"
        NF_EXIT_CODE=1
    else
        log_info "Result JSON generated: ${RESULT_JSON}"
        update_progress "RESULTS" "Result JSON generated"
    fi

    # ── Generate output TAR archive ──────────────────────────────────────
    # This is the second file daemon monitors: {ORDER_ID}.output.tar
    update_progress "RESULTS" "Creating output archive"

    TAR_FILE="${ORDER_OUTPUT_DIR}/${ORDER_ID}.output.tar"
    tar -cf "${TAR_FILE}" \
        -C "${ORDER_OUTPUT_DIR}" \
        --exclude="*.output.tar" \
        --exclude=".completion_status.json" \
        . 2>&1 | tee -a "${ORDER_LOG_DIR}/tar.log"

    if [[ ! -f "${TAR_FILE}" ]]; then
        log_error "Failed to create output TAR"
        update_progress "ERROR" "Failed to create output TAR"
        NF_EXIT_CODE=1
    else
        TAR_SIZE=$(du -h "${TAR_FILE}" | cut -f1)
        log_info "Output TAR created: ${TAR_FILE} (${TAR_SIZE})"
        update_progress "RESULTS" "Output archive created (${TAR_SIZE})"
    fi

    # ── Write completion marker ──────────────────────────────────────────
    if [[ $NF_EXIT_CODE -eq 0 ]]; then
        echo "{
    \"order_id\": \"${ORDER_ID}\",
    \"status\": \"COMPLETED\",
    \"service_code\": \"sgNIPT\",
    \"pipeline\": \"single_gene_nipt\",
    \"pipeline_version\": \"${PIPELINE_VERSION}\",
    \"panel\": \"Twist_UCL_SingleGeneNIPT_TE-96276661_hg38\",
    \"timestamp\": \"$(date -Iseconds)\",
    \"output_dir\": \"${ORDER_OUTPUT_DIR}\",
    \"samples\": \"${SAMPLE_NAMES}\",
    \"result_json\": \"${RESULT_JSON}\",
    \"output_tar\": \"${TAR_FILE}\",
    \"exit_code\": 0
}" > "${ORDER_OUTPUT_DIR}/.completion_status.json"

        log_info "Pipeline completed SUCCESSFULLY"
        update_progress "COMPLETE" "Pipeline completed successfully"
    fi

else
    log_error "============================================================"
    log_error "  Pipeline FAILED (exit code: ${NF_EXIT_CODE})"
    log_error "============================================================"
    log_error "Check logs: ${NF_LOG}"

    # Generate failure marker for daemon polling
    echo "{
    \"order_id\": \"${ORDER_ID}\",
    \"status\": \"FAILED\",
    \"service_code\": \"sgNIPT\",
    \"pipeline\": \"single_gene_nipt\",
    \"pipeline_version\": \"${PIPELINE_VERSION}\",
    \"panel\": \"Twist_UCL_SingleGeneNIPT_TE-96276661_hg38\",
    \"timestamp\": \"$(date -Iseconds)\",
    \"output_dir\": \"${ORDER_OUTPUT_DIR}\",
    \"samples\": \"${SAMPLE_NAMES}\",
    \"exit_code\": ${NF_EXIT_CODE}
}" > "${ORDER_OUTPUT_DIR}/.completion_status.json"

    update_progress "ERROR" "Pipeline FAILED (exit code: ${NF_EXIT_CODE})"
fi

# ── Cleanup Nextflow work directory (optional) ───────────────────────────────
if [[ $NF_EXIT_CODE -eq 0 ]] && [[ "${CLEANUP_WORK_DIR:-false}" == "true" ]]; then
    log_info "Cleaning up Nextflow work directory: ${WORK_DIR}"
    rm -rf "${WORK_DIR}"
    update_progress "CLEANUP" "Work directory cleaned up"
fi

log_info "Entrypoint finished with exit code: ${NF_EXIT_CODE}"
exit ${NF_EXIT_CODE}
