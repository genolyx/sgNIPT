#!/bin/bash
###############################################################################
#  Single Gene NIPT - Docker Run Wrapper Script
###############################################################################
#  This script is called by nipt-daemon (or manually) to launch a single
#  Order analysis inside a Docker container.
#
#  Designed to be compatible with the nipt-daemon runner.py pattern:
#    - Container name = ORDER_ID (daemon uses docker ps -f name=ORDER_ID)
#    - Output path follows YYMM work_dir pattern
#    - daemon polls for {ORDER_ID}.json + {ORDER_ID}.output.tar
#
#  Usage (minimal — FASTQ under fastq/{work_dir}/{order_id}/):
#    export SGNIPT_ROOT_DIR=/path/to/sgNIPT
#    ./run_sgnipt.sh --order_id NA12878_Twist_Exome --work_dir 2603
#
#  Explicit FASTQ paths (optional if auto-discovery works):
#    ./run_sgnipt.sh \
#        --order_id NA12878_Twist_Exome \
#        [--sample_name "NA12878_Twist_Exome"]   # default: same as order_id \
#        [--fastq_r1 "2603/NA12878_Twist_Exome/file_R1.fastq.gz"] \
#        [--fastq_r2 "2603/NA12878_Twist_Exome/file_R2.fastq.gz"] \
#        [--work_dir 2603] \
#        [--ref_dir ./data/refs]                # logged only; data mount unchanged \
#        [--target_bed panel_v1.bed] \
#        [--variant_caller bcftools] \
#        [--detached]
#
#  Directory structure: {type}/{YYMM}/{SAMPLE}/
#    fastq/2603/NA12878_Twist_Exome/   - input FASTQ
#    analysis/2603/NA12878_Twist_Exome/ - Nextflow work, intermediates
#    output/2603/NA12878_Twist_Exome/  - JSON, PNG for service-daemon
#    log/2603/NA12878_Twist_Exome/     - analysis logs
#
#  Environment Variables (set by daemon or .env file):
#    SGNIPT_ROOT_DIR     - Root directory for SgNIPT on host (required)
#    SGNIPT_DOCKER_IMAGE - Docker image name (default: sgnipt:latest)
#    SGNIPT_CPUS         - CPU limit for container (default: 16)
#    SGNIPT_MEMORY       - Memory limit for container (default: 32g)
#    SGNIPT_CLEANUP_WORK - Clean up Nextflow work dir after success (default: false)
###############################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "${SCRIPT_DIR}")"

# ── Load .env if present ─────────────────────────────────────────────────────
if [[ -f "${PROJECT_DIR}/.env" ]]; then
    # shellcheck disable=SC1091
    source "${PROJECT_DIR}/.env"
fi

# ── Configuration with defaults ──────────────────────────────────────────────
SGNIPT_ROOT_DIR="${SGNIPT_ROOT_DIR:?'SGNIPT_ROOT_DIR must be set'}"
DOCKER_IMAGE="${SGNIPT_DOCKER_IMAGE:-sgnipt:latest}"
CONTAINER_CPUS="${SGNIPT_CPUS:-16}"
CONTAINER_MEMORY="${SGNIPT_MEMORY:-32g}"
CLEANUP_WORK="${SGNIPT_CLEANUP_WORK:-false}"

# ── Host directory paths (derived from ROOT_DIR) ─────────────────────────────
HOST_FASTQ_DIR="${SGNIPT_ROOT_DIR}/fastq"
HOST_DATA_DIR="${SGNIPT_ROOT_DIR}/data"
HOST_CONFIG_DIR="${SGNIPT_ROOT_DIR}/config"
HOST_ANALYSIS_DIR="${SGNIPT_ROOT_DIR}/analysis"
HOST_OUTPUT_DIR="${SGNIPT_ROOT_DIR}/output"
HOST_LOG_DIR="${SGNIPT_ROOT_DIR}/log"

# ── Parse arguments ──────────────────────────────────────────────────────────
ORDER_ID=""
SAMPLE_NAMES=""
FASTQ_R1=""
FASTQ_R2=""
TARGET_BED=""
VARIANT_CALLER="bcftools"
FF_VARIANT_CALLER="bcftools"
EXTRA_ARGS=""
DETACHED=false
WORK_DIR_ARG=""
REF_DIR=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --order_id)          ORDER_ID="$2";          shift 2 ;;
        --sample_name)       SAMPLE_NAMES="$2";      shift 2 ;;
        --fastq_r1)          FASTQ_R1="$2";          shift 2 ;;
        --fastq_r2)          FASTQ_R2="$2";          shift 2 ;;
        --work_dir|-w)       WORK_DIR_ARG="$2";      shift 2 ;;
        --ref_dir|-r)        REF_DIR="$2";           shift 2 ;;
        --target_bed)        TARGET_BED="$2";        shift 2 ;;
        --variant_caller)    VARIANT_CALLER="$2";    shift 2 ;;
        --ff_variant_caller) FF_VARIANT_CALLER="$2"; shift 2 ;;
        --extra_args)        EXTRA_ARGS="$2";        shift 2 ;;
        --detached)          DETACHED=true;           shift   ;;
        *)
            echo "[ERROR] Unknown argument: $1"
            exit 1
            ;;
    esac
done

# ── work_dir: YYMM pattern ───────────────────────────────────────────────────
#  Structure: fastq|analysis|output|log / {YYMM} / {SAMPLE_NAME}/
#  Priority: --work_dir > SGNIPT_WORK_DIR env > current date
WORK_DIR="${WORK_DIR_ARG:-${SGNIPT_WORK_DIR:-$(date '+%y%m')}}"

# ── Defaults (carrier-style: one sample id drives order + sample name) ─────
if [[ -z "$SAMPLE_NAMES" ]]; then
    SAMPLE_NAMES="${ORDER_ID}"
fi

# ── Validate required arguments ──────────────────────────────────────────────
if [[ -z "$ORDER_ID" ]]; then
    echo "[ERROR] --order_id is required"
    exit 1
fi

# ── Auto-discover FASTQ under fastq/{work_dir}/{order_id}/ ─────────────────
#  Expects e.g. * _R1.fastq.gz / _R2.fastq.gz (Illumina-style)
if [[ -z "$FASTQ_R1" ]]; then
    FASTQ_HOST_DIR="${HOST_FASTQ_DIR}/${WORK_DIR}/${ORDER_ID}"
    if [[ ! -d "$FASTQ_HOST_DIR" ]]; then
        echo "[ERROR] FASTQ directory not found: ${FASTQ_HOST_DIR}"
        echo "[HINT]  Pass --work_dir and --order_id so that FASTQs live under fastq/{work_dir}/{order_id}/"
        echo "       or set --fastq_r1 / --fastq_r2 explicitly."
        exit 1
    fi
    mapfile -t _R1_CAND < <(find "$FASTQ_HOST_DIR" -maxdepth 1 -type f \( -name '*_R1.fastq.gz' -o -name '*_R1.fq.gz' \) | sort)
    if [[ ${#_R1_CAND[@]} -eq 0 ]]; then
        echo "[ERROR] No R1 FASTQ (*_R1.fastq.gz) found under: ${FASTQ_HOST_DIR}"
        exit 1
    fi
    if [[ ${#_R1_CAND[@]} -gt 1 ]]; then
        echo "[ERROR] Multiple R1 FASTQ files found; specify --fastq_r1 and --fastq_r2 explicitly:"
        printf '  %s\n' "${_R1_CAND[@]}"
        exit 1
    fi
    _R1_ABS="${_R1_CAND[0]}"
    _BASE="$(basename "$_R1_ABS")"
    _R2_ABS="${FASTQ_HOST_DIR}/${_BASE/_R1.fastq.gz/_R2.fastq.gz}"
    [[ -f "$_R2_ABS" ]] || _R2_ABS="${FASTQ_HOST_DIR}/${_BASE/_R1.fq.gz/_R2.fq.gz}"
    if [[ ! -f "$_R2_ABS" ]]; then
        echo "[ERROR] Paired R2 not found for: ${_R1_ABS}"
        echo "[ERROR] Expected: ${_R2_ABS}"
        exit 1
    fi
    # Paths relative to HOST_FASTQ_DIR (container: /Work/SgNIPT/fastq/...)
    FASTQ_R1="${WORK_DIR}/${ORDER_ID}/${_BASE}"
    FASTQ_R2="${WORK_DIR}/${ORDER_ID}/$(basename "$_R2_ABS")"
    echo "[INFO] Auto-selected FASTQ pair:"
    echo "       R1: ${FASTQ_R1}"
    echo "       R2: ${FASTQ_R2}"
fi

if [[ -z "$FASTQ_R1" ]]; then
    echo "[ERROR] --fastq_r1 is required (or use auto-discovery with fastq/${WORK_DIR}/${ORDER_ID}/)"
    exit 1
fi

# ── Ensure host directories exist ────────────────────────────────────────────
for DIR in "$HOST_FASTQ_DIR" "$HOST_DATA_DIR" "$HOST_CONFIG_DIR" \
           "$HOST_ANALYSIS_DIR" "$HOST_OUTPUT_DIR" "$HOST_LOG_DIR"; do
    if [[ ! -d "$DIR" ]]; then
        echo "[WARN] Creating directory: ${DIR}"
        mkdir -p "$DIR"
    fi
done

# Per-sample paths: {type}/{YYMM}/{ORDER_ID}/
#  - analysis: Nextflow work dir, intermediates, samplesheet
#  - output:   Final JSON, PNG for service-daemon
#  - log:      Analysis logs
ORDER_ANALYSIS_HOST="${HOST_ANALYSIS_DIR}/${WORK_DIR}/${ORDER_ID}"
ORDER_OUTPUT_HOST="${HOST_OUTPUT_DIR}/${WORK_DIR}/${ORDER_ID}"
ORDER_LOG_HOST="${HOST_LOG_DIR}/${WORK_DIR}/${ORDER_ID}"
mkdir -p "${ORDER_ANALYSIS_HOST}" "${ORDER_OUTPUT_HOST}" "${ORDER_LOG_HOST}"

# ── Container name = ORDER_ID (daemon uses docker ps -f name=ORDER_ID) ──────
CONTAINER_NAME="${ORDER_ID}"

# Remove existing container with same name (if any)
if docker ps -a --format '{{.Names}}' | grep -q "^${CONTAINER_NAME}$"; then
    echo "[WARN] Removing existing container: ${CONTAINER_NAME}"
    docker rm -f "${CONTAINER_NAME}" 2>/dev/null || true
fi

# ── Build entrypoint arguments ───────────────────────────────────────────────
ENTRYPOINT_ARGS=(
    --order_id "$ORDER_ID"
    --sample_name "$SAMPLE_NAMES"
    --fastq_r1 "$FASTQ_R1"
    --variant_caller "$VARIANT_CALLER"
    --ff_variant_caller "$FF_VARIANT_CALLER"
)

if [[ -n "$FASTQ_R2" ]]; then
    ENTRYPOINT_ARGS+=(--fastq_r2 "$FASTQ_R2")
fi

if [[ -n "$TARGET_BED" ]]; then
    ENTRYPOINT_ARGS+=(--target_bed "$TARGET_BED")
fi

if [[ -n "$EXTRA_ARGS" ]]; then
    ENTRYPOINT_ARGS+=(--extra_args "$EXTRA_ARGS")
fi

# ── Launch Docker container ──────────────────────────────────────────────────
echo "============================================================"
echo "  Single Gene NIPT - Launching Docker Container"
echo "============================================================"
echo "  Order ID       : ${ORDER_ID}"
echo "  Sample name    : ${SAMPLE_NAMES}"
if [[ -n "${REF_DIR}" ]]; then
    REF_DISPLAY="$(realpath "${REF_DIR}" 2>/dev/null || echo "${REF_DIR}")"
    echo "  Reference dir  : ${REF_DISPLAY}  (informational; pipeline uses ${HOST_DATA_DIR} → /Work/SgNIPT/data)"
fi
echo "  Container Name : ${CONTAINER_NAME}"
echo "  Work Dir       : ${WORK_DIR}"
echo "  Analysis       : ${ORDER_ANALYSIS_HOST}"
echo "  Output         : ${ORDER_OUTPUT_HOST}"
echo "  Log            : ${ORDER_LOG_HOST}"
echo "  Docker Image   : ${DOCKER_IMAGE}"
echo "  CPU / Memory   : ${CONTAINER_CPUS} / ${CONTAINER_MEMORY}"
echo "============================================================"

#  Volume mounts: {type}/{YYMM}/{ORDER_ID}/ structure
#    fastq/     → whole tree (fastq/2603/NA12878/...)
#    analysis/  → analysis/2603/ORDER_ID (Nextflow work, intermediates)
#    output/    → output/2603/ORDER_ID (JSON, PNG for service-daemon)
#    log/       → log/2603/ORDER_ID (analysis logs)

docker run \
    --user "$(id -u):$(id -g)" \
    --name "${CONTAINER_NAME}" \
    -d \
    -e TZ=Asia/Seoul \
    -e USER="$(whoami)" \
    -e USERNAME="$(whoami)" \
    -e HOME=/tmp \
    -e FONTCONFIG_PATH=/tmp \
    -e CLEANUP_WORK_DIR="${CLEANUP_WORK}" \
    -e NXF_HOME="/tmp/.nextflow" \
    -e NXF_TEMP="/Work/SgNIPT/analysis/tmp" \
    -e WORK_DIR="${WORK_DIR}" \
    -e BWA2="bwa-mem2" \
    -e SAMTOOLS="samtools" \
    -e BCFTOOLS="bcftools" \
    -e PYTHON3="python3" \
    -v "${HOST_FASTQ_DIR}:/Work/SgNIPT/fastq" \
    -v "${HOST_DATA_DIR}:/Work/SgNIPT/data" \
    -v "${HOST_CONFIG_DIR}:/Work/SgNIPT/config" \
    -v "${ORDER_ANALYSIS_HOST}:/Work/SgNIPT/analysis" \
    -v "${ORDER_OUTPUT_HOST}:/Work/SgNIPT/output/${ORDER_ID}" \
    -v "${ORDER_LOG_HOST}:/Work/SgNIPT/log/${ORDER_ID}" \
    --cpus "${CONTAINER_CPUS}" \
    --memory "${CONTAINER_MEMORY}" \
    "${DOCKER_IMAGE}" \
    "${ENTRYPOINT_ARGS[@]}"

DOCKER_EXIT=$?

if [[ $DOCKER_EXIT -eq 0 ]]; then
    echo "[INFO] Container '${CONTAINER_NAME}' started successfully (detached)"
    echo "[INFO] Monitor logs     : docker logs -f ${CONTAINER_NAME}"
    echo "[INFO] Analysis dir    : ${ORDER_ANALYSIS_HOST}"
    echo "[INFO] Output dir      : ${ORDER_OUTPUT_HOST}"
    echo "[INFO] Log dir         : ${ORDER_LOG_HOST}"
    echo "[INFO] Service-daemon  : ${ORDER_OUTPUT_HOST}/${ORDER_ID}.json"
    echo "[INFO]                    ${ORDER_OUTPUT_HOST}/${ORDER_ID}.output.tar"
else
    echo "[ERROR] Failed to start container '${CONTAINER_NAME}' (exit: ${DOCKER_EXIT})"
    exit $DOCKER_EXIT
fi
