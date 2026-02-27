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
#  Usage:
#    ./run_sgnipt.sh \
#        --order_id ORDER_001 \
#        --sample_name "SAMPLE_A,SAMPLE_B" \
#        --fastq_r1 "SAMPLE_A_R1.fastq.gz,SAMPLE_B_R1.fastq.gz" \
#        --fastq_r2 "SAMPLE_A_R2.fastq.gz,SAMPLE_B_R2.fastq.gz" \
#        [--target_bed panel_v1.bed] \
#        [--variant_caller bcftools] \
#        [--detached]
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

while [[ $# -gt 0 ]]; do
    case "$1" in
        --order_id)          ORDER_ID="$2";          shift 2 ;;
        --sample_name)       SAMPLE_NAMES="$2";      shift 2 ;;
        --fastq_r1)          FASTQ_R1="$2";          shift 2 ;;
        --fastq_r2)          FASTQ_R2="$2";          shift 2 ;;
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

# ── Validate required arguments ──────────────────────────────────────────────
if [[ -z "$ORDER_ID" ]]; then
    echo "[ERROR] --order_id is required"
    exit 1
fi
if [[ -z "$SAMPLE_NAMES" ]]; then
    echo "[ERROR] --sample_name is required"
    exit 1
fi
if [[ -z "$FASTQ_R1" ]]; then
    echo "[ERROR] --fastq_r1 is required"
    exit 1
fi

# ── work_dir: YYMM pattern (compatible with daemon's extract_work_dir) ───────
#  The daemon uses extract_work_dir(order_id) to derive YYMM from order_id.
#  Here we use the current date as fallback. The daemon can override this
#  by setting SGNIPT_WORK_DIR environment variable before calling this script.
WORK_DIR="${SGNIPT_WORK_DIR:-$(date '+%y%m')}"

# ── Ensure host directories exist ────────────────────────────────────────────
for DIR in "$HOST_FASTQ_DIR" "$HOST_DATA_DIR" "$HOST_CONFIG_DIR" \
           "$HOST_ANALYSIS_DIR" "$HOST_OUTPUT_DIR" "$HOST_LOG_DIR"; do
    if [[ ! -d "$DIR" ]]; then
        echo "[WARN] Creating directory: ${DIR}"
        mkdir -p "$DIR"
    fi
done

# Create work_dir-based output directory (daemon monitors this path)
# Pattern: {ROOT}/output/{YYMM}/{ORDER_ID}/{ORDER_ID}.json
ORDER_OUTPUT_HOST="${HOST_OUTPUT_DIR}/${WORK_DIR}/${ORDER_ID}"
ORDER_LOG_HOST="${HOST_LOG_DIR}/${WORK_DIR}/${ORDER_ID}"
mkdir -p "${ORDER_OUTPUT_HOST}" "${ORDER_LOG_HOST}"

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
echo "  Container Name : ${CONTAINER_NAME}"
echo "  Docker Image   : ${DOCKER_IMAGE}"
echo "  CPU / Memory   : ${CONTAINER_CPUS} / ${CONTAINER_MEMORY}"
echo "  Work Dir       : ${WORK_DIR}"
echo "  Root Dir       : ${SGNIPT_ROOT_DIR}"
echo "  Detached       : ${DETACHED}"
echo "============================================================"

#  Volume mounts follow the same pattern as sWGS NIPT:
#    HOST_DIR → /Work/SgNIPT/{subdir}
#
#  The container's entrypoint.sh uses these internal paths.
#  Output goes to /Work/SgNIPT/output/{ORDER_ID}/ which maps to:
#    {ROOT}/output/{YYMM}/{ORDER_ID}/ on the host
#
#  NOTE: We mount the YYMM-level output dir so the container writes
#  directly to the daemon-monitored path.

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
    -v "${HOST_ANALYSIS_DIR}:/Work/SgNIPT/analysis" \
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
    echo "[INFO] Daemon polls for : ${ORDER_OUTPUT_HOST}/${ORDER_ID}.json"
    echo "[INFO]                    ${ORDER_OUTPUT_HOST}/${ORDER_ID}.output.tar"
    echo "[INFO] Progress file    : ${ORDER_OUTPUT_HOST}/${ORDER_ID}_progress.txt"
else
    echo "[ERROR] Failed to start container '${CONTAINER_NAME}' (exit: ${DOCKER_EXIT})"
    exit $DOCKER_EXIT
fi
