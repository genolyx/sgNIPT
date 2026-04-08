#!/bin/bash
###############################################################################
#  Single Gene NIPT - Docker Run Wrapper Script
###############################################################################
#  This script is called by nipt-daemon (or manually) to launch a single
#  Order analysis inside a Docker container.
#
#  service-daemon / carrier_screening parity (default):
#    - Foreground ``docker run`` (no -d): the shell blocks until the pipeline exits, like
#      carrier_screening/src/run_analysis.sh — subprocess exit code == pipeline exit code.
#    - docker ps NAMES = --order_id (docker run --name <order_id>; label sgnipt.order_id)
#    - Output path follows YYMM work_dir pattern; results: {ORDER_ID}.json, .output.tar
#  Legacy: pass --detached for ``docker run -d`` + background permission repair (fire-and-forget).
#
#  Usage (minimal — FASTQ under fastq/{work_dir}/{order_id}/):
#    ./src/run_sgnipt.sh --order-id NA12878_Twist_Exome --work-id 2603
#    (SGNIPT_ROOT_DIR defaults to the repo root: parent of src/)
#    Override: export SGNIPT_ROOT_DIR=/other/root   # e.g. Portal/daemon layout
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
#        [--fresh]                            # remove analysis/<YYMM>/<id>/work + .nextflow, no Nextflow -resume
#        [--detached]                         # docker -d (only if you need non-blocking start)
#
#  Directory structure: {type}/{YYMM}/{SAMPLE}/
#    fastq/2603/NA12878_Twist_Exome/   - input FASTQ
#    analysis/2603/NA12878_Twist_Exome/ - Nextflow work, intermediates
#    output/2603/NA12878_Twist_Exome/  - JSON, PNG for service-daemon
#    log/2603/NA12878_Twist_Exome/     - analysis logs
#
#  Environment Variables (set by daemon, .env, or shell):
#    SGNIPT_ROOT_DIR     - Root directory for SgNIPT on host (default: repo root = dirname of src/)
#    SGNIPT_DOCKER_IMAGE - Docker image name (default: sgnipt:latest)
#    SGNIPT_CPUS         - CPU limit for container (default: 16)
#    SGNIPT_MEMORY       - Memory limit for container (default: 32g)
#    SGNIPT_CLEANUP_WORK - Clean up Nextflow work dir after success (default: false)
#    SGNIPT_NO_RESUME    - If 1/true, Nextflow runs without -resume (default: 0). Prefer --fresh to also delete work/.
#    CHOWN_SPEC          - uid:gid for output ownership (default: $(id -u):$(id -g))
#                          After each run, analysis/output/log/<WORK_DIR> gets
#                          chown+g+rwX+setgid so the group can add samples without
#                          755-only-owner permission issues.
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
# Default: clone/layout root (analysis|output|log|fastq|data live under here).
# Portal/daemon may set SGNIPT_ROOT_DIR=/data/sgnipt_work — that overrides this.
SGNIPT_ROOT_DIR="${SGNIPT_ROOT_DIR:-${PROJECT_DIR}}"
DOCKER_IMAGE="${SGNIPT_DOCKER_IMAGE:-sgnipt:latest}"
CONTAINER_CPUS="${SGNIPT_CPUS:-16}"
CONTAINER_MEMORY="${SGNIPT_MEMORY:-32g}"
CLEANUP_WORK="${SGNIPT_CLEANUP_WORK:-false}"

# 1-A: Output ownership spec (uid:gid). Override via env to match group write requirements.
CHOWN_SPEC="${CHOWN_SPEC:-$(id -u):$(id -g)}"

# ── Host directory paths (derived from ROOT_DIR) ─────────────────────────────
HOST_FASTQ_DIR="${SGNIPT_FASTQ_DIR:-${SGNIPT_ROOT_DIR}/fastq}"
# data/config: allow override so daemon (running in a container) can pass the
# real HOST-visible path for docker -v, which must be a HOST path not a
# container-internal bind-mount path.  Set SGNIPT_DATA_DIR / SGNIPT_CONFIG_DIR
# to the layout host path (e.g. /home/ken/sgNIPT/data) in .env.docker.
HOST_DATA_DIR="${SGNIPT_DATA_DIR:-${SGNIPT_ROOT_DIR}/data}"
HOST_CONFIG_DIR="${SGNIPT_CONFIG_DIR:-${SGNIPT_ROOT_DIR}/config}"
HOST_ANALYSIS_DIR="${SGNIPT_ROOT_DIR}/analysis"
HOST_OUTPUT_DIR="${SGNIPT_ROOT_DIR}/output"
HOST_LOG_DIR="${SGNIPT_ROOT_DIR}/log"

# 1-B: Repair permissions on analysis|output|log/<WORK_DIR> trees.
# Tries: docker alpine → passwordless sudo → chmod only.
# $1 = "quiet" suppresses progress messages (used at startup).
repair_order_tree_permissions() {
    local quiet="${1:-}" t any=0
    for base in "$HOST_ANALYSIS_DIR" "$HOST_OUTPUT_DIR" "$HOST_LOG_DIR"; do
        t="${base}/${WORK_DIR}"
        [[ -d "$t" ]] || continue
        any=1
    done
    [[ "$any" -eq 1 ]] || return 0

    if command -v docker >/dev/null 2>&1 && id -nG | grep -qw docker && docker info &>/dev/null; then
        if docker run --rm \
            -e CHOWN_SPEC="$CHOWN_SPEC" \
            -e WORK_DIR="$WORK_DIR" \
            -v "${HOST_ANALYSIS_DIR}:/fa" \
            -v "${HOST_OUTPUT_DIR}:/fo" \
            -v "${HOST_LOG_DIR}:/fl" \
            alpine:3.20 \
            sh -c 'for base in /fa /fo /fl; do
                t="$base/$WORK_DIR"
                [ -d "$t" ] || continue
                chown -R "$CHOWN_SPEC" "$t"
                chmod -R g+rwX "$t"
                find "$t" -type d -exec chmod g+s {} + 2>/dev/null || true
            done'; then
            [[ "$quiet" == quiet ]] || echo "  Permissions fixed (${WORK_DIR}): ${CHOWN_SPEC}, g+rwX, setgid (docker)"
            return 0
        fi
    fi
    if sudo -n true 2>/dev/null; then
        for base in "$HOST_ANALYSIS_DIR" "$HOST_OUTPUT_DIR" "$HOST_LOG_DIR"; do
            t="${base}/${WORK_DIR}"
            [[ -d "$t" ]] || continue
            sudo chown -R "${CHOWN_SPEC}" "$t" 2>/dev/null || true
            sudo chmod -R g+rwX "$t" 2>/dev/null || true
            sudo find "$t" -type d -exec chmod g+s {} + 2>/dev/null || true
        done
        [[ "$quiet" == quiet ]] || echo "  Permissions fixed (${WORK_DIR}): ${CHOWN_SPEC}, g+rwX, setgid (sudo)"
        return 0
    fi
    for base in "$HOST_ANALYSIS_DIR" "$HOST_OUTPUT_DIR" "$HOST_LOG_DIR"; do
        t="${base}/${WORK_DIR}"
        [[ -d "$t" ]] || continue
        chmod -R g+rwX "$t" 2>/dev/null || true
        find "$t" -type d -exec chmod g+s {} + 2>/dev/null || true
    done
    [[ "$quiet" == quiet ]] || echo "[WARN] Permissions: chmod only (chown requires docker group or passwordless sudo)"
}

# 1-C: Create per-sample directories robustly.
# Falls back to docker alpine if host mkdir fails (parent dir may be root-owned).
ensure_sample_output_dirs() {
    local a o l
    a="${ORDER_ANALYSIS_HOST}"
    o="${ORDER_OUTPUT_HOST}"
    l="${ORDER_LOG_HOST}"
    if mkdir -p "${a}" "$o" "$l" 2>/dev/null; then
        repair_order_tree_permissions quiet
        return 0
    fi
    if ! command -v docker >/dev/null 2>&1 || ! docker info &>/dev/null; then
        echo "[ERROR] Cannot create sample directories (permission denied) and docker unavailable."
        echo "  Fix with: sudo chown -R \"$(id -u):$(id -g)\" \"${HOST_ANALYSIS_DIR}/${WORK_DIR}\" \"${HOST_OUTPUT_DIR}/${WORK_DIR}\" \"${HOST_LOG_DIR}/${WORK_DIR}\""
        exit 1
    fi
    echo "[WARN] Creating output dirs via docker (host mkdir failed — parent may be root-owned)..."
    docker run --rm \
        -e WORK_DIR="$WORK_DIR" \
        -e ORDER_ID="$ORDER_ID" \
        -e CHOWN_SPEC="$CHOWN_SPEC" \
        -v "${HOST_ANALYSIS_DIR}:/fa" \
        -v "${HOST_OUTPUT_DIR}:/fo" \
        -v "${HOST_LOG_DIR}:/fl" \
        alpine:3.20 \
        sh -c 'mkdir -p "/fa/${WORK_DIR}/${ORDER_ID}" "/fo/${WORK_DIR}/${ORDER_ID}" "/fl/${WORK_DIR}/${ORDER_ID}" && \
            for base in /fa /fo /fl; do \
                t="$base/${WORK_DIR}"; \
                [ -d "$t" ] || continue; \
                chown -R "$CHOWN_SPEC" "$t"; \
                chmod -R g+rwX "$t"; \
                find "$t" -type d -exec chmod g+s {} + 2>/dev/null || true; \
            done'
}

# ── Parse arguments ──────────────────────────────────────────────────────────
require_arg() {
    [ $# -ge 2 ] && [ -n "$2" ] || { echo "[ERROR] $1 requires a value"; exit 1; }
}

ORDER_ID=""
SAMPLE_NAMES=""
FASTQ_R1=""
FASTQ_R2=""
TARGET_BED=""
EXTRA_ARGS=""
DETACHED=false
WORK_DIR_ARG=""
REF_DIR=""
FRESH=""
INPUT_BAM_HOST=""   # host path to bam_samplesheet.csv  (BAM mode; skips FASTQ steps)

while [[ $# -gt 0 ]]; do
    case "$1" in
        --order_id|--order-id)   require_arg "$1" "${2:-}"; ORDER_ID="$2";         shift 2 ;;
        --sample_name)           require_arg "$1" "${2:-}"; SAMPLE_NAMES="$2";     shift 2 ;;
        --fastq_r1)              require_arg "$1" "${2:-}"; FASTQ_R1="$2";         shift 2 ;;
        --fastq_r2)              require_arg "$1" "${2:-}"; FASTQ_R2="$2";         shift 2 ;;
        --work_dir|--work-id|-w) require_arg "$1" "${2:-}"; WORK_DIR_ARG="$2";     shift 2 ;;
        --ref_dir|-r)            require_arg "$1" "${2:-}"; REF_DIR="$2";          shift 2 ;;
        --target_bed)            require_arg "$1" "${2:-}"; TARGET_BED="$2";       shift 2 ;;
        --extra_args)            require_arg "$1" "${2:-}"; EXTRA_ARGS="$2";       shift 2 ;;
        --input-bam|--input_bam) require_arg "$1" "${2:-}"; INPUT_BAM_HOST="$2";   shift 2 ;;
        --fresh)                 FRESH="1";          shift   ;;
        --detached)              DETACHED=true;       shift   ;;
        *)
            echo "[ERROR] Unknown argument: $1"
            exit 1
            ;;
    esac
done

# Trim leading/trailing whitespace so docker --name matches `docker ps` / daemon filters
ORDER_ID="${ORDER_ID#"${ORDER_ID%%[![:space:]]*}"}"
ORDER_ID="${ORDER_ID%"${ORDER_ID##*[![:space:]]}"}"
if [[ -n "${SAMPLE_NAMES}" ]]; then
    SAMPLE_NAMES="${SAMPLE_NAMES#"${SAMPLE_NAMES%%[![:space:]]*}"}"
    SAMPLE_NAMES="${SAMPLE_NAMES%"${SAMPLE_NAMES##*[![:space:]]}"}"
fi

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

# Docker container NAME must be order_id so `docker ps -a` and nipt-daemon (`docker ps -f name=<order_id>`) work.
if [[ ! "${ORDER_ID}" =~ ^[a-zA-Z0-9][a-zA-Z0-9_.-]*$ ]]; then
    echo "[ERROR] --order_id must be a valid Docker container name: [a-zA-Z0-9][a-zA-Z0-9_.-]*"
    echo "  Got: ${ORDER_ID}"
    exit 1
fi

# ── BAM mode vs FASTQ mode ──────────────────────────────────────────────────
if [[ -n "$INPUT_BAM_HOST" ]]; then
    # ── BAM mode: skip all FASTQ discovery ──────────────────────────────────
    if [[ ! -f "$INPUT_BAM_HOST" ]]; then
        echo "[ERROR] BAM samplesheet not found: ${INPUT_BAM_HOST}"
        exit 1
    fi
    # Convert host path → container-internal path.
    # The BAM CSV lives under HOST_DATA_DIR (mounted at /Work/SgNIPT/data).
    INPUT_BAM_CONTAINER="/Work/SgNIPT/data/${INPUT_BAM_HOST#${HOST_DATA_DIR}/}"
    echo "[INFO] BAM mode: using samplesheet ${INPUT_BAM_HOST}"
    echo "       Container path: ${INPUT_BAM_CONTAINER}"
else
    # ── FASTQ mode: auto-discover or use explicit paths ──────────────────────
    if [[ -z "$FASTQ_R1" ]]; then
        FASTQ_HOST_DIR="${HOST_FASTQ_DIR}/${WORK_DIR}/${ORDER_ID}"
        if [[ ! -d "$FASTQ_HOST_DIR" ]]; then
            echo "[ERROR] FASTQ directory not found: ${FASTQ_HOST_DIR}"
            echo "[HINT]  Pass --work-id and --order-id so that FASTQs live under fastq/{work_id}/{order_id}/"
            echo "       or set --fastq_r1 / --fastq_r2 explicitly, or use --input-bam for BAM mode."
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
        FASTQ_R1="${WORK_DIR}/${ORDER_ID}/${_BASE}"
        FASTQ_R2="${WORK_DIR}/${ORDER_ID}/$(basename "$_R2_ABS")"
        echo "[INFO] Auto-selected FASTQ pair:"
        echo "       R1: ${FASTQ_R1}"
        echo "       R2: ${FASTQ_R2}"
    fi

    if [[ -z "$FASTQ_R1" ]]; then
        echo "[ERROR] --fastq_r1 is required (or use --input-bam for BAM mode)"
        exit 1
    fi
fi

# ── Ensure top-level directories exist ───────────────────────────────────────
for DIR in "$HOST_FASTQ_DIR" "$HOST_DATA_DIR" "$HOST_CONFIG_DIR" \
           "$HOST_ANALYSIS_DIR" "$HOST_OUTPUT_DIR" "$HOST_LOG_DIR"; do
    if [[ ! -d "$DIR" ]]; then
        echo "[WARN] Creating directory: ${DIR}"
        mkdir -p "$DIR"
    fi
done

# Per-sample paths: {type}/{YYMM}/{ORDER_ID}/
ORDER_ANALYSIS_HOST="${HOST_ANALYSIS_DIR}/${WORK_DIR}/${ORDER_ID}"
ORDER_OUTPUT_HOST="${HOST_OUTPUT_DIR}/${WORK_DIR}/${ORDER_ID}"
ORDER_LOG_HOST="${HOST_LOG_DIR}/${WORK_DIR}/${ORDER_ID}"

# 1-C: Create per-sample directories with permission repair
ensure_sample_output_dirs

# --fresh (carrier_screening run_analysis.sh parity): wipe Nextflow cache on host, disable -resume
if [[ -n "$FRESH" ]]; then
    echo "[WARN] --fresh: removing Nextflow work/ and .nextflow under ${ORDER_ANALYSIS_HOST} ..."
    rm -rf "${ORDER_ANALYSIS_HOST}/work" "${ORDER_ANALYSIS_HOST}/.nextflow" 2>/dev/null || true
    if [[ -e "${ORDER_ANALYSIS_HOST}/work" ]] || [[ -e "${ORDER_ANALYSIS_HOST}/.nextflow" ]]; then
        echo "[WARN] Paths still present (often root-owned); trying sudo..."
        if ! sudo rm -rf "${ORDER_ANALYSIS_HOST}/work" "${ORDER_ANALYSIS_HOST}/.nextflow"; then
            echo "[ERROR] Could not remove work/ or .nextflow/. Run:"
            echo "  sudo rm -rf \"${ORDER_ANALYSIS_HOST}/work\" \"${ORDER_ANALYSIS_HOST}/.nextflow\""
            exit 1
        fi
    fi
    echo "  Cleared: ${ORDER_ANALYSIS_HOST}/work"
    echo "  Cleared: ${ORDER_ANALYSIS_HOST}/.nextflow"
    echo ""
fi

# Nextflow -resume: off when --fresh (or SGNIPT_NO_RESUME=1)
SGNIPT_NO_RESUME="${SGNIPT_NO_RESUME:-0}"
[[ -n "$FRESH" ]] && SGNIPT_NO_RESUME="1"

# ── Container name = ORDER_ID (daemon uses docker ps -f name=ORDER_ID) ──────
CONTAINER_NAME="${ORDER_ID}"

# Remove existing container with same name (if any)
if docker ps -a --format '{{.Names}}' | grep -q "^${CONTAINER_NAME}$"; then
    echo "[WARN] Removing existing container: ${CONTAINER_NAME}"
    docker rm -f "${CONTAINER_NAME}" 2>/dev/null || true
fi

# ── Build entrypoint arguments ───────────────────────────────────────────────
ENTRYPOINT_ARGS=(
    --order_id     "$ORDER_ID"
    --sample_name  "$SAMPLE_NAMES"
)

if [[ -n "$INPUT_BAM_HOST" ]]; then
    ENTRYPOINT_ARGS+=(--input_bam "$INPUT_BAM_CONTAINER")
else
    ENTRYPOINT_ARGS+=(--fastq_r1 "$FASTQ_R1")
    [[ -n "$FASTQ_R2" ]] && ENTRYPOINT_ARGS+=(--fastq_r2 "$FASTQ_R2")
fi

[[ -n "$TARGET_BED"  ]] && ENTRYPOINT_ARGS+=(--target_bed "$TARGET_BED")
[[ -n "$EXTRA_ARGS"  ]] && ENTRYPOINT_ARGS+=(--extra_args "$EXTRA_ARGS")

# 1-D: docker.sock is root:docker → container needs the docker GID as supplementary group
DOCKER_GROUP_ARGS=()
if DOCKER_SOCK_GID="$(getent group docker 2>/dev/null | cut -d: -f3)" && [[ -n "${DOCKER_SOCK_GID}" ]]; then
    DOCKER_GROUP_ARGS=(--group-add "${DOCKER_SOCK_GID}")
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
echo "  SGNIPT_ROOT    : ${SGNIPT_ROOT_DIR}"
echo "  Work Dir       : ${WORK_DIR}"
echo "  Analysis       : ${ORDER_ANALYSIS_HOST}"
echo "  Output         : ${ORDER_OUTPUT_HOST}"
echo "  Log            : ${ORDER_LOG_HOST}"
echo "  Docker Image   : ${DOCKER_IMAGE}"
echo "  CPU / Memory   : ${CONTAINER_CPUS} / ${CONTAINER_MEMORY}"
if [[ -n "$FRESH" ]]; then
    echo "  Fresh run      : yes (work + .nextflow cleared, no Nextflow -resume)"
else
    echo "  Fresh run      : no (use --fresh for full re-run)"
fi
echo "============================================================"

#  Volume mounts: {type}/{YYMM}/{ORDER_ID}/ structure
#    fastq/     → whole tree (fastq/2603/NA12878/...)
#    analysis/  → analysis/2603/ORDER_ID (Nextflow work, intermediates)
#    output/    → output/2603/ORDER_ID (JSON, PNG for service-daemon)
#    log/       → log/2603/ORDER_ID (analysis logs)

#  Default: foreground (carrier_screening run_analysis.sh); optional --detached for -d.
COMMON_DOCKER_ARGS=(
    --user "${CHOWN_SPEC}"
    "${DOCKER_GROUP_ARGS[@]}"
    --name "${CONTAINER_NAME}"
    --label "sgnipt.order_id=${ORDER_ID}"
    -e TZ=Asia/Seoul
    -e USER="$(whoami)"
    -e USERNAME="$(whoami)"
    -e HOME=/tmp
    -e FONTCONFIG_PATH=/tmp
    -e CLEANUP_WORK_DIR="${CLEANUP_WORK}"
    -e SGNIPT_NO_RESUME="${SGNIPT_NO_RESUME}"
    -e NXF_HOME="/tmp/.nextflow"
    -e NXF_TEMP="/Work/SgNIPT/analysis/tmp"
    -e WORK_DIR="${WORK_DIR}"
    -e BWA2="bwa-mem2"
    -e SAMTOOLS="samtools"
    -e BCFTOOLS="bcftools"
    -e PYTHON3="python3"
    -v "${HOST_FASTQ_DIR}:/Work/SgNIPT/fastq"
    -v "${HOST_DATA_DIR}:/Work/SgNIPT/data"
    -v "${HOST_CONFIG_DIR}:/Work/SgNIPT/config"
    -v "${ORDER_ANALYSIS_HOST}:/Work/SgNIPT/analysis"
    -v "${ORDER_OUTPUT_HOST}:/Work/SgNIPT/output/${ORDER_ID}"
    -v "${ORDER_LOG_HOST}:/Work/SgNIPT/log/${ORDER_ID}"
    -v /data/reference:/data/reference:ro
    --cpus "${CONTAINER_CPUS}"
    --memory "${CONTAINER_MEMORY}"
    "${DOCKER_IMAGE}"
    "${ENTRYPOINT_ARGS[@]}"
)

set +e
if [[ "${DETACHED}" == "true" ]]; then
    docker run -d "${COMMON_DOCKER_ARGS[@]}"
else
    docker run "${COMMON_DOCKER_ARGS[@]}"
fi
DOCKER_EXIT=$?
set -e

if [[ "${DETACHED}" == "true" ]]; then
    if [[ $DOCKER_EXIT -eq 0 ]]; then
        echo "[INFO] Container '${CONTAINER_NAME}' started successfully (detached)"
        echo "[INFO] Monitor logs     : docker logs -f ${CONTAINER_NAME}"
        echo "[INFO] Analysis dir    : ${ORDER_ANALYSIS_HOST}"
        echo "[INFO] Output dir      : ${ORDER_OUTPUT_HOST}"
        echo "[INFO] Log dir         : ${ORDER_LOG_HOST}"
        echo "[INFO] Service-daemon  : ${ORDER_OUTPUT_HOST}/${ORDER_ID}.json"
        echo "[INFO]                    ${ORDER_OUTPUT_HOST}/${ORDER_ID}.output.tar"
        (
            docker wait "${CONTAINER_NAME}" >/dev/null 2>&1 || true
            repair_order_tree_permissions || true
        ) &
        disown $!
        exit 0
    else
        echo "[ERROR] Failed to start container '${CONTAINER_NAME}' (exit: ${DOCKER_EXIT})"
        exit $DOCKER_EXIT
    fi
else
    repair_order_tree_permissions || true
    if [[ $DOCKER_EXIT -eq 0 ]]; then
        echo "[INFO] Pipeline completed (exit 0)"
        echo "[INFO] Output: ${ORDER_OUTPUT_HOST}/${ORDER_ID}.json"
    else
        echo "[ERROR] Pipeline failed with exit code: ${DOCKER_EXIT}"
    fi
    exit "${DOCKER_EXIT}"
fi
