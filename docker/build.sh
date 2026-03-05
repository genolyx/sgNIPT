#!/bin/bash
###############################################################################
#  Single Gene NIPT - Docker Image Build Script
###############################################################################
#  Usage:
#    ./build.sh                    # Build with default tag (sgnipt:latest)
#    ./build.sh v1.0.0             # Build with specific version tag
#    ./build.sh v1.0.0 --no-cache  # Build without cache
###############################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "${SCRIPT_DIR}")"

VERSION="${1:-latest}"
BUILD_ARGS="${2:-}"

IMAGE_NAME="sgnipt"
FULL_TAG="${IMAGE_NAME}:${VERSION}"

echo "============================================================"
echo "  Building Single Gene NIPT Docker Image"
echo "============================================================"
echo "  Project Dir : ${PROJECT_DIR}"
echo "  Image Tag   : ${FULL_TAG}"
echo "  Build Args  : ${BUILD_ARGS:-none}"
echo "============================================================"

docker build \
    ${BUILD_ARGS} \
    -t "${FULL_TAG}" \
    -f "${SCRIPT_DIR}/Dockerfile" \
    "${PROJECT_DIR}"

# Also tag as latest if version was specified
if [[ "$VERSION" != "latest" ]]; then
    docker tag "${FULL_TAG}" "${IMAGE_NAME}:latest"
    echo "[INFO] Also tagged as ${IMAGE_NAME}:latest"
fi

echo "============================================================"
echo "  Build complete: ${FULL_TAG}"
echo "============================================================"
echo ""
echo "  To verify:  docker run --rm ${FULL_TAG} --help"
echo "  Image size: $(docker images ${FULL_TAG} --format '{{.Size}}')"
echo ""
