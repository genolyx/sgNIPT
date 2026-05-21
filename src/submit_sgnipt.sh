#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# submit_sgnipt.sh  —  service-daemon API를 통해 sgNIPT 분석 제출
#
# 사용법:
#   submit_sgnipt.sh --order-id <ID> --work-id <WID> --input-bam <CSV> [--fresh]
#   submit_sgnipt.sh --order-id <ID> --work-id <WID> --fastq-r1 <R1> --fastq-r2 <R2> [--fresh]
#
# run_sgnipt.sh와 달리 daemon API를 거치므로:
#   - process_results() 자동 실행 → result.json + variant_report 병합
#   - DB(orders.db) 자동 갱신 → 포털 즉시 반영
#   - 포털 진행상황 실시간 업데이트
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail

DAEMON_URL="${SGNIPT_DAEMON_URL:-http://localhost:8003}"
ORDER_ID=""
WORK_ID=""
INPUT_BAM_CSV=""
FASTQ_R1=""
FASTQ_R2=""
FRESH=false

usage() {
    echo "Usage:"
    echo "  $0 --order-id <ID> --work-id <WID> --input-bam <BAM_CSV> [--fresh]"
    echo "  $0 --order-id <ID> --work-id <WID> --fastq-r1 <R1> --fastq-r2 <R2> [--fresh]"
    echo ""
    echo "Options:"
    echo "  --order-id   Order ID (unique per run)"
    echo "  --work-id    Work directory ID (e.g. 2603)"
    echo "  --input-bam  Path to BAM samplesheet CSV (BAM mode)"
    echo "  --fastq-r1   R1 FASTQ path (FASTQ mode)"
    echo "  --fastq-r2   R2 FASTQ path (FASTQ mode)"
    echo "  --fresh      Force fresh run (discard Nextflow cache)"
    echo "  --daemon-url Daemon URL (default: http://localhost:8003)"
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --order-id)   ORDER_ID="$2";      shift 2 ;;
        --work-id)    WORK_ID="$2";       shift 2 ;;
        --input-bam)  INPUT_BAM_CSV="$2"; shift 2 ;;
        --fastq-r1)   FASTQ_R1="$2";      shift 2 ;;
        --fastq-r2)   FASTQ_R2="$2";      shift 2 ;;
        --fresh)      FRESH=true;         shift ;;
        --daemon-url) DAEMON_URL="$2";    shift 2 ;;
        -h|--help)    usage ;;
        *) echo "Unknown option: $1"; usage ;;
    esac
done

[[ -z "$ORDER_ID" ]]  && { echo "ERROR: --order-id required"; usage; }
[[ -z "$WORK_ID" ]]   && { echo "ERROR: --work-id required";  usage; }
[[ -z "$INPUT_BAM_CSV" && -z "$FASTQ_R1" ]] && { echo "ERROR: --input-bam or --fastq-r1/r2 required"; usage; }

# ── Daemon health check ───────────────────────────────────────────────────────
echo "[INFO] Checking daemon at ${DAEMON_URL} ..."
if ! curl -sf "${DAEMON_URL}/health" > /dev/null 2>&1; then
    echo "ERROR: Daemon not reachable at ${DAEMON_URL}"
    exit 1
fi
echo "[INFO] Daemon OK"

# ── Build params JSON ─────────────────────────────────────────────────────────
PY_FRESH="True" ; [[ "$FRESH" == "false" ]] && PY_FRESH="False"

# Build submit body
if [[ -n "$INPUT_BAM_CSV" ]]; then
    BODY=$(python3 -c "
import json
print(json.dumps({
    'order_id': '${ORDER_ID}',
    'service_code': 'sgnipt',
    'work_dir': '${WORK_ID}',
    'params': {'input_bam_csv': '${INPUT_BAM_CSV}', '_pipeline_fresh': ${PY_FRESH}},
}))
")
else
    BODY=$(python3 -c "
import json
print(json.dumps({
    'order_id': '${ORDER_ID}',
    'service_code': 'sgnipt',
    'work_dir': '${WORK_ID}',
    'fastq_r1_path': '${FASTQ_R1}',
    'fastq_r2_path': '${FASTQ_R2}',
    'params': {'_pipeline_fresh': ${PY_FRESH}},
}))
")
fi

# ── Check existing order status ───────────────────────────────────────────────
EXISTING=$(curl -sf "${DAEMON_URL}/order/${ORDER_ID}/status" 2>/dev/null || echo "null")
EXISTING_STATUS=$(echo "$EXISTING" | python3 -c "import json,sys; d=json.load(sys.stdin); print(d.get('status',''))" 2>/dev/null || echo "")

if [[ -n "$EXISTING_STATUS" && "$EXISTING_STATUS" != "FAILED" && "$EXISTING_STATUS" != "CANCELLED" && "$EXISTING_STATUS" != "COMPLETED" && "$EXISTING_STATUS" != "REPORT_READY" && "$EXISTING_STATUS" != "SAVED" ]]; then
    echo "[WARN] Order ${ORDER_ID} already exists with status: ${EXISTING_STATUS}"
    echo "       Use --fresh or wait for completion."
    exit 1
fi

# ── Save order (or re-use existing) ──────────────────────────────────────────
if [[ -z "$EXISTING_STATUS" ]]; then
    echo "[INFO] Saving order ${ORDER_ID} ..."
    SAVE_RESP=$(curl -sf -X POST \
        -H "Content-Type: application/json" \
        -d "$BODY" \
        "${DAEMON_URL}/order/sgnipt/save")
    echo "[INFO] Save response: $SAVE_RESP"
else
    echo "[INFO] Re-using existing order ${ORDER_ID} (status: ${EXISTING_STATUS})"
fi

# ── Start (queue) the order ───────────────────────────────────────────────────
START_BODY="{\"fresh\": $FRESH}"
echo "[INFO] Starting order ${ORDER_ID} (fresh=${FRESH}) ..."
START_RESP=$(curl -sf -X POST \
    -H "Content-Type: application/json" \
    -d "$START_BODY" \
    "${DAEMON_URL}/order/${ORDER_ID}/start")
echo "[INFO] Start response: $START_RESP"

QUEUE_POS=$(echo "$START_RESP" | python3 -c "import json,sys; d=json.load(sys.stdin); print(d.get('queue_position','?'))" 2>/dev/null || echo "?")
echo ""
echo "✓ Order submitted: ${ORDER_ID}"
echo "  Queue position : ${QUEUE_POS}"
echo "  Monitor status : curl -s ${DAEMON_URL}/order/${ORDER_ID}/status"
echo "  Or watch logs  : tail -f /home/ken/sgNIPT/log/${WORK_ID}/${ORDER_ID}/*.log"
