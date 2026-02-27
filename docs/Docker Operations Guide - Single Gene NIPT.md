# Docker Operations Guide - Single Gene NIPT

본 문서는 Single Gene NIPT 파이프라인을 Docker 컨테이너 기반으로 운영하는 방법을 설명합니다. 기존 sWGS NIPT의 nipt-daemon 아키텍처와 동일한 패턴으로, Order 단위로 Docker container를 실행하고 daemon이 동시 실행 수를 관리하는 구조입니다.

## Architecture Overview

```
┌─────────────────────────────────────────────────────────────────┐
│  Order Platform                                                 │
│  (LIMS / Order Management System)                               │
└──────────────────────┬──────────────────────────────────────────┘
                       │ REST API / Message Queue
                       ▼
┌─────────────────────────────────────────────────────────────────┐
│  nipt-daemon                                                    │
│  ┌───────────────────────────────────────────────────────────┐  │
│  │  • Order 큐 관리                                          │  │
│  │  • max_parallel 제어 (동시 컨테이너 수 제한)               │  │
│  │  • .completion_status.json 폴링으로 완료 감지             │  │
│  │  • Order Platform에 상태 보고                             │  │
│  └───────────────────────────────────────────────────────────┘  │
│       │           │           │                                  │
│       ▼           ▼           ▼                                  │
│  ┌─────────┐ ┌─────────┐ ┌─────────┐                           │
│  │Container│ │Container│ │Container│  ← max_parallel=3 예시     │
│  │Order_001│ │Order_002│ │Order_003│                            │
│  │         │ │         │ │         │                            │
│  │Nextflow │ │Nextflow │ │Nextflow │  ← 각 컨테이너 내부에서   │
│  │Pipeline │ │Pipeline │ │Pipeline │    Nextflow가 샘플 병렬 처리│
│  └─────────┘ └─────────┘ └─────────┘                           │
└─────────────────────────────────────────────────────────────────┘
```

## Host Directory Structure

`setup_host_dirs.sh` 스크립트로 초기 구조를 생성합니다.

```bash
./docker/setup_host_dirs.sh /data/SgNIPT
```

생성되는 구조는 다음과 같습니다.

```
/data/SgNIPT/                          ← SGNIPT_ROOT_DIR
├── fastq/                             ← FASTQ 파일 저장소
│   ├── SAMPLE_A_R1.fastq.gz
│   ├── SAMPLE_A_R2.fastq.gz
│   ├── SAMPLE_B_R1.fastq.gz
│   └── SAMPLE_B_R2.fastq.gz
├── data/                              ← 레퍼런스 데이터 (읽기 전용)
│   ├── reference/
│   │   ├── hg38.fa
│   │   ├── hg38.fa.fai
│   │   └── hg38.dict
│   ├── bwa_index/
│   │   ├── hg38.fa.0123
│   │   ├── hg38.fa.amb
│   │   ├── hg38.fa.ann
│   │   ├── hg38.fa.bwt.2bit.64
│   │   ├── hg38.fa.pac
│   │   └── hg38.fa.sa
│   ├── target_bed/
│   │   └── sg_nipt_panel_v1.bed
│   └── annotation/                    (optional)
│       └── clinvar.vcf.gz
├── config/                            ← 설정 파일 (optional)
│   ├── params.yaml                    (Nextflow 파라미터 오버라이드)
│   └── qc_thresholds.json             (QC 임계값 커스터마이징)
├── analysis/                          ← Nextflow work 디렉토리 (중간 파일)
│   ├── ORDER_001_samplesheet.csv
│   ├── ORDER_001_work/
│   └── ORDER_002_work/
├── output/                            ← 최종 결과물
│   ├── ORDER_001/
│   │   ├── SAMPLE_A/
│   │   ├── SAMPLE_B/
│   │   ├── multiqc/
│   │   └── .completion_status.json
│   └── ORDER_002/
│       └── ...
└── log/                               ← 실행 로그
    ├── ORDER_001/
    │   ├── nextflow_20260222_103000.log
    │   ├── trace.txt
    │   ├── timeline.html
    │   └── report.html
    └── entrypoint.log
```

## Volume Mount Mapping

Docker container 내부에서는 다음 경로로 접근합니다.

| Host Path | Container Path | 용도 | 접근 모드 |
|-----------|---------------|------|----------|
| `$ROOT/fastq` | `/Work/SgNIPT/fastq` | FASTQ 입력 파일 | Read-Only |
| `$ROOT/data` | `/Work/SgNIPT/data` | Reference, Index, BED | Read-Only |
| `$ROOT/config` | `/Work/SgNIPT/config` | 설정 파일 | Read-Only |
| `$ROOT/analysis` | `/Work/SgNIPT/analysis` | Nextflow work dir | Read-Write |
| `$ROOT/output` | `/Work/SgNIPT/output` | 최종 결과물 | Read-Write |
| `$ROOT/log` | `/Work/SgNIPT/log` | 실행 로그 | Read-Write |

## Docker Image Build

```bash
# 기본 빌드 (sgnipt:latest)
./docker/build.sh

# 버전 태그 지정
./docker/build.sh v1.0.0

# 캐시 없이 빌드
./docker/build.sh v1.0.0 --no-cache
```

## Running an Order

### 1. 환경 설정

```bash
# .env 파일 생성
cp docker/.env.example docker/.env

# 편집하여 경로 설정
vi docker/.env
```

### 2. 단일 Order 실행

```bash
./docker/run_sgnipt.sh \
    --order_id ORDER_001 \
    --sample_name "SAMPLE_A,SAMPLE_B" \
    --fastq_r1 "SAMPLE_A_R1.fastq.gz,SAMPLE_B_R1.fastq.gz" \
    --fastq_r2 "SAMPLE_A_R2.fastq.gz,SAMPLE_B_R2.fastq.gz"
```

### 3. 옵션 지정 실행

```bash
./docker/run_sgnipt.sh \
    --order_id ORDER_002 \
    --sample_name "SAMPLE_C" \
    --fastq_r1 "SAMPLE_C_R1.fastq.gz" \
    --fastq_r2 "SAMPLE_C_R2.fastq.gz" \
    --target_bed "panel_v2.bed" \
    --variant_caller freebayes \
    --extra_args "--vc_min_depth 100 --mpileup_max_depth 50000"
```

### 4. 모니터링

```bash
# 실시간 로그 확인
docker logs -f sgnipt_ORDER_001

# 컨테이너 상태 확인
docker ps --filter "name=sgnipt_"

# 완료 상태 확인 (daemon 폴링용)
cat /data/SgNIPT/output/ORDER_001/.completion_status.json
```

## Completion Status JSON

각 Order 실행이 완료되면 `.completion_status.json` 파일이 생성됩니다. nipt-daemon은 이 파일을 폴링하여 완료 여부를 감지합니다.

성공 시:
```json
{
    "order_id": "ORDER_001",
    "status": "COMPLETED",
    "timestamp": "2026-02-22T10:45:00+09:00",
    "output_dir": "/Work/SgNIPT/output/ORDER_001",
    "samples": "SAMPLE_A,SAMPLE_B",
    "exit_code": 0
}
```

실패 시:
```json
{
    "order_id": "ORDER_001",
    "status": "FAILED",
    "timestamp": "2026-02-22T10:45:00+09:00",
    "output_dir": "/Work/SgNIPT/output/ORDER_001",
    "samples": "SAMPLE_A,SAMPLE_B",
    "exit_code": 1
}
```

## nipt-daemon Integration

기존 nipt-daemon에서 Single Gene NIPT를 호출하려면, sWGS NIPT의 `docker run` 호출 부분을 `run_sgnipt.sh` 호출로 교체하면 됩니다. daemon이 관리해야 하는 핵심 사항은 다음과 같습니다.

**Order 큐잉**: Order Platform에서 수신한 분석 요청을 큐에 적재합니다.

**동시 실행 제어**: `max_parallel` 설정에 따라 동시에 실행 중인 `sgnipt_*` 컨테이너 수를 제한합니다. 현재 실행 중인 컨테이너 수는 `docker ps --filter "name=sgnipt_" -q | wc -l` 명령으로 확인할 수 있습니다.

**완료 감지**: `$ROOT/output/{ORDER_ID}/.completion_status.json` 파일의 존재 여부와 `status` 필드를 주기적으로 폴링합니다.

**상태 보고**: 완료 감지 후 Order Platform에 결과를 보고하고, 필요시 컨테이너를 정리(`docker rm`)합니다.

**리소스 관리**: 각 컨테이너의 CPU/메모리 제한은 `SGNIPT_CPUS`, `SGNIPT_MEMORY` 환경변수로 제어합니다. 서버 전체 리소스를 고려하여 `max_parallel × container_resources ≤ total_server_resources` 조건을 유지해야 합니다.

## Resource Planning

서버 리소스에 따른 권장 설정 예시입니다.

| Server Spec | max_parallel | SGNIPT_CPUS | SGNIPT_MEMORY | 비고 |
|-------------|-------------|-------------|---------------|------|
| 32 CPU / 64 GB | 2 | 16 | 32g | 표준 구성 |
| 64 CPU / 128 GB | 4 | 16 | 32g | 고성능 서버 |
| 16 CPU / 32 GB | 1 | 16 | 28g | 단일 실행 |

## Troubleshooting

**컨테이너가 즉시 종료되는 경우**: `docker logs sgnipt_ORDER_ID`로 entrypoint 로그를 확인합니다. 대부분 레퍼런스 파일 누락이나 FASTQ 경로 오류가 원인입니다.

**Nextflow 실패 후 재실행**: 동일한 Order ID로 다시 실행하면 기존 컨테이너가 자동 제거되고 새로 시작됩니다. Nextflow work 디렉토리가 남아있으면 `-resume` 기능이 작동합니다. 이를 활용하려면 `EXTRA_ARGS`에 `-resume`을 추가하세요.

**디스크 공간 관리**: `analysis/` 디렉토리에 Nextflow work 파일이 누적됩니다. `SGNIPT_CLEANUP_WORK=true`로 설정하면 성공 후 자동 삭제되지만, 이 경우 `-resume` 기능을 사용할 수 없습니다.
