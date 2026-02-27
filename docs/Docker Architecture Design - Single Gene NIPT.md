# Docker Architecture Design - Single Gene NIPT

## 기존 sWGS NIPT 디렉토리 매핑 (참고)

| Host | Container | 용도 |
|------|-----------|------|
| $HOST_FASTQ_DIR | /Work/NIPT/fastq | FASTQ 입력 |
| $HOST_ANALYSIS_DIR | /Work/NIPT/analysis | 중간 분석 파일 |
| $HOST_LOG_DIR | /Work/NIPT/log | 로그 |
| $HOST_DATA_DIR | /Work/NIPT/data | 레퍼런스, 인덱스 등 |
| $HOST_CONFIG_DIR | /Work/NIPT/config | 설정 파일 |
| $HOST_OUTPUT_DIR | /Work/NIPT/output | 최종 결과물 |

## Single Gene NIPT 디렉토리 매핑 설계

| Host | Container | 용도 |
|------|-----------|------|
| $HOST_FASTQ_DIR | /Work/SgNIPT/fastq | FASTQ 입력 파일 |
| $HOST_DATA_DIR | /Work/SgNIPT/data | hg38 reference, BWA index, target BED |
| $HOST_CONFIG_DIR | /Work/SgNIPT/config | nextflow params, QC thresholds |
| $HOST_ANALYSIS_DIR | /Work/SgNIPT/analysis | Nextflow work dir (중간 파일) |
| $HOST_OUTPUT_DIR | /Work/SgNIPT/output | 최종 결과물 (reports, VCFs) |
| $HOST_LOG_DIR | /Work/SgNIPT/log | Nextflow 로그, pipeline_info |

## Container 내부 고정 경로 (파이프라인 코드)

| Container Path | 내용 |
|----------------|------|
| /Work/SgNIPT/pipeline/ | main.nf, modules/, bin/, nextflow.config |
| /Work/SgNIPT/pipeline/bin/ | Python 분석 스크립트 |
| /Work/SgNIPT/pipeline/modules/ | Nextflow 프로세스 모듈 |

## data 디렉토리 구조 (host에서 준비)

```
$HOST_DATA_DIR/
├── reference/
│   ├── hg38.fa
│   ├── hg38.fa.fai
│   └── hg38.dict
├── bwa_index/
│   ├── hg38.fa.0123
│   ├── hg38.fa.amb
│   ├── hg38.fa.ann
│   ├── hg38.fa.bwt.2bit.64
│   ├── hg38.fa.pac
│   └── hg38.fa.sa  (BWA-MEM2 인덱스 파일들)
├── target_bed/
│   └── sg_nipt_panel.bed
└── annotation/  (optional)
    └── clinvar.vcf.gz
```
