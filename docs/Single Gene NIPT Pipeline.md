# Single Gene NIPT Pipeline

임산부 혈액 내 cell-free DNA(cfDNA)를 이용하여 태아의 단일 유전자 질환(Single Gene Disorder)을 비침습적으로 스크리닝하는 자동화 파이프라인입니다. Illumina 시퀀서에서 생성된 FASTQ 파일을 입력으로 받아 hg38 레퍼런스 게놈에 정렬한 뒤, 타겟 패널 영역에서 태아 특이적 변이를 검출하고 Fetal Fraction을 추정합니다.

## Pipeline Overview

파이프라인은 Nextflow DSL2로 구현되어 있으며, 각 분석 단계는 독립적인 모듈로 분리되어 유지보수와 확장이 용이합니다. 핵심 분석 로직은 Python으로 작성되었고, 생물정보학 도구들(fastp, BWA-MEM2, samtools, bcftools, mosdepth 등)을 조합하여 완전한 워크플로우를 구성합니다.

```
FASTQ Input
    │
    ▼
┌─────────────────────┐
│  1. Preprocessing   │  fastp: adapter trimming, quality filtering
│     (FASTP)         │
└────────┬────────────┘
         │
         ▼
┌─────────────────────┐
│  2. Alignment       │  BWA-MEM2: hg38 reference mapping
│     (BWA-MEM2)      │
└────────┬────────────┘
         │
         ▼
┌─────────────────────┐
│  3. Deduplication   │  samtools markdup / Picard MarkDuplicates
│     (MarkDup)       │
└────────┬────────────┘
         │
    ┌────┴─────────────────────────┐
    │                              │
    ▼                              ▼
┌──────────────┐          ┌──────────────────┐
│ 4. Target    │          │ 5. BAM QC        │
│    Extraction│          │    (flagstat,     │
│              │          │     mosdepth,     │
└──────┬───────┘          │     stats)        │
       │                  └──────────────────┘
       │
  ┌────┴──────────────────────┐
  │                           │
  ▼                           ▼
┌──────────────────┐  ┌──────────────────────┐
│ 6. Fetal         │  │ 7. Variant Calling   │
│    Fraction      │  │    (bcftools/        │
│    Estimation    │  │     freebayes/GATK)  │
│    (SNP-based)   │  │                      │
└────────┬─────────┘  └──────────┬───────────┘
         │                       │
         └───────────┬───────────┘
                     │
                     ▼
         ┌───────────────────────┐
         │ 8. Variant Analysis   │
         │    (Fetal genotype    │
         │     determination)    │
         └───────────┬───────────┘
                     │
                     ▼
         ┌───────────────────────┐
         │ 9. Report Generation  │
         │    (JSON + HTML)      │
         └───────────────────────┘
```

## Requirements

파이프라인 실행에 필요한 소프트웨어는 다음과 같습니다. Conda 환경 파일(`environment.yml`)을 통해 일괄 설치할 수 있습니다.

| Category | Tool | Version | Purpose |
|----------|------|---------|---------|
| Workflow | Nextflow | ≥23.04 | 파이프라인 오케스트레이션 |
| Language | Python | ≥3.9 | 분석 스크립트 실행 |
| Preprocessing | fastp | ≥0.23 | 어댑터 트리밍, 품질 필터링 |
| Alignment | BWA-MEM2 | ≥2.2 | 시퀀싱 리드 정렬 |
| BAM Processing | samtools | ≥1.17 | BAM 처리, 통계, 인덱싱 |
| BAM Processing | Picard | ≥3.0 | 중복 마킹 (선택) |
| Variant Calling | bcftools | ≥1.17 | 변이 호출 및 필터링 |
| Variant Calling | freebayes | ≥1.3 | 변이 호출 (선택) |
| Coverage | mosdepth | ≥0.3 | 타겟 커버리지 분석 |
| QC | MultiQC | ≥1.14 | QC 메트릭 통합 보고서 |

## Quick Start

### 1. 환경 설정

```bash
# Conda 환경 생성
conda env create -f environment.yml
conda activate single-gene-nipt

# 또는 개별 설치
# pip install numpy pandas
```

### 2. Reference 준비

```bash
# hg38 reference genome 다운로드 (이미 있다면 생략)
# BWA-MEM2 인덱스 생성
bwa-mem2 index /path/to/hg38.fa
```

### 3. Samplesheet 작성

CSV 형식으로 샘플 정보를 기입합니다. Paired-end와 Single-end 모두 지원합니다.

```csv
sample_id,fastq_1,fastq_2
NIPT_001,/data/NIPT_001_R1.fastq.gz,/data/NIPT_001_R2.fastq.gz
NIPT_002,/data/NIPT_002_R1.fastq.gz,/data/NIPT_002_R2.fastq.gz
```

### 4. 파이프라인 실행

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --target_bed /path/to/target_panel.bed \
    --reference_fasta /path/to/hg38.fa \
    --bwa_index_dir /path/to/bwa-mem2-index/ \
    --outdir ./results \
    -profile standard
```

파라미터 파일을 사용하는 경우:

```bash
nextflow run main.nf -params-file conf/params.example.yaml -profile standard
```

### 5. SLURM 클러스터 실행

```bash
nextflow run main.nf \
    -params-file conf/params.example.yaml \
    -profile slurm \
    -resume
```

## Input Files

### Samplesheet (CSV)

| Column | Required | Description |
|--------|----------|-------------|
| `sample_id` | Yes | 고유 샘플 식별자 |
| `fastq_1` | Yes | Read 1 FASTQ 파일 경로 |
| `fastq_2` | No | Read 2 FASTQ 파일 경로 (paired-end) |

### Target BED File

타겟 패널 BED 파일은 스크리닝 대상 유전자 변이 위치를 정의합니다. 확장 컬럼을 통해 유전자명, 질환명, 병원성 변이 정보를 포함할 수 있습니다.

```
#chrom  start       end         name            score  strand  gene   disease              pathogenic_variant
chr7    117559590   117559790   CFTR_F508del    .      +       CFTR   Cystic_Fibrosis      c.1521_1523delCTT
chr11   5226774     5226974     HBB_E6V         .      -       HBB    Sickle_Cell_Disease   c.20A>T
```

| Column | Index | Required | Description |
|--------|-------|----------|-------------|
| chrom | 1 | Yes | 염색체 |
| start | 2 | Yes | 시작 위치 (0-based) |
| end | 3 | Yes | 종료 위치 |
| name | 4 | Recommended | 타겟 이름 |
| score | 5 | No | 스코어 |
| strand | 6 | No | 가닥 방향 |
| gene | 7 | Recommended | 유전자명 |
| disease | 8 | Recommended | 관련 질환명 |
| pathogenic_variant | 9 | No | 병원성 변이 HGVS 표기 |

## Output Structure

```
results/
├── NIPT_001/
│   ├── fastp/
│   │   ├── NIPT_001.fastp.json          # fastp QC 원본 데이터
│   │   └── NIPT_001.fastp.html          # fastp QC HTML 보고서
│   ├── alignment/
│   │   ├── NIPT_001.dedup.bam           # 중복 제거된 BAM
│   │   ├── NIPT_001.dedup.bam.bai       # BAM 인덱스
│   │   ├── NIPT_001.target.bam          # 타겟 영역 BAM
│   │   └── NIPT_001.dedup_metrics.txt   # 중복 제거 메트릭
│   ├── qc/
│   │   ├── NIPT_001.fastq_qc.json       # FASTQ QC 평가 결과
│   │   ├── NIPT_001.bam_qc.json         # BAM QC 평가 결과
│   │   ├── NIPT_001.flagstat.txt         # samtools flagstat
│   │   ├── NIPT_001.stats.txt            # samtools stats
│   │   ├── NIPT_001.idxstats.txt         # samtools idxstats
│   │   └── NIPT_001.mosdepth.*           # mosdepth 커버리지
│   ├── fetal_fraction/
│   │   ├── NIPT_001.fetal_fraction.json  # FF 추정 결과
│   │   ├── NIPT_001.fragment_sizes.tsv   # Fragment size 분포
│   │   ├── NIPT_001.fragment_stats.json  # Fragment 통계
│   │   └── NIPT_001.ff_variants.vcf      # FF 추정용 VCF
│   ├── variants/
│   │   ├── NIPT_001.target_variants.vcf.gz     # 타겟 변이 VCF
│   │   ├── NIPT_001.filtered_variants.vcf      # 필터링된 변이
│   │   └── NIPT_001.variant_report.json        # 변이 분석 보고서
│   └── report/
│       ├── NIPT_001.final_report.json    # 최종 JSON 보고서
│       └── NIPT_001.final_report.html    # 최종 HTML 보고서
├── multiqc/
│   ├── multiqc_report.html               # MultiQC 통합 보고서
│   └── multiqc_data/                     # MultiQC 원본 데이터
└── pipeline_info/
    ├── timeline.html                     # 실행 타임라인
    ├── report.html                       # Nextflow 실행 보고서
    ├── trace.txt                         # 프로세스 추적 로그
    └── dag.html                          # 워크플로우 DAG
```

## Fetal Fraction Estimation

Fetal Fraction(FF)은 모체 혈장 내 cfDNA 중 태아 유래 DNA의 비율을 의미합니다. 본 파이프라인은 세 가지 방법을 조합하여 FF를 추정합니다.

### Primary Method: SNP-based Informative SNP

타겟 영역의 SNP 위치에서 모체가 homozygous reference(AA)이고 태아가 heterozygous(AB)인 **informative SNP**를 식별합니다. 이러한 위치에서 minor allele frequency(MAF)는 FF/2에 비례하므로, 다음 공식으로 FF를 계산합니다:

> **FF_i = 2 × (alt_reads / total_reads)** at each informative SNP
>
> **FF = median(FF_1, FF_2, ..., FF_n)**

IQR 기반 이상치 제거와 bootstrap 신뢰구간 계산을 통해 추정의 견고성을 확보합니다.

### Supplementary Methods

**Y-Chromosome Method**: Y 염색체 매핑 리드 비율로 FF를 추정합니다. 남아 태아에서만 적용 가능하며, 태아 성별 예측에도 활용됩니다.

**Fragment Size Method**: cfDNA fragment size 분포를 분석합니다. 태아 유래 cfDNA는 모체보다 짧은 경향(~143bp vs ~166bp)이 있어, 짧은 fragment 비율로 FF를 보조적으로 추정합니다.

## Variant Analysis

변이 분석 모듈은 추정된 Fetal Fraction을 기반으로 각 변이의 기원을 분류합니다.

| Classification | VAF Pattern | Description |
|----------------|-------------|-------------|
| Maternal Hom Alt | VAF ≈ 100% | 모체 homozygous alt, 태아도 최소 1 copy 보유 |
| Maternal Het | VAF ≈ 50% | 모체 heterozygous, 태아 유전 여부 불확실 |
| Fetal-Specific | VAF ≈ FF/2 | 부계 유래 태아 특이적 변이 |
| Maternal Het + Fetal Inherited | VAF ≈ 50% + FF/2 | 모체 보인자, 태아가 변이 유전 |
| Maternal Het + Fetal Not Inherited | VAF ≈ 50% - FF/2 | 모체 보인자, 태아가 변이 미유전 |

통계적 유의성은 binomial test를 통해 평가하며, 각 변이에 대해 confidence level(high/medium/low)을 부여합니다.

## Configuration

### Variant Caller 선택

파이프라인은 세 가지 variant caller를 지원합니다.

| Caller | Parameter Value | 특징 |
|--------|----------------|------|
| bcftools | `bcftools` (default) | 빠르고 가벼움, mpileup + call 방식 |
| freebayes | `freebayes` | Bayesian 기반, pooled-continuous 모드로 저빈도 변이 검출에 유리 |
| GATK Mutect2 | `gatk` | 체세포 변이 검출에 최적화, 가장 정교하지만 느림 |

### QC Threshold 커스터마이징

`conf/qc_thresholds.json` 파일을 수정하거나 별도 JSON 파일을 작성하여 QC 기준을 조정할 수 있습니다. 각 모듈의 기본 임계값은 코드 내에 정의되어 있으며, 외부 파일로 오버라이드됩니다.

### Execution Profiles

| Profile | Description |
|---------|-------------|
| `standard` | 로컬 실행 |
| `docker` | Docker 컨테이너 사용 |
| `singularity` | Singularity 컨테이너 사용 |
| `conda` | Conda 환경 사용 |
| `slurm` | SLURM 클러스터 제출 |
| `sge` | SGE 클러스터 제출 |
| `test` | 테스트 데이터로 파이프라인 검증 |
| `high_depth` | 고심도 시퀀싱 데이터 최적화 |
| `low_resource` | 저사양 환경 최적화 |

## Project Structure

```
single-gene-nipt/
├── main.nf                    # 메인 워크플로우
├── nextflow.config            # Nextflow 설정
├── environment.yml            # Conda 환경 정의
├── README.md                  # 이 문서
├── bin/                       # Python 분석 스크립트
│   ├── fastq_qc.py           #   FASTQ QC 평가
│   ├── bam_qc.py             #   BAM QC 평가
│   ├── fetal_fraction.py     #   Fetal Fraction 추정
│   ├── variant_analysis.py   #   변이 분석 및 분류
│   ├── collect_fragment_sizes.py  # Fragment size 수집
│   └── generate_report.py    #   최종 보고서 생성
├── modules/                   # Nextflow 프로세스 모듈
│   ├── preprocessing.nf      #   FASTQ 전처리
│   ├── alignment.nf          #   정렬 및 BAM 처리
│   ├── bam_qc.nf             #   BAM QC 및 커버리지
│   ├── fetal_fraction.nf     #   Fetal Fraction 추정
│   ├── variant_calling.nf    #   변이 호출
│   └── reporting.nf          #   보고서 생성
├── conf/                      # 설정 파일
│   ├── params.example.yaml   #   파라미터 예시
│   └── qc_thresholds.json    #   QC 임계값 설정
├── assets/                    # 예시 파일
│   ├── samplesheet.example.csv
│   └── target_panel.example.bed
├── docs/                      # 추가 문서
└── tests/                     # 테스트 데이터
```

## Disclaimer

> 본 파이프라인의 결과는 연구 목적으로만 사용되어야 합니다. 임상적 의사결정은 반드시 확진 검사(진단적 검사)를 통해 확인되어야 합니다.

## License

MIT License
