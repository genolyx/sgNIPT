# Data Preparation Guide - Twist Panel BED Files

본 문서는 Single Gene NIPT 파이프라인에서 사용하는 Twist Custom Panel (UCL_SingleGeneNIPT_TE-96276661_hg38) BED 파일들의 배치 방법과 용도를 설명합니다.

## Panel Overview

| 항목 | 값 |
|------|-----|
| Panel ID | TE-96276661 |
| Reference | hg38 |
| Target regions | 5,138개 |
| Total target size | ~1.07 Mb |
| Probe footprint | ~780 kb |
| Genes | 183개 (22 chromosomes) |
| Zero-probe regions | 29개 |
| Y-chromosome targets | 없음 |

## Host 디렉토리 구조

```
{SGNIPT_WORK_DIR}/data/
├── reference/
│   ├── hg38.fa                     ← Reference FASTA
│   ├── hg38.fa.fai                 ← samtools faidx index
│   └── hg38.dict                   ← Sequence dictionary (GATK용)
├── bwa_index/
│   ├── hg38.fa.0123                ← BWA-MEM2 index files
│   ├── hg38.fa.amb
│   ├── hg38.fa.ann
│   ├── hg38.fa.bwt.2bit.64
│   ├── hg38.fa.pac
│   └── hg38.fa.sa                  ← (optional, BWA-MEM2에서 불필요할 수 있음)
├── target_bed/
│   ├── All_target_regions_UCL_SingleGeneNIPT_TE-96276661_hg38_*.bed
│   ├── Probes_merged_ok_shareable_UCL_SingleGeneNIPT_TE-96276661_hg38_*.bed
│   ├── Target_regions_with_zero_probes_UCL_SingleGeneNIPT_TE-96276661_hg38_*.bed
│   ├── ff_snps_hg38.bed            ← 사용자가 생성해야 함 (아래 참조)
│   └── gene_list.tsv               ← singleGene_NIPT.xlsx에서 변환
├── annotation/
│   ├── clinvar.vcf.gz              ← (optional) ClinVar VCF
│   └── clinvar.vcf.gz.tbi
└── ...
```

## BED 파일별 용도

### 1. All_target_regions (필수)

Twist에서 제공하는 전체 타겟 영역 BED 파일입니다. 파이프라인의 **핵심 BED 파일**로, 다음 단계에서 사용됩니다:

- Variant Calling: 이 영역 내에서만 변이를 호출합니다.
- BAM QC: 타겟 영역 커버리지 및 on-target rate를 계산합니다.
- Target Read Extraction: 이 영역과 겹치는 리드만 추출합니다.

**파일 이름 패턴**: `All_target_regions*.bed`

**자동 감지**: entrypoint.sh가 `data/target_bed/` 디렉토리에서 자동으로 감지합니다.

### 2. ff_snps BED (권장)

Fetal Fraction 추정에 사용할 **Common SNP 위치 목록**입니다. Twist에서 제공하지 않으므로 **사용자가 직접 생성**해야 합니다.

**생성 방법**:

```bash
# 방법 1: dbSNP common SNPs + target regions 교차
# dbSNP common_all VCF에서 MAF > 5% SNP 추출 후 target regions과 교차
bcftools query -f '%CHROM\t%POS0\t%END\t%ID\n' \
    -i 'INFO/COMMON=1 && INFO/CAF[1] > 0.05 && INFO/CAF[1] < 0.95' \
    dbSNP_common_all.vcf.gz \
| bedtools intersect -a stdin -b All_target_regions*.bed \
> ff_snps_hg38.bed

# 방법 2: KRGDB (한국인 게놈 DB) 활용
# KRGDB에서 MAF > 5% SNP 추출 후 동일하게 교차
```

**예상 SNP 수**: 약 150~230개 (target 영역 ~1.07 Mb 기준)

**파일 이름 패턴**: `ff_snps*.bed`

**미제공 시**: target_bed를 대신 사용하지만, FF 추정 정확도가 떨어질 수 있습니다.

### 3. Probes_merged BED (권장)

Twist에서 제공하는 실제 프로브가 커버하는 영역입니다. QC에서 **프로브 커버리지**를 평가하는 데 사용됩니다.

**파일 이름 패턴**: `Probes_merged*.bed`

**미제공 시**: 프로브 커버리지 QC가 스킵됩니다.

### 4. Target_regions_with_zero_probes BED (권장)

프로브가 설계되지 않은 29개 타겟 영역입니다. 이 영역에서 발견된 변이에는 **경고 플래그**가 추가됩니다.

**파일 이름 패턴**: `Target_regions_with_zero_probes*.bed`

**미제공 시**: Zero-probe 영역 경고가 스킵됩니다.

### 5. Gene List (권장)

singleGene_NIPT.xlsx에서 변환한 유전자 목록입니다. 패널이 모든 타겟 유전자를 커버하는지 검증합니다.

**변환 방법**:

```bash
# Excel → TSV 변환 (Python)
python3 -c "
import pandas as pd
df = pd.read_excel('singleGene_NIPT.xlsx')
df.to_csv('gene_list.tsv', sep='\t', index=False)
"
```

**파일 이름 패턴**: `gene_list*` 또는 `singleGene*`

## 파이프라인 파라미터 매핑

| BED 파일 | Nextflow 파라미터 | 자동 감지 패턴 |
|----------|-------------------|----------------|
| All_target_regions | `--target_bed` | `All_target_regions*.bed` |
| ff_snps | `--ff_snps_bed` | `ff_snps*.bed` |
| Probes_merged | `--probes_bed` | `Probes_merged*.bed` |
| Zero-probe regions | `--zero_probe_bed` | `Target_regions_with_zero_probes*.bed` |
| Gene list | `--gene_list` | `gene_list*` 또는 `singleGene*` |

## Fetal Fraction 추정 방법

본 패널에는 Y 염색체 타겟이 없으므로, 다음 3가지 방법을 조합하여 FF를 추정합니다:

1. **SNP-based (Primary)**: ff_snps_bed 내 common SNP에서 informative SNP를 식별하고, minor allele frequency 패턴으로 FF를 계산합니다.
2. **ChrX-based (Supplementary)**: 남아 태아의 경우, X 염색체 타겟 영역(18개 유전자, ~97 kb)의 리드 비율로 FF를 추정합니다.
3. **Fragment size-based (Supplementary)**: cfDNA 단편 크기 분포에서 태아 유래 짧은 단편 비율로 FF를 추정합니다.

## 주의사항

- 모든 BED 파일은 **hg38 좌표계**여야 합니다.
- BED 파일에 `track` 또는 `browser` 헤더 라인이 있어도 파이프라인이 자동으로 무시합니다.
- ff_snps_bed가 없으면 target_bed 전체를 사용하므로, 비-SNP 영역에서 불필요한 계산이 발생합니다.
- Zero-probe 영역(29개)에서 검출된 변이는 리포트에 경고와 함께 표시되지만, 임상적 판단에서 제외하는 것을 권장합니다.
