# Portal sgNIPT Variant Table — 개발 요청 Prompt

---

## 배경 및 목적

sgNIPT(Single Gene NIPT) 파이프라인은 산모 cfDNA에서 태아의 단일 유전자 질환 변이를 검출합니다.  
파이프라인 최종 결과는 Daemon이 수집하는 **Order-level JSON** (`{ORDER_ID}.json`) 으로 출력됩니다.  
이 JSON 안의 `variant_analysis.fetal_specific_detail` 및 `pathogenic_details` 배열을 기반으로  
Portal에 **Variant 결과 테이블**을 구현해 주세요.

UI의 기준은 기존 **gx-exome Carrier test 결과 테이블**과 동일한 컬럼 구성을 따르되,  
sgNIPT 고유 컬럼(Origin, VAF, Confidence)을 추가합니다.

---

## 1. JSON 구조 (소스 데이터)

아래는 Daemon JSON의 관련 부분 스키마입니다 (`samples[n].variant_analysis`):

```json
{
  "variant_analysis": {
    "vep_annotated": true,
    "total_target_loci": 500,
    "covered_loci": 498,
    "total_variants": 45,
    "fetal_variants": 3,
    "pathogenic_variants": 2,

    "fetal_specific_detail": [
      {
        "chrom": "chr17",
        "pos": 43092918,
        "ref": "C",
        "alt": "T",
        "vaf": 0.077,
        "depth": 52,
        "origin": "fetal_specific",
        "gene": "MECP2",
        "disease": "Rett syndrome",
        "target_name": "MECP2_exon4",
        "pathogenic_variant": "NM_004992.4:c.916C>T",
        "fetal_genotype": "het",
        "maternal_genotype": "hom_ref",
        "confidence": "high",
        "details": "VAF 7.7% is consistent with fetal-specific het at FF=15.4%",
        "in_zero_probe_region": false,
        "warning": null,
        "vep": {
          "hgvsc": "NM_004992.4:c.916C>T",
          "hgvsp": "NM_004992.4:p.Arg306Cys",
          "consequence": "missense_variant",
          "impact": "MODERATE",
          "symbol": "MECP2",
          "mane": "NM_004992.4",
          "gnomad_af": null,
          "existing": "rs28935468",
          "sift": "deleterious(0)",
          "polyphen": "probably_damaging(0.998)",
          "clinvar": {
            "clnsig": "Pathogenic",
            "clndn": "Rett syndrome",
            "clnrevstat": "criteria_provided,_multiple_submitters,_no_conflicts",
            "clnhgvs": "NM_004992.4:c.916C>T"
          }
        }
      }
    ],

    "pathogenic_details": [
      /* fetal_specific_detail과 동일한 스키마. origin 필드 포함. */
      /* pathogenic_variant 필드가 비어 있지 않거나 vep.clinvar.clnsig 가 P/LP인 항목 */
    ]
  }
}
```

### 테이블 데이터 소스 결정 규칙

| 상황 | 사용할 배열 |
|------|------------|
| 전체 변이 보기 (`All` 탭) | `fetal_specific_detail` (최대 20개) |
| Pathogenic/LP 필터 | `pathogenic_details` (최대 10개) |
| Fetal only 필터 | `fetal_specific_detail` 중 `origin === "fetal_specific"` |
| Maternal 필터 | `fetal_specific_detail` 중 `origin` 이 `maternal` 또는 `maternal_het` |

---

## 2. 컬럼 정의

### Tier 1 — 항상 표시 (숨김 불가)

| 컬럼명 | JSON 경로 | 표시 형식 | 비고 |
|--------|-----------|-----------|------|
| **Gene** | `gene` (또는 `vep.symbol`) | Bold 텍스트 | vep.symbol 우선, 없으면 gene |
| **HGVSc** | `vep.hgvsc` | monospace, transcript prefix(NM_.../ENST...) 제거하고 `c.xxx` 부분만 표시 | vep 없으면 `pathogenic_variant` 필드 |
| **HGVSp** | `vep.hgvsp` | monospace, `p.xxx` 부분만 표시, HGVSc 아래 secondary text | missense/nonsense만 유의미. 없으면 `—` |
| **ClinVar** | `vep.clinvar.clnsig` | Pill: Pathogenic=빨강, Likely_pathogenic=주황, Uncertain_significance=회색, Benign=초록 | vep 없으면 `—` |
| **ACMG** | `vep.clinvar.clnsig`에서 파생 | Pill (P/LP/VUS/LB/B): 색상 ClinVar와 동일 | 매핑: `Pathogenic→P`, `Likely_pathogenic→LP`, `Uncertain_significance→VUS`, `Likely_benign→LB`, `Benign→B` |
| **Origin** | `origin` | Pill: `fetal_specific`=파랑, `maternal_het`=회색, `maternal`=회색, `ambiguous`=주황 | 표시 레이블: `fetal_specific→Fetal`, `maternal_het→Mat. Het`, `maternal→Maternal`, `ambiguous→Ambiguous` |
| **VAF** | `vaf` | `(vaf * 100).toFixed(1) + "%"` Bold | 오른쪽 정렬 |
| **Confidence** | `confidence` | Pill: `high`=초록, `medium`=주황, `low`=회색 | |
| **Action** | — | `Detail` 버튼 + `Coverage` 버튼 | 클릭 시 detail drawer 또는 coverage 페이지 이동 |

### Tier 2 — 기본 표시 (사용자 설정으로 숨김 가능)

| 컬럼명 | JSON 경로 | 표시 형식 | 비고 |
|--------|-----------|-----------|------|
| **Transcript (NM)** | `vep.mane` | monospace small text | 없으면 `—` |
| **Effect** | `vep.consequence` | `_variant` 접미사 제거, `_`→공백 | e.g. `missense variant` |
| **Impact** | `vep.impact` | Pill: `HIGH`=빨강, `MODERATE`=주황, `LOW`=회색 | Effect 컬럼 아래에 small Pill로 병기 권장 |
| **Allele Depth** | 계산값: `alt_count = Math.round(vaf * depth)`, `ref_count = depth - alt_count` | `REF(N)=ref_count / ALT(N)=alt_count` 형식 (gx-exome 동일) | depth 필드는 총 depth |
| **Zygosity** | `fetal_genotype` + `maternal_genotype` | F: `het`/`hom`/`at_least_het` / M: `hom_ref`/`het`/`hom_alt` | 2줄로 표시 |
| **gnomAD AF** | `vep.gnomad_af` | `null` → `Novel` (파란 pill), 숫자 → `(gnomad_af * 100).toFixed(3) + "%"` | 0에 가까울수록 임상적 의의 높음 |
| **Tags** | 파생 (아래 규칙 참고) | Pill 배열 | |
| **Disease** | `disease` | Small text | 없으면 `—` |

#### Tags 파생 규칙

```
origin === "fetal_specific"                     → "Fetal" 태그 (파랑)
vep.clinvar.clnsig === "Pathogenic"             → "Pathogenic" 태그 (빨강)
vep.clinvar.clnsig === "Likely_pathogenic"      → "Likely Path." 태그 (주황)
vep.clinvar.clnsig === "Uncertain_significance" → "VUS" 태그 (회색)
confidence === "low"                            → "Low conf." 태그 (회색)
in_zero_probe_region === true                   → "Low cov. region" 태그 (주황)
vep.gnomad_af === null AND vep.existing === ""  → "Novel" 태그 (파랑)
```

### Tier 3 — Detail expand 전용 (row 클릭 시 하단 확장)

| 컬럼/항목 | JSON 경로 | 비고 |
|----------|-----------|------|
| SIFT | `vep.sift` | `deleterious(0)` → `deleterious` + 점수 표시 |
| PolyPhen | `vep.polyphen` | `probably_damaging(0.998)` → 분리 표시 |
| dbSNP / COSV | `vep.existing` | `&` 구분자로 여러 ID 가능 |
| ClinVar HGVSc | `vep.clinvar.clnhgvs` | 전체 HGVSc (transcript 포함) |
| ClinVar Disease | `vep.clinvar.clndn` | `_` → 공백, `|` 구분자 존재 시 첫 번째만 |
| ClinVar Review Status | `vep.clinvar.clnrevstat` | 문자열 그대로 표시 |
| Genomic Position | `chrom:pos ref>alt` | e.g. `chr17:43092918 C>T` |
| Details (sgNIPT) | `details` | 파이프라인 분류 근거 텍스트 |
| Warning | `warning` | null 아닐 때만 노란 callout으로 표시 |

---

## 3. 행(Row) 하이라이트 색상

```
1순위 (빨강): origin === "fetal_specific" AND (clnsig === "Pathogenic" OR "Likely_pathogenic")
2순위 (파랑): origin === "fetal_specific" (pathogenic 아닌 fetal)
3순위 (주황): origin !== "fetal_specific" AND (clnsig === "Pathogenic" OR "Likely_pathogenic")
나머지: 기본 (하이라이트 없음)
```

---

## 4. 상단 Summary Stats

테이블 위에 다음 4개 stat 카드를 표시합니다:

```
Total variants       → variant_analysis.total_variants
Fetal-specific       → variant_analysis.fetal_variants
Pathogenic / LP      → variant_analysis.pathogenic_variants
Fetal + Pathogenic   → fetal_specific_detail 중 clnsig가 P/LP인 항목 수 (클라이언트 계산)
```

---

## 5. 필터 탭 (Table 상단)

```
All         → fetal_specific_detail 전체
Fetal       → fetal_specific_detail 중 origin === "fetal_specific"
P / LP      → pathogenic_details 배열 사용
Maternal    → fetal_specific_detail 중 origin이 "maternal" 또는 "maternal_het"
```

---

## 6. VEP 미주석 시 Fallback

`variant_analysis.vep_annotated === false` 이거나 특정 변이에 `vep` 필드가 없는 경우:

- **Gene**: `gene` 필드 사용
- **HGVSc**: `pathogenic_variant` 필드 사용 (있으면)
- **HGVSp**: `—`
- **ClinVar / ACMG / gnomAD AF / SIFT / PolyPhen**: `—`
- **Consequence / Impact**: `—`

---

## 7. 컬럼 순서 (좌 → 우)

```
Gene | HGVSc/HGVSp | Transcript | Effect/Impact | Origin | VAF | Allele depth | Zygosity | gnomAD AF | ClinVar | ACMG | Confidence | Tags | Disease | Action
```

sgNIPT 고유 컬럼(`Origin`, `VAF`, `Confidence`)에는 컬럼 헤더에 `[sgNIPT]` badge 또는 tooltip 표시 권장.

---

## 8. 기타 UX 요구사항

- **Sticky header**: 행이 많을 때 헤더 고정
- **Striped rows**: 짝수 행 연한 배경
- **Column 토글**: 사용자가 Tier 2 컬럼을 표시/숨김 설정 가능
- **Row expand**: 행 클릭 시 하단에 Tier 3 항목 표시
- **Coverage 버튼**: 해당 locus의 coverage 시각화 페이지로 이동 (기존 gx-exome Coverage 뷰 동일 방식)
- **Detail 버튼**: 해당 변이의 전체 VEP 주석 + pipeline 분류 근거(`details` 필드) 표시
- **Warning callout**: `warning` 필드가 null이 아닌 행에는 행 내 경고 아이콘 + tooltip 표시

---

## 9. 참고

- 기존 gx-exome Carrier test 테이블과 동일한 컬럼 헤더명 및 Pill 색상 체계 유지
- `vep.hgvsc`, `vep.hgvsp`의 transcript prefix(`NM_004992.4:`) 는 테이블 셀에서 제거하고 Detail expand에서 전체 표기 노출
- `gnomad_af`가 `null`인 경우 Novel로 처리 (Novel = gnomAD에 등재되지 않은 변이)
- ACMG 컬럼은 ClinVar `clnsig`에서 파생이므로 별도 API 호출 불필요
- `alt_count` = `Math.round(vaf * depth)`, `ref_count` = `depth - alt_count` 로 계산
