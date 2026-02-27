# nipt-daemon 통합 가이드

이 문서는 기존 nipt-daemon에 Single Gene NIPT(SgNIPT) 파이프라인을 통합하는 방법을 설명합니다.

---

## 1. 아키텍처 개요

기존 sWGS NIPT와 Single Gene NIPT는 동일한 daemon 인프라를 공유합니다.

```
Order Platform (AWS)
    │
    ▼
nipt-daemon (FastAPI)
    ├── POST /analysis/order/{order_id}/submit
    │       ├── type == "sWGS"    → run_nipt_v3.sh (기존)
    │       └── type == "SgNIPT"  → run_sgnipt.sh  (신규)
    │
    ├── Queue Manager (asyncio.Semaphore)
    │       └── max_parallel 공유 또는 분리 가능
    │
    └── Result Monitor
            ├── sWGS:   {output}/{YYMM}/{order_id}/{order_id}.json + .output.tar
            └── SgNIPT: {output}/{YYMM}/{order_id}/{order_id}.json + .output.tar
```

---

## 2. models.py 변경사항

### 2.1 SubmitOrderDto에 type 필드 활용

기존 `SubmitOrderDto.type` 필드를 사용하여 분석 유형을 구분합니다:

```python
# app/models.py - 기존 SubmitOrderDto에 추가 사항 없음
# type 필드가 이미 존재: "sWGS" 또는 "SgNIPT"

class SubmitOrderDto(BaseModel):
    patientBirthDate: str
    sequencingDataMethod: str  # "REMOTE", "LOCAL"
    labIdentifier: List[str]
    type: str                  # "sWGS" | "SgNIPT" ← 이 필드로 분기
    sampleBarcode: Optional[str] = None
    
    # SgNIPT 전용 필드 (Optional)
    targetPanel: Optional[str] = None      # Target BED 파일명
    variantCaller: Optional[str] = None    # bcftools, freebayes, gatk
```

### 2.2 FullOrder 확장

```python
class FullOrder(BaseModel):
    orderId: str
    clientId: str
    sequencingDataMethod: str
    lab: str
    age: int
    r1_id: str
    r2_id: str
    r1_path: str
    r2_path: str
    # SgNIPT 전용 필드
    analysis_type: str = "sWGS"          # "sWGS" | "SgNIPT"
    target_panel: Optional[str] = None
    variant_caller: Optional[str] = "bcftools"
```

---

## 3. runner.py 변경사항

### 3.1 파이프라인 실행 분기

`_run_pipeline` 메서드에서 `analysis_type`에 따라 분기합니다:

```python
async def _run_pipeline(self, full_order: FullOrder, r1_local: str, r2_local: str, 
                        detached: bool = True):
    """파이프라인 실행 - sWGS / SgNIPT 분기"""
    
    if full_order.analysis_type == "SgNIPT":
        await self._run_sgnipt_pipeline(full_order, r1_local, r2_local, detached)
    else:
        await self._run_swgs_pipeline(full_order, r1_local, r2_local, detached)


async def _run_sgnipt_pipeline(self, full_order: FullOrder, r1_local: str, r2_local: str,
                                detached: bool = True):
    """Single Gene NIPT 파이프라인 실행"""
    order_id = full_order.orderId
    work_dir = extract_work_dir(order_id)
    
    # SgNIPT root directory (별도 설정 가능)
    sgnipt_root = getattr(settings, 'sgnipt_root_dir', '/data/SgNIPT')
    
    # run_sgnipt.sh 호출
    script = Path(sgnipt_root) / "docker" / "run_sgnipt.sh"
    
    cmd = [
        "/usr/bin/bash", str(script),
        "--order_id", order_id,
        "--sample_name", order_id,  # 단일 샘플인 경우
        "--fastq_r1", os.path.basename(r1_local),
        "--fastq_r2", os.path.basename(r2_local),
    ]
    
    # Optional: target panel
    if full_order.target_panel:
        cmd.extend(["--target_bed", full_order.target_panel])
    
    # Optional: variant caller
    if full_order.variant_caller:
        cmd.extend(["--variant_caller", full_order.variant_caller])
    
    env = os.environ.copy()
    env["SGNIPT_ROOT_DIR"] = sgnipt_root
    env["SGNIPT_WORK_DIR"] = work_dir
    
    logger.info("Starting SgNIPT pipeline for %s", order_id)
    
    if detached:
        run_log = Path(sgnipt_root) / "log" / f"{order_id}.startup.log"
        run_log.parent.mkdir(parents=True, exist_ok=True)
        
        with open(run_log, "ab", buffering=0) as lf:
            proc = await asyncio.create_subprocess_exec(
                *cmd,
                stdout=lf, stderr=lf,
                env=env,
            )
            rc = await proc.wait()
        
        if rc != 0:
            raise RuntimeError(f"SgNIPT pipeline failed to start (rc={rc})")
        
        logger.info("SgNIPT pipeline started for %s", order_id)
    else:
        proc = await asyncio.create_subprocess_exec(
            *cmd,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE,
            env=env,
        )
        out, err = await proc.communicate()
        if proc.returncode != 0:
            raise RuntimeError(f"SgNIPT pipeline failed (rc={proc.returncode})")
```

### 3.2 모니터링 경로 분기

```python
async def _monitor_completion(self, order_id: str, full_order: FullOrder):
    """분석 완료 모니터링 - sWGS / SgNIPT 공통"""
    work_dir = extract_work_dir(order_id)
    
    if full_order.analysis_type == "SgNIPT":
        sgnipt_root = getattr(settings, 'sgnipt_root_dir', '/data/SgNIPT')
        json_file = f"{sgnipt_root}/output/{work_dir}/{order_id}/{order_id}.json"
        tar_file = f"{sgnipt_root}/output/{work_dir}/{order_id}/{order_id}.output.tar"
    else:
        json_file = f"/home/ken/ken-nipt/output/{work_dir}/{order_id}/{order_id}.json"
        tar_file = f"/home/ken/ken-nipt/output/{work_dir}/{order_id}/{order_id}.output.tar"
    
    logger.info("Monitoring %s (%s), waiting for: %s", 
                order_id, full_order.analysis_type, json_file)
    
    # 이하 기존 폴링 로직 동일 ...
```

### 3.3 컨테이너 이름 규칙

기존 sWGS: 컨테이너 이름 = `{order_id}`
SgNIPT:    컨테이너 이름 = `{order_id}` (run_sgnipt.sh에서 동일하게 설정)

```python
async def _check_container_failed(self, order_id: str, json_file: str) -> bool:
    """컨테이너 실패 여부 확인 - 동일 로직 사용 가능"""
    proc = await asyncio.create_subprocess_exec(
        "docker", "ps", "-q", "-f", f"name={order_id}",
        stdout=asyncio.subprocess.PIPE
    )
    out, _ = await proc.communicate()
    container_running = bool(out.decode().strip())
    return not container_running and not os.path.exists(json_file)
```

---

## 4. config.py 변경사항

```python
class Settings(BaseSettings):
    # ... 기존 설정 ...
    
    # === SgNIPT 전용 설정 ===
    sgnipt_root_dir: str = Field(
        default="/data/SgNIPT",
        env="SGNIPT_ROOT_DIR",
        description="Single Gene NIPT 루트 디렉토리"
    )
    
    sgnipt_docker_image: str = Field(
        default="sgnipt:latest",
        env="SGNIPT_DOCKER_IMAGE",
        description="SgNIPT Docker 이미지"
    )
    
    sgnipt_max_parallel: int = Field(
        default=2,
        env="SGNIPT_MAX_PARALLEL",
        description="SgNIPT 최대 동시 분석 수 (sWGS와 별도 관리 시)",
        gt=0, le=10
    )
```

---

## 5. main.py 엔드포인트 변경

### 5.1 submit 엔드포인트 분기

```python
@app.post("/analysis/order/{order_id}/submit")
async def submit_order(order_id: str, dto: SubmitOrderDto):
    """주문 분석 제출 - sWGS / SgNIPT 자동 분기"""
    
    # 기존 로직 동일 ...
    od = OrderDetailSubmit(
        id=order_id,
        clientId="...",
        # ...
    )
    
    # type에 따라 큐에 추가 (큐 매니저가 내부적으로 분기)
    await queue_manager.enqueue(order_id, od, dto)
    
    return {"status": "queued", "order_id": order_id, "type": dto.type}
```

### 5.2 진행 상황 조회 분기

```python
@app.get("/status/{order_id}/progress")
async def get_order_progress(order_id: str):
    """샘플 진행 상황 조회 - sWGS / SgNIPT 자동 감지"""
    work_dir = extract_work_dir(order_id)
    
    # SgNIPT progress 파일 먼저 확인
    sgnipt_progress = f"{settings.sgnipt_root_dir}/output/{work_dir}/{order_id}/{order_id}_progress.txt"
    swgs_progress = f"/home/ken/ken-nipt/output/{work_dir}/{order_id}_progress.txt"
    
    progress_file = sgnipt_progress if os.path.exists(sgnipt_progress) else swgs_progress
    # ... 이하 동일
```

---

## 6. 리포트 생성 통합

### 6.1 daemon의 report_functions.py에서 SgNIPT 리포트 생성

```python
# app/report_functions.py

async def generate_sgnipt_report(order_id: str, order_detail, template_dir: str):
    """SgNIPT 분석 결과를 기반으로 PDF 리포트 생성"""
    work_dir = extract_work_dir(order_id)
    
    # 1. 파이프라인 결과 JSON 로드
    result_json_path = f"{settings.sgnipt_root_dir}/output/{work_dir}/{order_id}/{order_id}.json"
    with open(result_json_path) as f:
        pipeline_result = json.load(f)
    
    # 2. 환자 정보 + 파이프라인 결과 병합 (기존 make_report_json 패턴)
    report_json = make_sgnipt_report_json(order_detail, pipeline_result)
    
    # 3. PPTX 템플릿 기반 PDF 생성
    from bin.generate_report import generate_from_json
    result = generate_from_json(
        report_json=report_json,
        template_dir=os.path.join(template_dir, "sgnipt"),  # SgNIPT 전용 템플릿
        output_dir=f"{settings.sgnipt_root_dir}/output/{work_dir}/{order_id}",
    )
    
    return result
```

### 6.2 SgNIPT 전용 PPTX 템플릿

`data/templates/sgnipt/` 디렉토리에 다음 템플릿 파일을 배치합니다:

```
data/templates/sgnipt/
├── SgNIPT_Basic_EN.pptx      # 영문 기본 리포트
├── SgNIPT_Basic_KR.pptx      # 한글 기본 리포트
├── SgNIPT_Detailed_EN.pptx   # 영문 상세 리포트
├── SgNIPT_Detailed_KR.pptx   # 한글 상세 리포트
├── No_Call_EN.pptx            # No Call 리포트
└── No_Call_KR.pptx            # No Call 리포트
```

템플릿 내 placeholder 목록:

| Placeholder | 설명 | 예시 값 |
|---|---|---|
| `<<Order ID>>` | 주문 ID | ORD-2602-001 |
| `<<Patient Name>>` | 환자명 | Jane Doe |
| `<<Maternal Age>>` | 산모 나이 | 35 |
| `<<Gestational Age>>` | 임신 주수 | 12w3d |
| `<<FF>>` | Fetal Fraction | 12.5% |
| `<<FF Method>>` | FF 추정 방법 | SNP-based |
| `<<Informative SNPs>>` | Informative SNP 수 | 45 |
| `<<Result>>` | 최종 결과 | NEGATIVE |
| `<<Total Loci>>` | 검사 대상 loci 수 | 250 |
| `<<Covered Loci>>` | 커버된 loci 수 | 248 |
| `<<Pathogenic Count>>` | 병원성 변이 수 | 0 |
| `<<Var1 Gene>>` ~ `<<Var10 Gene>>` | 변이 유전자명 | CFTR |
| `<<Var1 Classification>>` ~ | 변이 분류 | PATHOGENIC |
| `<<Report Date>>` | 리포트 생성일 | 2026-02-22 |

---

## 7. docker-compose 통합

### 7.1 nipt-daemon docker-compose에 SgNIPT 볼륨 추가

```yaml
# docker-compose.yml (production)
services:
  nipt-daemon:
    image: nipt-daemon:latest
    volumes:
      # 기존 sWGS 마운트
      - /home/ken/ken-nipt:/home/ken/ken-nipt
      - /var/run/docker.sock:/var/run/docker.sock
      # SgNIPT 마운트 추가
      - /data/SgNIPT:/data/SgNIPT
      # SgNIPT docker run 스크립트
      - /data/SgNIPT/docker/run_sgnipt.sh:/opt/sgnipt/run_sgnipt.sh:ro
    environment:
      # 기존 환경변수 ...
      # SgNIPT 추가
      - SGNIPT_ROOT_DIR=/data/SgNIPT
      - SGNIPT_DOCKER_IMAGE=sgnipt:latest
      - SGNIPT_MAX_PARALLEL=2
```

---

## 8. 배포 체크리스트

1. **SgNIPT Docker 이미지 빌드**
   ```bash
   cd /path/to/single-gene-nipt
   ./docker/build.sh v1.0.0
   ```

2. **Host 디렉토리 구조 생성**
   ```bash
   ./docker/setup_host_dirs.sh /data/SgNIPT
   ```

3. **Reference 데이터 배치**
   ```
   /data/SgNIPT/data/
   ├── reference/hg38.fa (+ .fai, .dict)
   ├── bwa_index/hg38.*
   ├── target_bed/sgnipt_panel_v1.bed
   ├── annotation/clinvar.vcf.gz (optional)
   └── templates/sgnipt/*.pptx
   ```

4. **nipt-daemon 코드 업데이트**
   - `models.py`: FullOrder에 `analysis_type` 필드 추가
   - `runner.py`: `_run_sgnipt_pipeline` 메서드 추가
   - `config.py`: SgNIPT 설정 추가
   - `docker-compose.yml`: SgNIPT 볼륨 마운트 추가

5. **nipt-daemon 재시작**
   ```bash
   docker compose down && docker compose up -d
   ```

6. **테스트 주문 제출**
   ```bash
   curl -X POST http://localhost:8000/analysis/order/TEST_001/submit \
     -H "Content-Type: application/json" \
     -d '{
       "patientBirthDate": "1991-01-01",
       "sequencingDataMethod": "LOCAL",
       "labIdentifier": ["lab1"],
       "type": "SgNIPT"
     }'
   ```

---

## 9. BAM 정리 정책

SgNIPT도 sWGS와 동일한 BAM 정리 정책을 적용합니다:

```python
def _remove_sgnipt_bams(self, order_id: str, analysis_dir: str):
    """SgNIPT BAM 파일 정리"""
    patterns = [
        f"{order_id}.sorted.bam", f"{order_id}.sorted.bam.bai",
        f"{order_id}.dedup.bam",  f"{order_id}.dedup.bam.bai",
        f"{order_id}.target.bam", f"{order_id}.target.bam.bai",
    ]
    # ... 기존 _remove_bams_internal과 동일한 로직
```

---

## 10. 모니터링 대시보드

기존 `dashboard_compact.html`에 SgNIPT 상태를 추가합니다:

```javascript
// /queue/summary 응답에 SgNIPT 통계 포함
{
    "swgs": {
        "running": "2/3",
        "completed_today": 15,
        "failed_today": 0
    },
    "sgnipt": {
        "running": "1/2",
        "completed_today": 8,
        "failed_today": 0
    }
}
```
