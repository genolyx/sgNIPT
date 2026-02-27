# nipt-daemon 소스 분석 요약

## 1. 아키텍처 개요
- **FastAPI 기반 비동기 웹서버** (uvicorn)
- **docker-compose**로 daemon 자체가 컨테이너로 실행됨
- daemon 컨테이너 내에서 **docker.sock 마운트**하여 호스트의 Docker를 제어
- prod/dev 환경 분리: docker-compose.yml (prod, port 8000) / docker-compose-dev.yml (dev, port 8001)

## 2. 핵심 모듈 구조
- `main.py` - FastAPI 앱, API 엔드포인트 정의
- `queue_manager.py` - asyncio.Semaphore 기반 동시 실행 제어
- `runner.py` - AnalysisTracker: Docker 실행, 모니터링, 결과 처리
- `config.py` - pydantic-settings 기반 설정 (환경변수 + .env 파일)
- `models.py` - Pydantic 모델 (FullOrder, OrderDetailSubmit, SubmitOrderDto 등)
- `notifier.py` - AWS API 통신 (결과 알림, 파일 업로드)
- `aws_client.py` - AWS API 클라이언트 (주문 조회, FASTQ 다운로드)
- `report_functions.py` - PPTX 템플릿 → PDF 리포트 생성
- `auth_client.py` - 인증 클라이언트

## 3. Order 처리 흐름
1. `POST /analysis/order/{order_id}/submit` → 큐에 추가
2. QueueManager.worker() → Semaphore로 max_parallel 제한
3. AnalysisTracker._run_analysis_workflow():
   - FASTQ 준비 (로컬 확인 또는 S3 다운로드)
   - run_nipt_v3.sh 실행 (bash 스크립트, detached 모드)
4. _monitor_completion(): 30초 폴링으로 JSON+TAR 파일 감시
5. _process_results(): TAR 업로드 → AWS 알림 → 컨테이너 정리 → BAM 삭제

## 4. 리포트 생성 흐름
1. `POST /analysis/order/{order_id}/report` → review_json 수신
2. make_report_json(): AWS에서 주문 상세 조회 + review 데이터 결합
3. generate_from_json(): 
   - 클라이언트별 PPTX 템플릿 디렉토리 탐색
   - 언어별(EN, ID, CN 등) 템플릿 파일 선택
   - placeholder 치환 (<<FF>>, <<Gender>>, <<T21 Result>> 등)
   - PPTX → PDF 변환 (LibreOffice)
   - 다국어 PDF 병합
4. upload_pdf_report(): AWS에 업로드

## 5. 디렉토리 구조 패턴
- work_dir = order_id에서 YYMM 추출 (예: GNCI25060001 → 2506)
- fastq: /home/ken/ken-nipt/fastq/{work_dir}/{order_id}/
- analysis: /home/ken/ken-nipt/analysis/{work_dir}/{order_id}/
- output: /home/ken/ken-nipt/output/{work_dir}/{order_id}/
- templates: /home/ken/nipt-daemon/data/templates/{ClientName}/
- signature: /home/ken/nipt-daemon/data/signature/

## 6. 설정 키
- MAX_PARALLEL: 동시 분석 수 (기본 2, 최대 10)
- WORKER_TIMEOUT: 2시간
- REMOVE_BAM: BAM 파일 자동 삭제 여부
- APP_ENV: dev/prod 환경 구분

## 7. Single Gene NIPT 파이프라인에 적용할 사항
- **동일한 daemon 구조 활용 가능**: submit → queue → run → monitor → notify
- **docker-compose 패턴 동일**: daemon이 docker.sock으로 분석 컨테이너 제어
- **리포트**: PPTX 템플릿 기반 placeholder 치환 → PDF 변환 패턴 재사용
- **디렉토리 매핑**: work_dir(YYMM) 패턴 동일 적용
- **API 엔드포인트**: /submit, /report, /sign-report 패턴 재사용
- **모니터링**: JSON+TAR 완료 파일 폴링 패턴 재사용
- **차이점**: 
  - sWGS는 bash 스크립트(run_nipt_v3.sh) 실행, SgNIPT는 Nextflow 실행
  - SgNIPT는 target BED 파일 필요
  - SgNIPT 결과 JSON 스키마가 다름 (variant 정보, fetal fraction 방법 등)
  - SgNIPT 리포트 템플릿 별도 필요
