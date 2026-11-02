# ⟦AI⟧ PSR DB 검수 기록 · 2026-10-30

## 독립 원문 검수

수집 담당자와 다른 검토자가 27편을 나누어 검수. 원문·수치 정의·운전점 연결·판본·모델 범위 대조. 전체 figure digitization이나 실험 재현은 수행하지 않음.

| 대상 | 원문 검수 기록 | 결과 |
|---|---|---|
| 실험 9편 | [review_scientific_experiments.json](review_scientific_experiments.json) | PASS |
| 공간 모드·전파·센싱 9편 | [review_scientific_modes.json](review_scientific_modes.json) | 지적 수정 후 PASS |
| 기초 이론·초기 실험 5편 | [review_scientific_foundations.json](review_scientific_foundations.json) | 표기·출처 수정 후 PASS |
| 잡음·전파·배열 문맥 4편 | [review_scientific_context.json](review_scientific_context.json) | 표기·출처 수정 후 PASS |

해결 사항:

- Zhang 2013: Gaussian 입력의 FWHM을 vortex 모드 크기로 복사한 항목 제거. vortex 빔폭은 미보고로 보존.
- Valente 2015: 약 40 mW를 `approximate`로 수정. 같은 William & Mary hot-vapor 장치 그룹 표기 통일.
- Matsko 2002: 원문의 `s = sqrt(Vvac/Vmin)`와 `10 log10(s)` 표기 주의 기록. 그 이론 dB 값을 통상적인 측정 noise-power dB로 직접 비교하지 않음.
- Lezama·Josse·Mikhailov·Wang: 실제 보관 판본에 맞춰 절·페이지·그림 locator 수정. Novikova의 enriched isotope 표기 보완.

원문 자체의 불일치도 보존:

- Agha 2010의 anti-squeezing 본문값과 figure caption 값은 차이 표시. 하나를 조용히 선택하지 않음.
- Lima 2024는 직교 seed가 있는 실험. 보고한 검출기 보정값과 단순 vacuum-loss 환산의 불일치는 AI note로 유지; 논문값을 임의 재작성하지 않음.
- Monsa 2026은 최종 출판본 PDF 기준. 이전 preprint의 더 낮은 squeezing 및 내부 표기 차이는 별도 출처 기록.

## 데이터·열람본 검증

[validation.json](validation.json): JSON Schema, ID, 원문 연결·SHA-256·PDF 페이지 수, 구간·미보고·환산식 검사.

[validation_detail.json](validation_detail.json): CSV 전체 필드 round-trip, 통합 JSON 무손실, 생성물·manifest 일치, 읽기 전용 `--check`, 잘못된 레코드 거부, abstract-only 구분 확인. 검증 성공은 구조·추적성 범위이며 과학적 예측 검증 아님.

[viewer_validation.json](viewer_validation.json): 실제 브라우저의 검색·연구 유형·주제·포함 범위 필터, 근거 펼치기, 로컬 출처 링크 및 화면 점검 기록.

레거시 FWM CSV 이동은 별도 [이동 검증](../../archive/fwm_squeezing_csv_migration_2026-10-30.validation.json) 참조. 원본 바이트 동일, 현재 reader·링크 수정, FWM DB 10편 legacy 보존 확인. 과거 시점의 source manifest·proofread snapshot·기존 PDF는 당시 경로를 유지하고 복구 manifest로 대응.

검증 범위 밖: GABES 예측값·미시적 noise 모델 구현·보정 학습·GC 캠페인. 이번 작업은 문헌 DB와 파일 경로 이전.
