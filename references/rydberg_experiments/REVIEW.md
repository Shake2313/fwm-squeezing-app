# ⟦AI⟧ 문헌 DB 검수 기록

검수일 2026-11-03. 같은 구축 세션에서 원문과 생성물을 교차 확인. 독립 검토자 검수나 물리 모델 검증으로 표시하지 않음.

## 원문·과학적 구분

| 확인 대상 | 반영 |
|---|---|
| Gold 지정 | 직접 사용자 요청·날짜·당분간 적용 상태 보존. Ju 논문만 Gold |
| Gold 실험 조건 | PDF pp. 2–5, 8–10 대조. 50 mm·약 20 °C·6 µW·30 mW·직경 0.15 mm·37.07104 GHz·IF 40 kHz |
| Gold 직경·주파수·OD | 1/e²/FWHM·절대 probe detuning·실측 OD 미지정. AutoOD 대화 계산을 실측값으로 입력하지 않음 |
| Gold 초미세 표기 | 본문 2절 F=3→F′=4와 그림 1(a) F=2→F′=3의 불일치 보존. 메인 전이는 본문 채택 |
| Gold 보충자료 | 원문이 별도 보충문서를 언급. 현재 main article 15쪽 보관·검수, 별도 supplement 미확보 |
| Gold 감도·온도 | 실측 12.5(8), PSN 11.2, 계산 최적 33 °C 각각 별도 quantity·실험 기록 |
| Sedlacek 2012 | 출판 제목과 preprint 제목 구분. 실제 셀 10 cm와 유효 길이 7.5 cm 별도 |
| Gordon 2019 | 출판 PDF의 79 µV/m/√Hz·최소장 46±2 µV/m 보존 |
| Prajapati 2021 | Repump on/off와 잡음 예산의 probe 파워를 별도 운전점으로 저장. 입력 560 µW와 검출 220 µW 구분 |
| Cai 2022 | 직접 변조 감도 35.75와 잡음 스펙트럼 추정 12.50 nV/cm/√Hz 분리. RF −6 MHz를 probe detuning으로 사용하지 않음 |
| Dixon 2023 | 공개 preprint 초록의 EIT 선폭 <5 MHz, 감도 2.3 µV/cm/√Hz·최소장 66 µV/cm 확인 |
| Su 2022 | 온도 스캔 27–65 °C·직경 정의·EIT peak contrast와 RF 감도 구분 |
| Meyer·Yang 이론 | 공명 small-signal OD·속도 평균·phase diffusion의 역할 확인. 실험 감도 표본으로 집계하지 않음 |
| Tu 2024 | 냉각 원자·100 Hz 반복·2.7 ms 검출창. 연속 열증기와 가동률 구분 |
| Artusio-Glimpse 2024 | 광 comb span 200 MHz, RF 대역폭 12±1 MHz, 감도 9 mV/m/√Hz 각각 별도 |
| Ma 2026 | 최소장 2.8±0.3 mV/cm에 Hz 정규화 단위 추가하지 않음 |
| Viray 2026 | Table 1의 carrier별 감도 사용. 40.68 MHz의 별도 Fig. 7 추정 352.2±13.7과 표 333±6 구분 |

미보고 불확도·상태 정의·빔 크기·기준면은 임의 보완하지 않음. 보고 범위·초록 최대값·공동 스캔을 개별 측정점으로 늘리지 않음. 출처를 그대로 찾을 수 있도록 PDF 페이지와 절·그림·표 병기.

## 구조·열람 검증

`build_db.py`, `build_db.py --check`, `verify_db.py` 실행 결과는 `validation.json`·`validation_detail.json`에 저장. 검사 범위는 schema, ID·근거 연결, PDF 형식·해시·페이지 수, CSV 전 필드 왕복, 생성물 stale 감지, 잘못된 자료 거부, Gold 정책 일치.

브라우저의 검색·주제/연구유형 필터·근거 펼치기·로컬 원문 링크·Gold 우선 표시 확인 결과는 `viewer_validation.json`에 기록.

변경 범위는 DB·루트 README의 탐색 링크·추적용 ignore 규칙. GABES 물리 계산의 검증은 이번 자료 구축 범위 밖. 수치 엔진 회귀·GC 캠페인 실행 생략.
