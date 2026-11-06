# ⟦AI⟧ 리드버그 뷰 3단계 검증 — 2026-11-05

직접 사용자 요청: 공용 지도·선택 체계로 EIT/AT 지도, 출력·온도 분석, 유한 IF 감도를 연결.

## 구현과 모델 범위

- EIT 기본 비교를 Probe power 지도와 대비·FWHM·입사 출력×대비(µW) 추세로 변경. 원래 ±5 MHz local-floor 특징 추출을 사용한다. Coupling power·Temperature·Residual width 비교, 선택 조건의 잔류 폭 보조 곡선과 가정 opacity 시나리오를 제공한다. 잔류 폭 0에서는 보조 곡선 중복을 피한다.
- AT 기본 RF Rabi 지도를 유지하고 기존 추출기의 투과 최대점·분리 간격·입력/추정장을 연결했다. 입력 RF와 광학 간격, RF off·미분해를 구별한다. RF detuning·probe/coupling power·temperature·IF frequency도 같은 공용 행 스캔으로 비교한다.
- 온도 입력의 실제 파라미터를 그대로 사용한다. Linked는 temp_c, Separated는 effective_temp_c 또는 cold_spot_temp_c이며 다른 온도와 heater_setpoint_c를 고정한다. 기존 추가 분석처럼 암묵적으로 Linked 또는 최소 LO를 강제하지 않는다. AT 감도는 현재 RF 조건을 사용하며 RF off에서는 빈 값이다.
- 클라이언트 Spectra/Sensitivity 전환은 같은 데이터와 선택을 사용한다. 기존 유한 IF Liouvillian 응답의 크기, current phasor 위상, PSN/총 잡음 등가 장 ASD를 연결했다. 정적 dT/dν 또는 기존 finite-difference proxy로 감도를 대체하지 않는다.
- 잡음 예산의 기존 per-point 후보에서 PSN 배열과 current phasor 배열을 노출했다. 새 감도 공식을 도입하지 않았다. 장 환산의 RMS/Peak와 전류/장 같은 진폭 규약의 응답도·ASD 불변성을 구별한다. 조건부·미보정 예측과 가정 opacity 범위를 화면·내보내기에 보존한다.
- 각 observable의 실제 optical x 배열을 지원한다. Doppler on IF의 static-slope 표본은 점으로만 표시하고 사이를 채우거나 연결하지 않는다. 위상은 ±180° 경계에서 선을 끊고 주기 색 범위를 사용한다. RF off·zero opacity·미계산을 0 감도로 채우지 않으며 전체 빈 지도에는 임의 숫자 색 범위를 표시하지 않는다.
- Inspect sampled optimum은 현재 조건의 최선 IF 표본으로 검사 커서만 이동한다. 조건 판독과 제어 변수 추세의 sampled optimum을 구별한다. 선택 조건은 빈 원, 계산한 비교 조건 중 최적점은 마름모로 표시한다. opacity 시나리오는 별도 곡선이며 신뢰구간이나 불완전 범위의 유한 상한으로 제시하지 않는다.
- 원자 솔버·물리 기본값·801점 optical 격자 무수정. 사용자 가이드 원본과 배포 HTML 갱신. 기존 추가 분석은 legacy 상세 화면에서 유지한다.

## 검증

- 수정 전 직렬화한 **10개 조건의 기존 raw/readout 필드 bit-identical**. EIT/AT 기본, probe 0.1/20 µW, IF 1/500 kHz, RF off, RF detuning, nonzero fixed detuning의 coupling scan, separated temperature를 포함한다. 추가된 필드는 비교에서 제외하고 기존 모든 필드·배열·metrics·noise/sensitivity·atomic response를 재귀 비교했다.
- 영향 범위 **21개 테스트 파일 267 passed (31.41 s)**. Rydberg 전기장/IF/opacity/상태·기존 EIT 및 실험/온도 helpers, 공용 지도/선택/제스처/좌표/캐시/내보내기, 조건·DEFAULTS·모드, Lambda와 polarization adapter를 포함한다. 기존 export의 datetime.utcnow 폐기 경고 2건은 범위 밖이다. 전체/GC 회귀 결과로 부르지 않는다.
- 마지막 빈 지도·잔류 폭 0 처리와 가이드 변경 후 관련 **9개 파일 80 passed (13.76 s)**. 실제 Doppler IF 표본·서로 다른 observable 좌표·위상 wrap·미계산/무한 값·서로 다른 최적점 마커·프로파일/opacity 검증을 포함한다. Impeccable source detector 검출 없음.
- 내장 브라우저: EIT 기본 7개 출력 조건, 잔류 폭/opacity 전환, Linked 온도 지도와 직접 입력 위치, AT Spectra/Sensitivity·최적점 검사, IF 주파수 지도→40에서100 kHz 입력 연동을 확인했다.
- IF Uniform samples=5에 정확한 현재 100 kHz 조건을 추가해 **6/6 complete, 3 cached**. 실제 JSON에서 6행의 801점 optical/IF 좌표, sensitivity profile·선택·완료 상태, 내부 raw/view 미포함을 확인했다.
- Chrome 390/320 px에서 가로 넘침 없음. 8개 SVG의 viewBox/표시 폭이 각각 356/286 px로 일치하며 응답·위상·log ASD 단일 열 표시를 확인했다. 임시 뷰포트 복구. 내장 브라우저 1280×720 완료 미리보기 보존.
- 최신 주 실행본 18582 유지, 작업용 18586 종료. 커밋·푸시 없음. 브라우저/서버 error 로그 없음. 최종 문서 검사 4 passed, 상대 링크·registry JSON·서버 health 검사 통과. 원래 배열 bytes/dtype와 유한 scalar float hex까지 기존 10조건 동일함을 재확인했다.

## 후속 범위

3단계 완료. Magneto 4단계와 세 뷰 통합 5단계는 대기한다. 새 RF/시간영역 물리, Doppler IF의 전체 detuning 최적화, 실험적 절대 감도 검증·opacity 보정은 이 UI 패치의 완료 범위가 아니다.
