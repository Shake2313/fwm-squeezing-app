# GABES 일일 리뷰 — 2026-11-04 Rydberg-EIT

## ⟦AI⟧ 선택·변경·판정

- 실제 시작 **2026-11-04 수요일 04:55:27 KST**. 설정상 예정 05:00보다 이른 전달. 실제 요일대로 `rydberg_eit`, `_SCHEMES` 등록 확인. `weekday-v1`; mod 5 미사용.
- 출발점: [10-28 같은 주제 보고서](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/daily_report/2026-10-28_scheme-3_rydberg-eit.md>)·[상태](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/daily_report/review_state.json>)·[공식 체크리스트](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json>). 10-19 취합 cutoff 10-18 유지.
- `main@676a7a61ee693dbde5ca4940c9500f24f965f8cb`. 이전 42개 source 중 **16개 변경**. 현재 관련 **56개 hash 종료까지 동일**. CP/새 광학축 구현은 기존 **미커밋 작업 트리**에 존재; HEAD만의 결과 아님. 사용자 편집 보존.
- **새 지적 1건, 기존 양성성 지적 해결 후보 1건. 리뷰 완료 / 검증 판정 failed.** 기존 시험 **107 passed**이나 새 배치 광학축 계약 probe 실패. scalar 모델의 조건부 실험 reference 가치는 유지; 비영점 detuning 배치 결과는 아래 수정 전 실험 판단에 사용 보류.

근거: [동반 evidence](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/daily_report/2026-11-04_scheme-3_rydberg-eit_evidence.json>)·[소스/환경](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/daily_report/2026-11-04_rydberg/snapshot.json>)·[git 한정 조회](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/daily_report/2026-11-04_rydberg/git.json>).

## ⟦AI⟧ 주요 관측·다음 행동

### 1. `review-2026-11-04-rydberg-batch-optical-axis` — new · 잠정 P1 · 체크리스트 미등록

새 probe/coupling 독립 스캔이 [본체 공통 Hamiltonian](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/rydberg.py:841>)에는 반영됐으나 [배치 adapter](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/analysis/rydberg_cell_heating/adapters.py:436>)는 기존 `diag(H)/2π = [0,−s,−s,−s−ΔMW]`를 재조립. 실제 식은 `[0,−Δp,−(Δp+Δc),−(Δp+Δc+ΔMW)]`. 양쪽 고정 detuning 0인 기본 fixture는 일치; 비영점 counterpart·coupling 스캔에서 다른 원자 응답을 올바른 static transmission과 조합함.

- **Probe 스캔, 고정 Δc=+12 MHz:** 중심 Δp=−12 MHz에서 configured atomic phasor `−1.49340e−12 + 2.05160e−14 i`, 본체와 같은 물리점 직접 응답 `6.14097e−11 − 9.57783e−11 i` (단위 `(rad/s)⁻¹`, 상대 복소오차 **100.73%**). RF/detector 환산 전 불일치; 그 차이로 감도 오차율을 직접 단정하지 않음.
- **Coupling 스캔, 고정 Δp=+12 MHz:** scheme 모드 배치 export가 최적 probe를 **−10.446 MHz**로 기록. 실제 probe는 **+12 MHz**, −10.446은 최적 coupling 좌표. [좌표 복사 위치](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/analysis/rydberg_cell_heating/adapters.py:367>).
- [Static 요약](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/analysis/rydberg_cell_heating/adapters.py:165>)도 resonance를 `argmin(abs(x))`로 찾음. 위 probe 스캔에서 `transmission_at_resonance` **0.787767** 대 본체 **0.980231**. 실제 두 광자 중심은 −12 MHz.

재현: [probe 소스](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/daily_report/2026-11-04_rydberg/integration_probe.py>)·[세 fixture 원시 결과](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/daily_report/2026-11-04_rydberg/integration_probe.json>)·[명시 계약 판정](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/daily_report/2026-11-04_rydberg/integration_assessment.json>). 5개 검색점만 사용한 identity 진단; 수렴된 감도 최적화 아님.

**다음/완료 조건:** 배치 static/configured/scheme 경로를 공통 Hamiltonian·resonance·물리 detuning export와 연결. 기본 parity 유지하며 양쪽 축·비영점 counterpart·필드 단위를 검증. 같은 원인 3증상을 1건으로 추적. 이번 리뷰에서 제품/시험/체크리스트 수정 안 함; 구현 승인은 별도.

### 2. `review-2026-10-22-rydberg_eit-selective-dephasing-positivity` — resolution_candidate

[현행 채널](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/rydberg.py:464>)은 `Lopt=√(2γg) diag(0,0,1,1)`, `LRF=√(2γRF) diag(0,0,0,1)`. 따라서 `γ02=γg`, `γ23=γRF`, `γ03=γg+γRF`. 독립 pair 감쇠와 **다른 물리 모델**. 새 [21개 모델 시험](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_rydberg_model_physics.py>)에서 conditional CP·유한 채널 양성/trace·이전 음의 정상상태 fixture·양쪽 optical axis/velocity·cache/export 통과. 과거 반례를 현행 결함으로 반복하지 않음. 실측 laser phase/transit rate와의 동일성은 별도 미검증. 공식 checklist 없는 기존 AI 관측의 해결 후보이며 자동 과제 폐쇄 아님.

### 3. `rydberg-readout-validity-guardrails` — reconfirmed · P1 needs_decision

현행 801점·20 °C·probe 6 µW·coupling 30 mW·직경 0.15 mm·셀 50 mm·opacity 0.001에서:

| 모드 | 현재 조건부 출력 | 범위 |
|---|---|---|
| EIT Doppler off / on | 선폭 **1.603494 / 2.878587 MHz** | optical-grid 수렴 인증 아님 |
| AT Doppler off / on, LO 3.7 MHz·IF 40 kHz | total=PSN **13.2413 / 2229.8020 nV/cm/√Hz** | RIN·electronics=0, QE 0.95·path 0.5·RMS dipole 규약 |
| AT Doppler on | static AT **미분해**, finite-IF 감도는 PREDICTED | static split과 finite-IF 유효성 별도 |

[원시 입력·출력·receiver 계측](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/daily_report/2026-11-04_rydberg/probe.json>). Doppler-off AT 기준광은 **2.519709–2.930825 µW**, 801점마다 signal과 정확히 balance. 고정 기준광 장치 검증 아님. `numerical_convergence`·`experimental_visibility` 모두 `NOT_ASSESSED`. 속도축 refinement PASS가 optical-grid/edge·측정 SNR을 대신하지 않음.

**다음:** 기존 acceptance대로 관측량별 수렴·edge 기준과 고정/단일점 balance reference 계약 검증. 임의 투과율 cutoff로 실험 불가 판정하거나 static AT 미분해만으로 finite-IF 감도를 일괄 기각하지 않음.

## ⟦AI⟧ 완료 유지·정정·잔여 항목

| ID | 이번 판정 |
|---|---|
| `rydberg-doppler-on-pole-residue` | **changed, done 유지**. 11-03 `e251527`에서 exact Maxwell 평균·finite-IF resolving grid 채택. 현행 affine/guard/fallback·velocity refinement 시험 PASS. 기존 완료를 미완료로 재지적 안 함. |
| `rydberg-effective-opacity-and-sweep-trust` | **reconfirmed, done**. 현재 CP 원자응답에서 opacity 1회 적용·zero response·시나리오·표시/export 시험 PASS. 옛 non-CP 원자 배열 불변 주장은 하지 않음. |
| `rydberg-extra-view-dependency-correctness` | **reconfirmed, done**. 입력 선언·cache·sweep 시험 통과. 잔여 계층 성능 과제와 구분. |
| `scheme-current-documentation-and-examples` | **changed, ready 유지**. 시험 docstring의 static-only 설명 이미 수정. Guide 그림도 정적 예시/finite-IF 검증 아님을 명시. CP·새 축·현행 source/config 예제 동기화는 별도 잔여 범위. |
| `rydberg-experimental-validation-dataset` | **reconfirmed, blocked_external**. [11-03 사용자 Gold 지정](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/references/rydberg_experiments/reference_policy.json>)과 [26편 문헌 DB](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/references/rydberg_experiments/README.md>)는 새 근거. 보유 측정 EIT/RF/PSD·receiver·온도 metadata를 갖춘 held-out 자료 완료와 다름. 예제 실측 입력 PENDING 유지. |
| `rydberg-csv-acquisition-order` | **reconfirmed, ready**. loader·numeric helper hash 동일; branch/방향 손실 구조 유지. 10-28 synthetic 숫자는 historical, 재실행 안 함. |
| `rydberg-layered-sweep-cache` / `rydberg-finite-if-affine-postprocess` | **reconfirmed, ready**. 기존 계층 재사용·조립/배열화 범위 유지. 새 benchmark·속도 후보 없음. |

중요 정정: exact-Maxwell 완료 기록의 **Doppler-on 1557.57 nV/cm/√Hz**, “Doppler-off bitwise”는 **CP 변경 전 소스 비교**에 한정. 현재 작업 트리는 `rydberg-eit-v8-cp-optical-detunings`; 오늘 **2229.80**을 해당 과거 수치와 동일시하지 않음. 이번 probe는 변경별 감도 기여를 분리한 실험이 아님. Gold 감도와 가까운 Doppler-off **13.24**도 zero technical noise·미보정 opacity·Rabi/transit 가정의 조건부 결과; 독립 실험 재현 근거로 승격하지 않음.

## ⟦AI⟧ 검증 범위·생략

- **Fresh 100 passed, 23.02 s, 2 deprecation warnings:** `python -m pytest -q tests/test_rydberg_model_physics.py tests/test_affine_doppler.py tests/test_rydberg_eit.py tests/test_rydberg_opacity.py --durations=8` (+ JUnit). 모듈별 21/18/37/24. [명령/실행시간](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/daily_report/2026-11-04_rydberg/checks.json>)·[로그](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/daily_report/2026-11-04_rydberg/pytest.txt>). 경고는 export의 `datetime.utcnow()`; 물리 실패로 분류 안 함.
- **Fresh 7 passed, 28.02 s:** adapter integration 전체 + workflow의 calibrated finite-IF 시험 1개. [명령](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/daily_report/2026-11-04_rydberg/integration_checks.json>)·[로그](<G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/daily_report/2026-11-04_rydberg/integration_pytest.txt>). **기존 시험은 새 비영점 광학축 probe 실패를 검출하지 못함.** 전체 판정 failed를 green으로 덮지 않음.
- 환경: 원본 G: 작업공간, Python 3.12.6·NumPy 2.2.0·SciPy 1.18.0·Numba 0.65.1, `MPLBACKEND=Agg`. threading/backend는 snapshot에 저장. 시간은 시험 수행시간이며 cold/warm 성능 benchmark 아님.
- **Reused:** source 동일한 CSV 방향 손실의 정적 구현 관측. **Historical:** 10-27 optical-grid/positivity/속도 수치, 10-28 receiver 숫자, exact-Maxwell commit의 92,532행·시간·모든 slider 검증. 날짜·소스 범위 유지; 오늘 broad PASS로 재사용 안 함.
- 전체 pytest·GC audit/campaign·sealed 연구·대형 sweep/fit·ARC 재생성·원격 issue/신규 문헌 검색·브라우저 검토 생략. 이번 세션 제품 변경 없음; 새 CP·광학축과 영향받는 배치 caller를 집중 확인. 전체 엔진/실험/주간 UI PASS 주장 없음. 한정 TODO/FIXME 검색 결과 없음은 원격 issue 부재 근거 아님.
- 재현 뒤 제품 수리는 수행하지 않음. **미해결 필수 확인:** 새 batch optical-axis 계약은 수정 후 실패 fixture·기존 parity 재검증 필요. 기존 관측량별 수렴·실측 검증 과제도 열림. 이번 리뷰의 링크/hash·상태 게시 외 추가 필수 검증 보류 없음.

다음 우선 행동: **배치 새 광학축 연결 수정**, 이후 기존 readout 수렴·receiver 계약. 미채택 AI 검토 제안; 이번 리뷰로 사용자 목표·우선순위·구현 승인 추가 안 함. 다음 실제 목요일은 `magneto`; 누락이 생겨도 배정 이동 없음.
