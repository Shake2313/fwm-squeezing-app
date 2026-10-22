# GABES 일일 Scheme 물리 검토 — 2026-09-22

## 선택·범위

Asia/Seoul 현지 일자는 22일. n = (22 mod 5) + 1 = **3**. 현재 드롭다운 순서는 gabes/schemes/__init__.py:19-25 기준 (1) OD/SAS, (2) Lambda coherence EIT/AT/CPT, (3) **Rydberg-EIT electrometry**, (4) Hanle/EIA/NMOR, (5) FWM. 같은 파일 27-39행의 Lambda/Magneto alias는 드롭다운 항목이 아니다.

현재 HEAD는 ae456ab. Rydberg scheme, electrometry·experiment·CSV helper, 해당 테스트와 cell-heating 분석 workflow는 HEAD 대비 작업 트리 변경이 없다. CLAUDE.md, README.md, docs/checklist.json, 기존 Scheme 3 일일 보고서, 소스의 TODO/FIXME, 로컬 issue-note 파일명, 테스트, 예제 설정을 먼저 검색했다. Rydberg 전용 TODO/FIXME 및 로컬 issue-note 파일은 발견하지 못했다. **개선안은 이미** docs/checklist.json:2936-3018,3357-3371,3540-3619 및 docs/daily_report/2026-09-17_scheme-3_rydberg-eit.md에 있다. 이 검토는 원격 이슈 현황을 새로 단정하지 않는다.

## 연구자 판정

**실제 물리 있음. 반정량적 실험 설계·교육 참조로 유용. 절대 감도나 가열 최적점의 독립 기준으로는 미검증.**

- ⁸⁵Rb 5S F=3 → 5P F′=4 → 40D → 39F 네 준위, 780/481 nm 광학 cascade와 약 37 GHz RF leg를 쓴다(gabes/schemes/rydberg.py:1-15,72-91,417-437). 스펙트럼은 detuning별 정상상태 Liouvillian을 푸는 4준위 OBE이며, 단순 피크 그림은 아니다(507-527,558-627).
- Probe·coupling Rabi는 광파워의 제곱근과 빔 지름 역수로 변하지만 기준 운전점에 맞춘 **유효 구동 anchor**다. 특히 probe의 bare D2 dipole 계산은 생략된 optical pumping/Zeeman 구조 때문에 맞지 않아 별도 조정됐다고 코드가 밝힌다(41-55,439-461). 1.6 MHz EIT 폭과의 일치는 독립 예측 검증이 아니다.
- Cold spot은 포화 증기압을, beam 위치 유효 온도는 국소 밀도·속도·transit/dephasing을 결정한다(463-505,558-590). 실제 가열 셀의 구분에 맞는다. Doppler-on은 5P의 probe-k와 Rydberg level의 잔여 2광자 k를 velocity class별로 적용한다(72-91,610-619). 기본 Doppler-off는 이 운동 효과를 모두 끈 보정 기준이다. 두 모드는 정확도만 다른 동일 calibration이 아니다.
- LO가 dressing한 정상상태 주위에서 유한 IF의 약한 SIG에 대한 Liouvillian 선형 응답을 풀고 Hermitian conjugate sideband를 복원한다(gabes/rydberg_electrometry.py:116-194). Doppler-on phasor는 복소수로 평균한다(gabes/schemes/rydberg.py:787-806). RF dipole/편광 인자·RMS/peak·balanced photodiode PSN/RIN/electronics를 명시해 field ASD로 변환한다(gabes/rydberg_electrometry.py:231-316,403-483). 문헌 감도 숫자를 결과에 주입하지 않는 것은 장점이다.
- 실험 파이프라인은 온도 상태, axial Beer–Lambert, 유효 원자 수, SAM 불확도, EIT/RF/PSD CSV 단위 및 SHA-256 provenance를 제공한다(gabes/rydberg_experiment.py:80-95,191-338,378-448,485-594; gabes/rydberg_experimental_csv.py:441-691). analysis/rydberg_cell_heating/README.md:22-39의 MEASURED/REFERENCE/PREDICTED/ASSUMED/PENDING 구분도 좋은 참조 습관이다.

### 오늘의 재계산

같은 worktree의 기본값으로 직접 계산했다. 첫 EIT 호출에는 Numba 초기화가 포함됐으므로 그 wall time을 steady-state 성능으로 쓰지 않는다.

| 조건 | 결과 | 뜻 |
|---|---:|---|
| 기본 EIT T(0), 대비, FWHM | 0.940016615, 0.135391492, 1.613058583 MHz | 기준 폭 1.6 MHz와 규모 일치 |
| 기본 AT T(0), 피크 간격 | 0.806029283, 3.496500 MHz | LO 3.7 MHz와 같은 규모 |
| 기본 유한-IF 최적 detuning, 총 ASD | −2.220 MHz, 6.082387 nV/cm/√Hz | 조건부 계산값 |
| 같은 AT 원자 solve, 유효 opacity 0.0005 / 0.001 / 0.002 | 11.781990 / 6.082387 / 3.238776 nV/cm/√Hz | 숨은 scalar의 강한 영향 |
| 80 °C EIT T(0), 대비, 출력 FWHM | 1.798923×10⁻⁶, 1.798923×10⁻⁶, 0.307015 MHz | 숫자 폭이 곧 관측 가능 폭은 아님 |
| 기본 AT, Doppler-on 총 ASD | 1007.281 nV/cm/√Hz | Doppler-off와 calibration regime 격차 큼 |

문헌 맥락 파일 references/rydberg_ju_experiment_context.md:20-22,94-100,140-155는 PSN 11.2(4), 실측 12.5(8) nV/cm/√Hz 및 문헌 온도 모델 최적점 약 33 °C를 기록한다. 기본 GABES 6.08은 이 값보다 낙관적이다. 직전 일일 보고서는 내장 가열 sweep 최적점 47 °C를 기록했다(docs/daily_report/2026-09-17_scheme-3_rydberg-eit.md:170-198). 이번 실행에서는 sweep 자체를 재측정하지 않았다. 문헌의 33 °C와 차이가 있으므로 온도 sweep은 운전점 권고가 아닌 **PREDICTED 시나리오**다.

## 현재 참조 한계

1. gabes/schemes/rydberg.py:643의 ls=0.001은 공개된 ⁸⁵Rb line strength가 아니라 여러 미모델 효과를 흡수한 유효 opacity다. 672-681행은 이를 Beer–Lambert 전달에 넣는다. 위 민감도처럼 절대 field ASD를 강하게 좌우하므로 provenance·불확도 없이는 절대값을 인용할 수 없다.
2. EIT FWHM은 대비가 양수이면 산출된다(979-1003). 80 °C 반례는 낮은 광량/대비에도 정밀한 폭이 남음을 보인다. 미해결 AT 조건에서도 Total sensitivity가 hero로 남는다(1127-1207; tests/test_rydberg_eit.py:134-141). 분해능, scan 가장자리, 검출 noise 기준이 함께 필요하다.
3. 감도 scan의 reference arm 광파워는 매 detuning의 signal 광파워에 맞춰 다시 설정된다(902-927). 고정 reference arm 또는 한 점에서 balance한 실험과 다르다. 이 scan은 이상적 점별 balance envelope로 읽어야 한다.
4. Doppler-on finite-IF는 full detuning 탐색 대신 최대 static slope 한 점을 선택한다(787-819). 기본 off/on 감도 격차는 현재 유효 opacity를 두 mode에 그대로 재사용할 수 없다는 경고다. Zeeman/편광 manifold, RF standing wave, 비선형 time-domain lock-in도 포함되지 않는다(1-15).
5. 예제는 workflow와 provenance의 시범 자료다. analysis/rydberg_cell_heating/example_config.json:8-27에 문헌 운전 조건과 가정한 detector/경로 수치가 섞여 있고, 측정 EIT/RF/온도 입력은 PENDING이다(53-78). RF dipole만 REFERENCE, electrometry block은 ASSUMED다(88-102). held-out 측정 검증을 대체하지 못한다.

## 기존 개선안의 계산 비용과 물리 보존

| 기존 제안·위치 | 비용 및 보존 조건 | 판단 |
|---|---|---|
| 유효 opacity 노출, sweep trust — docs/checklist.json:2955-2976 | 숨은 0.001을 이름·provenance·bounds가 있는 입력으로 옮기면 기본 solve 비용 **0**. 기본값을 고정하면 기존 spectrum·감도 보존. 불확도 sweep만 선택 시 추가 계산. | 최우선 P1. optical line strength와 이 participation을 중복 곱하지 않음. |
| 관측 가능성 gate와 고정 reference arm — 2978-3001 | 대비, 검출 광파워, 폭당 샘플 수, edge margin은 이미 얻은 801점 배열에 O(N) 후처리. 고정 arm은 원자 solve 그대로이고 검출기 모델만 현실화. 국소 refinement는 몇 점 추가 solve 비용. | P1. noise 기준이 없으면 “관측성 미평가”로 표시. |
| extra-view 의존성·계층 cache — 2936-2953 | cache key는 O(1). 원자 spectrum / IF phasor / Beer–Lambert / detector-SAM을 분리하면 detector knob 변경 때 원자 solve 재사용. 결과 동일성 gate 필요. | P1. 현재 ExtraView key·phasor LRU는 일부 구현. |
| 측정된 readout transfer function — 3567-3591 | 측정 complex H(IF)를 phasor에 곱하는 비용은 O(N_IF), 사실상 작음. 기존 선형 원자 응답 보존. | 측정 gain·phase 자료 확보 뒤 시행. full time-domain 전 단계. |
| Doppler-on pole–residue — 3357-3371 | 현재 2 m/s, 4σ **이산 Maxwell 측도**와 coherent phasor 규칙을 유지하면 동작 보존 가속. 연속 Gaussian 평균으로 바꾸면 결과도 바뀌므로 별도 물리 검토. | off가 기본이므로 P3 타당. 속도 개선은 opacity 유효성 문제를 해결하지 않음. |
| held-out 검증·Zeeman/편광·RF map — 3003-3018,3540-3564,3593-3619 | 측정 데이터와 상태/편광/공간 grid가 필요. 원자 차원과 공간 점 수가 늘어 무시 가능한 overhead가 아님. | 외부 데이터 전에는 opt-in 연구 경로로 보류. |

## 동작 불변 코드 최적화

1. gabes/schemes/rydberg.py:807-815는 Doppler-off finite-IF의 801개 detuning마다 Hamiltonian/Liouvillian을 Python에서 새로 만든다. detuning은 affine이므로 기준 행렬과 계수 행렬을 각 한 번 만들고 broadcast할 수 있다. 직전 보고서의 동일 행렬 벤치마크는 36.281 → 3.093 ms, np.array_equal 참이었다(docs/daily_report/2026-09-17_scheme-3_rydberg-eit.md:253-268). 행렬 조립만 바꾸며 물리·grid 불변.
2. 908-946행의 801개 detector dataclass/noise/sensitivity 생성은 식을 배열로 계산하고 최적점 하나만 객체화할 수 있다. 직전 보고서의 phasor-cache hit 후 readout 20.32 ms가 이 비용의 규모다(같은 보고서 270-278행). Schottky·RIN 식은 gabes/rydberg_electrometry.py:415-429와 scalar 결과를 교차검증해야 한다.
3. 가열 extra view는 10개 온도에서 EIT/AT solve를 반복한다(1384-1465). 네 cache 경계를 위 표처럼 나누면 detector-only 변경에는 원자 재계산이 없다. 1412-1423행은 표시용 slope 문자열을 split해서 숫자로 재파싱하므로 readout에 numeric field를 직접 보관하면 포맷 의존성도 없어진다. 두 변경 모두 값의 정의는 그대로다.

## 문서·검증 상태

README.md:115-137은 finite-IF, 온도 분리, SAM 및 한계를 적절히 설명한다. 반면 docs/Userguide/GABES_User_Guide_v2.html:617-625는 정적 EIT/AT 중심이고 현재 감도/온도/provenance 범위를 따라가지 못한다. tests/test_rydberg_eit.py:1-5의 “public UI shows the static spectrum only”도 오래됐다.

테스트는 4개 Scheme/helper 파일과 cell-heating workflow·adapter 통합 파일에서 **69 passed (31.19 s)**. 특히 저 IF와 정적 LO 미분 비교(tests/test_rydberg_eit.py:309-337), Hermitian sideband와 solve 횟수(tests/test_rydberg_electrometry.py:39-94), 온도 분리(tests/test_rydberg_eit.py:456-475), CSV 단위·provenance(tests/test_rydberg_experimental_csv.py:14-105)가 유효한 내부 검증이다. AGENTS.md 지시에 따라 실행한 python -m pytest -q는 63.31 s 뒤 references/AutoOD/test_autood_gui_adjustments.py 수집 중 선택 의존성 sympy 부재로 중단됐다. 전체 회귀 통과를 뜻하지 않는다.

**우선순위:** opacity provenance와 관측 가능성 gate → 고정 arm·cache 의존성 → 측정 EIT/RF/PSD held-out 검증. 성능 patch는 affine 조립 → detector 배열화 → Doppler-on pole 순서가 적절하다.
