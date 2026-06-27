# 2026-06-25 Scheme 1 Review: OD / SAS

## 오늘의 선택 스킴

- 오늘 현지 날짜는 `2026-06-25`이고, `n = (day mod 5) + 1 = (25 mod 5) + 1 = 1`이다.
- 현재 등록된 다섯 개 스킴 순서는 `OD / SAS -> Lambda coherence -> Rydberg-EIT electrometry -> Hanle / EIA / NMOR -> FWM`이다 (`gabes/schemes/__init__.py:19-24`, `README.md:12-16`).
- 따라서 오늘 검토 대상은 1번째 스킴 `SASScheme`이며, UI 제목은 `Absorption spectroscopy (OD / SAS)`이다 (`gabes/schemes/sas.py:53-62`).

## 이번 검토에 사용한 근거

- 스킴 등록 순서: `gabes/schemes/__init__.py:19-24`
- 현재 스킴 정의와 파라미터: `gabes/schemes/sas.py:53-106`
- species 기반 다중선 OD/SAS 계산 경로: `gabes/schemes/sas.py:156-246`
- generic fallback 계산 경로: `gabes/schemes/sas.py:248-312`
- 관측량/플롯/metric 구성: `gabes/schemes/sas.py:317-401`
- 원자 데이터와 manifold 구성: `gabes/species.py:192-370`
- buffer-gas 관련 현재 구현과 TODO: `gabes/constants.py:50-60`
- 스킴 설명과 현재 문서화된 주장: `README.md:12`, `README.md:177-198`, `docs/GABES_User_Guide_v2.html`
- 물리 테스트: `tests/test_sas.py:95-179`, 보조 OD 검증: `tests/test_absorption.py:49-106`

## 현재 정의와 구조

이 스킴은 “약한 probe + 역전파 pump” 흡수분광을 하나의 엔진으로 묶는다. `pump_power_mw = 0`이면 선형 Doppler-broadened OD이고, pump를 올리면 같은 스펙트럼 위에 Lamb dip과 crossover가 새겨지는 SAS로 연속적으로 넘어간다 (`gabes/schemes/sas.py:66-70`, `README.md:177-179`).

핵심 구현은 꽤 실험친화적이다.

- species 모드에서는 `species.build_manifold()`가 isotope/line별 전체 `Fg -> Fe` hyperfine manifold와 CG-branched decay, transit relaxation을 만든다 (`gabes/species.py:336-370`, `gabes/schemes/sas.py:169-178`).
- pump power는 beam waist와 `I_sat`을 써서 pump Rabi로 바뀐다. 즉, 사용자가 만지는 광파워가 계산 해밀토니안에 직접 연결된다 (`gabes/species.py:272-280`, `gabes/schemes/sas.py:174`).
- probe는 약하게 고정되어 population readout만 하고, pump에 의해 준비된 population difference를 probe Lorentzian에 실어 Doppler 평균한다 (`gabes/schemes/sas.py:214-245`).
- 85Rb D1의 pump-off limit은 기존 AutoOD 절대 스케일에 맞춰져 있다. 이는 단순 “모양만 비슷한” 수준이 아니라, 적어도 이 기준점에서는 흡수 절대 크기까지 맞추려는 의도가 분명하다 (`gabes/schemes/sas.py:124-145`, `gabes/species.py:192-205`, `gabes/species.py:248-270`, `README.md:179-183`).

정리하면 이 스킴은 “OD와 SAS를 서로 다른 toy model 두 개로 붙인 것”이 아니라, 같은 hyperfine medium 위에서 pump를 0으로 내리면 OD, 올리면 SAS가 되는 식으로 설계되어 있다. 실험물리 관점에서 이 점이 가장 좋다.

## 테스트와 오늘 확인한 실행 결과

오늘 직접 실행한 테스트:

- `python -m pytest tests/test_sas.py tests/test_absorption.py -q`

결과:

- 총 `23 passed`였다.

이 테스트들이 실제로 확인하는 물리는 다음과 같다.

- pump-off 85Rb D1이 AutoOD 기준과 적분값/peak 값에서 1% 이내로 맞는지 (`tests/test_sas.py:95-105`)
- pump-off 스펙트럼이 매끈한 Doppler 배경이고, 85Rb D1의 `F=3 / F=2` manifold 세기비가 검증된 `49/25`에 맞는지 (`tests/test_sas.py:108-113`, `gabes/species.py:370`)
- pump-on에서 sub-Doppler feature가 실제로 생기고, transit relaxation을 줄일수록 crossover transmission이 커지는지. 즉 hyperfine optical pumping signature가 살아 있는지 (`tests/test_sas.py:117-141`)
- natural Rb가 85Rb와 87Rb marker를 함께 올리는지 (`tests/test_sas.py:144-148`)
- generic fallback에서도 단일 Lamb dip과 two-line crossover가 되는지 (`tests/test_sas.py:169-179`)

추가로 오늘 기본적인 성능/동작을 직접 재보니, RB85 D2 SAS 기본형(`30 C`, `50 mm`, `1.5 mW`)에서 warm run 기준 `compute()`는 약 `0.46-0.59 s`, `observables()`는 약 `0.25 s`였다. 즉 데스크탑 상호작용용으로는 충분히 실용적이지만, 이미 “아주 가벼운 analytic toy”는 아니다.

## 실험물리 연구자의 관점에서 본 평가

### 무엇이 실제 물리에 유용한가

이 스킴은 실제 원자분광 실험에서 바로 감이 오는 축들을 잘 구현하고 있다.

- isotope/line 선택이 실제 D1/D2 hyperfine manifold와 연결된다 (`gabes/schemes/sas.py:71-79`, `gabes/species.py:336-370`).
- pump power, beam waist, cell temperature, cell length, transit relaxation이 모두 스펙트럼 형상에 물리적으로 납득 가능한 방향으로 연결된다 (`gabes/schemes/sas.py:66-98`, `gabes/species.py:192-205`, `gabes/species.py:272-280`).
- SAS의 핵심인 “단순 saturation hole”이 아니라, CG-branched decay를 통한 hyperfine optical pumping이 crossover의 강화/반전으로 나타나도록 모델링되어 있다 (`gabes/schemes/sas.py:127-145`, `README.md:187-198`, `tests/test_sas.py:127-141`).
- natural Rb overlay까지 제공하므로, 실제 랩에서 가장 흔한 “셀에 자연비 rubidium이 들어 있는 경우”를 빠르게 대조하기 좋다 (`tests/test_sas.py:144-148`).

실험 준비 단계에서 유용한 질문, 예를 들면 아래와 같은 것에는 꽤 쓸 만하다.

- 지금 보이는 dip/crossover가 어느 hyperfine 그룹인지
- 펌프를 올리면 crossover가 더 강해져야 하는지
- transit relaxation이나 waist를 바꾸면 feature contrast가 어떤 방향으로 바뀌는지
- natural Rb 셀과 isotope-enriched 셀의 qualitative 차이가 어떠한지

### 어디까지 레퍼런스로 쓸 수 있는가

내 평가는 다음과 같다.

- **정성적 reference**로는 충분히 좋다.
- **semi-quantitative lab reference**로도 꽤 유용하다.
- 다만 **정밀한 절대 계측 reference**로 쓰기에는 아직 조심해야 한다.

그 이유는 분명하다.

- 강점: pump-off 85Rb D1 절대 스케일이 AutoOD에 맞춰져 있고, 상대 line weight도 검증되어 있다 (`tests/test_sas.py:95-113`, `tests/test_absorption.py:80-99`).
- 강점: pump-on physics도 최소한 “실험실에서 가장 중요한 현상학”인 crossover enhancement/inversion을 직접 재현한다 (`tests/test_sas.py:127-141`, `README.md:187-198`).
- 한계: 현재 buffer gas는 단순 homogeneous broadening으로만 들어간다. pressure shift, Dicke narrowing, velocity-changing collision은 없다 (`gabes/constants.py:50-60`, `gabes/schemes/sas.py:83-86`, `162`, `257`).
- 한계: probe beam alignment mismatch, pump/probe polarization impurity, etalon/fringe, RAM 같은 실제 SAS 셋업의 흔한 비이상성은 모델 밖이다.
- 한계: transit relaxation은 효과적으로 한 개의 phenomenological rate로만 들어간다. 실제 beam profile, diffusion, wall collision hierarchy는 분해되지 않는다 (`gabes/schemes/sas.py:90-94`, `gabes/species.py:336-356`).

따라서 이 코드는 “레이저 잠금용 SAS 셋업을 설계하고 스펙트럼을 해석하는 연구자”에게는 이미 충분히 참고할 만하지만, pressure-buffered reference cell의 절대 line center shift나 세밀한 linewidth budget까지 바로 믿고 가져갈 정도의 최종 reference는 아직 아니다.

## 저장소에 이미 있던 개선안과 계산 부하

이번 검토에서 `OD / SAS` 전용으로 따로 적힌 개선안은 찾지 못했다. 기존 daily report, README, TODO성 메모, `docs/checklist.json`을 확인했지만, 현재 남아 있는 명시적 개선안 중 이 스킴에 직접 연결되는 것은 공통 항목 하나뿐이었다.

- `buffer-gas-pressure-shift` (`docs/checklist.json:18-23`)

이 항목은 현재 구현이 `neon_buffer_broadening()` 하나로 끝나기 때문에 생긴 TODO다 (`gabes/constants.py:50-60`). 지금은 pressure가 단지 `buffer_gamma`를 키워 homogeneous width에 더해지는 구조다 (`gabes/schemes/sas.py:162`, `175`, `257-258`).

이 개선안을 부하 관점에서 나누면 다음과 같다.

### 1. pressure shift 추가

- 구현 위치는 거의 `constants.py`의 coefficient table과 scan축 기준선 이동 쪽에 국한될 가능성이 크다.
- 계산량은 사실상 거의 늘지 않는다. 기존 선형대수 solve 수를 늘리지 않고, 각 transition 중심 또는 effective detuning에 상수항을 더하는 수준이기 때문이다.
- 실험물리적 가치는 크다. buffer cell을 쓰는 SAS/OD에서는 lock point 해석과 line assignment에서 pressure shift가 바로 문제로 들어오기 때문이다.

### 2. phenomenological Dicke narrowing 추가

- 이것도 우선은 `gamma_eff` 또는 Lorentzian/Voigt width를 pressure-dependent effective width로 바꾸는 식의 저차 보정으로 시작할 수 있다.
- 이 단계 역시 solve 수 증가 없이 scalar algebra만 조금 늘어나는 수준이라 계산 부하는 매우 작다.
- 실제로는 pressure를 올렸을 때 feature가 무조건 넓어지기만 하는 현재 모델의 편향을 완화해 줄 수 있다.

### 3. full velocity-changing collision 모델

- 여기부터는 얘기가 달라진다.
- velocity class 간 population/coherence 재분배까지 넣으려면 현재의 독립 Maxwell weight 합 구조 (`gabes/schemes/sas.py:211-245`)를 깨고, 속도 클래스 결합이나 추가 relaxation operator를 풀어야 할 가능성이 높다.
- 그러면 solve 크기나 solve 횟수가 증가해 interactive 성능이 꽤 희생될 수 있다.

정리하면, **저장소에 이미 적혀 있는 개선안 중 pressure shift + phenomenological Dicke narrowing 정도는 계산 부하를 거의 만들지 않으면서 물리를 유의미하게 끌어올릴 수 있는 좋은 후보**이고, full VCC까지 가면 그때부터는 무거워진다.

## 이전 개선안이 부족한 부분에서 조심스럽게 제안할 수 있는 물리 개선

OD/SAS 전용 TODO가 따로 없었기 때문에, 코드와 모델을 보고 조심스럽게 적을 만한 후보는 아래 정도다.

### 1. pump/probe polarization 선택지

현재 스킴 설명은 hyperfine pumping을 잘 담고 있지만, 편광에 따른 transition selectivity는 사용자가 직접 만지는 축으로 드러나지 않는다. 실제 SAS 실험에서는 linear/circular, pump-probe polarization mismatch가 crossover contrast를 크게 바꾼다. 이를 full Zeeman manifold까지 가지 않더라도 effective branching selector 정도로 부분 도입하면 실험 대응력이 더 좋아질 수 있다. 다만 이건 pressure shift보다 부하와 모델 복잡도가 더 크다.

### 2. isotope/line별 buffer-gas coefficient table

지금 TODO도 사실 이 방향을 암시한다 (`gabes/constants.py:52-53`). Ne만 아니라 buffer gas 종류, species, line마다 broadening/shift coefficient를 테이블화하면 실제 reference-cell 사용성과 문헌 대조력이 한 단계 올라간다. 이건 계산 부하가 거의 없다.

### 3. lock-point readout

현재는 transmission/OD/marker 중심의 스펙트럼 판독이다 (`gabes/schemes/sas.py:324-401`). 실험자는 종종 “어디에 락을 걸면 slope가 가장 큰가”를 바로 보고 싶어 한다. 따라서 `dT/dΔ` 최대 지점, crossover 중심 slope, lock discriminator proxy를 추가 metric으로 주면 물리 엔진을 바꾸지 않고도 실험적 유용성이 크게 오른다.

## 순수 코딩 측면의 속도 개선 후보

물리를 건드리지 않고 속도를 더 줄일 수 있는 부분도 보인다.

### 1. `compute()`와 figure 생성 분리

현재 `observables()`는 metric 계산과 Matplotlib figure 생성을 같이 한다 (`gabes/schemes/sas.py:324-401`). 자동 보고서, 테스트, 배치 분석에서는 metric/table만 필요한 경우가 많으므로, figure-less fast path를 두면 체감 속도를 줄일 수 있다. 물리 결과는 변하지 않는다.

### 2. pump population interpolation 재사용 강화

현재는 각 isotope component마다 `deff` table을 만든 뒤, 각 transition에 대해 `np.interp`를 반복한다 (`gabes/schemes/sas.py:218-245`). 동일한 manifold/temperature/grid 조합에서 재사용 가능한 중간 배열을 캐시하면 연속 파라미터 스캔에서 Python 레벨 오버헤드를 줄일 여지가 있다.

### 3. marker/derived table 구성의 경량화

`markers` 문자열과 관측용 table은 매 run마다 다시 조립된다 (`gabes/schemes/sas.py:187-202`, `324-401`). isotope/line 고정 상태에서는 거의 정적 데이터이므로, recompute가 필요 없는 knob와 연동해 캐시하면 UI responsiveness를 조금 더 다듬을 수 있다.

## 종합 판단

`OD / SAS` 스킴은 현재 GABES 안에서 “실험실 분광 감각”이 가장 직접적으로 살아 있는 축 중 하나다. pump-off OD와 pump-on SAS가 한 medium 안에서 자연스럽게 이어지고, hyperfine optical pumping까지 포함해 crossover의 질적 특징을 재현한다는 점이 크다. 그래서 **레이저 잠금용 알칼리 흡수분광을 다루는 원자광학 실험물리학자에게는 이미 꽤 유용한 semi-quantitative reference**라고 판단한다.

다만 pressure-buffered cell의 정밀 기준선, 세밀한 linewidth budget, collision-rich 환경까지 바로 신뢰할 정도의 정량 레퍼런스는 아직 아니다. 오늘 기준 최우선 개선 방향은 새롭고 무거운 모델을 한꺼번에 넣는 것보다, **이미 저장소에 적혀 있는 buffer-gas coefficient 확장과 pressure shift/Dicke narrowing의 저부하 도입**이다. 이 방향이 물리적 실익 대비 계산 비용이 가장 좋다.
