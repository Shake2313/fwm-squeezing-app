# 2026-07-01 Scheme 2 Review: Lambda coherence

## 오늘 선택된 스킴

- 오늘 현지 날짜는 `2026-07-01`이고 day-of-month는 `1`이다.
- 규칙 `n = (day mod 5) + 1`에 따라 `n = (1 mod 5) + 1 = 2`이므로 오늘 검토 대상은 2번째 스킴이다.
- 현재 등록된 다섯 개 스킴 순서는 `OD / SAS -> Lambda coherence -> Rydberg-EIT electrometry -> Hanle / EIA / NMOR -> FWM`이다 (`gabes/schemes/__init__.py:19-24`, `README.md:10-16`).
- 따라서 검토 대상은 `LambdaScheme`의 `Lambda coherence (EIT / AT / CPT)`이다 (`gabes/schemes/absorption.py:462-667`).

## 읽은 범위와 이전 개선안 확인

- 스킴 등록과 순서: `gabes/schemes/__init__.py:19-24`
- Lambda 스킴 파라미터/기본값/계산 본체: `gabes/schemes/absorption.py:462-667`
- Lambda metric/readout: `gabes/schemes/absorption.py:669-743`
- Doppler 평균과 affine kernel 경로: `gabes/schemes/absorption.py:60-140`
- 3준위 Lambda 원자 모델: `gabes/atoms.py:157-165`
- beam power/diameter 및 residual-k helper: `gabes/beam.py:9-18`, `gabes/beam.py:48-53`
- Lambda 물리 회귀 테스트: `tests/test_absorption.py:110-175`
- headless observables 테스트: `tests/test_headless_observables.py:38-72`
- 기존 Lambda 일일 보고서: `docs/daily_report/2026-06-26_scheme-2_lambda-coherence.md`
- 저장소 TODO/checklist: `docs/checklist.json`

2026-06-26 Lambda 보고서에서 새 후보로 제안했던 세 가지 중 두 가지는 이미 구현되어 있다. 첫째, `coupling_power_mw`와 `coupling_diameter_mm`가 anchor Rabi를 `sqrt(P)/d`로 스케일한다 (`gabes/schemes/absorption.py:529-547`, `583-591`; `gabes/beam.py:9-18`). 둘째, `beam_angle_mrad`가 residual two-photon Doppler를 `atoms.lambda3(two_photon_doppler_ratio=...)`로 전달한다 (`gabes/schemes/absorption.py:537-542`, `593-640`; `gabes/beam.py:48-53`). 셋째, figureless/headless readout도 현재 API와 테스트에 들어와 있다 (`gabes/schemes/base.py:111-120`, `gabes/schemes/absorption.py:466`, `tests/test_headless_observables.py:38-72`).

반면 `buffer-gas-pressure-shift`는 아직 `deferred`이고, `lambda-hyperfine-resolved-manifold`는 `GROUP C` heavy item으로 남아 있다 (`docs/checklist.json:21-27`, `84-90`).

## 현재 구현의 물리 모델

- 이 스킴은 hyperfine/Zeeman-resolved manifold가 아니라 `g0, g1, e`의 축약 3준위 Lambda 모델이다. excited state는 두 ground state로 대칭 붕괴하고, dark resonance linewidth floor는 `buffer_ground_relax_khz`가 대표한다 (`gabes/atoms.py:157-165`, `gabes/schemes/absorption.py:549-553`).
- 매질은 `Rb (natural)`, `85Rb`, `87Rb`, `133Cs`, `Generic`과 D1/D2 선택을 통해 linewidth, wavelength, dipole, density를 바꾼다 (`gabes/schemes/absorption.py:418-443`, `521-526`).
- coupling은 이제 실험자가 만지는 power/beam diameter를 앞단에 두고, `coupling_rabi_mhz`는 Advanced anchor로 남겨 둔다. 이는 절대 dipole/beam calibration을 완전히 예측하지는 않지만, 실험실 knob와 계산 knob 사이의 연결을 꽤 잘 만든다 (`gabes/schemes/absorption.py:529-547`, `583-591`).
- Doppler-on 경로는 scan point마다 Maxwell velocity classes를 평균하고, residual angle mismatch가 없으면 two-photon resonance는 Doppler-free로 유지된다. 비영각에서는 `AtomModel.S_v`에 residual two-photon k-ratio가 들어가 EIT/CPT가 실제로 broadening/washout된다 (`gabes/schemes/absorption.py:118-140`, `593-640`).
- `observables()`는 `include_figures=False`일 때 Matplotlib 없이 metric/table만 만든다. 일반 UI figure path는 유지되므로 기능 손상 없이 batch/report path가 가벼워졌다 (`gabes/schemes/absorption.py:669-713`).

## 직접 실행 결과

실행 명령:

```bash
python -m pytest tests/test_absorption.py tests/test_kernels.py tests/test_headless_observables.py -q
```

결과는 `29 passed in 11.80s`였다.

추가로 기본값을 직접 측정했다. warm run의 대표값은 다음과 같다.

| Regime | compute | headless readout | figure readout | 주요 metric |
|---|---:|---:|---:|---|
| EIT | 67.8 ms | 0.3 ms | 139.0 ms | `T(res)=0.014`, `FWHM=0.46 MHz`, `n_g=1.003e5` |
| AT | 3.8 ms | 0.2 ms | 154.9 ms | `AT splitting=46.0 MHz`, expected `Omega_c=46.0 MHz`, `T(center)=0.989` |
| CPT | 3.8 ms | 0.2 ms | 204.6 ms | `T(res)=0.706`, `FWHM=923.32 kHz`, `n_g=9.525e4` |

lab knob scaling은 기대대로 동작했다. AT에서 anchor `Omega_c = 8 Gamma`일 때 `P=1 mW, d=1 mm -> 8 Gamma`, `P=4 mW, d=1 mm -> 16 Gamma`, `P=1 mW, d=2 mm -> 4 Gamma`가 나왔다 (`tests/test_absorption.py:121-129`도 같은 불변량을 본다).

beam-angle mismatch도 물리적으로 유의미하게 작동했다. EIT, `50 C`, `L=1 mm`, Doppler-on 조건에서:

| angle | residual `|Delta k|/k` | `T(res)` | EIT window FWHM |
|---:|---:|---:|---:|
| 0 mrad | `0.0000e+00` | 0.750724 | 0.4597 MHz |
| 1 mrad | `1.0000e-03` | 0.326469 | 0.9194 MHz |
| 5 mrad | `5.0000e-03` | 0.159116 | 4.1371 MHz |
| 10 mrad | `1.0000e-02` | 0.138894 | 7.8146 MHz |

이는 저장소 테스트의 “10 mrad에서 warm EIT가 넓어지고 약해져야 한다”는 조건과 일치한다 (`tests/test_absorption.py:132-148`).

## 실험 원자광학 연구자 관점의 평가

현재 `LambdaScheme`은 더 이상 단순 교과서 데모만은 아니다. `Omega_c`를 power/diameter로 연결하고, beam-angle mismatch가 dark resonance를 씻어내는 효과까지 들어갔기 때문에 실제 tabletop EIT/AT/CPT 실험의 초기 설계와 정렬 sanity check에는 꽤 유용하다.

특히 다음 용도에는 실험물리학적 reference로 쓸 만하다.

- coupling power나 beam diameter를 바꿨을 때 AT splitting이 어느 방향과 크기로 변하는지 빠르게 확인한다.
- residual beam angle이 warm-vapor EIT contrast와 linewidth를 얼마나 망가뜨리는지 정렬 감도로 본다.
- `buffer_ground_relax_khz`를 linewidth floor로 보고 wall/buffer/collision 환경이 dark resonance를 얼마나 흐리는지 반정량적으로 스캔한다.
- 같은 susceptibility에서 `cell_mm`만 바꾸어 OD/transmission 감도를 보는 planning을 한다 (`gabes/schemes/absorption.py:555-556`, `675-676`).

다만 absolute reference 수준에는 여전히 한계가 있다.

- 실제 Rb/Cs Lambda 실험의 hyperfine optical pumping, Zeeman sublevel redistribution, polarization selection은 아직 3-level lumped model 밖에 있다 (`gabes/atoms.py:157-165`).
- buffer gas는 ground-coherence relaxation knob로는 표현되지만, pressure shift, pressure broadening table, Dicke narrowing은 Lambda 본체에 아직 직접 들어가지 않는다 (`docs/checklist.json:21-27`).
- beam-angle residual Doppler는 매우 좋은 저차 모델이지만, full velocity-changing collision이나 diffusion-mediated Ramsey narrowing까지 포함하지는 않는다 (`gabes/schemes/absorption.py:118-140`).

따라서 오늘 기준 평가는 **실험 조건 변화의 방향성과 민감도, AT/EIT/CPT 핵심 scale을 보는 semi-quantitative lab reference**다. 논문 그림의 절대 linewidth/contrast를 그대로 fit하는 최종 모델이라기보다는, 실험 전에 knob sensitivity를 잡고 데이터 해석의 1차 sanity check를 하는 도구로 보는 것이 안전하다.

## 기존 개선안의 부하와 남은 물리 개선

### 이미 반영된 개선: lab-facing coupling Rabi

`coupling_power_mw`와 `coupling_diameter_mm`는 solve 전에 scalar Rabi만 바꾸므로 계산 부하는 사실상 없다. 다만 이 값들은 recompute knob라 slider를 움직이면 solve 자체는 다시 돈다 (`gabes/schemes/absorption.py:529-547`, `583-591`). 물리 효용은 높다. 실험가가 “내가 넣은 mW와 beam size”에서 AT splitting/EIT width가 어떻게 움직이는지 바로 보기 때문이다.

### 이미 반영된 개선: residual two-photon Doppler

`beam_angle_mrad`는 residual k-ratio를 `AtomModel.S_v`에 넣는 방식이라 solver dimension을 키우지 않는다 (`gabes/schemes/absorption.py:593-640`; `gabes/beam.py:48-53`). Doppler-on Lambda 경로는 원래 velocity-class 평균을 수행하므로, 추가 부하는 스칼라 coefficient 변화 정도에 가깝다. 반면 물리 효과는 매우 크다. 10 mrad에서 EIT 창 폭이 `0.4597 MHz -> 7.8146 MHz`로 커지고 공명 투과가 `0.750724 -> 0.138894`로 낮아졌다.

### 아직 남은 개선: buffer-gas pressure shift / Dicke narrowing

`docs/checklist.json`의 `buffer-gas-pressure-shift`는 가장 비용 대비 효과가 좋은 다음 물리 개선으로 보인다. pressure shift는 detuning offset scalar, pressure broadening은 homogeneous gamma scalar, phenomenological Dicke narrowing은 Doppler width 또는 effective residual-k correction 수준으로 시작할 수 있다. 이 방식이면 matrix dimension과 velocity grid coupling을 바꾸지 않아 부하가 작다. full velocity-changing-collision kernel까지 가면 현재 separable Maxwell averaging 구조가 깨져 큰 설계 변경이 된다 (`docs/checklist.json:21-27`).

### 무거운 개선: hyperfine-resolved Lambda manifold

`lambda-hyperfine-resolved-manifold`는 물리적으로는 가장 가치가 크다. 실제 line assignment, optical pumping, polarization-dependent contrast를 설명할 수 있기 때문이다. 그러나 현재 3준위 density matrix의 유효 자유도는 작고, hyperfine/Zeeman manifold로 가면 Liouvillian 차원이 급격히 커진다. 저장소도 이를 `GROUP C`로 분류한다 (`docs/checklist.json:84-90`). interactive reference 성격을 지키려면 우선 저차 proxy나 opt-in heavy mode로 설계해야 한다.

## 순수 코딩 최적화 후보

- Lambda의 headless path는 이미 들어가 있어 가장 큰 plotting 병목을 batch/report에서 피할 수 있다. 오늘 측정에서도 EIT `headless readout`은 `0.3 ms`였고 figure readout은 `139.0 ms`였다.
- 다음 미세 최적화는 `_medium_from_params()` 캐시다. species/line/temp에서 density, dipole, mass, wavelength를 매번 만들지만 같은 파라미터 조합은 반복된다 (`gabes/schemes/absorption.py:418-443`). 배치 스캔에서는 작은 `lru_cache`가 안전하다.
- `atoms.lambda3(...)`도 `(gamma, gamma_gg, two_photon_doppler_ratio)` 키로 캐시할 수 있다 (`gabes/schemes/absorption.py:638-640`). 큰 병목은 아니지만 자동 보고서와 parameter sweeps에서는 누적 비용을 줄인다.
- `observables.absorption_coefficient()`에서 만든 `alpha/xphys`는 `cell_mm` 변경과 독립이다 (`gabes/schemes/absorption.py:671-676`). `cell_mm`가 navigate-only knob인 만큼 raw 또는 cached observable intermediate로 유지하면 cell-length 스캔이 더 가벼워진다.

## 종합 판단

오늘 기준 `LambdaScheme`은 **실험가가 EIT/AT/CPT의 기본 scale, coupling-power sensitivity, beam alignment sensitivity를 빠르게 판단하는 semi-quantitative reference**로 쓸 수 있다. 지난 Lambda 리뷰에서 제안됐던 저부하 개선 중 핵심 두 가지가 이미 반영되어 물리적 유용성이 올라갔다.

다음으로 가장 현실적인 개선은 `buffer-gas-pressure-shift`를 Lambda에도 낮은 차수의 scalar correction으로 연결하는 것이다. 계산 부하를 거의 늘리지 않으면서 warm-vapor 실험의 line center, homogeneous width, dark-feature contrast를 더 잘 따라갈 수 있다. hyperfine-resolved Lambda는 가치가 크지만 부하와 모델 해석 변화가 커서 별도 설계가 필요한 heavy mode로 남기는 편이 맞다.
