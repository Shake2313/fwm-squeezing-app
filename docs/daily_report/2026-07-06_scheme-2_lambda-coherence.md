# 2026-07-06 Scheme 2 Review: Lambda coherence

## 오늘 선택된 스킴

- 오늘 현지 날짜는 `2026-07-06`이고 day-of-month는 `6`이다.
- 규칙 `n = (day mod 5) + 1`에 따라 `n = (6 mod 5) + 1 = 2`이다.
- 현재 등록된 다섯 개 스킴 순서는 `OD / SAS -> Lambda coherence -> Rydberg-EIT electrometry -> Hanle / EIA / NMOR -> FWM`이다 (`gabes/schemes/__init__.py:19-24`, `README.md:12-16`).
- 따라서 오늘 검토 대상은 두 번째 스킴인 `LambdaScheme`, 즉 `Lambda coherence (EIT / AT / CPT)`이다 (`gabes/schemes/absorption.py:462-721`).

## 읽은 범위

- 스킴 등록과 공통 계약: `gabes/schemes/__init__.py:19-24`, `gabes/schemes/base.py:68-106`
- Lambda 본체: `gabes/schemes/absorption.py:462-721`
- Doppler 평균 경로: `gabes/schemes/absorption.py:118-140`
- D-line medium 구성: `gabes/schemes/absorption.py:418-443`
- 3준위 Lambda 원자 모델: `gabes/atoms.py:157-170`
- beam power/diameter 및 residual-k helper: `gabes/beam.py:9-20`, `gabes/beam.py:50-55`
- 물리 회귀 테스트: `tests/test_absorption.py:121-148`, `tests/test_absorption.py:180-191`
- headless readout 테스트: `tests/test_headless_observables.py:38-88`
- 기존 개선 메모: `docs/checklist.json:21-25`, `docs/checklist.json:56-60`, `docs/checklist.json:84-88`
- 직전 Lambda 리뷰: `docs/daily_report/2026-07-01_scheme-2_lambda-coherence.md`

## 현재 구현의 물리

`LambdaScheme`은 `g0, g1, e`의 축약 3준위 Lambda OBE다. 여기서 excited state는 두 ground state로 대칭 붕괴하고, dark resonance의 바닥 linewidth는 `buffer_ground_relax_khz`가 대표한다 (`gabes/atoms.py:157-170`, `gabes/schemes/absorption.py:549-553`). 이 구조는 실제 Rb/Cs hyperfine manifold 전체를 풀지는 않지만, EIT 투명창, Autler-Townes splitting, CPT dark resonance라는 세 가지 대표 실험 모드를 한 엔진에서 비교하기에는 잘 맞는다.

실험가 관점에서 유용한 점은 두 가지가 이미 코드에 들어와 있다는 점이다. 첫째, `coupling_power_mw`와 `coupling_diameter_mm`가 anchor Rabi를 `sqrt(P) / diameter`로 스케일한다 (`gabes/schemes/absorption.py:529-547`, `583-591`; `gabes/beam.py:9-20`). 둘째, `beam_angle_mrad`가 residual two-photon Doppler를 `atoms.lambda3(two_photon_doppler_ratio=...)`로 전달한다 (`gabes/schemes/absorption.py:537-542`, `593-646`; `gabes/beam.py:50-55`). 즉 "결합광 파워를 올리면 AT split이 얼마나 커지는가", "빔 정렬이 mrad 수준으로 틀어지면 warm-vapor EIT가 얼마나 망가지는가"를 실험실 knob와 거의 같은 언어로 물어볼 수 있다.

반대로 이 스킴은 아직 absolute spectroscopic reference는 아니다. hyperfine optical pumping, Zeeman sublevel population redistribution, polarization selection rule, buffer-gas pressure shift, Dicke narrowing은 compact 3-level model 밖에 있다 (`docs/checklist.json:21-25`, `docs/checklist.json:84-88`). 따라서 현재 사용 등급은 **실험 조건 변화의 방향성과 핵심 scale을 보는 semi-quantitative lab reference**가 적절하다. 논문 피팅용 최종 모델이라기보다, 실험 전 파라미터 감도와 데이터 해석의 1차 sanity check에 강하다.

## 오늘 실행 결과

검증 명령:

```bash
python -m pytest tests/test_absorption.py tests/test_headless_observables.py -q
```

결과는 `15 passed in 59.17s`였다.

추가 측정은 warm run으로 다시 수행했다.

| Regime | compute | headless readout | figure readout | 주요 metric |
|---|---:|---:|---:|---|
| EIT | 955.39 ms | 3.087 ms | 2413.49 ms | `T(res)=0.014`, `FWHM=0.46 MHz`, `n_g=1.003e5` |
| AT | 56.71 ms | 0.433 ms | 847.67 ms | `AT splitting=46.0 MHz`, expected `Omega_c=46.0 MHz`, `T(center)=0.989` |
| CPT | 48.90 ms | 0.228 ms | 1000.28 ms | `T(res)=0.706`, `FWHM=923.32 kHz`, `n_g=9.525e4` |

기본 EIT가 AT/CPT보다 무거운 이유는 기본값이 Doppler-on이고, scan point마다 velocity-class 평균을 수행하기 때문이다 (`gabes/schemes/absorption.py:118-140`, `gabes/schemes/absorption.py:482-489`). 같은 Lambda 엔진이라도 "warm vapor EIT reference"는 solve-bound이고, "cold/compact AT/CPT reference"는 훨씬 가볍다.

beam-angle mismatch sweep도 재확인했다. `50 C`, `L=1 mm`, Doppler-on EIT 조건에서:

| angle | residual `|Delta k|/k` | `T(res)` | EIT window FWHM |
|---:|---:|---:|---:|
| 0 mrad | `0.0000e+00` | 0.772036 | 0.9704 MHz |
| 1 mrad | `1.0000e-03` | 0.352611 | 0.9704 MHz |
| 5 mrad | `5.0000e-03` | 0.163975 | 4.3668 MHz |
| 10 mrad | `1.0000e-02` | 0.141192 | 7.7632 MHz |

이 결과는 residual two-photon Doppler knob가 단순 장식이 아니라 실제 warm-vapor alignment sensitivity를 강하게 반영한다는 뜻이다. 특히 1 mrad에서 contrast가 이미 크게 줄고, 5-10 mrad에서 linewidth가 MHz scale로 넓어지는 방향은 실험적으로 타당하다.

## 기존 개선안과 부하 평가

기존 개선안은 명확히 존재한다. `lambda-residual-two-photon-doppler`는 `done` 상태이고, 현재 코드에서 residual-k를 원자 모델의 Doppler ratio로 넣는다 (`docs/checklist.json:56-60`, `gabes/schemes/absorption.py:593-646`). 이 방식은 solver dimension을 늘리지 않는다. Doppler-on 경로에서는 원래 velocity 평균을 하므로 추가 부하는 거의 residual coefficient 변경 수준이고, Doppler-off 경로에서는 실질적으로 scalar 전처리다. 물리 효용 대비 비용이 매우 좋다.

`buffer-gas-pressure-shift`는 아직 `deferred`다 (`docs/checklist.json:21-25`). Lambda에서는 이미 `buffer_ground_relax_khz`로 ground coherence relaxation을 조절하지만, optical transition의 pressure shift, homogeneous pressure broadening table, Dicke narrowing은 별도 모델로 들어가지 않는다 (`gabes/schemes/absorption.py:549-553`). 이것을 낮은 차수로 넣는다면 계산 부하는 작다. pressure shift는 scan center 또는 detuning offset scalar, pressure broadening은 `gamma` scalar, phenomenological Dicke narrowing은 Doppler width 또는 residual-k effective factor를 조절하는 방식으로 시작할 수 있다. full velocity-changing collision kernel을 넣는 순간 현재의 separable Maxwell averaging 구조를 건드려야 하므로 훨씬 무거운 설계가 된다 (`gabes/schemes/absorption.py:118-140`).

`lambda-hyperfine-resolved-manifold`도 `deferred`이며 GROUP C heavy item이다 (`docs/checklist.json:84-88`). 물리적으로는 가장 큰 개선이다. 실제 Lambda 실험에서 line assignment, optical pumping, polarization-dependent contrast가 여기에 걸리기 때문이다. 그러나 현재 3-level density matrix에서 hyperfine/Zeeman manifold로 확장하면 Liouvillian 차원과 파라미터 해석이 크게 바뀐다. interactive reference로서의 속도를 유지하려면 기본 스킴에 바로 넣기보다 opt-in heavy mode 또는 저차 proxy부터 설계하는 편이 맞다.

## 순수 코딩 최적화 후보

- 자동 보고서, batch sweep, 테스트에서는 `headless_observables()`를 기본으로 쓰는 것이 가장 효과적이다. 오늘 측정에서도 Lambda figure readout은 약 `0.85-2.41 s`였고 headless readout은 `0.23-3.09 ms`였다 (`gabes/schemes/base.py:92-106`, `gabes/schemes/absorption.py:677-721`).
- `_medium_from_params()`는 species, line, temperature에서 density, dipole, mass, wave vector를 매번 계산한다 (`gabes/schemes/absorption.py:418-443`). 같은 조합을 반복하는 batch scan에서는 작은 `functools.lru_cache` wrapper가 기능 변경 없이 누적 비용을 줄일 수 있다.
- `atoms.lambda3(gamma_gg, gamma, two_photon_doppler_ratio)`도 동일 조합이 자주 반복된다 (`gabes/schemes/absorption.py:644-646`). 객체 생성 자체가 주 병목은 아니지만, parameter sweep에서 base dissipator 구성을 반복하는 비용은 줄일 수 있다.
- `cell_mm`는 `recompute=False`이고 Beer-Lambert 변환에서만 쓰인다 (`gabes/schemes/absorption.py:555-556`, `677-684`). 현재도 solve 재실행은 피하지만, alpha/xphys 중간값을 observable cache 내부에서 더 명시적으로 재사용하면 cell-length만 바꾸는 sweep이 더 가벼워질 수 있다.

## 종합 판단

오늘 기준 `LambdaScheme`은 실험물리학자가 EIT/AT/CPT 실험의 기본 scale, coupling-power dependence, beam-alignment sensitivity를 빠르게 확인하는 데 충분히 유용하다. 특히 mW/mm 단위 결합광 knob와 mrad 단위 residual Doppler knob가 들어온 뒤로는 단순 교과서 그림을 넘어 실험 준비용 reference에 가까워졌다.

다음 물리 개선 우선순위는 여전히 `buffer-gas-pressure-shift`의 낮은 차수 구현이다. matrix dimension을 늘리지 않고도 line center, homogeneous width, warm-vapor narrowing trend를 더 잘 따라갈 수 있다. hyperfine-resolved Lambda manifold는 가치가 크지만, 속도와 모델 해석을 크게 바꾸므로 별도 heavy mode로 분리하는 것이 안전하다.
