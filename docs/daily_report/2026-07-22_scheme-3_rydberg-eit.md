# 2026-07-22 Scheme 3 물리 검토: Rydberg-EIT electrometry

## 1. 오늘의 선택과 현재 다섯 스킴 순서

서울 현지 날짜의 일(day)은 22이므로

```text
n = (22 mod 5) + 1 = 3
```

이다. 현재 드롭다운 등록 순서는 `gabes/schemes/__init__.py:19-25`와
`README.md:8-16`에서 다음과 같이 확인된다.

| 순번 | 등록 객체 | 사용자 스킴 | 핵심 출력 |
|---:|---|---|---|
| 1 | `SASScheme()` | Absorption OD / SAS | pump-off OD, pump-on SAS |
| 2 | `LambdaScheme()` | Lambda coherence | EIT / AT / CPT |
| **3** | **`RydbergEITScheme()`** | **Rydberg-EIT electrometry** | **cascade EIT, microwave AT** |
| 4 | `MagnetoScheme()` | Hanle / EIA / NMOR | 영자기장 투과와 광회전 |
| 5 | `FWMScheme()` | FWM | seeded squeezing, generic biphoton |

따라서 오늘의 대상은 `gabes/schemes/rydberg.py`의 세 번째 스킴이다.

## 2. 선행 제안, TODO, 일일 보고서, issue 노트 선검색

먼저 다음 자료를 확인했다.

- 계획/TODO: `docs/checklist.json:1-142`
- 이전 Scheme 3 보고서:
  `docs/daily_report/2026-06-22_scheme-3_rydberg-eit.md`,
  `2026-06-27_scheme-3_rydberg-eit.md`,
  `2026-07-02_scheme-3_rydberg-eit.md`,
  `2026-07-12_scheme-3_rydberg-eit.md`
- 현재 설명: `README.md:8-16`, `README.md:44-58`,
  `docs/GABES_User_Guide_v2.html:595-605`, `gabes/schemes/rydberg.py:653-670`
- 전용/공통 검증: `tests/test_rydberg_eit.py:29-257`,
  `tests/test_kernels.py:102-136`, `tests/test_headless_observables.py:38-67`,
  `tests/test_schemes_render.py:109-116`
- 원격 issue: 2026-07-22 현재
  [GitHub Issues](https://github.com/Shake2313/fwm-squeezing-app/issues?q=is%3Aissue)는
  open/closed 모두 0건이며, 저장소 안에도 별도의 issue-note 파일은 발견되지 않았다.

기존 제안은 분명히 존재한다.

1. `rydberg-power-to-rabi`는 완료 상태다. probe/coupling power와 beam diameter를
   anchored Rabi로 연결하고, residual Doppler, 두 dephasing channel, AT center shift,
   finite-IF proxy, temperature dephasing을 추가했다
   (`docs/checklist.json:7-13`, `gabes/schemes/rydberg.py:222-253`, `306-392`).
2. figureless/headless readout도 완료됐다
   (`docs/checklist.json:29-34`, `gabes/schemes/rydberg.py:439-606`).
3. `low-order-polarization-zeeman-proxies`는 deferred다
   (`docs/checklist.json:87-93`).
4. full Zeeman/polarization/stray-field/ionization 모델과 full time-domain
   superheterodyne chain은 GROUP C deferred다
   (`docs/checklist.json:129-141`).

독립적인 Rydberg 예제 스크립트는 없다. 대신 `extra_views()`의 16점 probe-power
sweep가 실행 가능한 Fig. 2(b) 예제이고 (`gabes/schemes/rydberg.py:608-651`),
정적 사용자 가이드 그림 `docs/userguide_assets/rydberg.png`와 설명이 사용 예를
제공한다 (`docs/GABES_User_Guide_v2.html:595-605`).

## 3. 현재 구현이 담는 실제 물리

현재 모델은 85Rb
`5S1/2 F=3 → 5P3/2 F'=4 → 40D5/2 → 39F7/2`의 정적 4준위
cascade OBE다. 780 nm probe, 481 nm counter-propagating coupling, 37 GHz RF를
Hamiltonian의 세 결합으로 넣는다 (`gabes/schemes/rydberg.py:4-12`, `338-346`).
5P의 자연붕괴, 두 Rydberg 준위의 축약 붕괴, 5S–40D 및 40D–39F coherence
dephasing을 Lindblad/현상론 항으로 넣는다 (`gabes/schemes/rydberg.py:59-78`,
`200-219`; `gabes/atoms.py:56-82`).

실험 knob와 solve의 연결도 유용하다.

- probe와 coupling Rabi는 기준점에 고정한 `sqrt(P)/d` 스케일을 따른다
  (`gabes/beam.py:9-20`, `gabes/schemes/rydberg.py:222-244`).
- beam diameter와 온도는 thermal transit broadening을 정한다
  (`gabes/beam.py:23-28`, `gabes/schemes/rydberg.py:246-253`).
- 온도는 85Rb vapor density에도 들어가고, 선택적으로 추가 dephasing에 들어간다
  (`gabes/species.py:193-208`, `gabes/schemes/rydberg.py:311-329`).
- 5P는 full probe wavevector, 40D/39F는 `(k_probe-k_coupling)/k_probe`를 가지므로
  `doppler="on"`에서는 1광자 및 residual 2광자 Doppler를 함께 Maxwell 평균한다
  (`gabes/schemes/rydberg.py:59-78`, `348-357`; `gabes/atoms.py:84-98`).
- steady-state coherence는 density, dipole, cell length와 결합되어 Beer–Lambert
  transmission이 된다 (`gabes/schemes/rydberg.py:367-410`,
  `gabes/observables.py:393-418`).
- compensated/uncompensated EIT, EIT FWHM, RF AT splitting, detuned-AT center,
  static slope와 finite-IF discriminator를 읽는다
  (`gabes/schemes/rydberg.py:359-365`, `394-543`).

JIT warm 상태의 기준 출력은 이전 보고서와 재현됐다.

| 조건 | 주요 결과 |
|---|---|
| EIT reference, Doppler off | linewidth 1.61 MHz, resonance T 0.940, max slope/IF 0.122 MHz⁻¹ |
| AT reference, Doppler off | RF split 3.50 MHz, center +0.00 MHz, resonance T 0.806 |
| EIT reference, Doppler on | linewidth 2.92 MHz, resonance T 0.9966, contrast 약 0.00183 |

따라서 EIT/AT의 존재, power broadening, transit penalty, microwave dressed-state
splitting과 detuning 방향성은 모두 실제 원자물리다. 논문 기준점 주변의 정적 스펙트럼
설계와 knob 방향성 확인에는 유용하다. Ju et al.의 1.6 MHz 선폭과 37 GHz,
6 µW/30 mW/0.15 mm/50 mm 기준도 코드와 테스트에 고정되어 있다
(`gabes/schemes/rydberg.py:92-112`; `tests/test_rydberg_eit.py:29-61`,
`205-241`; [arXiv:2606.04354](https://arxiv.org/abs/2606.04354)).

## 4. 실험물리 reference로서의 한계와 이번 검토의 새 관찰

### 4.1 가장 중요한 경계: Doppler-off와 숨은 `line_strength=0.001`

기본 `doppler="off"`는 residual 2광자 Doppler만 끄는 것이 아니라 `kv=[0]`을
사용해 중간 5P의 1광자 Doppler도 함께 끈다 (`gabes/schemes/rydberg.py:348-357`).
그 상태에서 absolute susceptibility 변환에는 설명 없는
`line_strength=0.001`이 들어간다 (`gabes/schemes/rydberg.py:367-378`, `401-406`).

같은 저장소의 hyperfine 함수로 85Rb D2 `F=3→F'=4`의 `C_F²`를 계산하면 1.0이다
(`gabes/species.py:226-237`). 즉 0.001은 알려진 cycling-line strength가 아니라,
미해상 Doppler velocity fraction, ground/mF participation, optical pumping, mode/polarization,
fit amplitude를 한꺼번에 흡수한 **유효 opacity anchor**로 해석해야 한다. 현재 derived
table과 help에는 이 의미가 드러나지 않는다.

진단적으로 같은 coherence에 factor 1을 넣으면 Doppler-off 기준 공진 투과는
`1.37×10⁻²⁷`까지 내려간다. 반대로 Doppler-on에서 현재 0.001을 유지하면 공진 투과가
0.9966으로 올라가 EIT contrast가 약 0.00183에 그친다. 이 극단적인 차이는
`doppler=off + 0.001` 조합이 논문 기준 모양을 맞춘 empirical reference이지,
temperature/cell/absolute OD를 독립 예측하는 모델은 아니라는 증거다.

권고는 full manifold보다 먼저 다음처럼 calibration boundary를 공개하는 것이다.

- `line_strength`를 `effective_opacity_participation` 같은 명시적 이름과 derived-table
  항목으로 노출하고, 0.001이 atomic `C_F²`가 아님을 help/info에 적는다.
- `calibrated stationary reference`와 `velocity-resolved` 모드를 구분하고, 후자의
  opacity anchor를 별도로 검증한다.
- default transmission/contrast 및 Doppler-grid convergence test를 추가한다. 현재
  Doppler test는 단지 `width_on > width_off`만 요구한다
  (`tests/test_rydberg_eit.py:174-184`).

스칼라 표시와 테스트는 solve 비용을 늘리지 않는다. velocity-resolved 경로를 기본으로
바꾸는 것은 별도 문제다. 현재 679개 velocity class(`dv=2 m/s`) 때문에 warm compute가
약 0.22 s로, Doppler-off 약 3.7 ms보다 약 50–60배 무겁다. 다만 `dv=2 → 1 → 0.5 m/s`에서
현재 0.001-anchor linewidth가 `2.921 → 2.894 → 2.894 MHz`로 수렴해, 기존 `dv=2`는
현재 모델의 what-if 용도로는 약 1% 수준의 합리적인 절충이다.

### 4.2 보이지 않는 신호도 정상 linewidth/IF hero로 보고된다

`_readout()`은 `_eit_features()`가 계산한 contrast를 버리고 width만 사용하며
(`gabes/schemes/rydberg.py:443-446`), slope와 IF optimum을 전체 scan에서 찾는다
(`gabes/schemes/rydberg.py:446-459`). UI 범위 안의 80 °C에서 다음 결과가 나왔다.

- resonance transmission: `1.80×10⁻⁶`
- central contrast: `1.80×10⁻⁶`
- 보고 linewidth: `0.307 MHz`
- IF optimum detuning: `-8.96 MHz` (scan 경계 -9 MHz 바로 안쪽)

즉 사실상 불투명한 조건에서도 좁은 linewidth가 hero로 보이고, electrometry
discriminator는 중앙 Rydberg feature가 아니라 scan edge를 고른다. 2026-07-20에 추가된
`EIT status`는 contrast가 정확히 0 이하일 때만 작동한다
(`gabes/schemes/rydberg.py:421-437`, `488-498`;
`tests/test_rydberg_eit.py:123-133`).

이미 계산한 contrast를 metric으로 내고, IF/slope search를 central EIT/AT feature window로
제한하며, 최적점이 경계에 붙으면 `clipped/background-dominated` 상태를 내는 개선이
필요하다. 이는 O(Nscan) 후처리뿐이라 solver 차원과 solve 횟수를 전혀 늘리지 않는다.
shot-noise threshold를 임의로 고정하기 어렵다면 우선 contrast와 transmitted power를
정직하게 표시하고 hero 승격만 보류해도 충분하다.

### 4.3 AT peak 위치의 표시 정밀도가 grid보다 높다

AT peak는 local maximum의 grid 좌표를 그대로 쓴다
(`gabes/schemes/rydberg.py:394-399`, `463-471`). scan은 항상 801점이고 Rabi가 커질수록
window가 넓어져 점 간격이 커진다 (`gabes/schemes/rydberg.py:292-304`).
`Ω_RF/2π=20 MHz`에서는 801점 간격이 0.15 MHz인데 split을 소수 둘째 자리까지
`19.20 MHz`로 표시한다. 1601점에서는 `19.35 MHz`가 되어 한 grid step 차이가 났다.

두세 점의 quadratic peak interpolation 또는 grid-step에 맞춘 유효숫자/status 표시는
추가 solve 없이 가능하다. 특히 Doppler-on에서 scan 점수를 두 배로 늘리는 것보다
후처리 interpolation이 훨씬 싸다. EIT FWHM crossing에는 이미 선형 interpolation을
쓰므로 (`gabes/schemes/rydberg.py:412-437`) 같은 수치정책을 AT에도 적용하기 쉽다.

### 4.4 절대 electrometry calibration은 아직 아니다

UI와 사용자 가이드는 AT split이 RF electric-field ruler라고 설명한다
(`docs/GABES_User_Guide_v2.html:595-605`). 그러나 현재 입력은 이미 MHz 단위인
`lo_rabi_mhz`이며 (`gabes/schemes/rydberg.py:146-153`), 40D–39F dipole matrix element,
편광/각운동량 계수, RMS convention을 이용한 `E = ħΩ/μ34` 변환은 없다.
reference sensitivity/PSN 숫자는 의도적으로 internal constants와 테스트에만 있다
(`gabes/schemes/rydberg.py:92-94`; `tests/test_rydberg_eit.py:64-76`).

또한 논문의 40 kHz LO–SIG beat, balanced detector, quantum efficiency, photon shot noise,
technical noise와 time-domain transfer function은 정적 finite-difference slope가 대신한다
(`gabes/schemes/rydberg.py:8-12`, `447-506`). 그러므로 현재 적절한 등급은
**정적 spectrum/knob exploration용 semi-quantitative reference**이며, SI-traceable field
calibration이나 sensitivity budget reference가 아니다.

### 4.5 모델 범위와 검증 깊이

Zeeman sublevel, orthogonal linear polarization selection, optical pumping, stray E/B mixing,
Rydberg ionization/charge noise는 축약 4준위에 없다. probe/coupling Rabi와 transit factor는
논문 기준점에 맞춘 anchor다 (`gabes/schemes/rydberg.py:29-43`, `222-253`). 논문 원자료도
공개되어 있지 않아 현재 테스트는 Fig. 2 전체 curve를 독립 검증하지 않고, 기준 linewidth,
AT split ratio, compensated/uncompensated 방향성, probe-power monotonicity를 검증한다
(`tests/test_rydberg_eit.py:43-61`, `214-257`). `info()`의 “reproduces Fig. 2” 표현은
“reference point와 qualitative trend를 재현”으로 좁히는 편이 실험 reference로 더 정직하다
(`gabes/schemes/rydberg.py:653-669`).

## 5. 기존 개선안의 계산비용 평가

| 기존 항목 | 현재 상태 | 비용 | 물리 보존 판단 |
|---|---|---:|---|
| power/diameter → Rabi, transit, 두 dephasing | 완료 | 스칼라 전처리 | 차원 증가 없이 유용한 실제 knob 방향성을 보존한다. |
| AT center shift, finite-IF, temperature term | 완료 | O(Nscan) 후처리 또는 스칼라 | 거의 무시 가능한 비용이다. 단 IF는 detector chain이 아닌 proxy라고 유지해야 한다. |
| headless observables | 완료 | 물리비용 0, figure 회피 | 현재 측정에서 headless 약 0.52 ms, figure 약 82.6 ms였다. batch/report의 우선 경로가 맞다. |
| low-order polarization/Zeeman proxy | deferred GROUP B | 스칼라/convolution이면 거의 0 | full manifold보다 먼저 uncertainty/participation proxy로 넣을 수 있다. 다만 0.001 opacity anchor의 의미부터 공개해야 중복 보정이 없다. |
| full Zeeman/polarization/field/ionization | deferred GROUP C | 매우 큼 | n-level density matrix의 dense solve는 대략 n⁶로 증가한다. 절대 calibration에는 필요하지만 negligible-overhead 개선이 아니다. |
| full time-domain superheterodyne | deferred GROUP C | time axis와 noise model 추가 | 정적 spectrum 물리를 보존한 단순 patch가 아니라 별도 engine layer다. |

따라서 기존 deferred 제안 중 물리를 거의 공짜로 보강할 수 있는 것은 low-order proxy와
표시/불확도 계층이다. full Zeeman 및 time-domain chain은 현재 4준위 정적 엔진에 억지로
넣지 말고 별도 범위 합의가 필요하다는 기존 checklist 판단이 타당하다.

## 6. 동작을 바꾸지 않는 순수 코드 최적화 후보

1. **hand-written LU 주위의 BLAS thread context 제거/상시 controller 재사용**

   Rydberg affine path는 `with core.blas_single_thread()` 안에서 Numba의 hand-written
   real LU를 호출한다 (`gabes/schemes/rydberg.py:260-269`,
   `gabes/kernels.py:413-451`). inner solve는 BLAS를 호출하지 않는데 Windows의
   `threadpoolctl` library discovery가 Doppler-off profile에서 10회 중 약 44 ms를 썼다.
   null context 진단은 output max difference 0으로, median compute를
   `3.70 → 1.76 ms`로 줄였다. Doppler-on도 `221 → 206 ms`였다. 작은 basis transform의
   BLAS-thread 안정성을 여러 플랫폼에서 확인한 뒤 이 affine path만 context 밖으로 빼는
   것이 안전하다.

2. **velocity loop 밖으로 `base + s*A_coef` hoist**

   `_affine_scan_chi_real()`은 같은 scan point에서 velocity class마다
   `base + s*A_coef + kv*B_coef` 전체를 다시 조립한다
   (`gabes/kernels.py:426-436`). `base + s*A_coef`를 `j` loop 밖에서 한 번 만들고,
   각 velocity에서 `kv*B_coef`만 더한 진단 구현은 bit-identical이었고 Doppler-on median을
   약 `190 → 166 ms`로 줄였다(측정 변동이 있으므로 약 10%대 후보로 해석). Lambda의
   affine Doppler path에도 같이 이득이 난다.

3. **bounded `_atom(ground_deph, rf_deph)` cache**

   `_atom()`은 호출마다 4-level Lindblad와 `S_v`를 재구성한다
   (`gabes/schemes/rydberg.py:200-219`; `gabes/atoms.py:45-98`). 오늘 측정은 약
   0.30 ms/call이었다. 동일 dephasing으로 probe power만 바꾸는 16점 extra view에서는
   같은 atom을 반복 생성하므로 bounded cache 또는 cached immutable dissipator가
   동작을 바꾸지 않고 수 ms를 줄인다. 현재 `AtomModel` 배열이 mutable이므로 객체 자체를
   cache한다면 read-only 보장을 함께 두는 것이 안전하다.

4. **이미 있는 headless 경로 유지**

   `_readout()` median은 약 0.26 ms, `headless_observables()`는 약 0.52 ms,
   figure 포함은 약 82.6 ms였다 (`gabes/schemes/rydberg.py:439-606`). Fig. 2(b) extra-view는
   compute 약 181 ms, render 약 124 ms였다. 자동 sweep/report가 figure를 만들지 않게 하는
   현재 API가 여전히 가장 큰 무위험 체감 최적화다.

velocity-grid 배열 자체의 cache는 우선순위가 낮다. Doppler-on profile의 대부분은
`affine_scan_chi`의 velocity-class LU solve였고 (`gabes/schemes/rydberg.py:352-365`),
grid 생성은 병목이 아니었다.

부수적으로 `gabes/kernels.py:482`는 아직 “Rydberg: `B_coef=S_v=0`”이라고 쓰지만,
현재 Rydberg atom은 Doppler-on을 위해 nonzero `S_v`를 가진다
(`gabes/schemes/rydberg.py:68-78`; `tests/test_rydberg_eit.py:187-203`). 실행에는 영향 없는
stale comment이므로 다음 문서 정리 때 고치는 편이 좋다.

## 7. 검증

실행:

```bash
python -m pytest tests/test_rydberg_eit.py tests/test_kernels.py \
  tests/test_headless_observables.py tests/test_schemes_render.py -q
python -m pytest -q
```

결과:

- 관련 테스트: **56 passed in 18.26s**
- 전체 테스트: **189 passed in 58.72s**

이번 작업에서는 production code를 변경하지 않았고 이 보고서만 추가했다. 기존 untracked
일일 보고서, 첨부물, `tmp/` 등 사용자 작업은 건드리지 않았다.

## 8. 결론

Scheme 3은 실제 4준위 cascade OBE, power/beam/transit 연결, microwave AT dressing,
residual Doppler와 정적 readout을 갖춰 **기준점 주변의 Rydberg-EIT/AT 실험 설계용으로
유용한 semi-quantitative reference**다. 그러나 default의 `doppler=off`와 숨은
`line_strength=0.001`이 절대 opacity와 온도 scaling을 empirical fit으로 만들고,
현재 metric selection은 불투명/edge-dominated 조건에서도 false precision을 낸다.

가장 비용 대비 좋은 다음 순서는 다음과 같다.

1. 0.001 opacity anchor의 의미와 mode별 calibration boundary를 공개한다.
2. EIT contrast 및 edge/clipping status를 표시하고 중앙 feature에서만 IF/slope를 고른다.
3. AT peak를 sub-grid interpolation하고 표시 정밀도를 grid에 맞춘다.
4. 그 다음 low-order polarization/Zeeman uncertainty proxy를 추가한다.

이 네 항목은 solver 차원을 늘리지 않거나 후처리만 바꾸므로 물리를 더 정직하게 만들면서
계산 오버헤드를 사실상 무시할 수 있다. full Zeeman과 full superheterodyne은 절대
electrometry를 위해 중요하지만 별도 고비용 모델로 남겨 두는 것이 맞다.
