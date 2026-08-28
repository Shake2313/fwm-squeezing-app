# 2026-07-16 Scheme 2 물리 검토: Lambda coherence

## 오늘 선택과 현재 다섯 스킴

- 서울 현지 날짜의 day-of-month는 `16`이고, `n = (16 mod 5) + 1 = 2`이다.
- 드롭다운 등록 순서는 다음과 같다 (`gabes/schemes/__init__.py:19-24`, `README.md:8-16`).

  1. **OD / SAS** — `SASScheme()`
  2. **Lambda coherence (EIT / AT / CPT)** — `LambdaScheme()`
  3. **Rydberg-EIT electrometry** — `RydbergEITScheme()`
  4. **Hanle / EIA / NMOR** — `MagnetoScheme()`
  5. **FWM** — `FWMScheme()`

- `eit`, `at`, `cpt` 단일-regime 인스턴스는 테스트와 직접 호출용 alias일 뿐 드롭다운 순서를 바꾸지 않는다 (`gabes/schemes/__init__.py:27-39`). 따라서 오늘 대상은 두 번째 `LambdaScheme`이다.

## 검토 범위와 기존 제안 탐색

- 코드: `gabes/schemes/absorption.py:41-141, 418-752`, `gabes/atoms.py:21-98, 157-179`, `gabes/beam.py:9-55`, `gabes/kernels.py:413-490`, `gabes/observables.py:388-426`
- 문서: `README.md:8-20, 42-47, 67-80, 171-186`, `docs/GABES_User_Guide_v2.html:575-593, 789-803`, `static/GABES_User_Guide.html`의 동일 섹션
- 테스트: `tests/test_absorption.py:110-192`, `tests/test_kernels.py:102-136`, `tests/test_headless_observables.py:38-89`, `tests/test_schemes_render.py:31-42`
- 기존 개선안: `docs/checklist.json:21-33, 80-85, 102-113`과 Scheme 2 일일 보고서 `2026-06-26`, `2026-07-01`, `2026-07-06`
- 별도 Lambda TODO/issue/proposal 파일은 발견하지 못했다. 실제 제안의 기준은 `docs/checklist.json`과 위 세 일일 보고서다. 코드 TODO는 공통 Ne buffer-gas 계수 일반화 한 건이다 (`gabes/constants.py:79-89`).
- Lambda 전용 standalone 예제 스크립트는 없다. 현재 실행 가능한 예는 `tests/test_absorption.py:110-192`의 직접 API 호출과 앱/사용자 가이드의 EIT·AT 그림뿐이다. 즉 테스트가 사실상 유일한 코드 예제다.

외부 물리 기준도 최소한으로 대조했다. Rb D1의 weak-probe/strong-control Lambda EIT와 Doppler-free 구성이 실제 실험에서 질적으로 재현된다는 근거는 Li와 Xiao의 Rb 실험에 있고, warm-vapor 정량 모델에는 Raman dephasing과 transit/population exchange, velocity-changing collisions가 서로 다른 채널로 필요하다는 근거는 Ghosh 등의 연구에 있다.

- Y.-Q. Li and M. Xiao, *Phys. Rev. A* **51**, R2703 (1995), <https://doi.org/10.1103/PhysRevA.51.R2703>
- J. Ghosh et al., *Realistic theory of electromagnetically-induced transparency and slow light in a hot vapor of atoms undergoing collisions* (2009), <https://arxiv.org/abs/0901.3790>

## 결론 요약

`LambdaScheme`은 **실제 3준위 정상상태 OBE, Maxwell 속도 평균, 실험 단위의 coupling power/beam size, Beer-Lambert 전파, EIT 분산과 AT splitting**을 계산한다. 따라서 단순한 선 모양 장난감은 아니며, EIT/AT의 동작점 탐색과 정렬 민감도에 유용한 **semi-quantitative 실험 설계 reference**다.

다만 현재 구현을 Rb/Cs 특정 hyperfine Lambda 실험의 절대 reference로 쓰면 안 된다. 대칭 branching의 scalar 3-level 모델이고, 편광·Zeeman·hyperfine optical pumping·population exchange·transit/diffusion·VCC가 없다. `CPT`는 별도 bichromatic clock/CPT 모델이 아니라 같은 weak-probe Lambda 엔진의 좁은-scan preset이다. 또한 오늘 검토에서 **beam-angle residual Doppler가 `Doppler=off`에서도 스펙트럼을 바꾸는 저비용 수정 가능한 정확성 문제**를 확인했다.

| 용도 | 현재 적합성 | 판단 |
|---|---|---|
| EIT 창, control-power/waist 경향 | 좋음 | 실제 OBE와 warm-vapor 평균을 사용 |
| AT splitting 대 `Omega_c` | 좋음 | analytic scale과 테스트가 직접 고정 |
| mrad 정렬 민감도 | 보정 후 유용 | warm 경향은 타당하지만 cold-limit 위반이 있음 |
| CPT clock/자기장·편광·light-shift reference | 낮음 | weak probe, scalar 3-level, sublevel 없음 |
| 절대 transmission, linewidth, group delay | 제한적 | effective line strength/dephasing이며 외부 실험 anchor 없음 |

## 구현된 물리

### 1. 3준위 Lambda OBE

`atoms.lambda3()`는 `g1, g2, e` 세 준위를 만들고 excited state가 두 ground state로 `Gamma/2`씩 붕괴하도록 한다. `buffer_ground_relax_khz`는 실제 population relaxation이 아니라 `rho_g1g2`와 그 켤레에 직접 넣는 **Raman-coherence dephasing**이다 (`gabes/atoms.py:157-179`, `gabes/atoms.py:79-82`). 따라서 UI의 “ground relaxation”은 실험적으로는 `T2`형 유효값으로 읽어야 한다.

Hamiltonian은 weak probe와 control을 직접 포함한다. probe는 항상 `1e-3 Gamma`로 고정되고, `H[1,1] = Delta_c - s`, `H[0,2] = Omega_p/2`, `H[1,2] = Omega_c/2`이다 (`gabes/schemes/absorption.py:29-31, 649-656`). 이 구조는 선형 weak-probe EIT/AT에는 적합하지만, 두 광장의 세기가 비슷한 CPT optical-pumping 실험을 일반적으로 대표하지는 않는다.

### 2. 매질과 lab-facing control

`Rb (natural)`, `85Rb`, `87Rb`, `133Cs`, `Generic`과 D1/D2 선택으로 온도 의존 밀도, 자연 linewidth, 파장, 질량, reduced dipole을 바꾼다 (`gabes/schemes/absorption.py:418-443, 521-526`). 그러나 특정 `Fg -> Fe` leg, CG coefficient, Zeeman population은 선택하지 않으므로 이는 **scalar D-line medium**이다. `line_strength=1`의 절대 susceptibility와 group index는 특정 실험 transition에 자동 보정된 값이 아니다.

Coupling Rabi는 Advanced anchor를 기준으로

`Omega_c proportional to sqrt(P) / diameter`

로 변한다 (`gabes/schemes/absorption.py:529-548, 583-591`, `gabes/beam.py:9-20`). 1/e^2 diameter라는 정의도 UI help에 명시돼 있다. 절대 dipole/모드 overlap을 억지로 예측하지 않고 실험자가 측정한 Rabi를 anchor로 쓰는 방식이라, 실험 planning에는 정직하고 실용적이다.

### 3. Doppler 평균과 readout

Doppler-on Lambda는 4-sigma Maxwell grid를 `dv=2 m/s`로 만들고, 601 scan point와 속도 class를 affine Liouvillian kernel에서 푼다 (`gabes/schemes/absorption.py:118-140, 609-631`). Numba 경로는 real Hermitian basis, scan 병렬화, velocity-weighted coherence contraction을 사용한다 (`gabes/kernels.py:413-490`).

`chi_bar`는 density·dipole과 결합해 `alpha = k Im(chi)`가 되고, cell length는 solve 없이 `T=exp(-alpha L)`에만 들어간다 (`gabes/schemes/absorption.py:677-685`, `gabes/observables.py:388-410`). EIT/CPT의 group index는 `Re(chi)`의 수치 미분이다 (`gabes/schemes/absorption.py:740-751`, `gabes/observables.py:418-426`). 이는 느린 빛의 분산 지표이지, finite pulse 전파·왜곡·검출을 계산한 group-delay 실험 모델은 아니다.

## 직접 재계산 결과

현재 HEAD, Windows, `MPLBACKEND=Agg`에서 측정했다. 첫 EIT 호출은 Numba/parallel 초기화까지 포함해 `1.40 s`, 같은 process의 warm median은 `84 ms`였다. 수치는 hardware와 process 상태에 의존한다.

| Regime 기본값 | warm compute | headless | figure readout | 주요 결과 |
|---|---:|---:|---:|---|
| EIT | 84 ms | 0.353 ms | 203 ms | `T(res)=0.014`, FWHM `0.46 MHz`, `n_g=1.003e5` |
| AT | 11.8 ms | 0.147 ms | 206 ms | split `46.0 MHz`, `Omega_c=46.0 MHz`, `T(center)=0.989` |
| CPT | 10.2 ms | 0.315 ms | 176 ms | `T(res)=0.706`, FWHM `923.32 kHz`, `n_g=9.525e4` |

`temp=50 C`, `L=1 mm`, warm EIT에서 기존 각도 sweep도 재현됐다.

| angle | residual `|Delta k|/k` | `T(res)` | window FWHM |
|---:|---:|---:|---:|
| 0 mrad | `0` | `0.750724` | `0.45968 MHz` |
| 1 mrad | `9.9999996e-4` | `0.326469` | `0.91936 MHz` |
| 5 mrad | `4.9999948e-3` | `0.159116` | `4.13712 MHz` |
| 10 mrad | `9.9999583e-3` | `0.138894` | `7.81456 MHz` |

## 가장 중요한 발견: residual Doppler의 cold-limit 위반

기존 `lambda-residual-two-photon-doppler` 제안은 `done`으로 표시돼 있다 (`docs/checklist.json:80-85`). 아이디어 자체는 맞고 solver dimension도 늘리지 않는 좋은 저비용 물리 개선이다. 하지만 현재 전달 방식에는 작은 구조적 오류가 있다.

1. `beam_angle_mrad`는 `r = |Delta k|/k`를 만든다 (`gabes/schemes/absorption.py:593-597`, `gabes/beam.py:42-55`).
2. `atoms.lambda3()`는 이 `r`을 ground level 1의 `S_v` coefficient로 넣는다 (`gabes/atoms.py:165-178`).
3. 그런데 Lambda affine 경로는 `S_v`를 velocity-only `kv`뿐 아니라 scan coefficient에도 `A_coef = dL/ds - S_v`로 넣는다 (`gabes/schemes/absorption.py:45-57`). fallback도 `steady_state_batched(..., s-kv, S_v)`를 사용한다 (`gabes/schemes/absorption.py:136-140`).

따라서 원하는 residual term `Delta k dot v` 외에 `r*s`가 scan detuning에 섞인다. 물리적으로 beam-angle residual Doppler는 `v=0`에서 사라져야 하지만, `Doppler=off`에서 0 mrad와 10 mrad를 비교하자 다음 차이가 남았다.

- `max |Delta chi| / max |chi| = 2.94e-2`
- 최대 transmission 차이 `0.00564` (`-2.069 MHz`에서 `0.27903 -> 0.27339`)
- 공명점 transmission은 같았기 때문에 현재 테스트의 center-only 검증으로는 잡히지 않는다.

### 저비용 수정 권고

- scan detuning operator와 velocity operator를 분리해 `r`이 `kv`에만 곱해지게 한다. 예를 들어 affine kernel에 독립 `A_scan`, `B_velocity`를 전달하거나 Lambda `build_H0`에서 현재 중복 scan term을 정확히 상쇄할 수 있다.
- `doppler="off"`에서 angle 0/10 mrad의 전체 `chi_bar`가 동일해야 한다는 회귀 테스트를 추가한다.
- solver 차원과 velocity grid는 그대로이므로 runtime 증가는 사실상 없다. 오히려 coefficient 조립만 바로잡는 변경이다.
- 이 수정 전까지 warm angle sweep은 정성적 alignment sensitivity로만 사용하고, 1% 수준의 scan-axis/linewidth 정량값은 절대 기준으로 인용하지 않는 편이 안전하다.

## 기존 개선안과 계산 비용 평가

### 이미 반영된 항목

1. **Power/diameter-driven Rabi**: scalar 전처리만 추가하며 solver 차원은 같다. 물리를 잃지 않고 실험 knob 대응을 크게 개선했다 (`gabes/schemes/absorption.py:583-591`).
2. **Residual two-photon Doppler**: 속도 class를 원래 풀고 있으므로 올바르게 coefficient만 전달하면 추가 solve가 없다. 위 cold-limit 오류를 고친 뒤에는 negligible-overhead 개선이다.
3. **Headless observables**: `include_figures=False`로 Matplotlib을 건너뛴다 (`docs/checklist.json:29-33`, `tests/test_headless_observables.py:46-67`). 오늘 측정에서 `0.15-0.35 ms`로, figure `176-206 ms`보다 수백 배 가볍다. 물리는 전혀 바뀌지 않는다.

### 남아 있는 저차 buffer-gas 제안

`buffer-gas-pressure-shift`는 gas/species/line coefficient table, pressure shift, phenomenological Dicke narrowing을 제안한다 (`docs/checklist.json:21-26`, `gabes/constants.py:79-89`).

- pressure shift는 detuning scalar offset, homogeneous broadening은 optical `gamma` scalar라 matrix 차원 증가가 없다.
- phenomenological Dicke narrowing도 effective Doppler/residual-k coefficient로 제한하면 기존 affine 구조를 보존하며 추가 비용이 거의 없다.
- 단, coefficient table과 적용 convention은 실제 species/line별 문헌값과 uncertainty를 함께 기록해야 한다. 임의 fit factor로 숨기면 reference 신뢰도가 오히려 떨어진다.
- full VCC는 별개다. 속도 class 사이를 결합하므로 현재 독립 Maxwell 평균을 깨고 iterative/expanded solve가 필요하다. checklist의 GROUP C 분류가 타당하다 (`docs/checklist.json:102-106`). Ghosh 등의 연구도 VCC, transit, Raman dephasing, population exchange가 서로 대체 가능한 하나의 gamma가 아님을 보여 준다.

### Hyperfine/Zeeman-resolved Lambda

`lambda-hyperfine-resolved-manifold`는 실제 line assignment, branching, polarization, optical pumping을 설명하는 가장 큰 물리 개선이지만 GROUP C가 맞다 (`docs/checklist.json:109-113`). 현재 3-level의 9개 density-matrix 성분에서 실제 hyperfine/Zeeman manifold로 가면 Liouvillian 차원과 속도 평균 비용이 크게 증가하고, 기존 scalar preset의 의미도 바뀐다. 기본 interactive 모드에 무조건 넣기보다 별도 opt-in reference mode로 설계해야 한다.

## 순수 코드 최적화 후보

현재 최적화의 기준점은 이미 좋다. default warm EIT에서 affine Numba 경로 `0.090 s`, NumPy fallback `1.667 s`로 `18.47x`였고 `max relative chi difference = 3.63e-15`였다. 공식 parity test도 모든 EIT/AT/CPT Doppler on/off case에서 `<1e-11`을 요구한다 (`tests/test_kernels.py:116-136`).

1. **affine inner-loop matrix assembly 축소**: `_affine_scan_chi_real()`은 각 scan point에서 모든 velocity마다 `base + s*A + kv*B`의 전체 9x9를 다시 조립한 뒤 trace row를 즉시 덮어쓴다 (`gabes/kernels.py:426-449`). `base + s*A`를 scan당 한 번 만들고, 덮어쓸 trace row는 애초에 조립하지 않으면 LU 결과와 물리는 그대로이며 assembly 연산·메모리 쓰기를 줄일 수 있다. 실제 speedup은 별도 benchmark로 확인해야 한다.
2. **3-level AtomModel/template cache**: `_medium_from_params()`는 약 `52 us`, `atoms.lambda3()` 구성은 약 `569 us/call`이었다. warm AT/CPT가 약 `10-12 ms`이므로 immutable template/Lindblad cache는 반복 batch에서 수%를 줄일 여지가 있다 (`gabes/schemes/absorption.py:418-443, 633-646`, `gabes/atoms.py:45-98`). 배열 mutation을 막거나 복사 정책을 명시해야 한다.
3. **headless를 batch 기본 경로로 유지**: 이미 가장 큰 behavior-preserving 이득이다. 자동 보고서·parameter sweep에서는 `observables()` 대신 `headless_observables()`를 사용해야 한다.
4. **낮은 우선순위의 observable intermediate cache**: `alpha`, `chi_phys`, group-index gradient는 cell length와 독립이므로 raw에 보관할 수 있다 (`gabes/schemes/absorption.py:677-685, 740-742`). 하지만 현재 headless 전체가 `0.35 ms` 이하이므로 복잡성을 늘릴 만큼 큰 병목은 아니다.

과거 보고서의 `_medium_from_params()`/`atoms.lambda3()` cache 제안은 여전히 동작 보존형이지만, 오늘 측정상 EIT의 주병목은 affine solve이고 figure path의 주병목은 Matplotlib이다. 따라서 kernel assembly와 headless 사용이 우선이고, small-object cache는 batch 미세 최적화로 보는 것이 맞다.

## 문서·테스트·예제의 reference 품질

- 사용자 가이드는 scalar Lambda와 full polarization/Zeeman 부재를 정직하게 밝힌다 (`docs/GABES_User_Guide_v2.html:789-803`).
- 반면 Lambda 설명은 아직 coupling Rabi와 기본 EIT/AT 그림 위주이고, 새 power/diameter 및 beam-angle knob, `buffer_ground_relax_khz`가 순수 Raman T2라는 점, CPT가 weak-probe preset이라는 점은 충분히 설명하지 않는다 (`docs/GABES_User_Guide_v2.html:575-593`). 문서 보강은 runtime 0이며 실험 reference 오용을 줄인다.
- 테스트는 AT split, Rabi scaling, warm-angle broadening, cold EIT transparency, sub-natural CPT, kernel parity를 잘 고정한다 (`tests/test_absorption.py:110-192`, `tests/test_kernels.py:116-136`). 그러나 외부 실험 linewidth/contrast, species/D-line 변경, group index/group delay, cold-angle invariance는 고정하지 않는다.
- standalone example이 없으므로 `analysis/lambda_lab_sweep.py` 같은 짧은 headless 예제로 power/diameter, angle, ground-coherence dephasing sweep과 CSV 출력을 보여 주면 실험가에게 유용하다. 이는 물리·runtime 변경이 없는 문서성 개선이다.

## 검증

- `python -m pytest tests/test_absorption.py tests/test_kernels.py tests/test_headless_observables.py tests/test_schemes_render.py -q` -> `30 passed in 35.57 s`
- `python -m pytest -q` (`MPLBACKEND=Agg`) -> `118 passed in 25.60 s` (이번 재검증)
- default EIT/AT/CPT, warm angle sweep, cold angle invariance, Numba/NumPy parity 및 microbenchmark를 별도로 실행했다.

## 최종 판단과 우선순위

Scheme 2는 **EIT/AT의 핵심 간섭·dressing·Doppler 평균과 실험 knob scaling을 빠르게 보는 데 실제로 유용**하다. 특히 AT splitting, coupling `sqrt(P)/d` scaling, warm alignment sensitivity, headless batch 속도는 좋은 실험 준비 도구다.

우선순위는 다음이 적절하다.

1. **Residual Doppler에서 scan/velocity coefficient를 분리하고 cold-limit 회귀 테스트 추가** — 물리 정확성 개선, 사실상 무부하.
2. **가이드와 예제에 scalar/weak-probe/T2/absolute-scale 한계를 명시** — 무부하, reference 신뢰도 개선.
3. **문헌 계수 기반 pressure shift + homogeneous broadening + 제한적 Dicke proxy** — solver 차원 증가 없이 높은 warm-vapor 효용.
4. **affine kernel assembly 축소 및 immutable AtomModel template cache benchmark** — 동작 불변의 순수 성능 개선.
5. **full VCC와 hyperfine/Zeeman manifold** — 가치가 크지만 별도 설계와 opt-in heavy mode가 필요한 GROUP C.

따라서 현재 코드는 **실험 경향과 scale을 보는 semi-quantitative reference**로는 추천할 수 있지만, cold-limit residual Doppler 오류를 고치기 전 angle 정량값과, hyperfine/Zeeman·transport calibration이 없는 절대 linewidth/contrast/group-index를 논문급 reference로 인용해서는 안 된다.
