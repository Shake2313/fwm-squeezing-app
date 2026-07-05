# 2026-07-03 Scheme 4 Review: Hanle / EIA / NMOR

## 선택 결과

- 현재 등록 순서는 `OD / SAS -> Lambda coherence -> Rydberg-EIT electrometry -> Hanle / EIA / NMOR -> FWM`이다. 코드상 `_SCHEMES`는 `SASScheme()`, `LambdaScheme()`, `RydbergEITScheme()`, `MagnetoScheme()`, `FWMScheme()` 순서로 등록된다 (`gabes/schemes/__init__.py:19-24`, `README.md:15`).
- 오늘 현지 날짜는 `2026-07-03`이므로 `n = (3 mod 5) + 1 = 4`이고, 네 번째 스킴 `MagnetoScheme`을 검토했다 (`gabes/schemes/magneto.py:127`).

## 읽은 근거

- 스킴 등록과 순서: `gabes/schemes/__init__.py:14-24`
- Magneto 스킴 정의: `gabes/schemes/magneto.py:127-187`
- 파라미터 및 regime: `gabes/schemes/magneto.py:187-327`
- 본 계산 경로: `gabes/schemes/magneto.py:378-488`
- paraffin two-region / buffer solver: `gabes/schemes/magneto.py:521-587`
- 관측량 및 headless figure 분리: `gabes/schemes/magneto.py:608-705`
- Zeeman manifold, CG coupling, grouped emission: `gabes/zeeman.py:70-139`
- numba magneto kernels: `gabes/kernels.py:199-238`
- 테스트: `tests/test_magneto.py:123-234`, `tests/test_headless_observables.py:38-87`, `tests/test_kernels.py:63-96`
- 기존 개선안: `docs/checklist.json:21-39`, `docs/checklist.json:77`
- FWM-Hanle 연계 분석: `analysis/resonant_hanle_squeezing_reference.py:396`, `analysis/resonant_hanle_squeezing_reference.py:829`, `analysis/resonant_hanle_squeezing_reference.py:900`, `analysis/resonant_hanle_squeezing_reference.md:60`
- 이전 Magneto 리뷰: `docs/daily_report/2026-06-23_scheme-4_hanle-eia-nmor.md`, `docs/daily_report/2026-06-28_scheme-4_hanle-eia-nmor.md`

## 현재 구현의 물리적 내용

`MagnetoScheme`은 단순한 Lorentzian toy model이 아니라 `87Rb D1`의 근영자기장 Zeeman manifold를 직접 세우고, `sigma+`, `sigma-`, `pi` 전이의 Clebsch-Gordan coupling과 spontaneous-emission branching을 포함한다 (`gabes/zeeman.py:70-139`). 특히 emission을 polarization별 grouped jump operator로 묶어 transfer of coherence 경로를 보존하기 때문에, `Fe = Fg + 1` 전이에서 intrinsic EIA가 생기는지 테스트로 확인된다 (`tests/test_magneto.py:179-192`).

실험물리학 관점에서 좋은 점은 세 가지다.

1. `qwp_deg`, residual transverse field, transition 선택이 Hanle dip, EIA/MIA peak, circular-light LCA의 부호와 형태를 실제 line-shape 결과로 바꾼다 (`gabes/schemes/magneto.py:113-124`, `gabes/schemes/magneto.py:426-438`, `tests/test_magneto.py:123-138`, `tests/test_magneto.py:195-211`).
2. paraffin cell은 light region과 dark region을 분리한 two-region exchange model로 구현되어, wall coherence가 좁은 Ramsey feature를 만드는 방향성을 표현한다 (`gabes/schemes/magneto.py:401-416`, `gabes/schemes/magneto.py:521-565`, `tests/test_magneto.py:141-152`).
3. NMOR는 transmission을 재포장한 값이 아니라 `chi_+ - chi_-`의 실수부 차이에서 rotation을 만들고, zero-crossing slope를 readout metric으로 준다 (`gabes/schemes/magneto.py:653-679`, `tests/test_magneto.py:214-225`).

따라서 이 스킴은 Hanle/EIA/NMOR 실험을 준비하는 연구자가 "이 조건이면 dip인가 peak인가", "transverse residual field가 LCA를 살리는가", "wall lifetime을 늘리면 중앙 feature가 좁아지는가", "최적 bias field는 대략 어디인가"를 빠르게 판단하는 semi-quantitative reference로는 꽤 유용하다.

## 오늘 재현 실험

짧은 로컬 실험은 `scan_points=121`, `velocity_classes=5`, `doppler=on`으로 돌렸다. `numba` kernel은 사용 가능했다 (`gabes/kernels.py:43-45`, `gabes/kernels.py:199-238`).

| 조건 | compute | headless observables | figure observables | 핵심 결과 |
|---|---:|---:|---:|---|
| paraffin, linear QWP | 0.465 s | 0.0005 s | 0.851 s | EIT-like dip, feature amp `-2.01e-2`, FWHM `0.197 uT` |
| paraffin, circular QWP | 0.129 s | 0.0004 s | 0.126 s | EIA/MIA-like peak, feature amp `+7.26e-2`, FWHM `0.216 uT` |
| NMOR | 0.125 s | 0.0006 s | 0.110 s | `B=0` rotation `0`, slope `-1.02 mrad/uT` |
| buffer 20 Torr | 0.064 s | 0.0003 s | 0.108 s | broad EIT-like dip, FWHM `7.405 uT`, `T(0)=0.761` |

첫 paraffin 계산은 kernel/cache warm-up 성격이 섞여 있고, 이후 같은 규모의 계산은 대략 0.06-0.13초 범위였다. Headless 관측량 경로는 이미 매우 가벼워서 batch/report 작업에서 Matplotlib 비용을 거의 제거한다 (`gabes/schemes/magneto.py:132`, `gabes/schemes/magneto.py:608-681`, `tests/test_headless_observables.py:38-87`).

buffer pressure만 바꾸고 `buffer_ground_relax_khz=20`, `collisional_depol_khz=2`를 고정하면 다음처럼 optical broadening만 직접 바뀐다.

| Ne pressure | optical broadening | feature amp |
|---:|---:|---:|
| 0 Torr | 0.0 MHz | `-1.99e1` |
| 20 Torr | 78.2 MHz | `-4.41e0` |
| 80 Torr | 312.8 MHz | `-8.29e-1` |

이 결과는 현재 코드의 한계를 잘 보여준다. `ne_pressure_torr`는 `constants.neon_buffer_broadening()`을 통해 optical homogeneous width에 들어가지만 (`gabes/constants.py:81-89`, `gabes/schemes/magneto.py:396`), diffusion-limited ground relaxation이나 spin depolarization은 여전히 별도 override knob이다 (`gabes/schemes/magneto.py:233-255`).

## 기존 개선안 반영 여부와 부하 평가

이전 Magneto 리뷰에서 제안된 `b_offset_ut`는 이미 구현되어 있다. UI 파라미터로 노출되고 (`gabes/schemes/magneto.py:220-223`), 실제 계산 축은 `b_physical_ut = b_ut + b_offset_ut`로 이동한다 (`gabes/schemes/magneto.py:421-423`). 연산 부하는 사실상 0에 가깝고, 실험적으로는 shielding/coil offset을 모델에 직접 맞출 수 있어 유용하다.

`figureless-observables-paths`도 완료 상태다 (`docs/checklist.json:28-32`). Magneto는 `supports_headless_observables=True`이고 `observables(..., include_figures=False)` 경로를 지원한다 (`gabes/schemes/magneto.py:132`, `gabes/schemes/magneto.py:608`). 물리를 바꾸지 않으면서 batch sweep과 자동 보고서 생성의 wall-clock을 크게 줄이는 좋은 코딩 개선이다.

아직 남은 주요 개선안은 두 가지다.

1. `buffer-gas-pressure-shift`: Ne broadening을 gas/species/line table로 확장하고 pressure shift와 phenomenological Dicke narrowing을 추가하자는 항목이다 (`docs/checklist.json:21-25`, `gabes/constants.py:81-89`). pressure shift는 detuning offset에 scalar를 더하는 수준이라 부하가 거의 없다. 저차 Dicke narrowing도 `gamma_opt` 또는 Doppler scalar correction에 넣으면 solver 차원을 키우지 않는다. 반면 full velocity-changing collision은 velocity class를 서로 결합하므로 현재의 separable `(B x velocity)` batched solve 구조를 깨고, interactive 성능을 크게 해칠 가능성이 높다 (`docs/checklist.json:77`).
2. `magneto-buffer-relaxation-map`: `ne_pressure_torr`에서 기본 `buffer_ground_relax_khz`, `collisional_depol_khz`를 경험식 또는 table로 채우자는 항목이다 (`docs/checklist.json:35-39`). 이는 heavy solve 앞단의 scalar default mapping이라 연산 부하는 거의 없다. 실험물리학적으로는 "20 Torr를 넣었는데 왜 ground relaxation은 사용자가 따로 골라야 하는가"라는 현재의 해석 부담을 줄인다.

## 실험물리학 레퍼런스로서의 등급

현재 구현은 "현상 분류와 파라미터 감도 탐색용 레퍼런스"로는 충분히 쓸 만하다. 특히 paraffin Hanle/NMOR, polarization sign switch, TOC 기반 EIA, circular-light LCA는 테스트와 코드 구조 모두에서 물리적 의도가 선명하다.

다만 "절대 magnetometer calibration reference"로는 아직 조심해야 한다. buffer-gas diffusion, pressure shift, Dicke narrowing, pressure-dependent ground relaxation이 아직 calibration 없이는 phenomenological하다는 점은 새 FWM-Hanle 연계 분석 문서도 명시한다 (`analysis/resonant_hanle_squeezing_reference.md:60`). 또한 paraffin cell 쪽도 `wall_coherence_ms`, `transit_relax_khz`, `dark_return_khz`가 material/geometry에서 자동 산출되는 예측식이 아니라 사용자가 조정하는 effective parameter다 (`gabes/schemes/magneto.py:238-260`, `gabes/schemes/magneto.py:401-416`).

따라서 연구자가 이 코드를 논문 수치의 직접 대체물로 쓰기보다는, 실험 설정을 설계하고 measured Hanle trace에 맞춰 offset/scale/relaxation을 보정하는 전단계 모델로 쓰는 것이 적절하다. `analysis/resonant_hanle_squeezing_reference.py`의 Hanle calibration fit 경로는 이 사용법과 잘 맞는다 (`analysis/resonant_hanle_squeezing_reference.py:396`, `analysis/resonant_hanle_squeezing_reference.py:829`, `analysis/resonant_hanle_squeezing_reference.py:900`).

## 순수 코딩 최적화 후보

- 이미 반영된 headless observables는 가장 효과적인 개선이다. 오늘 측정에서 headless readout은 0.3-0.6 ms 수준이고 figure 생성은 0.1초 이상 걸렸다.
- `_hamiltonian()`은 호출 때마다 angular momentum matrix를 만든다 (`gabes/schemes/magneto.py:493-509`). `(Fg, Fe)` 또는 block dimension별로 `Fx, Fy, Fz`를 캐시하면 물리를 건드리지 않고 작은 반복 비용을 줄일 수 있다.
- `zeeman.zeeman_manifold()`의 coupling topology와 emission skeleton은 relaxation rate와 부분적으로 분리 가능하다 (`gabes/zeeman.py:70-139`). template cache를 두면 repeated parameter sweep에서 Python object construction 비용을 줄일 수 있다.
- 그러나 큰 비용은 결국 paraffin two-region의 `2M x 2M` batched solve다 (`gabes/schemes/magneto.py:521-565`). 따라서 미세 캐시보다 `scan_points`, `velocity_classes`, figure 생성 여부, extra-view/batch 경로 분리가 실제 체감 성능에 더 중요하다.

## 검증

`python -m pytest tests/test_magneto.py tests/test_headless_observables.py tests/test_kernels.py -q`

결과: `32 passed in 9.92s`

## 최종 판단

`MagnetoScheme`은 오늘 기준으로도 GABES 안에서 물리 구현 밀도가 높은 스킴이다. 실험가가 Hanle/EIA/NMOR의 부호, 폭, offset, polarization sensitivity를 빠르게 확인하는 데는 실제로 쓸 수 있다. 남은 가장 값싼 물리 개선은 pressure shift, low-order Dicke narrowing, pressure-to-relaxation default mapping이며, 이들은 solver 차원을 늘리지 않고도 buffer-cell 해석력을 크게 올릴 수 있다. full velocity-changing collision은 물리적으로 매력적이지만 현재 interactive 목표에는 과한 GROUP C 작업으로 남겨두는 편이 맞다.
