# 2026-07-12 Scheme 3 물리 검토: Rydberg-EIT electrometry

## 선택과 검토 범위

- 현지 날짜는 `2026-07-12`이며 `n = (12 mod 5) + 1 = 3`이다.
- 등록 순서는 `OD / SAS → Lambda coherence → Rydberg-EIT electrometry → Hanle / EIA / NMOR → FWM`이므로 세 번째 `RydbergEITScheme`을 선택했다 (`gabes/schemes/__init__.py:19-24`, `README.md:8-16`).
- 코드 본체와 설명은 `gabes/schemes/rydberg.py:1-620`, 전용 검증은 `tests/test_rydberg_eit.py:29-215`, affine solver 동등성 검증은 `tests/test_kernels.py:102-136`, 기존 개선안은 `docs/checklist.json` 및 2026-06-22/27, 2026-07-02의 scheme-3 보고서를 확인했다. 별도 Rydberg 사용 예제 파일은 없고, `extra_views()`의 Ju et al. Fig. 2(b) sweep가 실행 가능한 예제 역할을 한다 (`gabes/schemes/rydberg.py:558-601`).

## 현재 정의와 실제 물리

이 스킴은 85Rb `5S1/2 F=3 → 5P3/2 F'=4 → 40D5/2 → 39F7/2`를 축약한 정적 4준위 cascade OBE이다. 780 nm probe, 481 nm counter-propagating coupling, 37 GHz microwave를 Hamiltonian의 세 결합으로 넣고 (`gabes/schemes/rydberg.py:339-347`), steady-state optical coherence를 밀도·dipole·cell length와 결합하여 Beer–Lambert transmission으로 바꾼다 (`gabes/schemes/rydberg.py:307-393`, `402-411`).

실험 계획 도구로 유용한 부분은 다음과 같다.

- probe/coupling power와 beam diameter가 anchored `sqrt(P)/d` Rabi로 실제 solve에 들어가며, diameter는 동시에 thermal transit broadening을 정한다 (`gabes/schemes/rydberg.py:222-253`, `332-346`). 따라서 장비 knob 변화의 방향성을 직접 볼 수 있다.
- compensated EIT와 residual-Zeeman inhomogeneous convolution을 비교하고 (`gabes/schemes/rydberg.py:360-366`), RF AT splitting과 detuned microwave의 center shift를 추출한다 (`gabes/schemes/rydberg.py:462-480`).
- residual two-photon Doppler를 선택적으로 Maxwell 평균할 수 있고, 기본값은 논문 기준의 transit-limited 보정 모델을 유지한다 (`gabes/schemes/rydberg.py:349-358`).
- 유한 IF transmission 차분은 실제 lock-in 전체 모델은 아니지만 정적 discriminator의 operating detuning을 고르는 데 쓸 수 있다 (`gabes/schemes/rydberg.py:448-491`).

직접 재계산한 기본 EIT 결과는 linewidth `1.61 MHz`, 공진 transmission `0.940`, 최대 slope 및 IF discriminator `0.122 /MHz`였고, AT 기본값은 splitting `3.50 MHz`, center shift `+0.00 MHz`, 공진 transmission `0.806`이었다. probe `1 → 6 → 10 µW`에서 linewidth는 `1.417 → 1.613 → 1.754 MHz`; diameter `0.10 → 0.15 → 0.30 mm`에서 `2.674 → 1.613 → 0.680 MHz`로 변했다. power broadening과 작은 beam의 transit penalty라는 실험적 방향성이 일관된다.

따라서 현재 수준은 **정적 spectrum, AT splitting, beam/power/dephasing 민감도 및 operating-point 탐색을 위한 semi-quantitative 실험물리 reference**로 평가한다. 그러나 절대 전기장 감도나 장비 calibration reference로 쓰기에는 부족하다. Zeeman sublevel·편광 선택규칙·optical pumping·stray E/B mixing·ionization/charge noise가 없는 lumped ladder이고, detector noise와 time-domain lock-in/superheterodyne 전달함수도 없다 (`gabes/schemes/rydberg.py:4-12`, `603-619`; `docs/checklist.json`의 `rydberg-full-zeeman-polarization-field-model`, `full-time-domain-superhet-demodulation`).

## 기존 개선안과 계산 부하

이전 개선안은 존재한다. `rydberg-power-to-rabi`는 이미 완료되어 power/waist coupling, 두 dephasing channel, residual Doppler, AT center shift, finite-IF proxy와 temperature-linked dephasing까지 반영됐다 (`docs/checklist.json`, `gabes/schemes/rydberg.py:96-184`, `222-253`, `306-393`, `440-494`). 2026-07-02 보고서가 제안한 headless readout도 `headless_observables()`와 `include_figures=False` 경로로 반영되어 있다 (`gabes/schemes/rydberg.py:440-556`, `tests/test_headless_observables.py:32-64`).

부하는 다음처럼 구분된다.

- power-to-Rabi, transit/temperature dephasing, center shift, IF discriminator는 스칼라 전처리 또는 이미 구한 spectrum의 후처리여서 solver 차원과 scan 수를 늘리지 않는다. 물리적 유용성 대비 사실상 무시 가능한 비용이다.
- Doppler-on은 같은 4준위를 유지하지만 velocity class를 추가한다. 이번 JIT warm 상태 3회 측정에서 off는 `3.8–4.3 ms`, on은 `448–519 ms`로 약 100배 이상 느렸다. 절대 시간은 CPU 부하와 BLAS/JIT 상태에 민감하지만, velocity-class solve가 지배적이라는 결론은 분명하다. calibrated default가 아니라 opt-in what-if로 두는 판단은 여전히 맞다.
- full Zeeman/polarization/field 및 full time-domain detector 모델은 체크리스트의 GROUP C다. 전자는 상태수 증가로 density-matrix/Liouvillian 비용을 크게 키우고, 후자는 시간축·noise bandwidth·검출기 모델을 새로 요구한다. 작은 수정으로 가장할 수 없으며 목적과 검증 데이터 합의 후 별도 모델로 구현해야 한다.
- 저부하 대안은 체크리스트의 `low-order-polarization-zeeman-proxies`처럼 polarization purity, effective participating fraction, asymmetric broadening을 scalar proxy로 두는 것이다. solver 차원을 유지해 비용은 거의 없지만, 이를 원자상태 분해 결과가 아닌 phenomenological sensitivity band로 명시해야 한다.

## 순수 코딩 최적화 검토

- 현재 가장 큰 최적화는 이미 들어간 affine scan kernel이다. scan마다 Liouvillian을 Python에서 재조립하는 대신 affine coefficient와 연속 배열을 넘긴다 (`gabes/schemes/rydberg.py:255-275`, `tests/test_kernels.py:102-136`). 기능 보존 관점에서 우선 유지해야 한다.
- `_probe_line()`과 `_cascade_skeleton()`은 불변 원자 상수를 `lru_cache`한다 (`gabes/schemes/rydberg.py:46-78`). 추가로 Doppler velocity grid를 `(T, dv, cutoff, mass)` 키로 캐시하면 repeated sweep에서 grid 생성 비용을 줄일 수 있지만, 주비용은 batched solve라 효과는 제한적이다.
- `extra_views()`는 16개 power마다 compensated solve 하나만 수행한 뒤 uncompensated curve를 convolution으로 재사용한다 (`gabes/schemes/rydberg.py:561-574`). 합리적이다. batch/report는 반드시 headless 경로를 써 Matplotlib 생성을 피하는 것이 가장 안전한 체감 최적화다.
- 더 공격적인 후보는 여러 probe power의 affine matrices를 한 번에 묶는 batch axis 확장이나 Liouvillian dissipator template 재사용이다. 다만 4준위 warm solve가 이미 수 ms이고 복잡성·메모리 증가 대비 이득이 작다. 우선순위는 낮다.

## 검증 및 결론

이번 실행에서 수치 민감도 실험을 수행했고, 저장소 정책에 따른 전체 `python -m pytest -q`는 `118 passed in 47.68s`로 통과했다. 현재 스킴은 논문 기준 조건 주변의 정적 Rydberg-EIT/AT 실험 설계에는 충분히 유용하지만, 절대 electrometry 감도 예측은 내부 reference 상수와 정적 slope를 넘어서는 검출·noise·sublevel 물리가 필요하다. 다음 저비용 물리 개선으로는 full manifold보다 **명시적으로 phenomenological인 polarization/Zeeman proxy와 uncertainty band**가 비용 대비 가장 현실적이다.
