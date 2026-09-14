# Ultra 성능 패치 — 정확한 Floquet 대칭의 재검증

2026-09-11. 검토한 [외부 제안](https://claude.ai/code/artifact/bf263a2e-0d04-4705-90b6-a9d14eca192a)의 항목 중 수학적으로 동일한 반복 계산을 제거했다. 변경 경로는 `gabes/kernels.py`의 `floquet_chi_grid`이며, 실제 앱의 `FWMScheme.compute`로 속도와 출력 배열을 비교한다. Ultra의 physics나 검증 등급을 높였다는 의미는 아니다.

## 동일성의 조건과 증명

Row-major vec에서 J(X)=P X* P라 정의하자. P는 density matrix의 두 Hilbert-space index를 교환하는 permutation이고, *는 complex conjugate다. 입력이

\[
J(L_0)=L_0,\quad J(C_\delta)=C_\delta,\quad J(S_v)=S_v,
\qquad C_-=J(C_+)
\]

를 만족하고 δ, Δ_eff, Ω_hf와 branch가 실수이면 모든 격자점의 L=L₀+δCδ−Δ_eff S_v가 Hermiticity를 보존한다. 새로운 wrapper는 이 coefficient identities를 **raw vec 기저에서 정확한 배열 비교**로 검사한다. 허용 오차 안의 작은 물리 항을 버리는 검사가 아니다. 한 항이라도 실패하면 기존 두 사슬 general kernel로 돌아간다. Complex frequency도 general 경로를 사용한다.

Hermitian operator basis의 열로 이루어진 unitary U를 사용하면 Lʳ=U†LU는 실수이고 C₋ʳ=(C₊ʳ)*다. 여기의 *는 transpose가 없는 원소별 켤레다. Raw vec 기저에서 P를 생략하는 잘못된 식은 테스트가 검출한다. 원자 상태의 h≠0 Fourier 성분은 여전히 복소수다.

같은 symmetric finite cutoff −N,…,+N에서 positive recurrence를

\[
R_h=-\left[L^r+ih\Omega+C_-^r R_{h+1}\right]^{-1}C_+^r,
\qquad R_{N+1}=0
\]

로 정의한다. 기존 negative recurrence는 Q₋ₙ₋₁=0에서 시작한다. Ω가 실수이므로 위 식에 켤레를 취한 결과가 Q₋ₕ의 recurrence와 같다. 따라서 **boundary에서의 귀납법으로 Q₋ₕ=Rₕ***다. 이는 N=1뿐 아니라 **모든 유한 N**에서 성립한다. 기존과 동일하게 recurrence block의 역행렬과 trace-normalized 해가 존재하는 영역을 전제한다.

Zero harmonic의 Schur complement는

\[
L_{\rm eff}^r=L^r+C_+^r Q_{-1}+C_-^r R_1
             =L^r+2\operatorname{Re}(C_-^rR_1)
\]

이므로 실수 LU로 풀 수 있다. U의 첫 n열이 population projectors이고 첫 열이 |0⟩⟨0|이므로 row 0를 population 합으로 바꾸는 trace constraint도 기존 raw solve와 정확히 대응한다. Polarization weights는 wʳ=wU로 바꾸고, ρ₊₁ʳ=R₁ρ₀ʳ를 수축한다.

**결론:** 동일한 N, 동일한 입력/격자, 동일한 finite-Floquet 연립방정식을 기저만 바꾸어 푼다. Weak-seed, small pump, adiabatic elimination 또는 N→∞ 수렴 가정이 아니다. 부동소수점 연산 순서와 pivot 기저의 차이로 bitwise equality를 요구하지는 않는다.

## 이후 수정자가 지켜야 할 경계

Hot loop에 위 귀납 근거와 이 문서 링크를 남겼다. 조건이 계속 성립한다면 음의 사슬을 다시 LU로 푸는 것은 동일한 해를 중복 계산하는 것이다. 엄밀성 검증에는 유지한 general two-chain kernel과 독립 dense-block reference를 사용한다. 모델이 변하면 먼저 wrapper의 가정이 여전히 성립하는지 판정한다.

- `_floquet_chi_grid_general`은 패치 전 `floquet_chi_grid`의 본체를 보존한다. Benchmark는 commit `c0c46f0`의 함수와 AST를 비교해 이름 외에 변한 것이 없는지 확인한다.
- `tests/test_floquet_symmetry.py`는 2/3/4-level, N=1/2/3/4/5, 양 branch, 복소 drive와 임의 readout weights를 독립 dense solver와 대조한다.
- 계수의 한 ULP 대칭 위반과 complex frequency가 fast path로 들어가지 않는지 검사한다.
- LU 호출 수를 계측해 각 격자점에서 **N개의 complex LU + 1개의 real LU**만 수행하는 것을 검증한다. 기존 N=3은 7개의 complex LU였다.
- 기존 FWM regression, compiled/NumPy parity 및 전체 order-gate 검사를 유지한다. 유한 probe/conjugate reference states 두 개는 서로 다른 물리 입력이므로 계속 각각 푼다.

## 외부 제안 중 다른 항목의 판정

| 제안 | 이번 판정 |
|---|---|
| ±harmonic 대칭과 real zero block | 위 조건·증명·독립 검산을 거쳐 적용 |
| adjacent-order 감사를 δ 행 10%에만 수행 | 전체 스캔 certificate와 동일하지 않음. 생략한 행의 불일치를 놓치는 반례를 보고서에 저장. 전체 검사를 유지 |
| N=3을 N=2로 변경 | 서로 다른 finite-boundary 문제. 특정 점에서 작은 차이는 대수적 동일성이 아님. N 유지 |
| 128-node/5σ 또는 composite velocity quadrature | 기존 grid·cutoff와 다른 근사. 128/256 및 5/6σ를 직접 재계산한 결과를 기록. Composite panels와 영역 전체 오차 검증이 없는 상태에서 기본값을 교체하지 않음 |
| Δ ladder를 k·dv의 정수배로 정렬 | 정렬된 새 입력에서 열 재사용은 가능하지만 기존 Δ값이 약 0.6% 이동. 현재 앱의 사용자가 정한 detuning을 바꾸는 속도 패치로 적용하지 않음 |
| operating_point의 cubic readout | 별도 보간 정확도/overshoot 검증이 필요한 출력 변경. 이번 계산량 절감과 분리 |
| probe/conjugate 두 seed solve 통합 | 두 finite-drive Hamiltonian이 다르므로 이 대칭으로 합칠 수 없음. Ω_ref 독립성을 잘못 주장한 주석만 수정 |

Quadrature 진단은 nonuniform nodes에서 χ를 **직접** 계산한다. 기존 uniform-grid interpolation formula에 Legendre nodes를 넣지 않는다. 평가한 S_dB는 현재 Ultra의 gain-referred indicator이며 microscopic intensity-noise spectrum이 아니다. 진단의 특정 지점 PASS를 전체 domain의 증명으로 확대하지 않는다.

## 재현 방법

```powershell
python -m analysis.ultra_performance.verify --output NEW.json
python -m pytest -q tests/test_floquet_symmetry.py tests/test_fwm_floquet.py tests/test_kernels.py tests/test_regression.py
python -m pytest -q
```

보고서는 새 파일에만 저장한다. `analysis/ultra_performance/ultra_app_before.json`과 `.npz`는 실제 패치 전 앱의 parameters, timings, 전체 array outputs와 kernel hash를 담는다. 추가 `ultra_before.*`는 넓은 기본 API scan의 별도 기록이며 앱 시간으로 혼용하지 않는다.

성능 측정은 두 compiled path를 예열하고 동일한 앱 설정으로 general→fast / fast→general 순서를 번갈아 실행한다. Numba 16 threads, BLAS 1 thread에서 JIT 시간은 제외하고 측정 중 다른 test/benchmark를 동시에 실행하지 않는다. 결과에는 pre-patch arrays와의 비교, 복소 transfer를 포함한 전체 array 오차, full-scan order-gate 및 stress points를 남긴다.

## 측정 결과

[검증 기록](../analysis/ultra_performance/verification.json)의 대조 실험에서 실제 앱 Ultra는 **13.3155초 → 6.5186초, 2.043배** 빨라졌다. 각 경로 2회 중앙값이며 JIT 예열 이후의 시간이다. 이 수치는 이 환경의 측정값으로, 다른 장치나 설정에서 같은 배율을 보장하지 않는다.

| 항목 | 기존 general | 대칭 적용 |
|---|---:|---:|
| 실행 시간, 1차 | 13.3028 s | 6.6231 s |
| 실행 시간, 2차 | 13.3282 s | 6.4142 s |
| detuning / velocity 점 수 | 401 / 1,573 | 401 / 1,573 |
| Floquet 차수 / 비교 차수 | 3 / 2 | 3 / 2 |
| 차수 검사 범위 | 전체 401점 | 전체 401점 |
| 차수 검사 판정 | CONVERGED | CONVERGED |

미리 저장한 패치 전 앱 배열과 비교하면 G_s/G_c의 최대 절대 오차는 각각 9.87×10⁻¹¹ / 8.94×10⁻¹¹이다. 각 전체 배열의 최대 절대값으로 정규화한 오차는 2.65×10⁻¹⁴ / 2.37×10⁻¹⁴이며, small-signal gain/complex transfer도 3.25×10⁻¹⁴ 이하이다. S_dB 최대 절대 차이는 4.03×10⁻¹² dB이다. 이는 부동소수점 기저 변환에 따른 수치적 차이의 측정이며 수학적 동일성 자체의 근거는 위 증명이다.

추가로 반대 branch, 2 W pump, Δ=−2.1 GHz 및 near-resonant Δ=−0.2 GHz를 두 구현으로 다시 계산했다. 네 경우 모두 출력 비교 기준과 동일한 차수 검사 판정을 통과했다. 반대 branch 검사는 기존 구현과의 일치만 검증하며 그 branch의 물리적 타당성을 높이지 않는다.

별도 quadrature 재검증에서는 Δ=−0.2 GHz, 100 °C, 200 mW의 δ=−8 MHz 점에서 128→256 nodes(5σ)가 현재 gain-referred indicator를 1.878 dB 바꿨고, 256 nodes에서 5→6σ도 2.363 dB 바꿨다. 세 detuning 사례 모두 전체 81점 스캔에서는 인접 quadrature 간 최대 차이가 0.01 dB를 넘었다. 높은 node 수를 정답으로 선언하지 않으며, 이 결과는 plain Legendre 기본값을 바로 채택할 근거가 부족함을 보여 준다. FWM susceptibility 적분만 교체해 진단했고 별도 absorption 계산은 원래 설정을 유지했다.

제안의 약 21배 추정에는 sampled audit와 다른 격자 근사가 함께 포함되어 있다. 이번 **동일 문제의 정확한 계산량 절감**에 대해 검증된 실측치는 위 2.043배다.

## 테스트 결과

- 신규 대칭/guard/dense-reference/LU-count 검사: **14 passed**.
- 기존 FWM Floquet, compiled kernel 및 기준 배열 regression 검사: **36 passed**.
- `GABES_DISABLE_NUMBA=1`의 신규 검사: **13 passed, 1 skipped**. Skip은 Numba 전용 LU 호출 계측이다.
- 최종 `python -m pytest -q`: **778 passed, 1 failed**, 169.59초. 유일한 실패는 `tests/test_docs_consistency.py::test_repository_visibility_wording_is_consistently_public`이며, 작업 전부터 누락된 `docs/FWM physics and analytic reconstruction/FWM_physics.tex`를 읽으려다 발생한 `FileNotFoundError`다. 이번 변경으로 복원하거나 검사를 생략하지 않았다.
- 검증 보고서의 현재 source hashes와 패치 전 snapshot hashes, checklist JSON 및 근거 파일 링크를 확인했다. 수치 비교와 제안 판정의 `expected_controls_passed`는 true다.
