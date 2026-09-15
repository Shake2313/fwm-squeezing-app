# FWM Squeezing Fast/Balanced 리마스터 — 개발 기록

2026-09-11. 범위는 FWM scheme **Squeezing 모드의 Fast/Balanced** 두 tier다. Ultra는 손대지 않았고
(출력 bit-identical을 게이트 G3.7로 확인), Grand Challenge 작업과 독립이다. 완료 기준은
[GATE.md](GATE.md)에 구현 전 동결했고, 판정은 `python -m analysis.fwm_lite.gate`의 보고서로만 한다.

**측정 환경.** 탐색·게이트 측정은 커밋 `c0c46f0` 위의 공유 작업트리에서 실행했다. 그 트리에는 다른 세션의
미커밋 변경(예: `kernels.py`의 Hermitian 대칭 Floquet 커널)이 함께 있었다. 커밋에는 이 패치만 담았고, 커밋 단독
트리의 검증은 맨 아래 "커밋 단독 검증"에 적는다.

**한 줄 요약.** Fast/Balanced는 이제 Ultra와 *같은 모델·같은 Maxwell 측도·같은 readout*을 계산한다.
속도 클래스마다 Floquet 방정식을 푸는 대신, 고정 δ에서 응답이 Δ_eff의 **정확한 유리함수(부분분수)** 라는
사실을 이용해 δ당 고유분해 한 번으로 Maxwell 합을 끝낸다. Balanced는 401개 δ를 모두 이렇게 풀고,
Fast는 δ를 적응적으로 골라 풀고 나머지는 스플라인으로 채운다.

---

## 0. 출발점 (변경 전, 이 PC, numba 16 threads, warm)

| tier | 시간 | δ 점 / 속도 노드 | S @ δ=−8 MHz | G_s peak | Ultra 대비 S 최대 오차 |
|---|---:|---|---:|---:|---:|
| Fast | 0.58 s | 181 / 297 (4 m/s, 3σ) | −8.654 dB | 4016 | 3.05 dB |
| Balanced | 1.87 s | 301 / 591 (2 m/s, 3σ) | −8.660 dB | 4016 | 1.50 dB |
| Ultra | 6.2–6.7 s | 401 / 1573 (1 m/s, 4σ) | −7.800 dB | 3724 | — |

증거: `before/tier_timings_before.json`, `before/ultra_default_before.npz`, `before/pytest_before_summary.txt`
(856 passed, 1 failed — 작업 전부터 누락된 `FWM_physics.tex` 때문).

관찰 세 가지가 설계를 정했다.

1. **개형 차이의 주원인은 격자가 아니라 readout 물리.** 기존 Fast/Balanced는 `PHASE_BALANCED/FINE`
   (1/16 segment, pump scatter·depletion·overlap 없음)이라 dip 전체가 Ultra보다 0.85 dB 깊었다.
   → 두 tier 모두 Ultra readout(`PHASE_ULTRA`)을 그대로 쓴다.
2. **Ultra 시간의 96%는 χ̄ 표**(δ × 속도 클래스 × 2 seed × N_F 3·2). readout 0.52 s 중 0.36 s는
   진단용 흡수표(`_hyperfine_alpha` 3회), 0.045 s는 매번 다시 만드는 24-level Zeeman 상수.
3. **3σ 속도 절단 자체가 2–20 % χ̄ 오차**(1 m/s 표 오프라인 분석). 기존 Fast가 빨랐던 이유 일부는 틀렸기 때문.

## 1. 탐색 순서와 근거

모든 실험은 12개 동작점(`gate.py`의 `CASES`: default, v6 최적점, frontier, 근공명, 청색 근/원, −3 GHz,
60 °C/50 mW, 150 °C, 135 mW, 200 µW seed, 50 mm cell)을 Ultra 기준으로 쟀다. 개형 지표는 |δ|≤500 MHz 창의
S_dB tube(가로 ±1 MHz 여유), p95, log10 G_s tube, 동작점 차이다.

### 1.1 속도 구적: Gauss 계열은 실패, 균일 사다리꼴이 최적

1 m/s·±6σ χ̄ 표에서 후보 구적을 비교했다(최대 상대오차, 4σ 기준 측도).

| 후보 | 노드 | default | near_res |
|---|---:|---:|---:|
| 균일 4 m/s | 393 | 1.2e-3 | 3.0e-2 |
| 균일 8 m/s | 197 | 2.5e-2 | 2.6e-1 |
| Gauss–Legendre 128 @4σ | 128 | 3.6e-1 | 9.8e-1 |
| Gauss–Hermite 16–96 | ≤96 | 0.8–2.4 | 3.6–8.3 |
| composite GL 0.25σ×4 | 128 | 2.0e-1 | 7.7e-1 |

χ̄(Δ_eff)에는 FWHM 15–60 MHz(근공명에서 ~Γ)의 공명이 있고, Gauss 규칙의 중앙 노드 간격(~24 MHz)이
그보다 넓다. 균일 사다리꼴은 Lorentzian 폭 γ에 대해 오차가 **e^{−2πγ/h}** 로 줄어드는 초수렴을 한다
(Trefethen–Weideman). 페르미식으로 γ ≥ Γ/2 ≈ 3 MHz를 넣으면 h ≈ 2πγ/ln(1/ε) → ε=1e-3에서 h≈2.7 MHz(2 m/s):
Balanced의 옛 2 m/s가 바로 이 값이다. **결론: 속도 노드 수는 줄일 여지가 거의 없다.**

### 1.2 δ 축: 균일 솎아내기는 실패, 적응 이분법은 성공

Ultra의 401점 평균 χ̄를 솎아 보간한 뒤 **Ultra readout을 그대로 통과**시켜 비교했다.

- 균일 stride 2(5.5 MHz 간격) cubic: v6_opt 3.84 dB, low_pump 7.42 dB 오차.
  원인은 폭 ~5 MHz의 Raman 공명(G_s가 한 격자에서 1→30)으로, Ultra 401점도 겨우 표본화한다.
- **적응 이분법**(stride 32에서 시작, 구간 중점을 예측→정확히 풀기→관측량으로 비교, 틀린 구간만 분할):
  27–121 노드로 전 케이스 tube ≤ 0.04 dB.
- 보간 대상은 **선형응답량 χ̄**(매끄러움)이고, 판정은 **사용자가 보는 S_dB·log G_s**(지수화된 양)로 한다.
- 이분법이 좁은 공명을 먼 곳에서 찾는 이유: 분산형 성분의 **1/δ 꼬리**(흡수형 1/δ²보다 느리게 감쇠).
- not-a-knot cubic이 PCHIP/Akima보다 좋다(PCHIP는 극값을 평탄화, long_cell에서 3.26 dB).
- 곡선 자체의 범위가 작을 때(cold_weak: 0.04 dB)는 절대 허용치가 개형을 뭉갠다 →
  허용치 = min(절대값, 곡선 범위의 5 %), 하한 2e-3 dB / 2e-4 dex.

### 1.3 커널 미세 최적화의 한계

Cp는 Hermitian 기저에서 rank 6(probe seed)/8(conjugate seed)이다. 조화마다 16개 대신 r개 RHS만 푸는
저계수 커널을 만들었고 1e-14까지 일치했지만 **1.36×** 에 그쳤다. M=16에서는 점당 ~36k flop이 바닥이다
(Ultra 성능 진단에서 얻은 "M=16에서는 LU가 행렬곱보다 싸다"는 교훈과 같다).
**점의 수를 줄이는 구조 변경만 이긴다.**

### 1.4 핵심 전환: 극점–잔류(pole–residue) 형태

고정 δ에서 확장 Floquet 계는 Δ_eff에 대해 아핀이다: (𝔸₀(δ) − Δ_eff·𝕊) x = e, 𝕊 = I ⊗ S_v.
S_v는 광학 coherence(조화 블록당 8좌표)에만 작용하므로 좌표를 S의 치역 R과 영공간 N으로 나누면
Δ_eff는 R–R 블록에만 나타나고, N의 소거(Schur 보수)는 Δ_eff와 무관하며 **정확**하다.

```
K = A_RR − A_RN A_NN⁻¹ A_NR,   Z = S_RR⁻¹ K,   g = S_RR⁻¹(A_RN A_NN⁻¹ e_N − e_R)
χ(Δ_eff) = c₀ − Σ_k res_k / (λ_k − Δ_eff),   λ = eig(Z),   res = (w_eff V) ∘ (V⁻¹ g)
```

즉 χ는 N_F=3에서 56개, N_F=2에서 40개 극점을 갖는 **정확한 유리함수**다. 그러면

- 임의의 이산 측도(Ultra의 1 m/s·4σ 격자 + 선형보간 가중치)는 Σ_j w_j/(λ_k − D_j)라는 **값싼 산술 합**이 되고,
- 연속 Gaussian 평균은 유수 정리로 **닫힌 형태**: ⟨1/(λ−X)⟩ = −i√π w(z)/(σ√2)(Im z>0),
  +i√π w(−z)/(σ√2)(Im z<0), z=(λ−Δ)/(σ√2), w는 Faddeeva 함수.

검증: 직접 해 대비 1e-11–4e-9, Faddeeva 평균 대비 0.5 m/s·8σ 격자 1e-10–1e-12, cond(V) 7.5e2–8.9e3,
실축에서 가장 가까운 극점 3.2–4.6 MHz(=Γ/2 부근, 1.1절의 페르미 추정과 일치).

### 1.5 "더 정확한 답"이 게이트를 깨뜨린 사건 → 모델 정의를 맞춘다

처음에는 비절단 Maxwell(Faddeeva)을 기본으로 했다. 11개 케이스는 Ultra와 ≤0.03 dB였지만
**hot(150 °C, Δ=0.9 GHz)에서 tube 0.112 dB**, frontier p95 0.035 dB로 동결 기준을 넘었다.
원인은 Ultra의 **4σ 절단이 F′ 공명 무리(Δ_eff≈−0.4…0 GHz, 이 조건에서 3.4–5σ)를 자르기** 때문이다.

기준을 완화하지 않았다. 4σ 절단은 Ultra tier가 정의하는 **모델의 일부**이므로, Fast/Balanced가
"같은 모델의 싼 계산"이 되려면 같은 측도를 써야 한다. 그래서 극점 형태를 **Ultra의 이산 측도 그대로**
합산하도록 바꿨고 hot은 0.0005 dB가 됐다. 비절단 평균은 `PoleResidues.gaussian_mean`으로 남겼다.
Ultra 절단의 물리적 타당성은 이 작업 범위 밖이라 [checklist](../../docs/checklist.json)에 후속 항목으로 기록했다.

### 1.6 실수형(real form): 𝕁 대칭

Hermitian 기저에서 L, C_δ, S_v는 실수이고 C₋=conj(C₊)이므로, 켤레 + 조화 반전(h↔−h)이 해를 해로 보낸다:
x₋ₕ = conj(xₕ). xₕ = aₕ + i·bₕ로 쓰면 같은 크기의 **실수 연립계**가 된다
(`pole_doppler.floquet_real_form`, 0차 방정식은 L a₀ + 2P a₁ + 2Q b₁, P+iQ=C₊). 고유분해가
zgeev 1.22 ms → dgeev **0.44 ms**(56×56).

### 1.7 구현 제약에서 나온 결정

- **numba에서 ctypes로 LAPACK dgeev 호출**은 동작하지만(401×56² 48 ms) `cache=True`가 불가("dynamic globals")
  → 프로세스마다 재컴파일. **분업**으로 해결: 구조적 산술(Schur 보수·잔류·극점 합)은 캐시되는 numba 커널,
  고유분해만 NumPy LAPACK을 스레드 풀에서(LAPACK은 GIL 해제).
- **scipy는 핵심 requirements에 없다** → Faddeeva는 Weideman(1994) 32항 유리 전개(Horner 평가, scipy 대비
  3e-13), 스플라인은 Thomas 알고리즘 기반 not-a-knot(scipy CubicSpline 대비 4e-14)로 직접 구현.
- **진단용 흡수**: 약한 probe 2-level OBE의 수치 Doppler 평균은 Voigt이므로 해석형(`voigt_profile`)으로
  바꿨다(0.1 ms vs 150 ms). pump OD 차이 ≤1.5e-3 → S 영향 ≤0.001 dB. 먼 꼬리의 arm OD 배열(표시·S에 미사용)은
  Ultra의 4σ 절단 때문에 최대 24 % 다르며, provenance에 "analytic Voigt"로 명시했다.
  선 가중치(C_F²만 또는 p_F·C_F²)는 해석형이 다시 적지 않는다. 프로세스당 한 번 `absorption._hyperfine_alpha`를
  네 선 중심에서 평가해 바닥 manifold마다 1 또는 p_F로 **정확히 분류**하고(`_reference_population_factors`),
  2 % 안에 맞는 규약이 없으면 Ultra의 수치 표로 되돌아간다. 이 규약은 2026-09-06에 교정 중이며(열적 부준위
  인구를 한 번만 셈), 커밋 단독 검증에서 HEAD와 작업트리의 규약이 달라 Fast/Balanced와 Ultra의 S_dB가 최대
  0.68 dB 벌어진 것으로 필요가 드러났다. 작업트리(교정된 규약)에서는 이 변경 전후 출력이 비트 단위로 같다.
- **가드**: 풀이한 모든 δ 행에서 극점 형태(Δ_eff=Δ)를 기존 compiled kernel과 대조(허용 1e-6).
  실패 행은 그 행만 격자 해로 교체하고, 적응 모드에서는 교정된 값으로 세분을 다시 판정한다.
  12 케이스 실측 최대 5e-8, 교체 0행.

### 1.8 프로파일 기반 다듬기 (Fast default)

| 단계 | 시간 | 병목 → 조치 |
|---|---:|---|
| 첫 통합 | 148 ms | 스코어러 55 ms(64-segment Python 루프, 라운드당 2회), eig 풀 285회 제출, 라운드마다 가드 12회 |
| 배치·풀링·지연 가드 | 109 ms | 예측/정확 점을 한 번에 채점, 네 계의 eig를 크기별로 모아 최소 청크 4, 가드는 끝에서 한 번 |
| 컴파일 스코어러 | 94 ms | Ultra의 segmented depletion 재귀를 같은 식(닫힌형 2×2 exp, 지수 클램프, pump 차감)으로 numba화 — 세분 판정 전용, 표시 readout은 Ultra 공유 코드 |

컴파일 스코어러는 Python 판과 rtol 1e-8 이내(Manley–Rowe 한계를 넘는 행은 pump 차감 피드백이
반올림을 ~1e-10으로 증폭)이며, 세분 결정은 이전과 동일한 노드 집합을 냈다.

## 기법

요청된 네 관점별로 실제로 쓴 것만 적는다.

**프로그래머의 최적화**
- 측정 우선(Amdahl): 매 단계 cProfile로 병목을 확인하고 나서 고쳤다(1.8절).
- 작업량 제거 > 미세 최적화: 커널 1.36×(1.3절)보다 "δ당 solve 1573→eig 1" 구조 변경이 결정적.
- 배치·풀링·지연 검증: 라운드별 호출을 합치고, 검증은 끝에서 한 번.
- GIL-free·캐시 가능한 커널과 벤더 LAPACK의 분업, 결정론적 청크 경계(같은 입력 → 같은 비트).

**옛날 프로그래머의 경량화**
- 룩업+보간 테이블 사고방식: 비싼 물리는 표본점에서만, 나머지는 Horner 평가 스플라인.
- 적응 세분(adaptive plotting/Bezier flattening식 이분법): 틀린 구간만 쪼갠다, 푼 점은 버리지 않는다.
- Thomas(삼중대각) 알고리즘, Weideman 유리 전개처럼 의존성 없는 짧은 수치 루틴.
- 보간 가중치를 노드 가중치로 미리 접어 두기(선형보간 + 가중합 → 단일 가중합).
- 상수 메모이제이션(Zeeman 일관성 진단), 표시하지 않는 것을 비싸게 계산하지 않기(진단 흡수의 해석형).

**물리학자의 추정**
- 페르미 추정으로 격자 한계를 먼저 계산: 사다리꼴 오차 e^{−2πγ/h}, γ ≥ Γ/2 → h ≲ 2–4 m/s(1.1절),
  이후 실측 극점 거리 3.2 MHz로 확인.
- 유수 정리/부분분수: 속도 적분을 극점 합과 Faddeeva 닫힌 형태로(1.4절).
- 대칭으로 차원·산술 줄이기: 𝕁(켤레∘조화 반전) → 실수형(1.6절), D 의존 블록만 남기는 Schur 보수.
- 선형인 변수에서 보간하고 관측량으로 판정, 분산 꼬리로 공명 탐지(1.2절).
- "더 정확한 답"과 "같은 모델"을 구분: 비교 대상의 모델 정의(4σ)를 먼저 확인(1.5절).

## 기각

| 시도 | 결과 | 기각 이유 |
|---|---|---|
| Gauss–Legendre/Hermite 16–128 노드 | 25–740 % 오차 | 공명 폭 < 노드 간격 |
| composite GL 패널 | 20–77 % | 같은 이유 |
| sinh 사상 사다리꼴(꼬리에 노드 성기게) | 같은 노드 수에서 균일보다 나쁨 | 꼬리의 공명도 중심과 같은 해상도 필요 |
| 균일 δ 솎아내기(stride 2/4/8) | 최대 3.8–8.4 dB | 폭 ~5 MHz Raman 공명 누락 |
| PCHIP/Akima 보간 | cubic 대비 수배–수십 배 오차 | 극값 평탄화 |
| 저계수 RHS Floquet 커널 | 1.36× | M=16에서 점당 비용 바닥 |
| AAA 유리근사(속도축) | 적합 비용 ≥ 격자 풀이 | 행마다 SVD 반복 |
| 윤곽 이동(복소 Δ_eff) | 불가 | 𝕁 대칭으로 극점이 실축 양쪽 |
| 비절단 Maxwell을 기본값으로 | hot tube 0.112 dB | Ultra 모델 정의(4σ)와 다른 적분 — 옵션으로만 유지 |
| numba+ctypes dgeev | 빠르지만 캐시 불가 | 매 프로세스 재컴파일(G1.4 위반) |
| 복소 eig(zgeev) | 실수형 대비 2.8× 느림 | 𝕁 대칭 미사용 |
| 라운드마다 커널 가드 | Fast 시간 8 % | 결과 동일, 끝에서 한 번이면 충분 |
| 1-segment 근사 스코어러 | 미채택 | 사용자가 보는 depletion 물리와 달라 과소세분 위험 |

## 측정

게이트 최종 판정: **PASS** ([gate_report.json](gate_report.json), 2026-09-11, 이 PC: Python 3.14.2 · NumPy 2.4.1 ·
numba 0.65.1 · 16 cores). 12개 동작점, 각 tier 5회 중앙값. Ultra는 같은 실행에서 다시 계산한 기준이다(케이스당
6.2–7.5 s).

| 항목 | 기준 (Fast / Balanced) | Fast | Balanced |
|---|---|---:|---:|
| default 시간 | ≤0.12 / ≤0.25 s | 0.101 s | 0.236 s |
| 12 케이스 최악 시간 | ≤0.25 / ≤0.40 s | 0.123 s | 0.257 s |
| 변경 전 대비 (default) | ≥4.5× / ≥7× | 5.8× (0.584→0.101 s) | 8.0× (1.898→0.236 s) |
| Ultra 대비 (대략) | — | ~65× | ~28× |
| 풀이한 δ 수 | — | 33–115 / 401 | 401 / 401 |
| S_dB tube 최대 | ≤0.10 / ≤0.05 dB | 0.018 dB | 0.0023 dB |
| S_dB p95 최대 | ≤0.03 / ≤0.015 dB | 0.0055 dB | 0.0010 dB |
| log10 G_s tube 최대 | ≤0.010 / ≤0.005 | 0.0009 | <1e-4 |
| 동작점 \|ΔS\| 최대 | ≤0.05 / ≤0.02 dB | 0.0018 dB | 0.0008 dB |
| 동작점 \|ΔG_s\|/G_s 최대 | ≤1 / ≤0.5 % | 0.004 % | <0.001 % |
| S 최솟값 차 최대 | ≤0.05 / ≤0.03 dB | 0.0036 dB | 0.0023 dB |
| G_s 최댓값 상대차 최대 | ≤2 / ≤1 % | 0.006 % | <0.001 % |
| Pearson(S_dB) 최소 | ≥0.998 / ≥0.9995 | 0.999992 | 1.000000 |
| 극점 가드 최대 / 격자 교체 행 | ≤1e-6 | 2.3e-11 / 0 | 2.3e-11 / 0 |
| Floquet N_F=3 vs 2 감사 | CONVERGED | 12/12 | 12/12 |

그 밖의 항목: G1.4 새 프로세스의 첫 Fast 호출 0.66 s(import 0.87 s 별도), G3.7 Ultra default 출력 bit-identical,
G3.8 numba 비활성 경로의 최대 상대차 1.7e-11, G4.1 전체 pytest 921 passed / 1 failed(작업 전과 같은
`test_docs_consistency` 누락 파일, 신규 실패 0), G4.2 신규 테스트 numba 32 passed · numba 없이 27 passed 5 skipped
(numba 전용 검사), G4.3 regression baseline 불변, G5 기록.

**여유가 얇은 항목.** Balanced default 0.236 s는 기준 0.25 s에 약 5 % 여유다. 첫 게이트 실행에서는 같은 항목이
0.202 s였다(그 실행은 게이트 스크립트의 G3.8 단계 문법 오류로 중단되어 수정 후 처음부터 다시 돌렸다). 실행 간
±15 % 흔들림이 있으므로 다른 작업이 CPU를 쓰는 상태에서는 넘을 수 있다. Fast default는 두 실행에서 0.091 / 0.101 s.

## 재사용

이 방법은 **속도 의존이 고정 행렬 S에 곱해진 아핀 이동(−D·S)** 이고 관측량이 상태에 선형인 모든 정상상태에
그대로 쓸 수 있다. 확인 순서:

1. 스캔 한 행을 `(A00 + t·A01 − D·S) x = rhs`로 쓸 수 있는가? `t`는 스캔 변수(δ, B, …), D는 Doppler 이동.
   - Λ EIT/AT/CPT: `absorption._affine_scan_coeffs`가 이미 `base + s·A + kv·B` 형태 → A00=base(+trace 행),
     A01=A_coef, S=−B_coef.
   - magneto buffer (B, v) 격자: Zeeman이 B에 아핀이면 t=B.
   - finite-Floquet: `pole_doppler.floquet_real_form`(C₋=conj(C₊), 실수 L·S일 때).
2. `AffineShiftSystem(...).residues(t)`로 극점 행을 얻고, 비교 대상과 **같은 측도**를 쓴다
   (`discrete_mean(nodes, weights)` — 기존 격자 재현, `gaussian_mean` — 비절단 Maxwell).
3. 가드: 몇 개 Δ에서 기존 solver와 대조(여기서는 모든 행, 1e-6).
4. 표시 축이 촘촘하면 `adaptive_scan.refine_scan_nodes` + `cubic_spline`: 선형응답량을 보간하고 관측량으로 판정,
   허용치는 곡선 범위에 비례시킨다.
5. provenance에 방법·풀이 점 수·보간 여부·감사 범위를 남기고, 보간이 있으면 claim gate 사유에 적는다.

주의: Doppler가 δ(두 광자)에도 들어가는 비공선 기하, 속도에 따라 행렬 구조가 바뀌는 모델(예: 속도 의존
충돌률)은 이 가정 밖이다. S_RR이 특이하면(`AffineShiftSystem`이 거부) 적용할 수 없다.

## 후속 (체크리스트에 기록)

- Ultra의 4σ 절단이 공명 무리를 자르는 조건(hot 등)에서 ~0.1 dB 영향 — tier 정의 재검토.
- `analysis/squeezing` 스캐너와 다른 scheme(OD/SAS, Λ, magneto)에 적용.

## Full probe scan 전환 (2026-09-14)

두 branch를 −8…12 GHz로 그리는 extra view(브랜치당 약 701점: 66.7 MHz 간격 301점 + 브랜치 중심 ±80 MHz의
0.4 MHz 간격 401점)도 Fast/Balanced에서 극점 엔진으로 옮겼다. 이전에는 옛 격자 설정(3σ, `PHASE_BALANCED/FINE`)을
썼고 default에서 Fast 4.07 s, Balanced 8.13 s, Ultra 22–24 s였다.

**기준(구현 전 고정).** 같은 모델(`PHASE_ULTRA`, 1 m/s·4σ)의 Ultra full scan 대비, 네 동작점(default, frontier,
near_res, hot)의 두 branch 전 구간에서
- Balanced: S_dB 점별 최대 ≤ 0.05 dB, p95 ≤ 0.015 dB, log10 G_s 점별 최대 ≤ 0.005
- Fast: S_dB tube 최대 ≤ 0.10 dB, p95 ≤ 0.03 dB, log10 G_s tube 최대 ≤ 0.010
- 공통: 모든 풀이 행 가드 ≤ 1e-6, Floquet 판정이 Ultra와 같음, 이전보다 빠름

**경과.**
1. Fast에 적응 세분을 그대로 쓰니 인덱스 stride 32가 66.7 MHz 구간에서 약 2 GHz가 되어, 조밀 창 바로 밖의
   이득 절벽과 먼 공명을 건너뛰었다. 보간한 점에서 S 최대 7.1 dB, 풀이한 점의 오차는 0이었다.
2. 초기 노드 간격을 "stride × 네 칸 연속 유지되는 최소 표시 간격" 이하로 묶었다(`adaptive_scan.initial_scan_nodes`).
   거친 구간은 전부 풀리고, 균일 축인 메인 곡선의 노드와 출력은 비트 단위로 그대로다. 단순 최솟값 규칙은 조밀 창
   점과 1 kHz 떨어진 거친 점 하나 때문에 모든 점을 푸는 것으로 무너져서 네 칸 규칙으로 바꿨다.
3. 그래도 150 °C에서 조밀 창 안의 폭 약 2 MHz 이득 구조(G_s가 2 MHz 안에서 3–4배)를 놓쳤다(log10 G_s tube
   0.156). 그 구간의 S_dB는 평평해 판정에 드러나지 않았다. 적응 세분은 표시 간격 2.75 MHz에서 검증했는데 이
   창은 7배 촘촘하다.
4. 기준을 완화하지 않았다. full scan에서는 **두 극점 tier 모두 모든 표시 점을 푼다.** 적응 세분은 검증된 메인
   401점 곡선에만 쓴다.

**결과(PASS).** 네 동작점, 두 branch 모두 S_dB 점별 최대 ≤ 0.0008 dB, log10 G_s 차 < 1e-4, 가드 ≤ 2.5e-10,
격자 교체 0행, Floquet CONVERGED(Ultra와 동일). Default view 시간은 Fast 0.72 s(5.6×), Balanced 0.67–1.27 s
(측정 실행 간 흔들림, 약 7–12×)이며 full scan에서 두 tier는 같은 계산이다. 메인 Fast/Balanced 출력은 이 변경
전후로 비트 단위로 같다.

**교훈.** 적응 보간을 믿을 수 있는 범위는 검증한 표시 해상도까지다. 더 촘촘한 표시 창에서는 좁은 이득·분산
구조가 보이므로 모든 점을 푸는 편이 정확하고, 극점 엔진 덕분에 비용도 작다. 기준·검증 스크립트는 탐색
스크립트처럼 세션 scratch에서 실행했고, 운영 경로의 계약은 `tests/test_fwm_fast_tiers.py`(view가 넘기는 방법·모델,
넓은 스캔 두 branch의 극점=격자)와 `tests/test_pole_doppler.py`(초기 노드 규칙)가 고정한다.

## 재현

```powershell
python -m pytest -q tests/test_pole_doppler.py tests/test_fwm_fast_tiers.py
python -m analysis.fwm_lite.gate --output analysis/fwm_lite/gate_report.json
```

탐색 스크립트(오프라인 구적 연구, δ 솎아내기, 속도 간격 민감도, 극점 검증)는 세션 scratch에서 실행했고,
그 결론과 수치는 위 표에 옮겼다. 제품 코드의 정확성 주장은 `tests/`와 게이트 보고서만 근거로 한다.

## 커밋 단독 검증

커밋은 부모 `c0c46f0`에 이 패치만 더한 트리다. 그 트리를 분리된 worktree에 체크아웃해 다시 검증했다.

- **첫 시도에서 5개 테스트 실패.** HEAD의 `absorption._hyperfine_alpha`는 선을 p_F·C_F²로 가중하는데, 해석형
  twin은 작업트리의 교정된 C_F² 규약을 적어 두었다. pump OD가 12/5·12/7배 달랐고 Ultra-phase S_dB가 최대
  0.68 dB 벌어졌다. 규약을 분류해 따르도록 고쳤다(1.7절).
- **수정 후 분리 트리**(규약 p_F·C_F² 자동 감지):
  - 전체 pytest 622 passed, 1 skipped, 실패 0 (FWM·극점 관련 125 passed 포함).
  - 같은 트리에서 다시 계산한 Ultra 대비 6개 동작점(default, frontier, hot, near_res, low_pump, blue_far)이
    모두 GATE.md G2 한계를 통과했다. Fast tube ≤ 0.015 dB·p95 ≤ 0.0055 dB, Balanced tube ≤ 0.0012 dB,
    가드 ≤ 5e-11, Floquet 감사 CONVERGED.
- **작업트리**(교정된 규약)에서는 수정 전후 Fast/Balanced 출력이 비트 단위로 같으므로 위 게이트 보고서가 유효하다.
- **재현 범위.** `gate_report.json`과 `before/`의 증거는 교정된 흡수 규약이 들어 있던 작업트리에서 만들었다.
  분리 트리에서 Ultra default의 G_s/G_c는 그 캡처와 1e-13까지 같지만 S_dB는 규약에 따라 달라지므로,
  다른 세션의 흡수 규약 교정이 반영되기 전에는 `gate.py`의 G3.7(Ultra 불변 비교)이 이 커밋 단독으로는
  재현되지 않는다. 이는 Ultra 코드 변경이 아니라 기준 캡처의 전제 차이다.

## 후속 gain 핫픽스 (2026-09-15)

Fast/Balanced에 별도 [유효 결합 교정](../fwm_gain_hotfix/DEVLOG.md) 추가.
원자 응답·극점 가드·Maxwell 측도는 위 기록과 동일. 기본값에서는 전파 전 비선형 결합만 교정.
`gain_closure_enabled=False`이면 기존 Ultra와 같은 물리의 결과 재현.
이 문서의 solver acceleration gate는 명시적으로 correction off를 사용.
Ultra production physics 유지. Gold calibration과 held-out 한계는 새 개발 기록 참고.
