# Fast/Balanced absolute-gain 핫픽스

2026-09-15. 범위: seeded ⁸⁵Rb D1 Fast/Balanced. 현재 `main` 공유 작업트리에서 작업.
기존 미커밋 변경 보존. Ultra·Grand Challenge의 원자/양자 물리 경로 유지.

## 결과·사용법

- Fast/Balanced 기본값: **단일 유효 비선형 결합 계수** 적용. Gold 표시 이득 약 **394.87 → 15.50**.
- `Advanced → Gain estimate → Semi-empirical gain correction`에서 on/off 비교.
- API: `gain_closure_enabled=False`로 변경 전 재현. 명시적 Fast/Balanced fidelity에서만 적용.
  fidelity 없는 `compute_spectrum`, `full_spectrum` 및 Ultra에는 미적용. 기존 fidelity 별칭도 지원.
- 원시 출력 `gain_closure`: 적용 여부, 계수, 정의, 출처, 교정점, 미확인 불확도, 알려진 실패 기록.
- **Gold 부근의 임시 교정 모델. 광범위한 절대 이득 검증 미완료.**
  Liu 비교 개선, McCormick 비교 악화. physical squeezing 예측·검증 불가 상태 유지.

## 1. Gain convention 고정

정본: [reference_points.json](reference_points.json).
저장소 literature CSV와 [Sim 원문](https://www.nature.com/articles/s41598-025-86479-w)의
실험 설정·gain spectra·squeezing 조건 대조. 모든 입력에 개별 출처 기록.

| 항목 | Gold |
|---|---:|
| OPD: F=2→F′=3 기준 pump detuning | +0.9 GHz |
| TPD: ω_seed = ω_pump − ν_HF + δ | −8 MHz |
| 온도 / cell length | 121 °C / 12.5 mm |
| Pump / seed input | 600 mW / 8 µW |
| Pump / seed waist | 530 / 330 µm, 1/e² **반경** |
| Crossing angle | 0.32° |
| 문헌 표의 대표 probe gain | **15.5**, 보고 범위 15–16의 대표값 |
| 별도 반올림 probe/conjugate output | 111 / 109 µW |
| 별도 raw probe ratio | 111/8 = **13.875** |
| 별도 raw conjugate/seed ratio | 109/8 = **13.625** |

교정 목표 `G_s = P_probe,out / P_seed,in = 15.5`.
15.5는 정밀 측정의 중앙 추정량 아님. 15–16도 오차막대·신뢰구간 아님.
Raw power 비와 대표 gain은 반올림·측정 조건/서술 차이로 수치 불일치.
동일한 정밀 교정값으로 합치지 않음. 임의 검출효율·광손실 역보정 없음.
검출 η=0.8694 등은 기존 noise readout 입력으로만 유지.

## 2. Gain-length 예산과 채택 모델

이상적 parametric 비교용 `x = acosh(sqrt(G_s))`. 일반 손실·위상불일치 전파의
엄밀한 단일 q를 측정한 값은 아님.

| 정확한 δ=−8 MHz 평가 | 값 |
|---|---:|
| 교정 전 G_s / G_c | 394.299058 / 396.520206 |
| x_model | 3.681067 |
| x_target, G_s=15.5 | 2.047033 |
| x_target / x_model | **0.556098** |
| raw power ratio 13.875의 x / x_model | 1.989666 / 3.681067 = 0.540513 |

401점 화면에서 보간한 −8 MHz 값은 교정 전 394.87. 위 표의 394.30은
정확한 −8 MHz를 포함한 3점 평가. 교정도 정확한 detuning에서 시행.
둘 사이 차이는 display interpolation이며 calibration target 교체가 아님.

**채택:**

```text
χ_sc,eff = C_mix χ_sc
χ_cs,eff = C_mix χ_cs
C_mix = 0.5594938027
```

이 계수로 Maxwell 행렬·64-segment 전파·기존 pump budget을 다시 계산.
`G_corrected = constant × G_raw` 형태 아님. χ_ss/χ_cc, 원자 정상상태,
밀도, 자연폭, transit, 기존 1/12 정규화, inherited 0.74 및 lab penalties 유지.

`C_mix`는 **Gold 하나에 교정한 유효 참여 계수**. 원자·기하 효과를 독립적으로 도출한
숫자로 제시하지 않음. Zeeman/polarization/angular/Raman/spatial 항의 개별 기여 식별 불가.
불확도 미확인. 계수는 기존 0.74·100 kHz transit·현재 finite-seed 모델에 조건부.
교정값의 해석은 qL 감소와 일치하지만 x 비와 정확히 동일할 필요 없음.
실제 전파에는 대각 손실·분산·위상불일치·pump budget이 포함되기 때문.

코드: [fwm_gain_closure.py](../../gabes/fwm_gain_closure.py).
교정 재현: [candidates.py](exploration_physics/candidates.py),
[candidate_results.json](exploration_physics/candidate_results.json).

## 3. 후보 비교·정확도/비용 선택

| 후보 | Gold G_s | 추가 원자 solve | 판정 |
|---|---:|---:|---|
| 기존 모델 | 394.30 | 0 | 기준 |
| off-diagonal C=0.5 | 10.13 | 0 | qL 규모 추정 지지 |
| Gaussian participation만 | 50.32 | 0 | 충분한 closure 아님 |
| 모든 χ에 C=0.5 | 8.46 | 0 | 대각 흡수·분산까지 변경; 불필요 |
| Gaussian × fitted residual | 15.50 | 0 | waist 추세 반전 때문에 기각 |
| **고정 C_mix, off-diagonal만** | **15.50** | **0** | **채택** |
| full Zeeman / 2D / self-consistent radial propagation | 미측정 | 다수 | 이번 interactive 경로 제외 |

Gaussian 후보: `C_spatial = w_p²/(w_p²+w_s²) = 0.720626`.
`q(r) ∝ I_p(r)` 가정에서 Gaussian 모드 적분으로 도출 가능.
하지만 Gold Ω_p/2π≈707.55 MHz로 강한 구동. 실제 χ는 pump intensity에 비선형.
Gold residual 0.776400을 붙이면 아래와 같이 미확인 waist 의존성이 크게 추가됨.

| Pump waist | 기존 G_s | Gaussian 후보 | 채택한 고정 C_mix |
|---|---:|---:|---:|
| 477 µm | 472.34 | 13.47 | 17.35 |
| 530 µm | 394.30 | 15.50 | 15.50 |
| 583 µm | 318.05 | 16.44 | 13.45 |

추가 radial 상태 해석 없이 반전을 물리 효과로 주장 불가. Gaussian 항은 **미적용 진단값**으로만
남김. 고정 C_mix는 이 추가 기울기를 만들지 않음.
비교 재현: [alternative.py](exploration_physics/alternative.py),
[alternative_results.json](exploration_physics/alternative_results.json).

다른 후보:

- **Zeeman/polarization:** 현재 drive/readout에 3·C_F², macroscopic coupling에 1/12 포함.
  24-level CG 합 진단값은 정의상 1. 추가 CG/평형 인구 보정은 중복. 별도 참여율 근거 부족.
- **Angular Doppler:** Gold transverse Raman RMS≈1.380 MHz.
  bare transit 폭 0.1 MHz를 곧바로 suppression 분모로 사용 불가.
  근방 χ를 Gaussian convolution한 단순 후보의 off-diagonal 변화는 약 +0.016%.
  별도 2D reference의 fixed-lab-beat 적분과도 다름. 독립 suppression 추가하지 않음.
- **Raman:** transit 50/100/200 kHz에서 기존 Gold G_s≈617/394/184.
  민감하지만 독립 측정 없음. gain 하나로 transit·C_mix 동시 fit하지 않음.
- **Spatial:** 기존 normalized axial crossing profile 유지. 새 Gaussian radial 법칙 미적용.

[sensitivity_results.json](exploration_physics/sensitivity_results.json)에 근거 수치.
채택 경로는 원자 solve/velocity node 수를 늘리지 않고 두 복소 배열만 곱함.

## 4. Calibration과 held-out 비교 분리

교정 이후 C_mix 고정. 다른 문헌에 재교정 없음.

| 동작점 | 문헌 probe gain | 기존 | 채택 모델 | 판정 |
|---|---:|---:|---:|---|
| Sim Gold | 대표 15.5 | 394.299 | 15.500 | calibration |
| Liu 2011 | 8 | 127.577 | 6.445 | 조건부 비교 개선; 약 19% 낮음 |
| McCormick 2008 | 9 | 7.133 | 1.926 | **비교 악화; 약 4.67배 낮음** |

Liu conjugate: 5.610 vs 문헌 7. McCormick은 baseline이 더 가까움.
두 held-out 문헌 모두 OPD의 정확한 excited F′ 기준 등 입력 모호성 존재.
그 모호성이 관측된 실패를 없애지는 않음. **범용 absolute gain model 검증 실패/미완료**로 기록.
예측이 모두 한 자릿수 배 범위라는 사실만으로 정밀 정확도 인정 불가.

Gold 주변 7개 변수의 3점 검사: pump power, 온도, pump/seed waist, angle, OPD, TPD.
고정 C_mix는 원래 모델의 국소 변화 방향 유지. 공명 이동 때문에 모든 변수가 전역 단조일 이유 없음.
이 검사는 수치적 추세 검사이며 실험적 parameter-dependence 검증 아님.
추가 문헌 7개는 isotope·missing gain·불완전 입력 등 이유로 제외; 사유는 reference JSON에 기록.

더 넓은 6변수×5점 화면 스캔도 시행. OPD 0.7–1.1 GHz, pump 300–900 mW,
101–141 °C, pump waist 400–660 µm, seed waist 230–430 µm, angle 0–0.64°.
**고온 외삽 미해결:** 131→141 °C에서 기존 이득 18,670→16,820은 pump cap 영향으로 감소하지만,
보정 이득은 약 320→12,368로 증가. 고온에서는 여전히 비현실적 증폭 가능.
나머지 비교 구간의 변화 방향은 기존 모델과 같음. 이 고온 실패와 McCormick 실패 때문에
넓은 운전 범위에서 정량적으로 신뢰할 수 있다는 acceptance claim은 하지 않음.

## 5. 회귀·성능 계약

- Pole guard는 원래 atomic response 검사. C_mix가 수치 오차를 숨기지 않도록 유지.
- Fast adaptive refinement는 **교정된** 최종 gain/noise를 채점.
- N_F=3/2 양쪽 transfer에 같은 교정 적용. 원자 χ 수렴 감사는 교정 전 값 유지.
- Ultra의 모든 top-level numeric **배열** 및 두 Fast/Balanced off fixture를 변경 전 캡처와 비교.
- Detector efficiency, post-cell loss, excess-noise 계수는 fitting/closure 입력 아님.
- 기존 `analysis/fwm_lite/gate.py`는 명시적으로 correction off를 사용해
  기존 같은-물리 solver acceleration 계약 유지. 신규 hotfix 검증은 별도 runner 사용.

### 성능

[validation_report.json](validation_report.json). 동일 PC·동일 실행에서 on/off 교대로 7회,
JIT/cache 예열 후 중앙값. Python 3.14.2, NumPy 2.4.1, Numba 0.65.1, 16 threads.
벤치마크 시작/종료 production source hash 동일.

| Tier | 수정 전 캡처 | 같은 실행 off | on | on/off 변화 | 기존 interactive 예산 |
|---|---:|---:|---:|---:|---:|
| Fast | 101.46 ms | 89.59 ms | **86.60 ms** | −3.34% | 120 ms |
| Balanced | 230.09 ms | 203.08 ms | **202.53 ms** | −0.27% | 250 ms |

수정 전 캡처 대비 −14.65% / −11.98%. 실행 시점 차이 포함하므로 실제 가속 주장에 사용하지 않음.
같은 실행 on/off도 작은 차이로, **runtime regression 관측 없음**이 결론.
Fast는 교정된 관측량 기준으로 적응 노드가 바뀔 수 있음. Balanced는 401점 유지.
Ultra default 변경 전 중앙값 6.638 s. 이번 핫픽스로 Ultra의 계산 방법 변경 없음.

화면상 Gold: Fast 15.504216 / Balanced 15.504247.
화면상 held-out: Liu 6.4944 vs 8, McCormick 1.9312 vs 9.
정확한 detuning 표와의 작은 차이는 401점 화면 보간 때문.
변경 전 캡처와 비교한 Fast/Balanced off 및 Ultra 두 조건의 모든 저장 배열 정확히 일치.

### 테스트

- `python -m pytest -q`: **1,699 passed, 3 skipped, 1 failed**, 261.90 s.
  실패는 `test_repository_visibility_wording_is_consistently_public`:
  작업 전부터 삭제된 `docs/FWM physics and analytic reconstruction/FWM_physics.tex` 누락.
  기존 삭제 복구·다른 작업의 문서 변경은 하지 않음. 핫픽스 신규 실패 없음.
- 이후 portable snapshot 허용오차만 정리하고 신규 검사 재실행:
  `python -m pytest -q tests/test_fwm_gain_hotfix.py` → **28 passed**, 16.86 s.
  플랫폼별 LAPACK/Numba 반올림 때문에 저장 fixture 비교는 relative 1e−7,
  absolute 1e−10, dB는 absolute 1e−5 사용. 동일 실행 on/off와 동일 환경
  before/after benchmark의 정확한 배열 일치 검사는 그대로 유지.
- 변경 전 pytest는 1,670 passed, 3 skipped, 2 failed.
  위 문서 누락 외 1건은 실행 중 production 파일 수정으로 발생한 quantum source-hash 검사 오염.
  기존 physics 실패로 분류하지 않음. 소스 고정 후 전체 실행에서 해당 검사 통과.
- production source hash가 벤치마크 이후 바뀌지 않았음을 마지막에 다시 확인.

로그: [pytest_after.log](pytest_after.log).
요약: [pytest_summary.json](pytest_summary.json).

## 6. 재현

```powershell
python -m pytest -q tests/test_fwm_gain_hotfix.py
python -m analysis.fwm_gain_hotfix.validate
python -m pytest -q
```

변경 전 증거: [before/capture.json](before/capture.json), numeric NPZ 4개.
공유 작업트리 source hash·Python/NumPy/Numba·thread 수 포함.
`.gitignore` 예외는 이 테스트 fixture와 실행 로그만 보존.

## 7. 남은 물리

Grand Challenge의 full Zeeman/polarization, non-collinear finite-seed Doppler,
독립 Raman decoherence, spatial pump/seed 상태, self-consistent depletion,
microscopic atomic noise, 충분한 held-out gain 검증 필요.
그 결과로 같은 출력 계약의 principled reduced model이 생기면 이 계수 제거/교체 가능.
현재 핫픽스가 그 closure를 완료했다는 주장 없음.

## 8. 핫픽스만 분리한 커밋 검증

커밋 전 `6061f0f` + 이 핫픽스만 담은 별도 소스 사본 생성. 기존 staged 변경 57개 경로와
미커밋 흡수 정규화·극점 재사용·Grand Challenge 작업은 제외.

첫 검증에서 기존 `before/` fixture 4건의 잡음 배열 불일치 확인.
기존 증거는 공유 작업트리의 `C_F²` 규약, parent commit은 `p_F·C_F²` 규약이기 때문.
원자 이득 교정이나 허용오차를 바꾸지 않고, **변경 전 parent FWM 코드와 parent 의존성**에서
새 기준 4건을 [before_parent/capture.json](before_parent/capture.json)에 별도 보존.
`snapshot_convention.py`가 실제 규약에 맞는 기준을 선택하고 알 수 없는 규약은 거부.
두 규약 모두 검사 가능하며 기존 증거 삭제·재교정 없음.

[commit_validation_report.json](commit_validation_report.json): 분리된 코드에서도
Gold 표시 이득 약 15.504, Fast 86.30 ms / Balanced 205.56 ms.
같은 실행 off는 90.86 / 205.42 ms. Ultra 두 조건 및 Fast/Balanced off의 저장 배열 모두
해당 parent 기준과 **정확히 일치**. 소스 hash도 검증 실행 중 동일.

분리된 전체 테스트에는 기존 양자 연구 보고서의 immutable source-hash 불일치가 있음.
보고서가 원래 작업트리의 혼합 CRLF/LF 및 이전 source byte를 고정해, clean parent에서도
`parent source changed: rerun the primary/reference audit` 발생.
해당 연구 기록 재생성이나 다른 작업의 테스트 수정을 이 핫픽스에 포함하지 않음.
최종 분리 테스트 요약은 [commit_pytest_summary.json](commit_pytest_summary.json)에 기록.

SHA-256로 연결된 이 폴더의 기록은 `.gitattributes`의 `-text`로 원래 byte 보존.

최종 분리 실행: **1,438 passed, 1 skipped, 11 failed, 25 errors**, 215.38 s.
핫픽스 **34개 검사 모두 통과**. 남은 36건은 모두 위 기존 quantum source-hash 검사/설정 오류.
다른 미커밋 연구 수정 없이 커밋하므로 해당 기존 오류는 별도 후속으로 유지.
