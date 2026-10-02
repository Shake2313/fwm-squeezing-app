# Fast/Balanced gain 핫픽스 → Grand Challenge 개발 인계

2026-09-16. 근거: 핫픽스 커밋 `08268ab`, [개발 기록](../../analysis/fwm_gain_hotfix/DEVLOG.md).
목적: 확인한 반례·수치 계약·재현 방법 전달. Grand Challenge와 milestone 상태는 진행 중 유지.
후속 항목 정본: [checklist.json](../checklist.json)의 `gain_hotfix_development_handoff`.

## 1. 교정값의 사용 범위

Fast/Balanced는 `χ_sc`, `χ_cs`에 **C_mix=0.5594938027** 적용 후 Maxwell 전파 재계산.
Gold 한 점의 대표 gain 15.5에 맞춘 값. 기존 residual 0.74·transit 100 kHz 등에도 조건부.
독립 Zeeman 참여율·공간 중첩·충돌 계수로 해석 불가. Grand Challenge no-fit 입력으로 승격 금지.

- Grand Challenge 비교에 production Fast/Balanced를 호출하면 `gain_closure_enabled=False` 명시.
- Ultra와 fidelity 없는 `compute_spectrum`·`full_spectrum`은 이번 교정 미적용.
  이 사실만으로 기존 경로 전체가 no-fit 또는 실험 검증됐다는 뜻은 아님.
- 계수는 평균장 교정. 미시적 diffusion·완전한 covariance·같은 조건 SQL을 제공하지 않음.
  gain 일치로 physical squeezing 상태 변경 금지.
- 후속 원리 기반 reduced model은 같은 출력 계약으로 교체 가능.
  각 physics 효과·수렴·held-out 성능을 먼저 검증하고 임시 계수 제거 여부 판단.

## 2. 수치 비교 전에 observable과 입력 규약 고정

정본: [reference_points.json](../../analysis/fwm_gain_hotfix/reference_points.json).

| 항목 | 이번 작업에서 확인한 구분 | 다음 개발에서 유지할 계약 |
|---|---|---|
| Probe power gain | `P_probe,out/P_seed,in`; 대표 15.5와 raw 111/8=13.875는 다른 기록 | 측정면·동작점·원시 power·대표값 여부 각각 저장 |
| Conjugate ratio | 109/8=13.625; probe gain과 별개 | 두 출력을 같은 gain 열로 합치지 않음 |
| 보고 범위 15–16 | 대표값의 근거; 불확도·신뢰구간 아님 | 오차 정보 없으면 unknown 유지 |
| 빔 크기 | pump/seed 530/330 µm는 1/e² 반경 | 반경/직경, amplitude/intensity 규약 저장 |
| 주파수 | OPD는 F=2→F′=3 기준 +0.9 GHz; TPD −8 MHz | excited-state 기준·부호·Hz/rad/s 명시; RF 축 별도 |
| 평가 격자 | 정확한 −8 MHz에서 기존 G=394.299; 화면 보간은 약 394.87 | 지정점 직접 평가와 화면 보간 오차 따로 기록 |

임의 η·광손실 역보정으로 숫자를 일치시키지 않음. Detector/noise 입력은 source coupling과 분리.
Grand Challenge의 독립 입력 ledger와 기존 `fwm-reference-parameter-provenance`에 연결.

## 3. Gain 오차는 결합·전파 단계에서 분해

진단량 `x=acosh(sqrt(G_s))`: 기존 3.6811, 목표 2.0470, 비율 약 0.5561.
이상적 parametric 비교에만 사용. 손실·위상불일치·depletion이 있는 일반 전파의 정확한 qL 아님.

- 최종 gain을 상수배해 수정하면 전파 구조·두 출력 관계를 놓침.
  원자 상태 → χ → transfer → pre-cap power → 최종 표시값을 함께 비교.
- χ 대각 흡수·분산과 비대각 mixing을 분리 진단. 이번 교정은 비대각만 변경.
  Grand Challenge의 physical variant는 Hamiltonian·reservoir·평균장·잡음을 일관되게 재계산.
- CG/line strength·ground population·1/12·기존 residual·새 계수의 사용 위치를 ledger로 추적.
  추가 Zeeman/인구 보정이 이미 평균된 요소를 다시 곱하지 않는지 확인.
- Transit 50/100/200 kHz에서 기존 Gold gain 약 617/394/184.
  Gain 하나로 transit·overlap·residual을 동시에 식별할 수 없음.
  독립 입력이 없으면 민감도와 식별 불가능한 조합 기록.

## 4. 실패 후보를 후속 모델의 반례로 재사용

근거: [후보 비교](../../analysis/fwm_gain_hotfix/exploration_physics/alternative_results.json),
[민감도](../../analysis/fwm_gain_hotfix/exploration_physics/sensitivity_results.json).

**Gaussian mode overlap:** `w_p²/(w_p²+w_s²)=0.720626`은 `q(r)∝I_p(r)` 가정의 진단값.
Gold Ω_p/2π 약 707.55 MHz에서는 국소 원자 응답의 포화·비선형을 따로 확인해야 함.

| Pump waist | 기존 gain | Gaussian × fitted residual | 채택한 고정 C_mix |
|---|---:|---:|---:|
| 477 µm | 472.34 | 13.47 | 17.35 |
| 530 µm | 394.30 | 15.50 | 15.50 |
| 583 µm | 318.05 | 16.44 | 13.45 |

Gaussian 후보는 Gold를 맞추지만 waist 추세 반전. 이번 핫픽스에서 기각.
후속 공간 모델에서는 pump power 고정/peak intensity 고정을 구분하고,
국소 atomic state·수집 모드·radial/axial 전파 수렴으로 기울기를 설명할 것.
기존 기울기를 무조건 정답으로 강제하는 회귀 fixture로 쓰지는 않음.

**Angular Doppler:** Gold transverse Raman RMS 약 1.380 MHz는 속도 분포 척도.
Gain 폭·RF bandwidth·추가 scalar decoherence와 동일하지 않음.
Bare transit 0.1 MHz를 suppression 분모로 대입하지 않음.
단순 χ convolution 진단에서 변화가 약 +0.016%였어도 실제 2D finite-seed 효과 부재의 증거 아님.
Fixed lab beat와 각 Doppler shift를 유지한 독립 2D reference·θ→0 한계 필요.

## 5. 한 점의 성공과 영역 검증 분리

근거: [최종 고정 계수의 held-out 비교](../../analysis/fwm_gain_hotfix/exploration_physics/alternative_results.json),
[넓은 스캔](../../analysis/fwm_gain_hotfix/validation_report.json).
아래는 **정확한 지정 detuning** 비교. 화면 보간 수치는 별도.

| 조건 | 문헌 probe gain | 기존 모델 | 교정 후 | 해석 |
|---|---:|---:|---:|---|
| Sim Gold | 대표 15.5 | 394.299 | 15.500 | 교정점 |
| Liu 2011 | 8 | 127.577 | 6.445 | 약 19% 낮음 |
| McCormick 2008 | 9 | 7.133 | 1.926 | 약 4.67배 낮음; 기존보다 악화 |

Held-out은 이번 계수 fit에서 제외한 문헌 조건. 이미 열람한 자료이며 새 blind holdout 아님.
문헌의 excited-state detuning 기준 등 모호성 보존. 실패를 지우는 근거로 사용하지 않음.

고온 화면 스캔: 131→141 °C에서 기존 gain 약 18,670→16,820, 교정 후 약 320→12,368.
기존 cap 영향과 교정 후 큰 증폭 모두 기록. 고정 계수로 고온 외삽 문제 해결 안 됨.

다음 비교에는 두 gain·transfer·cap 전후 power·원인 mask를 함께 저장.
Pump, T, waist, angle, OPD, TPD 변화와 fixed/retuned 운전 조건을 구분.
수치 수렴, 국소 추세, 기존 모델 일치, 실험 검증을 각각 판정.
McCormick·고온 실패는 후보 모델 진단에 유지하되 합격 수치를 새로 fit하지 않음.

## 6. 소스·물리 규약·수치 허용오차를 별도 관리

핫픽스 분리 커밋 과정에서 공유 작업트리 fixture는 `C_F²`, parent `6061f0f`는
`p_F·C_F²` 흡수 규약임을 확인. Source revision 하나만으로 dirty worktree를 재현할 수 없었음.

- 변경 전 parent 코드 **및 의존성**에서 [before_parent](../../analysis/fwm_gain_hotfix/before_parent/capture.json) 별도 캡처.
  원래 [before](../../analysis/fwm_gain_hotfix/before/capture.json) 증거 유지.
- [규약 선택기](../../analysis/fwm_gain_hotfix/snapshot_convention.py)는 알려진 두 규약만 허용.
  혼합·미등록 규약 거부. 물리 출력이나 허용오차를 바꾸어 fixture에 맞추지 않음.
  두 fixture 지원은 두 물리 규약의 타당성 인증과 별개.
- Portable snapshot: 일반 수치 배열 `rtol=1e-7, atol=1e-10`; dB `rtol=0, atol=1e-5`.
  Shape·정수/bool·증거 hash는 정확 비교. 이는 핫픽스 회귀 계약이며 GC 오차 기준 대체 불가.
- 같은 실행의 on/off 검사는 별도. [커밋 검증](../../analysis/fwm_gain_hotfix/commit_validation_report.json)에서
  Ultra 두 조건·Fast/Balanced off의 모든 저장 배열은 parent 기준과 **정확히 일치** 관측.
  Runner는 정확 일치와 허용오차 통과를 별도 기록하며, 성공 판정은 허용오차·소스 안정성 기준.
- 실행 중 파일 수정은 source-stability 실패. 물리 회귀와 구분하고 고정 소스에서 재실행.
  공유 작업트리 검사와 실제 배포할 커밋 내용의 검사를 구분.

기존 quantum 보고서의 raw hash·혼합 CRLF/LF 문제도 분리 커밋에서 재현됐음.
해당 실행은 1,438 passed / 1 skipped / 11 failed / 25 errors; 핫픽스 34개 통과.
36건은 그 분리 소스의 기존 provenance 검사/설정 오류. 현재 GC 상태를 뜻하는 수치 아님.

**GC에는 이미 대응 구현 존재:** [고정 소스 캠페인](portable_thermal_campaign.md).
상대 경로·Python universal newline identity와 raw hash를 함께 보존하고 고정 ZIP에서 실행.
검증 단위 테스트는 synthetic 임시 증거, 실제 연구 실행은 엄격한 source/reference 계약 사용.
새 구현 과제로 중복 등록하지 말고 다음 캠페인에서도 재사용.
소스 의미·물리 규약 변경까지 줄바꿈 정규화로 승인하거나 과거 report를 재해시하지 않음.
SHA-256로 연결된 NPZ·JSON·원본 로그는 Git checkout 전후 byte도 확인.
핫픽스는 해당 증거 폴더에 `.gitattributes`의 `-text` 적용. Python 소스의 newline 동일성 정책과 구분.

## 7. 속도 비교에도 같은 모델·같은 실행 조건 필요

커밋 검증의 예열 후 7회 중앙값: Fast on/off **86.30/90.86 ms**,
Balanced **205.56/205.42 ms**. 같은 실행에서 교대 측정, 소스 시작/종료 hash 확인.
이 결과는 runtime regression 미관측. 새 원자 알고리즘의 가속률 증거 아님.

- Fidelity·물리 variant·격자·adaptive 채점 관측량·Python/NumPy/Numba·thread 환경 기록.
- Fast adaptive refinement는 교정된 최종 관측량 기준이므로 on/off node 수가 달라질 수 있음.
  같은 물리의 solver 가속 비교는 correction off 등으로 모델을 고정.
- Atomic pole guard·교정 전 χ 수렴과 최종 관측량 수렴을 둘 다 유지.
  작은 최종 gain이 원자 응답 오차를 숨기지 않도록 함.
- GC 고정 캠페인의 실제 BLAS thread·solver work·cache hit/miss 계약 재사용.
  JIT 준비, 새 solve, cache 집계 시간을 구분. 이전 실행과의 벽시계 차이를 가속률로 주장하지 않음.

## 8. 다음 개발 체크리스트

- [ ] **GCH-1 입력/observable:** Gold 대표값·raw 두 power 비·평가 격자 분리. Unknown 필드를 유지한 비교표 작성.
- [ ] **GCH-2 결합/정규화:** C_mix를 target-fitted 비교 경로로 격리. 각 CG·population·residual의 한 번 사용 및 drift/noise 일관성 검증.
- [ ] **GCH-3 공간/각도 반례:** 477/530/583 µm와 θ→0 사례 재평가. Waist 추세는 독립 공간 계산, 2D는 fixed-lab-beat reference로 검증.
- [ ] **GCH-4 적용 영역:** Sim/Liu/McCormick 및 고온 스캔을 계수 재조정 없이 비교. Cap 전후 출력·실패·미확인 입력 함께 저장.
- [ ] **GCH-5 재현/성능:** 해당 단계의 고정 소스·규약에서 회귀와 독립 수렴, 같은 실행 timing 재검증. 과거 증거 보존.

위 체크는 **다음 GC 확장에 적용할 미완료 검사**. 핫픽스의 완료 기록·현재 thermal path 인증과 구분.
독립 실측·새 blind holdout 확보는 기존 `sabes-calibration-campaign` 및
`fwm-squeezing-spectrum-heldout-validation`에서 계속 추적.
