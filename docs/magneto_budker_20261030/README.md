# Budker 기준 자기 측정 디튜닝 구현

## 사용자 요청과 출처

2026-10-30, 현재 대화. 주 참고문헌: 사용자가 제공한
[PhysRevA.62.043403.pdf](<G:/내 드라이브/Docs/PhysRevA.62.043403.pdf>).

> Budker 논문을 메인 레퍼런스 삼고, 이론 논문들을 몇 개 더 찾아서, 자기 측정의 SNR이 완벽한 레조넌스가 아니라 약간의 디튜닝을 가질 때 더 측정이 잘 되는 현상을 자기 스킴에서 구현할 겁니다.

최초 요청은 구현 전 레퍼런스·주요 코드·패치 계획·영향·다른 스킴 재사용 요소 보고.
사전 보고 뒤 같은 날짜 사용자 승인: **“좋습니다. 이대로 한번 진행해주세요.”**
아래 설계·검증 기준은 AI 작성. 논문 서술은 사용자 지시가 아닌 과학적 근거.

## ⟦AI⟧ 2026-10-30 사용자 정정과 수정 범위

최초 구현 보고 뒤 사용자 직접 정정:

> 그냥 Magneto-optics 스킴 자체가 Budker 모델을 표현하도록 패치해달라고 한 거였는데, 혹시 디폴트에서 Budker를 선택했을 때만 저게 구현되는 건가요? 그러면 안 되는데...
> 또, 구현 화면 캡처에서 끊어진 부분이 나온 이유가 뭔가요?

최초 구현은 기존 단일 전이 B 스캔을 유지하고 Budker 계산을 선택형 광학 디튜닝
모드로 추가함. 사용자 의도와 불일치. 아래 최초 구현·검증 기록은 그 상태의
이력이며, 스킴 전체 교체 완료의 근거로 사용하지 않음.

수정 완료: 등록된 Magneto-optics의 B 스캔과 광학 디튜닝 스캔을 같은 다중 초미세
전이·공통 dark reservoir 물리 엔진으로 통합. 프리셋은 실험 조건만 변경하고 모델을
선택하지 않음. 직접 호출 별칭 Hanle/EIA/NMOR도 동일 엔진 사용.
단일 전이 구현은 사용자 선택지와 분리한 `SinglePairMagnetoReference` 수치 참조로 한정.
최초 기본 화면부터 새 엔진 사용. `Calculation` 및 중복 `budker_` 설정 제거.
**Scan axis → Magnetic field / Optical detuning**으로 스캔 축만 변경.
기본값은 ⁸⁵Rb D1의 B 스캔. D1/D2 프리셋은 조건만 변경.

### 캡처 곡선 단절의 원인

최초 D1 기준 스캔의 저장 결과에서 Δ=0, +300, +400 MHz의 회전 기울기는
±h와 ±h/2 유한 차분 간 차이가 각각 **5.667%, 5.324%, 5.159%**.
설정한 5% 문턱을 넘어 slope와 자기장 ASD를 null/NaN으로 치환했고,
그래프가 그 점에서 끊김. 물리적 불연속 또는 신호가 사라지는 현상으로 확인된 것이 아님.

세분 스캔에서는 h를 절반으로 줄였고 속도 적분 차수도 높임. 같은 세 Δ의 차이는
**1.593%, 1.459%, 1.440%**, 기울기는 모두 유한한 값으로 회복.
이 결과는 앞 캡처의 누락이 유한 차분 검사의 분해능 문제임을 뒷받침함.
저장 근거: [최초 D1 검사 JSON](../../analysis/magneto_budker_20261030/results/d1_validation.json)의
`cases.reference.readout.slope_valid`와 `slope_relative_change`.

수정: 공통 reservoir 방정식을 자기장으로 미분한 DC 선형 응답으로 기울기 계산.
임의 보간으로 누락을 연결하지 않음. 실제 기울기 영점에서는 등가 자기장 잡음이
발산하므로 해당 물리적 특이점과 계산 실패 상태를 구별해 표시.
수정본 49개 광학 표본 전체에서 해석적 기울기·ASD 유효. 이전 누락 세 지점도 포함.
미분 방정식 잔차 최대 5.42×10⁻¹². 49→97Δ·속도차수6→10 비교에서
공통 표본 기울기 차이/peak 5.56×10⁻⁷, 최소 감도 상대 변화 3.11×10⁻⁷.
이 검사는 선택한 유효 모델의 수치 확인이며 실험 정량 인증 아님.

### 현재 통합본 검증

새 기본 B 스캔과 Δ=−300 MHz B 스캔, 광학 49점·세분 97점 스캔을
실제 등록 스킴으로 실행. 동일 물리 작동점에서 두 축의 투과율·기울기·감도 일치.
SNR용 ±1 pT 검사와 해석적 감도는 별도. SNR 선형성 검사 실패로 기울기 곡선을
숨기지 않음. 횡방향 B와 원편광 투과 응답도 같은 다중 초미세 엔진에서 계산.
buffer 옵션은 단일 영역 reset 근사이며 압력 광학 넓어짐·확산 생략을 명시.

[현재 재현 스크립트](../../analysis/magneto_budker_20261030/validate_unified.py),
[현재 수치 결과](../../analysis/magneto_budker_20261030/results/unified_validation.json),
[기본 B 곡선](../../analysis/magneto_budker_20261030/results/unified_default_field.png),
[수정 광학 곡선](../../analysis/magneto_budker_20261030/results/unified_optical.png).

`python -m pytest -q tests/test_magneto.py tests/test_magneto_doppler.py tests/test_magneto_multilevel.py tests/test_magneto_full_obe_reference.py tests/test_magneto_unified.py tests/test_magneto_budker.py tests/test_magneto_noise.py tests/test_magneto_readout.py tests/test_ui_readout.py tests/test_defaults_store.py tests/test_ui_defaults.py tests/test_resonant_hanle_reference.py --durations=15 --junitxml=docs/magneto_budker_20261030/unified_tests.xml`

**249 passed**, 18.82 s. [기록](unified_tests.xml).
새 등록 모델·해석적 미분·UI/검출/기본값과 미등록 구형 참조의 보존 범위를 포함.
`pytest tests/test_schemes_render.py -k Magneto`: **8 passed**, 2.98 s
([기록](unified_render_tests.xml)).
`pytest tests/test_kernels.py -k magneto`: **6 passed**, 3.02 s
([구형 참조 kernel 기록](reference_kernel_tests.xml)). 총 **263개 고유 검사 통과**.
이후 SNR 도움말만 현행 검사 내용으로 정정; 관련 3개 검사 재확인.
과거 squeezing/Hanle 분석과 pole 검산은 구형 참조를 명시해 기존 물리를 보존.
공유 core/FWM 수치 변경·전체 suite·GC 캠페인 없음.

실제 브라우저에서 **프리셋을 조작하지 않은 최초 기본 화면**의 JSON을 내려받아
`model=Budker multihyperfine`, B 스캔 83/83 응답 유효 확인.
Scan axis만 Optical detuning으로 바꿔 동일 모델·동일 조건, 49/49 유한 기울기·감도 확인.
0, +300, +400 MHz의 감도도 각각 8.776×10⁻¹⁵, 2.105×10⁻¹⁵,
1.405×10⁻¹⁵ T/√Hz로 정상. 모델/Calculation 선택지 없음.
[브라우저 기록](unified_browser_smoke.json),
[수정 기본 화면](../../output/playwright/magneto_unified_default.png),
[연속 광학 곡선 화면](../../output/playwright/magneto_unified_optical.png).
오류 0건. 검사용 서버·브라우저 종료.

## ⟦AI⟧ 현재 Defaults — Hanle / NMOR

사용자 후속 요청에 따라 기본 프리셋을 **Hanle**, **NMOR** 두 개로 정리.
이 자기 스킴에서는 EIT/EIA/LCA를 Hanle 배치의 투과·흡수 공명 특징으로 묶고,
NMOR은 편광 회전 읽기로 구분. 독립 모델이 아닌 공통 엔진의 시작 설정.
두 프리셋의 기존 물리 입력은 유지. 상단 선택지는 추가하지 않음.
편광·셀·완화·광학선은 사이드바에서 조절 가능.
Defaults·통합 스킴·UI 검사 **78 passed**, Magneto 표시 검사 **8 passed**.
프리셋 정리 범위만 검증; 수치 엔진·전체 suite·GC 재검사 없음.

### 다섯 프리셋 복원 — 후속 단순화 전 이력

사용자 추가 요청: EIT dip / EIA peak / Hanle / Buffer LCA / NMOR을 상단에
노출하지 않고 Defaults 선택지로 제공. 통합 스킴 `presets()`에 다섯 이름 추가.
`recommended_defaults=None`, schema의 `applies_defaults=False` 유지하여 상단 노출 없음.
전체 스냅샷으로 B 축·신호·셀·편광을 복원. D1/D2 reference 항목도 유지.
모든 선택은 같은 다중 초미세 엔진 사용.

이전 kHz pure-dephasing 수치와 새 s⁻¹ isotropic reset은 의미가 다르므로 무조건 복사하지 않음.
⁸⁷Rb D1의 전이 위치는 성분 centroid 기준 −407.24/+407.24 MHz로 변환.
기타 값은 새 모델의 시작 조건이며 프리셋 이름 자체를 계산된 dip/peak 부호의 검증으로
사용하지 않음. 기본값·상단 비노출·공통 엔진·완전 설정 적용 집중 검사 포함.
`pytest tests/test_magneto_unified.py tests/test_ui_defaults.py tests/test_defaults_store.py`:
**78 passed**. `pytest tests/test_schemes_render.py -k Magneto`: **13 passed**.
새 다섯 조건의 Doppler-on 중앙/양 날개 계산도 정상상태 검사 통과.
수치 엔진 변경 없음. 전체 suite·GC 검사 생략.

## ⟦AI⟧ 문헌 근거

| 문헌 | 확인 위치 | 구현에서의 역할 |
|---|---|---|
| Budker et al., *PRA* **62**, 043403 (2000), [DOI](https://doi.org/10.1103/PhysRevA.62.043403) | 제공 PDF 원문, 식 (1)–(2), III–IV절, 그림 4–7 | 주 근거. 광학 디튜닝별 회전 기울기·투과광·광자 샷잡음 감도 비교 |
| Budker et al., *RMP* **74**, 1153 (2002), [원문](https://arxiv.org/pdf/physics/0203077) | V.C, VIII.E.1–2, XI.B 식 (31), XII.A 식 (33) | 겹친 초미세 전이의 펌프/프로브 교차항, 벽 충돌 속도 재분배, AOC, 광자·원자 잡음 구분 |
| Malakyan et al., *PRA* **69**, 013817 (2004), [원문](https://arxiv.org/pdf/physics/0309061) | 식 (39)–(42), (52), (55)–(58) | 약한 광펌핑·완전 속도 혼합 한계. 펌핑된 코히런스를 속도 평균한 뒤 속도별 프로브 응답 계산 |
| Rochester & Budker, *J. Mod. Opt.* **49**, 2543 (2002), [원문](https://arxiv.org/pdf/physics/0202071) | 식 (7), (10)–(17), IV절 | 흡수·회전·전파·검출 광자 수의 일관성. 고립된 전이 모델이며 속도 혼합은 제외 |
| Fleischhauer, Matsko & Scully, *PRA* **62**, 013808 (2000), [원문](https://arxiv.org/pdf/quant-ph/0001072) | 식 (1), (11)–(14), (60)–(63) | 신호 기울기로 환산한 감도와 광자 잡음/AC Stark 잡음 구분. 논문의 유효 비공명 결합 디튜닝 Δ₀를 UI 광학 디튜닝과 동일시하지 않음 |
| Reiter & Sørensen, *PRA* **85**, 032111 (2012), [원문](https://arxiv.org/pdf/1112.2806) | 식 (3)–(5), III절 | 약한 광학 여기에서 들뜬 상태를 소거하는 유효 Hamiltonian·Lindblad 연산자 근거 |

RMP VIII.E: ⁸⁵Rb D1의 $F_g=3\to F'_e=2,3$ 겹침 영역에서,
서로 다른 전이로 펌핑·검출한 원자의 회전 기여가 같은 전이의 기여와 반대 부호.
벽 충돌로 속도가 바뀌면 흡수선 중심에서 상쇄, 양 날개에서 회복 가능.
같은 속도를 유지하는 transit 응답과 구별. VIII.E.2는 밝은/어두운 영역·속도 혼합·
스핀 교환을 포함한 밀도행렬 확장을 요구하지만 완전한 수치 방정식은 제공하지 않음.
공통 어두운 영역 저장소와 블록 소거는 이를 구현하는 본 작업의 선택.

추가 확인 후보 PRL **85**, 2088 (2000)은 원문 확인 미완료.
웹 실패 후 요청된 브라우저 경로도 시도: APS 권한 거절, OSTI 빈 화면,
UNT CAPTCHA, 저자 저장소 링크 404. AOC 근거는 확인된 RMP V.C 사용.

### 주 논문의 실제 조건과 수치

| 항목 | 원문 조건 |
|---|---|
| 대상 | ⁸⁵Rb D1·D2, $F_g=3$ 성분; 선형 입사 편광, Faraday 배치 |
| 셀 | 진공·파라핀 코팅, 직경 10 cm |
| 온도·밀도 | 20 °C, 약 $4.5\times10^9\ \mathrm{cm^{-3}}$; 밀도 탐색은 20–30 °C |
| 광속 직경 | 약 3 mm, 그림 4–6 |
| 그림 4 | D1 wall Ramsey: 1.2 mW/cm²; transit: 0.6 mW/cm² |
| 그림 5·D1 최적 | 약 1 mW/cm², Doppler $F_g=3$ 흡수 성분 중심보다 약 −0.6 GHz; 약 $4\times10^{-12}\ \mathrm{G/\sqrt{Hz}}=0.4\ \mathrm{fT/\sqrt{Hz}}$ |
| 그림 6·D2 최적 | 약 4.5 mW/cm² 이상, 같은 방식의 중심보다 약 +0.4 GHz; 약 $3\times10^{-12}\ \mathrm{G/\sqrt{Hz}}=0.3\ \mathrm{fT/\sqrt{Hz}}$ |
| 기울기 측정 | $B_z=\pm1\ \mathrm{\mu G}$ 회전 차; 잔류장은 셀 평균 0.1 µG 미만 |
| 검출 | 셀 뒤 편광 변조 1 kHz·5 mrad, 첫 고조파 lock-in 검출 |

원문 수치는 측정한 회전 기울기·투과광으로 추정한 **샷잡음 한계**.
장치가 실제 샷잡음 한계에 도달했다는 결과 아님. 검출 변조 1 kHz는 ENBW 또는
원자 자기장 응답 대역폭과 다름. 원문의 약 ±0.4–0.6 GHz는 광학 디튜닝이며
Zeeman/Larmor 레조넌스로부터의 디튜닝 아님. 기존 단일 전이의 Δ=0과도 다름.

## ⟦AI⟧ 최초 구현 설계와 근사 범위 — 통합 전 이력

상태: 조건부 모델과 선택형 Budker 디튜닝 모드 구현. 정량 실험 재현은 미완료.
단일 고정 $F_g,F'_e$ 모델과 별도로 전체 초미세·Zeeman 결합을 사용.
현재 구현 파일: [magneto_multilevel.py](../../gabes/magneto_multilevel.py).

- 광학 들뜬 상태는 약한 여기에서 단열 소거. 각 바닥 $F_g$ 내부의 모든 $m,m'$
  코히런스·광펌핑·AC Stark Hamiltonian 유지; 서로 다른 바닥 $F_g$ 코히런스는 secular 근사.
- 자연선폭으로 절대 쌍극자 크기 결정. 부호를 보존한 Wigner-6j·CG 결합 사용.
  들뜬 각 상태의 전체 자발 방출률은 Γ, 바닥 초미세 분기와 누설 포함.
- 같은 입사 바닥 F·최종 바닥 F·방출 편광에 도달하는 가상 $F'_e$ 경로는
  진폭으로 합산. 서로 구별되는 최종 광자 채널을 무차별 합산하지 않음.
- 밝은 속도군과 공통 어두운 영역을 결합. 재진입 시 Maxwell 속도 분포로 재분배.
  독립 속도군 계산을 대조군으로 유지. Schur 소거는 이 유효 모델 내부의 수치 소거.
- 회전 기울기와 투과광을 동일한 원자 상태·광학 조건에서 계산.
  최적 디튜닝·상쇄점·부호는 결과로 탐색; 논문 최적점을 강제로 삽입하지 않음.

모델 Δ=0은 선택한 바닥 F 광학 성분의 선세기 가중 중심. 논문의 Doppler 흡수선
성분 중심과 동일함을 가정하지 않음. 정량 비교 전 주파수 기준 정렬 필요.
밀도는 선택 동위원소의 국소 밀도이며 존재비를 두 번 곱하지 않음.
원문에서 동위원소 농축 조건은 명시적으로 확인하지 못함.

유효 연산자는 $H_{\rm NH}=H_e-\frac i2\sum_qL_q^\dagger L_q$,

\[
H_{\rm eff}=H_g-\tfrac12V_-\left(H_{\rm NH}^{-1}+H_{\rm NH}^{-\dagger}\right)V_+,
\qquad L_q^{\rm eff}=L_qH_{\rm NH}^{-1}V_+.
\]

약한 **광학 여기**와 약한 **바닥상태 펌핑**은 다른 조건.
Γ보다 훨씬 느린 바닥상태 이완 동안 여러 번 펌핑·AOC가 일어날 수 있음.
따라서 AOC 기작은 유지할 수 있으나 광학 포화·들뜬상태 인구가 커지는 영역의
정량 정확도까지 보장하지 않음. 포화 진단과 모델 상태를 결과에 보존.
자기장 의존 원자 스핀 잡음, 기술잡음, 전자잡음, 셀 궤적·스핀 교환의 완전한
미시 모델, 실험의 전체 lock-in 검출 모델은 별도 범위.
광학 분모의 Zeeman shift와 선형 Faraday 배경도 생략. B≈0의 비선형 회전 모델임.
최초 결과는 조건부 광자 샷잡음 감도이며 논문 절대 수치 재현으로 표시하지 않음.
현재 전파는 균일 국소 광세기에서 얻은 복소 응답의 경로 적분·Beer–Lambert.
광세기 감쇠를 원자 상태에 재반영하는 종방향 자기일관 계산은 포함하지 않음.
코팅 셀의 종방향 확장은 속도와 위치를 함께 연결하는 공통 어두운 저장소가 필요;
z 구간마다 독립 저장소를 두는 방식도 별도 근사임.

## ⟦AI⟧ 검증 프로토콜

감도는 물리적 바이어스 $B_0$에서

\[
\eta_B(\Delta)=\frac{\sqrt{S_\theta(\Delta)}}{|\partial\theta/\partial B|_{B_0}},
\qquad
\sqrt{S_\theta}=\sqrt{\frac{h\nu}{2\eta_{\rm det}P_{\rm in}T(\Delta)}}
\]

사용. 단위 T/√Hz, 한쪽 PSD·이상적인 균형 polarimeter·선형 입력 기준.
ENBW W에서 RMS 최소장은 \(\eta_B\sqrt W\), 선언한 작은 신호
\(\delta B\)의 SNR은 \(|\delta B|/(\eta_B\sqrt W)\).
기존 B 스캔의 최대 회전 진폭/RMS와 구분.

| 검증 | 통과의 의미 |
|---|---|
| CG·6j·자발 방출 합규칙, 열적 비편극 약광 흡수 | 전이 세기·인구·절대 단위의 중복 적용 없음 |
| 유효 생성자의 trace 보존, Hermiticity, 정상상태 양성·잔차 | 선택한 유효 모델 수치 일관성 |
| 무광·무밀도·큰 디튜닝 한계, 균형 편광 $B\to-B$ 대칭 | 기본 한계와 부호 확인 |
| 유한 차분 $B$ 간격 축소, 속도 및 Δ 격자 세분화 | 기울기·최적점의 수치 안정성. 스캔 경계 최적은 내부 최적으로 인증하지 않음 |
| 공통 저장소/독립 속도군 비교 | 조건을 고정하고 속도 재분배에 의한 초미세 교차항 효과 분리 |
| ⁸⁵Rb D1의 중심 억제·날개 회복·비영점 감도 최적 | Budker 메커니즘의 정성적 구현. 모든 조건에서 비영점 최적을 강제하지 않음 |
| 검출 효율·광자 수·ENBW 변경 | 감도 ASD와 RMS/SNR의 서로 다른 스케일링 확인 |
| 기존 자기 스킴·공통 검출기 호출자·캐시/내보내기 | 선택 모드 도입에 따른 기존 동작 및 통합 회귀 범위 |

정량 논문 재현 판정은 별도: 주파수 영점·동위원소·셀·이완·실제 빔·검출 효율·
포화 조건을 맞추고 그림 데이터·측정 오차와 비교해야 함.
상대 경향 통과를 −0.6 GHz·0.4 fT/√Hz 또는 +0.4 GHz·0.3 fT/√Hz의 재현으로 바꾸지 않음.

### 실행 증거

`python -m pytest -q tests/test_magneto.py tests/test_magneto_doppler.py tests/test_magneto_multilevel.py tests/test_magneto_noise.py tests/test_magneto_budker.py tests/test_magneto_readout.py tests/test_ui_readout.py tests/test_defaults_store.py tests/test_ui_defaults.py --durations=10 --junitxml=docs/magneto_budker_20261030/focused_tests.xml`

**217 passed**, 14.39 s. [검사 기록](focused_tests.xml).
기존 export의 UTC 날짜 API deprecation warning 2건. 실패 없음.
새 모델·기존 자기 스킴·검출기·readout·기본값/캐시·export 집중 범위.
공유 core/도플러/FWM 수치 엔진 변경 없음. 전체 suite·GC 캠페인은 생략.

`python -m pytest -q tests/test_magneto_full_obe_reference.py --durations=5 --junitxml=docs/magneto_budker_20261030/obe_reference_tests.xml`

**2 passed**, 1.05 s. [별도 전체 OBE 비교](obe_reference_tests.xml).
⁸⁷Rb D1·단일 속도·I=10⁻⁴ mW/cm²에서 모든 16준위 OBE를 독립 조립.
정확한 선형 배경을 뺀 비선형 회전 비교: +100 MHz **0.124%**, −250 MHz **0.00068%** 차이.
정규화 보정 없음. 선택한 약광 두 점의 부호·크기 검산이며 코팅/강광 인증 아님.

독립 검토 뒤 정상상태 실패/탐색 경계/근사 범위의 compact status를 연결.
작은 자기장 신호의 sinusoidal peak=√2·RMS가 검사한 국소 장 범위를 넘으면
SNR을 보류하고 ASD는 유지. 검출 설정 변경으로 원자 재풀이하지 않는 계약 검사 포함.

### 디튜닝 수치 실험

재현 명령:

```powershell
python -m analysis.magneto_budker_20261030.validate --line D1
python -m analysis.magneto_budker_20261030.validate --line D2
```

[D1 JSON](../../analysis/magneto_budker_20261030/results/d1_validation.json),
[D2 JSON](../../analysis/magneto_budker_20261030/results/d2_validation.json).
소스 SHA-256·전체 입력·상태·기울기·투과율·ASD·소요시간 보존.
기준 25Δ·속도 적분 차수6·h=10⁻⁵ µT와 세분 49Δ·차수10·h=5×10⁻⁶ µT 비교.
조건: ⁸⁵Rb·20 °C·n=4.5×10⁹ cm⁻³·100 mm·d=3 mm·η=0.8·ENBW=1 Hz.
가정 이완율: exit 50000/s, return 20/s, wall reset 3/s, light reset 1/s.
벽 완화·spin exchange를 독립 실측으로 분해한 값 아님.

| 조건부 결과 | D1·I=1 mW/cm² | D2·I=4.5 mW/cm² |
|---|---:|---:|
| 세분 격자 최소점 | −300 MHz | +300 MHz |
| PSN 최소 감도 | 0.567 fT/√Hz | 0.1675 fT/√Hz |
| 중심 대비 SNR 개선 | 15.39배 | 2.089배 |
| 25→49점 최소 감도 변화 | −0.073% | −0.286% |
| 공통 유효 Δ의 기울기 최대 차이/peak | 0.410% | 0.333% |
| 속도 보존 대조군의 최소점 | −100 MHz | +100 MHz |

이는 **선언한 centroid 기준** 결과. 논문의 흡수 중심 기준과 미정렬.
비공진 이득과 속도 혼합에 따른 스펙트럼 변화 구현 확인; 논문 정량 재현 아님.
세분 격자 minimum 이동은 두 선 모두 0 MHz이나 50 MHz 사이 연속 최적은 미평가.

입력 광세기를 스캔 최소 투과율만큼 낮춘 균일 세기 시나리오에서도
D1 최소점의 감도 변화는 약 −1.1%, D2는 약 −9.4%.
원래 입력 광자 budget을 고정하고 원자 광세기 가정만 바꾼 민감도 확인.
공간 전파 해·오차 상계 아님. 특히 D2 최소 T≈0.683이므로 균일 세기 근사는
정량 검증 필요. 보수적 가상 여기 bound로 광학 소거 범위도 별도 표시.

### 최초 구현 사용법 — 통합 전 이력

자기 스킴의 **Defaults → Budker D1 conditions (model)** 또는 D2 선택.
혹은 **Calculation → Budker NMOR detuning scan**. 세 패널과 JSON/CSV 내보내기 사용.
단일 전이 B 스캔은 기본 선택으로 유지. 정량 후속 범위는
[체크리스트](../checklist.json)의 `magneto-budker-quantitative-regime`에 기록.

실제 Chromium 브라우저에서 D1 기본값 선택·3개 패널·상태 표시·JSON 다운로드/파싱 확인.
브라우저 오류 0건. [기록](browser_smoke.json),
[화면](../../output/playwright/magneto_budker.png).
최초 AppTest 전체 앱 검사에서는 기존 Guide v2 컴포넌트 등록의 harness 오류로 중단
([기록](ui_smoke.json)); 실제 브라우저에서는 재현되지 않음. 새 모드의 canvas와
export는 별도 집중 검사에서도 통과. 검사용 서버·브라우저는 종료.
