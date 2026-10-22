# GABES 일일 Scheme 물리 검토 — 2026-09-20

## 선택과 판정

- Asia/Seoul 현지 날짜 2026-09-20. n = (20 mod 5) + 1 = **1**.
- 대상: **Absorption spectroscopy (OD / SAS)**, registry key sas. 검토 HEAD 2998dc1과 현재 작업 트리. 기존 미완료 변경은 보존.
- **판정:** 실제 원자 물리 구현. 약한 probe의 Rb D선 OD, hyperfine 위치·상대 세기, Doppler 폭, 정성적 SAS crossover 및 실험 조건 탐색에 유용. 절대 pumped-SAS contrast·lock slope, 고온 Cs/Ne 셀, 코팅 셀, 자기장·편광, 유한 probe의 동시 pump/probe 동역학을 그대로 계측 표준값으로 쓰기에는 부족.

2026-09-15의 [Scheme 1 보고서][prior] 이후 확인한 핵심 소스는 같은 상태다. 이번 실행은 새 기능 제안보다 기존 미해결 신뢰 계약과 비용을 재평가했다.

## 다섯 scheme의 현재 정의와 순서

단일 등록 순서: [scheme registry][registry]. 내부 ODScheme은 등록되지 않은 검증 primitive다.

| 번호 | key | 정의 |
|---:|---|---|
| 1 | sas | OD / SAS: pump-off Doppler 흡수, pump-on 포화 흡수·hyperfine pumping |
| 2 | lambda | Λ coherence: EIT / AT / CPT |
| 3 | rydberg_eit | Rydberg-EIT 전기장 계측 |
| 4 | magneto | Hanle / EIA / NMOR |
| 5 | fwm | FWM: seeded gain 지표 / SFWM biphoton |

등록 객체의 실제 name/title을 Python에서 다시 조회했다. [README scheme 표][readme-scheme]와 일치한다.

## 선행 개선안 검색

- [2026-09-15 Scheme 1 보고서][prior]와 이전 일일 보고서, [통합 checklist][check-collision], [코드의 Ne TODO][ne], [pole 이식 기록][dense]을 먼저 확인했다.
- 별도 이름의 로컬 issue/TODO/proposal 파일은 검색되지 않았다. 원격 origin의 [공개 GitHub Issues](https://github.com/Shake2313/fwm-squeezing-app/issues)는 2026-09-20 API 조회 시 전체 0건. **개선안 부재가 아니다.** 실질적인 issue 기록은 checklist에 있다.
- 완료된 항목: [sub-Doppler resolution gate][check-resolution], [paraffin population memory][check-coating], [약한 probe 분산][check-dispersion], [밀도·Doppler 온도 분리][check-temp]. 새 제안으로 재등록할 이유 없다.
- 미완료 핵심: [충돌 계수·pressure shift][check-collision] P1 ready, [species/line별 Rabi][check-rabi] P1 ready, [dense velocity 연산][dense] P3 ready. [실험 spectrum fitting][check-fit]과 [segmented cell][check-segment]은 큰 작업, [Zeeman/polarization reference][check-zeeman]과 코팅 셀 확장은 연구 과제다.

## 구현된 물리와 실험 참고 범위

1. [species 자료·manifold][species-manifold]는 85Rb/87Rb/133Cs D1/D2, hyperfine A/B, 허용 전이, Wigner 세기, 분기 붕괴를 구성한다. [적분 흡수][species-area]는 C_F²에 이미 포함된 mF 합을 두 번 세지 않아 85Rb D1의 바닥 manifold 면적비 7/5를 준다. [6개 isotope/line 항등식 테스트][test-sum]는 수명에서 얻는 oscillator sum과 축퇴도 비를 검증한다. Steck의 [Rb85 원자 자료](https://steck.us/alkalidata/rubidium85numbers.pdf)가 적절한 외부 원자 상수 근거다.
2. [SAS 응답][sas-response]는 pump detuning Δ+kv와 역행 probe detuning Δ−kv를 같은 속도군에서 결합한다. CG 분기와 transit reset으로 hyperfine optical pumping과 crossover가 나온다. 이는 [Smith–Hughes 원 논문](https://home.sandiego.edu/~severn/p480w/rb_absorption_crossover_resonances_amjp.pdf)의 실험 물리와 방향이 맞다. 다만 현재 F별 lumped 상태는 mF·편광을 풀지 않는다.
3. pump-off 경로는 펌프 population solve를 생략하고 Maxwell 속도 평균·Beer–Lambert 전파를 쓴다([compute][sas-compute], [T와 OD][obs]). 약한 probe 절대 흡수의 활용 조건은 엄격하다. [Siddons 등 실험·이론 비교](https://arxiv.org/abs/0805.1139)는 충분히 약한 probe에서 높은 일치를 보인다. 현재 GABES의 AutoOD 일치는 **같은 계열 계산기와의 소프트웨어 비교**이고 독립 보류 실험 spectrum 검증은 아니다.
4. 약한 probe의 흡수와 분산은 [같은 복소 line profile][sas-response] 및 같은 pump population에서 나온다. [Voigt 기준·부호 테스트][test-dispersion]가 있다. finite probe에서는 [2-level 포화와 Gaussian z/r 전파][probe-prop]를 쓰며 pump-prepared population을 고정한다. 이 출력은 특히 D1/Cs와 동시 유한 pump/probe에서 정성적. [현재 범위 표시][sas-scope]는 이 제한을 설명한다.
5. 두 온도, 코팅 population memory, CSV의 상대/절대 교정을 구분한 점은 실험 설계에 유용하다([두 온도][sas-temps], [코팅][sas-coated], [CSV][csv]). 코팅의 25.1 ms는 한 셀의 참조값이지 범용 T1이 아니다.

## 현재 신뢰 한계

| 항목 | 확인된 현상 | 계측 판단 |
|---|---|---|
| OD convention | [공통 함수][obs]는 −log10(T)=αL/ln10. 원자 흡수 문헌의 optical depth αL=−ln(T)와 정확히 ln10배 차이 | α/T 계산 오류는 아니다. export에 자연로그 optical depth와 decadic absorbance를 명시해야 한다. |
| Pump power | [power→Rabi][rabi]는 모든 isotope/D선에 같은 85Rb D1 I_sat=4.484 mW/cm²를 사용 | P 축은 경향 탐색용. 절대 saturation/contrast 보정 근거 부족. |
| 충돌 | [self broadening][self]은 iso 인자를 받지만 단일 Rb 계수, [Ne][ne]는 3.91 MHz/Torr broadening만 추가. [probe Γ_eff와 pump Liouvillian][sas-compute]의 optical dephasing 계약도 다름 | 고온 Cs·buffer 셀의 정량 선폭과 crossover 불신. [Rb D1 self broadening 측정](https://arxiv.org/abs/1107.3092)은 적어도 Rb D1 계수의 물리적 근거다. |
| 수치 해상도 | 기본 1401점의 sub-Doppler feature는 FWHM당 2.342점. [UI gate][test-resolution]는 올바르게 resolution-limited를 출력 | 선폭·lock slope를 정량 인용하려면 고정 창 grid 수렴 필요. |
| 공간·상태 | pump는 입사 power로 한 번 계산. [uniform slab][sas-compute]과 [T=exp(−αL)][obs]로 pump depletion/gradient 없음. Zeeman·편광도 없음 | 고 OD pumped SAS, field-sensitive lock point에 직접 사용 불가. |
| cache·출처 | [cache_version 7][sas-class]은 2026-09 thermal normalization 수정 뒤 그대로. [Streamlit cache key][cache]에 사용 | 프로세스 내 이전 결과·export provenance 구별 실패 위험. |

## 이번 실행의 직접 확인

기본값: natural Rb D1, 40 °C, 75 mm, pump 0.5 mW, probe 0 µW, 1401 scan점. 외부 테스트 병렬 실행 중 측정된 wall time은 하드웨어 부하 영향이 있으므로 성능 보증치가 아니다.

| 경로 | 직접 재계산 결과 |
|---|---:|
| 기본 SAS | Tmin 0.1180495075; peak αL 2.136651188; peak −log10(T) 0.927935821 |
| 기본 SAS feature | FWHM 18.86836 MHz; 2.34181 samples/FWHM; resolution-limited |
| pump-off OD | Tmin 0.105095655; peak αL 2.252884348; peak −log10(T) 0.978415241 |
| 기본 SAS 프로파일 | cProfile 0.715 s; 두 isotope의 응답 0.695 s 누계. 그 안에서 np.interp 0.147 s, velocity average 0.117 s, pump populations 0.071 s |

[이전 보고서][prior]의 4×/8× fixed-window 수렴 검사는 [현 테스트][test-refine]에 남아 있다. 기본 grid의 resolution gate는 수렴 증명과 다르다.

## 기존 제안의 비용과 물리 보존성

| 기존 제안 | 계산 비용 | 물리·수치 의미 |
|---|---|---|
| OD/absorbance 명칭 분리, cache version 갱신, 기준 source hash 보존 | 추가 solve 0; O(Nscan) 변환과 O(1) metadata | **물리 exact.** 우선 trust 수정. [obs][obs], [cache][cache] |
| isotope/line/polarization Rabi, 출처 있는 collision table·shift·dephasing | 계수 조회 O(1), 같은 matrix 차원·solve 수. 자료 조사·검증 비용은 큼 | **실제 물리 교정은 거의 실행비용 없이 가능.** pumped 수치는 올바르게 바뀐다. zero-pressure·기존 reference anchor를 별도 gate. [checklist][check-collision] |
| dense scan×velocity fused loop | O(Nscan Nv) 유지, 큰 임시 배열·memory traffic 감소 | 현재 4σ/1 m/s measure와 np.interp 규칙을 지키면 결과를 부동소수점 허용오차 내 보존. 3× 목표는 아직 측정 전. [제안][dense] |
| analytic Voigt / pole-residue 평균 | velocity 차원 제거 가능 | 연속 Maxwell 적분은 현 이산 4σ measure와 정확히 같지 않다. pump population의 Γ/8 interpolation도 바뀔 수 있어 순수 코드 최적화로 분류하지 않는다. [응답][sas-response] |
| 자동 fitting / segmented pump depletion | fitting 약 N_eval×compute; segment는 적어도 N_z배 | 물리 활용도 있지만 무비용 불가. 독립 spectrum, 파라미터 식별성, z 수렴 필요. [fitting][check-fit], [segment][check-segment] |
| full VCC·mF/편광·코팅 transport | velocity·state·geometry 차원 증가 | 무시할 만한 비용으로 기존 physics를 보존하는 변경 아님. 측정 조건과 보류 데이터가 먼저. [Zeeman][check-zeeman] |

## 동작을 바꾸지 않는 코드 최적화

1. **pump-off scalar fast path.** [현재 helper][sas-average]는 population=1에도 큰 probe profile 배열을 1과 곱한 뒤 @weights. 이 경우 곧바로 profile@weights 가능. 현재 PC의 1401×1403 배열 미세 벤치: helper median 6.17 ms → 0.70 ms, 출력 array_equal. 전체 compute 이득은 이보다 작다.
2. **pump-off finite-probe 중복 matvec 제거.** [transition_alpha/thermal_alpha][sas-response]는 pump-off에서 같은 A_t×(Lp@wt)를 두 번 계산. thermal_alpha=transition_alpha 재사용 가능. 추가 physics 가정 없음.
3. **같은 보간 규칙의 compiled fused loop.** [probe_base/deff_grid/pop_at/Lp/Dp][sas-response]는 isotope마다 약 1401×1400 float64 임시 배열을 반복 생성한다(배열당 약 15 MiB). 직접 프로파일의 병목도 matrix solve보다 메모리 연산이다. 두 quadrature, coated/finite probe, Numba on/off 결과를 모두 비교해야 한다.
4. [Gaussian 전파][probe-prop]의 32점 Legendre node/weight 재생성은 상수 cache 가능. headless 재표시 최적화이나 기본 물리 solve에는 영향 작다.

## 문서·예제·테스트

- [README][readme-scheme]는 현재 보정된 thermal strength와 finite probe 범위를 상당 부분 설명한다. 반면 [in-app About][sas-about]의 “lab AutoOD <0.1%”는 independent experimental validation처럼 읽힐 수 있다.
- [User Guide OD/SAS 설명][guide]은 finite-probe 범위를 natural Rb D2로 좁혀 적고 기본값을 0 mW로 표기한다. 실제 [schema][sas-schema]는 모든 기본 alkali D선에 0–100 µW를 노출한다. [OD 그림][od-image]과 [SAS 그림][sas-image]은 2026-06-08 생성된 개념 예제다. 그림의 깊이와 base-10 OD peak는 오늘 기본값(각각 OD Tmin 0.1051, SAS 0.1180; OD10 peak 0.978/0.928)과 맞지 않는다. 독립 실행 가능한 Scheme 1 예제 script/config도 없다.
- [테스트][test-sum]는 oscillator sum, 같은 속도군, pumping crossover, 코팅 불변성, 분산 부호, 해상도, 두 온도, [finite probe software anchor][test-probe]를 다룬다. 다만 [pump-off AutoOD 비교][test-autood]는 같은 저장소의 ODScheme과 비교하고, 별도 AutoOD reference는 [gitignore][ignore] 아래라 일반 clone의 독립 재현 자산이 아니다.
- 실행: 관련 4개 모듈 **132 passed / 75.41 s**; 별도 AutoOD-NatRbD2 reference **29 passed / 4.62 s**. 지침의 전체 **python -m pytest -q**는 AutoOD GUI 참조 테스트 수집 중 optional sympy 부재로 **1 collection error / 36.11 s**. Scheme 1 실패로 해석하지 않는다.

## 권고

1. OD 자연로그/상용로그를 명시하고 cache_version·export provenance를 갱신한다. 예제 그림과 Guide의 수치·probe 범위를 현재 코드로 재생성한다.
2. 기존 P1 Rabi와 collision 계약을 sourced coefficient, 단위, polarization, pressure-shift 부호, pump optical dephasing까지 맞춰 구현한다. 새 실험 비교 없이는 “검증됨”보다 “소프트웨어 기준 일치”로 표기한다.
3. 이후 같은 속도 quadrature를 보존하는 fused arithmetic을 benchmark한다. 기본 1401점의 SAS lock 수치는 고정 창 grid 수렴 전까지 reference 수치로 내보내지 않는다.

### 근거 파일

[registry]: <G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/__init__.py:19>
[readme-scheme]: <G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/README.md:8>
[prior]: <G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/daily_report/2026-09-15_scheme-1_od-sas.md:1>
[sas-class]: <G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:175>
[sas-schema]: <G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:188>
[sas-about]: <G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:296>
[sas-scope]: <G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:113>
[sas-temps]: <G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:56>
[sas-compute]: <G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:365>
[sas-response]: <G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:539>
[sas-coated]: <G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:1090>
[sas-average]: <G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/schemes/sas.py:996>
[species-manifold]: <G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/species.py:337>
[species-area]: <G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/species.py:248>
[rabi]: <G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/species.py:272>
[self]: <G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/species.py:210>
[ne]: <G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/constants.py:79>
[obs]: <G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/observables.py:579>
[probe-prop]: <G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/probe_saturation.py:166>
[csv]: <G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/gabes/experimental_csv.py:214>
[cache]: <G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/streamlit_app.py:640>
[test-sum]: <G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_sas.py:91>
[test-autood]: <G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_sas.py:260>
[test-resolution]: <G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_sas.py:406>
[test-refine]: <G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_sas.py:442>
[test-dispersion]: <G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_sas.py:594>
[test-probe]: <G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/tests/test_probe_saturation.py:212>
[check-collision]: <G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:2630>
[check-rabi]: <G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:2670>
[check-resolution]: <G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:2691>
[check-coating]: <G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:2712>
[dense]: <G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:3320>
[check-dispersion]: <G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:3876>
[check-temp]: <G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:3904>
[check-fit]: <G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:3932>
[check-zeeman]: <G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:3998>
[check-segment]: <G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/checklist.json:4032>
[guide]: <G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/Userguide/GABES_User_Guide_v2.html:553>
[od-image]: <G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/Userguide/userguide_assets/od.png>
[sas-image]: <G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/docs/Userguide/userguide_assets/sas.png>
[ignore]: <G:/내 드라이브/Coding/GABES - Generic Atomic Bloch Equation Solver/.gitignore:19>
