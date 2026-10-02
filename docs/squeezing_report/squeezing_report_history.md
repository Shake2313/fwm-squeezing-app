# 스퀴징 보고서 변경·분석 기록

기록일: 2026-09-08. 이 파일은 버전 이력, 코드 검토, 과거 계산의 재평가와 출처 정정 사항을 보존한다. 최신 이론·최적점·실험 비교·변수 허용범위·코드 한계는 `squeezing_report_v7.tex`와 그 PDF에 독립적으로 서술한다. 본 기록을 읽어야 최신 보고서를 이해할 수 있도록 구성하지 않는다.

## 문서 편집 원칙

첫 v7 초안은 7월 이후 커밋 감사와 v6 결론의 재판정에 많은 지면을 사용했다. 사용자가 원한 것은 GABES와 스퀴징 이론으로 무엇을 계산했고, 실험과 어느 정도 부합하며, 변수별 허용범위와 결합 관계가 무엇인지 설명하는 완결된 과학 보고서였다. 이에 따라 같은 v7의 본문을 다음 원칙으로 다시 작성한다.

- `v1`부터 `v6`까지의 파일은 과거 버전 기록으로 보존한다.
- 최신 보고서에는 필요한 정의·수식·계산 조건·실험 비교·한계·개발 과제를 모두 담는다. 과거 보고서를 먼저 읽도록 요구하지 않는다.
- GitHub 커밋 목록, 이전 주장에 대한 유지·철회 판정, 출처 오류의 상세 분석은 이 Markdown 기록으로 분리한다.
- 문헌 실측값, 실측 이득을 입력한 조건부 이론값, GABES 평균장 결과, 역추정한 유효 파라미터를 구분한다. 실측 이득을 넣어 실측 스퀴징과 가까워진 것을 독립적인 GABES 예측 성공으로 표현하지 않는다.
- 논문들 사이의 설정 분포를 한 장치의 tolerance로 사용하지 않는다. 한 장치에서 다른 조건을 고정한 스캔과 명시한 성능 기준으로 허용범위를 정한다.
- 최신 보고서는 변경 이력 설명으로 시작하지 않고, 물리적 문제와 계산·실험 결과로 시작한다. 코드 미완료 사항은 그것이 제한하는 예측량 및 향후 개발 과제와 연결한다.

이 기록을 작성하면서 원본 reference CSV, 기존 TeX·README, 물리 엔진은 수정하지 않았다. 아래의 CSV 정정 목록은 보고서에서 채택한 해석과 후속 원자료 정리의 근거다.

## GitHub 검토 기준

첫 감사의 날짜 경계는 **2026-07-22 00:00 KST 이후**다. 보존된 `v7_assets/commits.json`의 전체 26개 항목을 아래에 수록한다. 이 목록은 최초 감사 시점의 snapshot이며 이후 GitHub가 진전해도 자동으로 최신화된 기록으로 해석하지 않는다.

- 경계 직전 revision: `bd44273089459217a5114e9897665f0499c1ae74`.
- 최초 감사의 GitHub `main` 및 로컬 HEAD: `c0c46f0a8b9db9c023f4c5b6ffe5185769f57018`.
- [전체 GitHub 비교](https://github.com/Shake2313/fwm-squeezing-app/compare/bd44273089459217a5114e9897665f0499c1ae74...c0c46f0a8b9db9c023f4c5b6ffe5185769f57018).
- v6 표지 날짜는 2026-07-06, tolerance 보충은 2026-07-08이었다. 7월 23일 변경은 분석 경로 재배치이며 새로운 최적화 계산이 아니었다.
- 첫 대상 커밋은 2026-07-23 02:44 KST, UTC로는 7월 22일이다. 시간대에 따른 표기 차이를 별도 물리 변경으로 보지 않는다.

| 날짜·시간 (KST) | 커밋 | 원문 subject |
|---|---|---|
| 2026-07-23 02:44:05 | [a82bbf7](https://github.com/Shake2313/fwm-squeezing-app/commit/a82bbf71b4f9554fb3bc54684afb6f970fa7df0f) | chore!: reorganize docs and analysis artifacts |
| 2026-08-04 18:04:39 | [f3898a5](https://github.com/Shake2313/fwm-squeezing-app/commit/f3898a57629cb3b7d38f42e5cbbe0e819d9edaaf) | feat(sabes): add a lab-facing simulator for the Sim et al. squeezing setup |
| 2026-08-06 18:23:11 | [3cb6b14](https://github.com/Shake2313/fwm-squeezing-app/commit/3cb6b1492a2ed33f9055ee8ab9be2943a1d3c3aa) | feat(sabes): plan the optical-table UI, and land its layout model and canvas |
| 2026-08-06 18:31:59 | [130706f](https://github.com/Shake2313/fwm-squeezing-app/commit/130706f89e559a0a42318561030d1621a6b2d5be) | fix(sabes): retitle the real browser tab, not the Community Cloud wrapper |
| 2026-08-06 18:32:47 | [df1837e](https://github.com/Shake2313/fwm-squeezing-app/commit/df1837e98056b3eb8fcc852018a6fb2b64e94960) | docs(sabes): record the canvas spike as confirmed on the deployment |
| 2026-08-06 18:52:13 | [5eadd08](https://github.com/Shake2313/fwm-squeezing-app/commit/5eadd088c1495e8516013abb8f2345986d049649) | feat(sabes): draw the optical table, and make it clickable |
| 2026-08-06 20:19:23 | [942f206](https://github.com/Shake2313/fwm-squeezing-app/commit/942f206626e17149b9436b58079a5f608b07aa4d) | feat(sabes): move every parameter onto its optic |
| 2026-08-07 01:17:32 | [e1e3cdc](https://github.com/Shake2313/fwm-squeezing-app/commit/e1e3cdc196d619cbe8df22d5059d4b170b50b3dc) | feat(sabes): add virtual instruments driven by the real model |
| 2026-08-07 11:24:22 | [0d16b69](https://github.com/Shake2313/fwm-squeezing-app/commit/0d16b69f294c13378f4c2f40bd3191b38a35cddd) | feat(sabes): give every post-amplifier optic its own transmission |
| 2026-08-18 16:00:36 | [c964a72](https://github.com/Shake2313/fwm-squeezing-app/commit/c964a724b1ad370e53af1a89c78f8df04fd37983) | fix(ui): remove duplicate top-bar User's Guide button |
| 2026-08-28 13:35:39 | [5cc7eee](https://github.com/Shake2313/fwm-squeezing-app/commit/5cc7eeef18818fa421cd528fee1da385c3cba315) | feat(sas): add paraffin population memory |
| 2026-08-28 18:09:13 | [7244f25](https://github.com/Shake2313/fwm-squeezing-app/commit/7244f25f249d92b8b07fe905a38b5cff1010c06b) | feat: harden model claims and validation |
| 2026-08-28 18:11:33 | [d03e024](https://github.com/Shake2313/fwm-squeezing-app/commit/d03e0244c0224885f9b268606269303121e29f9b) | fix(analysis): keep manifest paths portable |
| 2026-08-28 18:32:08 | [153329a](https://github.com/Shake2313/fwm-squeezing-app/commit/153329ac7d754ebca36cba5752f45be5b2e823f1) | fix(analysis): gate low-pump claims |
| 2026-08-28 18:47:58 | [a3d44d4](https://github.com/Shake2313/fwm-squeezing-app/commit/a3d44d431dae7e787cf6b660f6be2964100cb856) | fix(analysis): adapt transit diagnostic |
| 2026-08-28 19:04:37 | [ce60b5c](https://github.com/Shake2313/fwm-squeezing-app/commit/ce60b5c029495400e014c0f1ce7d281b17046d60) | docs: add reproducible analysis artifacts |
| 2026-08-31 08:55:10 | [567d6c5](https://github.com/Shake2313/fwm-squeezing-app/commit/567d6c5e971d53f78523c8731a1c9756df742b5b) | docs: add daily physics reviews |
| 2026-08-31 13:47:26 | [c922ca2](https://github.com/Shake2313/fwm-squeezing-app/commit/c922ca2703bfa3a3eda598be33de04985325130f) | fix(app): recover stale CSV imports |
| 2026-09-02 20:35:54 | [2595740](https://github.com/Shake2313/fwm-squeezing-app/commit/2595740f2b1a40e27c49c840befd47a2a9bf2fd8) | feat(od-sas): add the dispersive quadrature and split the two temperatures |
| 2026-09-02 20:36:15 | [39684ba](https://github.com/Shake2313/fwm-squeezing-app/commit/39684baa5e2633dde7338cdd0051e15d1efb95e8) | fix(plots): pin a font stack per figure so isotope glyphs survive |
| 2026-09-02 21:24:28 | [c5527fe](https://github.com/Shake2313/fwm-squeezing-app/commit/c5527fef82212b25736d6057f7e63ae29b063f81) | feat(od-sas): add probe saturation, hide Generic |
| 2026-09-04 14:33:54 | [f33b3e3](https://github.com/Shake2313/fwm-squeezing-app/commit/f33b3e3ccc92c4cf1603d3b437c88e0b5b9d1e6f) | feat(fwm): simplify controls and readouts |
| 2026-09-05 00:48:03 | [d98bdf4](https://github.com/Shake2313/fwm-squeezing-app/commit/d98bdf47d5a1f57663b97463c88cce68aa0c8264) | perf(core): streamline solves and readouts |
| 2026-09-07 10:33:17 | [8d8bbe7](https://github.com/Shake2313/fwm-squeezing-app/commit/8d8bbe77eeefee9aacbda73132784fa04b02aa08) | fix(fwm): state gain-only trust boundary |
| 2026-09-07 10:33:32 | [8f97f05](https://github.com/Shake2313/fwm-squeezing-app/commit/8f97f0577a8da989affd8afbf10e0d317a141fa2) | feat(app): export results with provenance |
| 2026-09-07 10:36:47 | [c0c46f0](https://github.com/Shake2313/fwm-squeezing-app/commit/c0c46f0a8b9db9c023f4c5b6ffe5185769f57018) | docs(checklist): record ready-item progress |

## 스퀴징 관련 코드 변경의 의미

### 원자 정규화·열적 이완·전파

주요 수정은 [7244f25](https://github.com/Shake2313/fwm-squeezing-app/commit/7244f25f249d92b8b07fe905a38b5cff1010c06b)에 집중되어 있다. Trace-one 원자 밀도행렬이 공급하는 pump-modified population에 외부 열평형 ground population을 다시 곱하던 중복을 제거했다. 85Rb D1 hyperfine 선세기는 구동, 분극 readout, spontaneous-emission branching에 반영되고, 거시 구조 계수는 `1/[2(2I+1)] = 1/12`로 정리되었다.

상속된 reference residual `0.74`는 이 정규화 수정 이후 독립 실험으로 다시 피팅된 상수가 아니다. 추가 mode-overlap, polarization, Zeeman-participation penalty의 기본값은 각각 1이다. 이 세 항은 미해결 결합 손실을 구분하는 입력이며 `0.74`를 세 가지 실제 측정 효율로 분해한 결과가 아니다.

Ground relaxation은 coherence-only 감쇠와 thermal transit reset을 구분하도록 바뀌었다. 현재 thermal reset은 `L_t rho = gamma_t (rho_th Tr(rho) - rho)`, `rho_th = diag(5/12, 7/12, 0, 0)` 형태다. Maxwell 감쇠 부호와 phase-mismatch 부호를 정리하고, Option-A 전파에서 bare wavevector를 사용해 대각 susceptibility에 이미 포함된 분산을 mismatch에 중복 반영하지 않도록 했다.

Ultra의 segment 전파와 pump-budget 감소는 자기일관적인 depleted three-field Maxwell–Bloch 해법과 구분해야 한다. 표시 이득의 post-hoc saturation cap도 원자 상태와 세 광장을 동시에 재계산하는 구현이 아니다.

### Floquet와 별도 응답 reference

일반 finite-order Floquet 계산, continued-fraction과 독립 dense-block reference, 기본 `N_F=3` 및 `N_F=2`와의 전체 표시 scan 비교가 도입되었다. Complex susceptibility·transfer, gain, wrapped phase, gain 최대 위치의 기준을 통과해야 `CONVERGED`를 반환한다. 밀도행렬 trace와 positivity만으로 차수 절단을 인증하지 않는다.

표준 red-seed Raman minus branch에 대해 self-consistent pump-only 정상상태와 trace-zero Nambu 약한 응답을 별도 reference로 계산한다. Production의 finite-seed Floquet 경로를 자동 대체하는 기능은 아니다. 기존 plus branch는 static-pump gauge parity 검증을 통과하지 않아 이 reference의 인증 범위에서 제외된다.

독립 `(v_z,v_x)` Maxwell 평균을 사용하는 비공선 2D Doppler reference도 추가되었다. Gold 온도·각도에서 Raman RMS 폭 약 1.380 MHz를 재현하는 별도 검사가 있지만 production 평균은 1D다. 이 두 reference 모두 microscopic Langevin diffusion을 제공하지 않는다.

### 이득과 스퀴징의 해석

Field transfer를 photon-flux amplitude로 변환해 정준성·commutator defect를 검사한다. Power gain 차이와 photon-flux gain 차이는 엄밀히 다른 양이다. 열린 소산계에서 drift만의 commutator defect는 reservoir 보완이 필요하다는 뜻이며, 실제 스퀴징 부재의 증명도 스퀴징의 인증도 아니다.

[f33b3e3](https://github.com/Shake2313/fwm-squeezing-app/commit/f33b3e3ccc92c4cf1603d3b437c88e0b5b9d1e6f)와 [8d8bbe7](https://github.com/Shake2313/fwm-squeezing-app/commit/8d8bbe77eeefee9aacbda73132784fa04b02aa08)은 표시·제어와 주장 범위를 정리했다. dB 출력은 `Squeezing indicator` 또는 gain-only 진단이며, frequency-dependent microscopic atomic covariance와 독립 실험 검증이 갖춰진 physical squeezing spectrum으로 표시하지 않는다.

### SABES 장치·검출 경로

[f3898a5](https://github.com/Shake2313/fwm-squeezing-app/commit/f3898a57629cb3b7d38f42e5cbbe0e819d9edaaf) 이후 optical-table UI, 소자별 파라미터, post-amplifier 투과율과 가상 계측기가 추가되었다. [d98bdf4](https://github.com/Shake2313/fwm-squeezing-app/commit/d98bdf47d5a1f57663b97463c88cce68aa0c8264)에는 계산 재사용·readout 정리와 함께 unwanted-mode shot noise 및 `I² RIN` 진단이 포함된다. EOM/etalon mode ledger와 detector headroom을 추적하는 개선이며, 원자 공분산과 측정된 RF transfer가 자동으로 완성된 것은 아니다.

[8f97f05](https://github.com/Shake2313/fwm-squeezing-app/commit/8f97f0577a8da989affd8afbf10e0d317a141fa2)은 결과·파라미터·revision export를 추가했고, [c0c46f0](https://github.com/Shake2313/fwm-squeezing-app/commit/c0c46f0a8b9db9c023f4c5b6ffe5185769f57018)은 checklist 진행 상태를 기록했다. Absorption/SAS/Rydberg/biphoton 기능의 개선을 seeded FWM의 양자잡음 검증으로 전용하지 않는다.

## 첫 감사에서 보존한 v6 판정

이 절은 최초 v7 초안에 있었던 역사적 검토를 보관한다. 아래 값과 권고를 최신 과학 보고서의 새 최적화 결과로 사용하지 않는다.

| v6 내용 | 당시 감사 판정 | 기록상 해석 |
|---|---|---|
| 표준 branch의 beat 부호 | 유지 | `Omega_beat = Omega_HF - delta`; RF analysis frequency와 구분 |
| OPD −1.50 GHz, 110°C, TPD −280 MHz에서 −8.102 dB | 보관된 모델 결과 | 당시 `N_F=1`과 hardened 목적함수의 optimum |
| 양의 OPD lobe보다 0.32 dB 우세 | 역사적 순위 | 현 모델 또는 실험의 음의 OPD 선호를 증명하지 않음 |
| eta=1에서 −15.621 dB source-limit | 물리적 source-limit 해석 철회 | 검출 효율만 1로 둔 과거 지표 |
| Low-angle trusted −20.70 dB | trusted 해석 철회 | gap gate가 정준성·수치수렴·미시적 확산 검증을 대체하지 못함 |
| TPD/OPD/온도/loss tolerance 창 | 역사적 민감도 | 과거 목적함수의 1D 허용구간 |
| Seed 1–64 µW 무감 | 모델 범위 내 결과 | 실제 seed RIN·불균형·포화에 대한 무감도는 아님 |
| 온도 1–2°C 높게, pump lock 약 10 MHz | 현 실험 권고에서 제외 | 검증되지 않은 optimum·tolerance에 의존 |
| Loss → TPD → OPD → 온도 → seed 점검 순서 | 재작성 | SQL·전자잡음·mode/RIN과 gain map을 함께 측정할 필요 |

여기서 양의 OPD는 pump의 optical detuning 부호이고 Raman plus branch를 뜻하지 않는다. v6와 첫 감사의 대표점 계산은 모두 표준 Raman minus branch였다.

8월의 별도 보관 분석은 v6의 과거 전파·감쇠 조건을 유지하고 Floquet 차수만 변경했다. 이는 현재 Ultra의 실행 결과가 아니다.

| 차수 | 과거 probe gain | 과거 conjugate 계수 | 과거 gain 차이 |
|---:|---:|---:|---:|
| 1 | 22.5913 | 21.9687 | +0.622609 |
| 2 | 19.3415 | 19.4472 | −0.105684 |
| 3 | 19.3415 | 19.4472 | −0.105684 |

차수 2 대비 1의 probe gain 차이는 약 16.8%였다. 같은 과거 목적함수의 재최적화는 TPD −280→−275 MHz, 지표 −8.10236→−7.94704 dB로 이동했다. 과거 cutoff 3→5 sigma 검사에서는 probe gain 약 2.25%, 작은 gain 차이 약 21.5%의 변화가 있었다. 차이의 제곱을 최소화하는 지표는 개별 gain보다 더 민감할 수 있다.

v6의 NPZ와 그림에는 full Git SHA·전체 소스 hash가 모두 포함되어 있지 않고 일부 그림 생성기는 현재 `gabes.observables`를 import한다. 현재 코드를 사용해 그림만 재생성하면 과거 gain과 현재 후처리가 혼합될 수 있으므로, 역사적 재현에는 당시 revision과 후처리의 동시 고정이 필요하다. 이 작업에서 v6 파일과 그림은 재생성하지 않았다.

## 첫 v7 감사의 계산 결과와 검증 범위

첫 감사 조건은 pump 600 mW, seed 8 µW, pump/probe/conjugate 반경 530/330/330 µm, cell 12.5 mm, angle 0.32°, eta 0.8694, reference residual 0.74, 추가 penalty 각각 1, Raman minus branch, Ultra, pump-scatter kappa 0.1, transit `2π × 100 kHz`였다. 각 중심 TPD ±40 MHz를 1 MHz 간격 81점, velocity step 5 m/s·cutoff 3 sigma, `N_F=3`으로 계산했다.

| 고정 좌표 | Probe gain | Conjugate gain | Ultra dB 지표 |
|---|---:|---:|---:|
| 문헌 좌표: +0.90 GHz / 121°C / −8 MHz | 395.058 | 397.279 | −7.801 |
| v6 적색 좌표: −1.50 GHz / 110°C / −280 MHz | 1.09359 | 0.17221 | −1.429 |
| v6 양의 OPD 좌표: +1.40 GHz / 125°C / 0 MHz | 2506.752 | 2526.893 | −6.792 |

이 표의 dB 값은 microscopic quantum-noise prediction이 아니다. 특히 문헌 좌표의 약 −7.8 dB 지표는 실측 gain 약15와 계산 gain 약395의 차이를 해결하지 않는다. 첫 감사의 국소 창에서 문헌 좌표의 gain 최대는 오른쪽 경계였으므로 그 스캔만으로 gain optimum을 확정하지 않았다.

세 국소 81점 scan 모두 `N_F=3/2` 비교 gate를 통과했다. 각 중심 주변 3점에서 velocity step 5→2.5 m/s와 cutoff 3→4 sigma를 동시에 변경했을 때 gain 최대 상대변화는 각각 약 0.193%, 1.453%, 0.0153%였다. 두 오차원을 분리한 최종 수렴 인증은 아니며, 2D 원자 운동·공간모드·양자잡음·실험 재현 검증도 아니다.

최초 감사 시점에는 `python -m pytest -q`가 652개 통과했다. 그 실행 기록은 `v7_assets/validation.json`에 보존되어 있다. 최신 scientific 재작성의 검증 결과는 최신 생성기·manifest와 최신 README의 기록을 우선한다. 최초 실행의 주요 source hash, HEAD와의 raw/줄바꿈 정규화 비교, 기존 staged/unstaged 상태 및 실행 전후 source 불변성은 `v7_assets/audit_manifest.json`에 있다.

## Reference CSV의 정정·출처 쟁점

검토 대상은 `references/fwm_squeezing_paper_parameters.csv`의 10개 논문이다. Gold와 Allen은 저장소의 로컬 최종 PDF를 읽었고, 나머지는 공개 원문 또는 저자 공개본으로 대조했다. 원자료의 미확인 항목을 0 또는 독립적으로 측정된 파라미터로 채우지 않는다.

1. **Gold beam-size convention.** 최종 논문은 530/330 µm를 명시적으로 **1/e² radius**라 부른다. CSV의 “README treats ... despite paper wording”이라는 유보는 최종본 문구와 맞지 않는다. 보고서는 반경으로 사용한다. [Sim 최종 원문, Experimental schematic](https://doi.org/10.1038/s41598-025-86479-w).

2. **Gold gain의 유효숫자와 측정 조건.** Gain spectrum은 약16, IDS 운전점 서술은 약15이며 동일 문단의 출력111 µW/input8 µW는 13.875다. CSV15.5는 정밀한 측정값으로 쓰지 않고 보고된 약15–16 범위로 취급한다. 출력비·gain scan·별도 noise 측정의 반올림 또는 측정면 차이가 원문에서 완전히 해소되지 않았으므로 임의로 원인을 확정하지 않는다. [Sim, Gain spectra / Squeezing in 85Rb](https://doi.org/10.1038/s41598-025-86479-w).

3. **Gold loss와 slope의 의미.** 논문은 optical loss5.5(2)% 및 system loss8.0(2)%를 합쳐 total13.5(2)%를 보고한다. 코드의 `0.945 × 0.92 = 0.8694`와 논문 총효율0.865는 같은 계산이 아니다. 원문의 system loss를 독립 photodiode QE 측정으로 단정하지 않는다. Slope ratio0.134, 약−8.7 dB는 상수 배경을 줄인 검출면 측정값이며 loss-corrected intrinsic squeezing이 아니다. 0.134가 명목 loss floor0.135보다 작다는 차이는 반올림·불확도 및 서로 다른 측정 정의를 고려해야 한다. 원문의 multimode bandwidth 변화 “0.8 GHz”는 전체 MHz 대역폭과 모순되어 정량 근거에서 제외한다. [Sim, Squeezing / Multimode effect](https://doi.org/10.1038/s41598-025-86479-w).

4. **Wu의 출력 파워가 서로 다른 run을 가리킴.** CSV의 probe129 µW와 conjugate17 µW는 twin-beam 한 쌍의 출력이 아니다. Fig.2 caption에서 seed6.4 µW와0.9 µW 조건의 같은 출력 채널이 각각129 µW와약17 µW다. 이 두 수로 power imbalance를 만들지 않는다. 또한 Fig.2의 gain20 조건과 최종 dual-seeding 저주파 결과를 구분한다. [Wu 원문, Fig.2](https://pmc.ncbi.nlm.nih.gov/articles/PMC8658049/).

5. **Vogl의 branch와 대표 gain.** 원문은 probe를 pump보다 약3 GHz 높은 blue side에 둔다. 표준 red-seed GABES reference의 직접 검증 대상이 아니다. Gain은4–10 범위이며 CSV7은 임의의 중간값이다. 본문의 “800 MHz” analysis frequency는 Fig.2 caption의800 kHz와 모순되므로800 kHz를 채택한다. 원문은 optimum OPD 폭 약200 MHz·TPD 수MHz라고 서술하지만 허용 dB threshold를 명시하지 않아 정밀 tolerance로 사용하지 않는다. [Vogl 원문, pp.1–3](https://arxiv.org/pdf/1209.2464).

6. **McCormick 2007 cell length의 판본 차이.** CSV에는19 mm이나 접근 가능한 arXiv v1에는 Brewster-tilted 유효 경로 약8 mm가 적혀 있다. 최종 출판 PDF의 해당 문구를 직접 확보하지 못했으므로 어느 쪽을 최종 사실로 강제 수정하지 않는다. 절대 gain의 길이 기반 재현에서는 이 입력을 미확인으로 분류한다. 총 검출효율0.8×0.82와 직접−3.5 dB/손실보정−8.1 dB는 원문에서 확인된다. [arXiv 원문](https://arxiv.org/pdf/physics/0607254), [최종 출판 정보](https://doi.org/10.1364/OL.32.000178).

7. **Liu pump 범위와 efficiency 불완전성.** 대표 IDS 운전점은400 mW이지만 Fig.4는 pump100–700 mW에서 SQL 교차 대역폭5.5–16.5 MHz를 측정했다. CSV pump min/max400/400을 논문 전체 탐색 범위로 해석하지 않는다. QE96%만으로 총효율을 확정할 수 없고, pump-off probe transmission90%는 내부 흡수 진단이라 후단 대칭 loss에 그대로 합치지 않는다. 논문 TPD는 pump–seed spacing minus HF로 정의되어 +4 MHz가 GABES −4 MHz에 해당한다. [Liu 저자 공개 원문, pp.2–3](https://zhifanzhou.com/files/Realization%20of%20low%20frequency%20and%20controllable%20bandwidth%20squeezing.pdf).

8. **Allen slope·SQL와 실제 tolerance.** 최고 직접−8.2 dB, slope−8.8 dB, 안정 run 평균−7.905 dB는 다른 결과다. 특히 저자는 보수적인 SQL 기준을 위해 낮은 conjugate power에 맞춰 비교하므로 일반적인 총출력 SQL과 측정면을 확인해야 한다. Slope ratio0.131을 손실보정 intrinsic squeezing으로 부르지 않는다. 원문 조건부 tolerance는 cell±1°C, etalon±2.5 m°C, seed plateau12–20 µW, pump 약450 mW 이상이며 다른 조건을 고정한 결과다. 재조정 종료기준 약11 MHz는 제어 알고리즘 설정으로서 IDS tolerance가 아니다. [Allen 최종 원문](https://doi.org/10.1364/OE.600911), 로컬 `references/2026 OE Long-term stabilization ids.pdf`, pp.9–12.

9. **Jain의 10% loss 산술.** 원문은 fiber 전−7.2 dB, 후−4.4 dB와 coupling>90%를 보고하고, 10% 대칭 optical loss만 예상하면−6.3 dB라고 적는다. 그러나 명시한 숫자에 `V_after = 0.9 V_before + 0.1`을 적용하면 **−5.66244 dB**이고 추가 차이는 **1.26244 dB**다. 보고서는 반올림된 입력의 자체 재계산값을 사용하며, 원문의−6.3 dB/1.9 dB를 산술 검증 없이 재사용하지 않는다. 저자의 미상관 공간모드 설명은 단순 vacuum-loss 모형을 넘는 제한을 뒷받침한다. [Jain 원문, Results](https://arxiv.org/pdf/2507.03755), [최종 DOI](https://doi.org/10.1364/OL.564233).

10. **Monsa preprint slope·주파수·효율의 미해결 항목.** 원문 slope0.135는 정확히−8.6967 dB인데 약−8 dB로 기재하고 직접−8 dB와 비교해 system noise 약1 dB라고 서술한다. Abstract의0.8 MHz, 본문의40–500 kHz와 Fig.6의약1 MHz도 각 측정 문맥을 구분해야 한다. 보고된300 mW·gain약15·직접IDS약−8 dB는 참고하되, 미상 총효율·beam waist를 추정값으로 채워 calibration에 사용하지 않는다. [Monsa v1 원문](https://arxiv.org/pdf/2601.13939v1).

11. **de Araujo의 관측량과 출처 경로.** 이 논문은 vacuum-seeded homodyne quadrature squeezing으로 직접 IDS식의 같은 행에 회귀할 대상이 아니다. CSV의 저주파5 Hz보다 아래인 sub-1 Hz까지 원문에서 보고하며, pump800 mW를 두 영역에 나눠 약400 mW씩 사용한다. CSV의 다른8개 논문 `source_file` 경로는 현재 저장소 `references/` 안에 존재하지 않아 공개 원문으로 확인했다. 이는 논문 미확인이 아니라 로컬 경로 provenance의 공백이다. [de Araujo 원문](https://pdfs.semanticscholar.org/3c89/e3090c114bfd10282b2189aa0326f8347475.pdf), [저자 기관 출판 정보](https://jqi.umd.edu/publications/characterizing-two-mode-squeezed-light-four-wave-mixing-rubidium-vapor-quantum-sensing), [DOI](https://doi.org/10.1364/OE.507727).

McCormick 2008은 별도로 measured effective gain과 intrinsic gain을 구분하고 distributed gain/loss 모형을 제시한다. 보고서의 간단한 `1 - eta + eta/(2G - 1)` 비교에서는 실측 이득을 쓴 조건부 계산이라는 한계를 밝힌다. [원문](https://arxiv.org/pdf/quant-ph/0703111), [최종 DOI](https://doi.org/10.1103/PhysRevA.78.043816).

## 산출물의 역할과 우선순위

- `squeezing_report_v7.tex` / `squeezing_report_v7.pdf`: 최신 독립형 과학 보고서.
- `squeezing_report_history.md`: 현재 파일. 편집 이력·전체 커밋·과거 판정·출처 정정 기록.
- `build_v7_audit.py`와 **`v7_assets/`는 첫 감사계산**을 보존한다. `commits.json`, `point_audit.csv`, `audit_manifest.json`, `validation.json` 및 그림은 당시 조건의 산출물이다.
- **`v7_science/`는 최신 scientific calculations**와 그림·표·provenance의 위치다. 최신 보고서의 정량 논의는 해당 계산 조건과 데이터 정의를 따른다.

두 디렉터리의 결과가 다른 경우 먼저 model fidelity, 수치격자, 고정 변수, 후처리 및 noise 가정이 동일한지 확인한다. 동일하지 않은 실행의 숫자를 단순 버전 차이 또는 실험 검증 여부로 해석하지 않는다.
