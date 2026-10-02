# Analytic reconstruction v3 보존·검증 기록

산출물: `docs/FWM physics and analytic reconstruction/squeezing_analytic_reconstruction_v3.tex`, 동명 PDF. 기존 v1/v2 TeX·PDF 수정 없음. 과학 본문 영어 유지. 개발 경과·버전 차이는 이 기록에만 보관.

## 읽기와 보존 범위

- `CLAUDE.md`, `README.md`, publication policy 확인. v2 3,107줄 전체 읽기 완료. v1과 v2 차이 대조, v1 고유 수치·반례 별도 확인.
- v2 label **133/133 보존**. v3 총 **188개**, 추가 55개. 중복 label·누락 label 없음. v3 단일 TeX에 모든 본문 포함. 이전 PDF/TeX를 읽어야 이해되는 유도 없음.
- v1 label 62개 중 57개 직접 유지. 나머지 5개는 아래 의미 매핑. label 보존만으로 의미 보존을 증명한다고 주장하지 않음. 아래 주제별 대조와 수치 복원이 별도 근거.

| 기존 내용 | v3 처리 |
|---|---|
| 원자 준위·carrier/beat/RF 주파수·Fourier 부호·SI 차원·Rabi/전기장/광자유속 변환 | Definitions 전부 유지. 새 reciprocal normalization과 구별 명시 |
| pump-only nullspace·gauge transformation·thermal reload·rate/Sylvester 근사·강한 pump 반례 | 유도·조건·수치 유지. 개발 경과 대신 모델 가정으로 표현 |
| trace-zero 응답·2D Raman Doppler·Floquet recurrence와 cutoff gate | 식·표·복소 응답/위상/최적점 gate 유지 |
| 고유벡터 pole sum·Jordan block·zero-mode residue·Lorentzian 부호·Voigt/Faddeeva prefactor | 전부 유지. GKSL 반평면 명제의 조건 강화 |
| canonical Maxwell 변환·geometric mismatch·Cayley–Hamilton·2×2 지수·gain/Manley–Rowe gap | 전부 유지. double-dispersion 예시는 명시적 반례로 유지 |
| Langevin commutator·covariance·Lyapunov entire function·고유기저/Van Loan 평가·분할 전파 | 전부 유지. Lyapunov 일반 고유합 개수 오류 수정 |
| 임의 independent-vacuum completion과 Caves 한계 구분 | 전부 유지. 미시적 atomic diffusion과 혼동 금지 |
| TMSV marginals·binomial thinning·factorial moments·unequal loss·weighted difference·displaced thermal·thermal attenuator | 모든 식·조건·예시 유지 |
| detector loss·photocurrent·SQL·gain-only diagnostic·passive counterexample | 전부 유지. finite RF mode/one-sided PSD/finite-seed 검출 조건 확장 |
| 해석 항등식 수치 검증·pole spectrum·Floquet N=1/2/3·finite-seed limit·state positivity 반례·속도 step/cutoff 표 | 원 수치 전부 유지. 물리 검증과 수치 일치 구분 |
| phenom drift·electric transfer·선택 진공 completion·raw experimental gain discrepancy | 원 식·수치 유지. microscopic 모델에 phenom gain을 사전 배정하던 표 수정 |
| Torrey cubic·차원 표·가정 표·parameter provenance·ell/transit 민감도·normalization ledger·누락 물리 | 전부 유지. 별도 microscopic/transport 모델의 존재 반영 |
| 기존 12개 문헌·재현 경로 | 유지. Gaussian channel 및 CF4 primary references 2개, 현재 미시적 수치 데이터 경로 추가 |

## v1 고유 내용 복원

Appendix `app:gc-countermodels`에 부정확한 모델을 **명시적 반례**로 복원. 현재 예측으로 승격하지 않음.

- static drift `(-1.330328+242.349727i, …)`와 transfer `(-0.376876+2.339296i, …)` 전체 행렬.
- gain 5.61434/4.66816, gap 0.946182, K 고유값 0.189391/3.153969.
- reservoir 생략 commutator 잔차 0.0551; vacuum completion local/global 잔차 4.44e-16/1.04e-15.
- source/detected ratio 0.09230/0.21084; DC weight 1.20269에서 0.09774/0.21558; ideal matched model 0.09776/0.21560. 해당 dB 전부 유지.
- 0.1/0.15/2/4 MHz flat-spectrum 예시. 파장·면적 같다는 근사 명시.
- ell 로그민감도 3.53/4.29; 90→110 kHz damping 변화의 gain 5.6138→5.6033, 4.6631→4.6573.
- 근사 비교 gain 15.5/13.6의 -63.8%/-65.7% 구분. 원시 power-ratio 비교와 혼동 금지.
- gain-only covariance ansatz 식 복원. -8.10/-15.62 dB는 해당 ansatz의 점수일 뿐임을 명시.

| 직접 재사용하지 않은 v1 label | 대체·보존 위치 |
|---|---|
| `tab:frequency-audit` | `tab:frequency-definitions`, Fourier·carrier 유도, 새 `sec:gc-normalization`·`sec:gc-readout` |
| `tab:spectrum-comparison` | `tab:model-layer-comparison`, `tab:gc-countermodel-noise` 및 4개 RF의 flat-spectrum 반례 |
| `tab:acceptance` | `tab:validation-hierarchy`, `tab:gc-thermal-convergence`, 명시적 경계·앙상블·장치 claim gate |
| `app:changelog`, `tab:changelog` | 개발 이력 자체는 과학 본문에서 제외. 해당 내용의 과학적 모델 차이·식·반례는 위 대조표 및 countermodel appendix에 보존 |

## 명시적 이론 수정

1. 2×2 Lyapunov 고유합은 일반적으로 **4개**. 두 대각 합과 서로 켤레인 두 교차 합. 특수 조건에서만 중복. 기존 '3개' 수정.
2. `TJT†=J` 군은 일반적으로 **U(1,1)**. 전체 위상 제거 후 SU(1,1). 이름 수정, 항등식 유지.
3. `ell=.74/12`은 한 개의 외부 population 중복을 없애는 phenom 정의. pumped F=2/F=3 manifold 각각의 정확한 약흡수 normalization 인증 아님. 필요한 `3 C_F² p_F/(2F+1)=p_F S_FF'/3` 유도 및 5/12·7/12 불일치 추가.
4. Lindblad 반평면 정리는 GKSL임이 입증된 generator에 조건부. spectral stability나 한 stationary state의 positivity만으로 generator CP 입증 불가. 별도 Choi 반례 추가.
5. perturbation coupling operator를 0으로 보내는 극한과 약한 probe 진폭→0을 구별. 후자는 finite susceptibility 가능.
6. 닫힌 canonical field의 B=D=0 한계와 atomic jumps=0을 구별. 후자는 initial atomic noise·비감쇠 모드 제거를 뜻하지 않음.
7. microscopic M/D가 전혀 없다는 낡은 상태 문장 제거. phenom drift의 미시적 완성과 이미 존재하는 reciprocal reduced-model M/D를 구별. microscopic 행의 403.535/405.854 사전 배정 제거.
8. 기존 sign/population/transit 변형 결과를 개발 순서 대신 가정이 다른 countermodel로 제시. 기존 scalar·gain-only 값은 물리 squeezing 예측으로 재해석하지 않음.

## Grand Challenge 추가 이론

- GKSL에서 real Hermitian-basis drift·Einstein product rule·source Gram positivity·비정상 covariance 미분·stationary commutator 항등식 유도.
- greater/lesser ordered spectra와 독립 full-density QRT. 비혼합·peripheral/elastic mode 한계. random GKSL 30개×RF 5개 독립 probe 수치.
- reciprocal SI drive/readout dipole·linewidth·manifold population normalization·carrier 기준.
- line density·원자 slice-noise normalization·field M/D 유도·면적 상쇄·속도군 covariance 가중·bath 중복 금지.
- 네 physical sideband의 여덟 real quadrature·Bogoliubov 변환·Gaussian CP/uncertainty·RF 유한모드와 버려진 vacuum-mode noise.
- finite coherent seed·periodic diffusion Fourier blocks·독립 one-period QRT tail·mean/tangent 차이·photocounting Gaussian closure·finite optical-band tail 한계.
- 두 moving phase·kinetic lattice·Q=0 coherent quotient·moving-grating 반례·유한 aperture 및 공간-mode projection.
- ballistic 두 나이 covariance·initial internal noise와 jump noise·적분 출력의 augmented solve·독립 adjoint 및 full-density QRT.
- Poisson arrival-number noise·공통 entry phase 평균·Maxwell flux·occupancy identity·경계 상태와 upstream pumping 민감도.
- nonlocal response/noise kernel과 optical boundary problem. integrated path response를 local drift로 오인하지 않음.
- exact constant-step divided exponential·complex density eigenoperator 취급·condition/residual guard·block fallback·CF4·primary/reference path gate.
- p2 세 scramble의 72 trajectories/360 solves 및 전체 6개 방향의 gate 실패. source-resolved 최대 43.6110%, occupancy ratios, path convergence와 ensemble convergence 구분.
- 동일 detector current의 one-sided SQL·Hz/angular PSD 변환·RF 평균 후 ratio/dB·parameter covariance와 model-conditioned uncertainty.

## 검증

- XeLaTeX 반복 빌드. 최종 PDF **53쪽 A4**, 595.28×841.89 pt.
- unresolved reference/citation 없음. overfull/underfull box 없음. missing character 없음. TeX 비정상 control character 없음.
- 모든 페이지 PyMuPDF 1.27.1로 1.3배 raster. 53 PNG·9 contact sheets. 전체 contact sheets 육안 확인. 수식 밀집 12/22/25/26쪽, 표 42/44/46쪽, 데이터 경로 51쪽 추가 확대 확인.
- 22쪽 RF-filter 식의 literal `frac12` 발견 후 수정, 재빌드·전체 재렌더·해당 식 확대 재확인 완료. 최종 PDF에는 정상 1/2.
- PDF text bounding boxes 검사: 페이지 밖 text **0개**. U+FFFD replacement character **0개**.
- 기존 longtable 44–47쪽에서 엔진의 `Infinite glue shrinkage ... being split` 진단 2개 남음. CJK italic/smallcaps font substitution 경고 2개 남음. 최종 해당 표·본문 렌더에서 겹침·잘림·누락 없음. 경고를 무검증으로 무시하지 않음.
- `python -m pytest -q tests/test_fwm_docs_consistency.py`: **7 passed**, 최종 source 기준 1.41s. 전체 pytest는 main agent가 병렬 실행·기록 담당. 새 ODE 대규모 계산 없음.
- `analytic_source_validation.json`: label inventory·TeX hash. `analytic_pdf_validation.json`: PDF hash·page/bounds 결과. `analytic_compile/`: aux/log. `analytic_visual_qa/`: raster·추출 text. `analytic_build.txt`, `analytic_gc_chapter.tex`, `analytic_v1_fixture.tex`, `analytic_render.txt`: 이 문서 전용 빌드·QA 자료. 최종 TeX는 이 자료를 include하지 않는 단일 소스.

최종 TeX SHA-256: `20b052e8fdbdb10fd3015a4f9725ff9e6171861a36d6b35ff2cd470ea9592527`.

최종 PDF SHA-256: `20d621c00a0f65b41582aee8aac412a6460a23d142f3c16b8e289834ef52c926`.

## 과학적 미해결 범위

thermal ensemble 수렴 없음. 주입 state의 물리적 준비 경로 미확정. full Zeeman/polarization strong-drive 모델 없음. 비국소 atomic kernel을 푸는 cell Maxwell boundary problem 미완성. 측정 collection/detector/SQL·reservoir 입력으로 닫힌 독립 장치 예측 없음. QRT·CP·commutator 일치는 해당 수학적 모델의 증거이며 실제 실험 squeezing 인증 아님.
