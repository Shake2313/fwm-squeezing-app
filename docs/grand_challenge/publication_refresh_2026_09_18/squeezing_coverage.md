# 스퀴징 독립 보고서의 내용 인벤토리와 편집 근거

작성 범위: 새 `docs/squeezing_report/squeezing_report_v8.tex`·PDF 및 전용 자산. 원본 v1–v7, `squeezing_paper`, 기존 수치·그림은 수정하지 않는다. 보고서는 실험 정합성·변수 탐색·frontier의 과학적 해석을 담당한다. 개발 경위와 버전 비교는 이 문서에만 둔다.

## 원본별 실질 내용과 반영

| 검토 원본 | 과학적 내용 | 독립 보고서에서의 처리 |
|---|---|---|
| v1 | 4준위 감수율·Maxwell 전파, pump cap, 검출 floor, 온도·detuning 지도, 효율 레버, frontier·고이득 점근식, 차분 초과잡음 | 유효한 floor·목적함수·cap 한계 보존. `5/(8G²)`는 IDS 공식으로 틀렸으므로 올바른 `1/(2G−1)`와 반례를 함께 설명. 원시 격자가 소실된 지도 및 잘못된 식의 optimum은 현재 결과로 채택하지 않음. |
| v2 | 셀 내부 흡수와 검출손실 분리, 유한 온도 optimum, 낮은 온도·OD frontier | 분포 Gaussian 공분산 모형과 정의된 Pareto/연결 허용영역으로 보존. 평균장 감쇠에 별도 OD를 중복 적용하는 처리는 정당화하지 않음. 잘못된 잡음 환산을 쓴 좌표는 현재 운전 권고에서 제외. |
| v3 | 현실 검출·η=1 탐색의 구분, 광역 19,389점 scan, −81.95 dB, 초미세 기준선·문헌 비교 | 이상검출과 실제 측정량의 구분 보존. −81.95 dB는 공식 불일치의 반례로만 보존. Pooser/Jasperse 참조는 과학적 문헌 맥락으로 복원. |
| v4 | 올바른 이상 IDS 점근, gain-gap 검사, formula cross-check, 충돌/흡수와 고온 경계 | 같은 gain만으로 일반 공분산을 정할 수 없음을 명시. `0.5≤Gs−Gc≤1.5`는 휴리스틱 선별이며 CP/잡음 인증 아님. −31.35 dB 및 frontier 좌표는 일반 squeezing 예측으로 승격하지 않음. |
| v5 | 양 팔 OD·scatter 목적함수, red/blue 근접성, TPD/angle 상세 및 low-angle ridge, 광펌핑 해석 | 손실의 기준면·중복 계상·geometry 한계를 복원. 높은 gain/낮은 OD 자체가 모순이라는 주장, 4준위 정상상태에는 광펌핑이 없다는 주장은 폐기하고 원자 상태·CG·잡음의 별도 검증으로 대체. |
| v6 | finite/ideal frontier, red/blue 좌표, 저각도 ridge, 다섯 변수 tolerance | 대표 저장 좌표·조건부 dB와 ± 공차를 데이터 진단 표에 보존. −8.10/−15.62/−20.70 dB는 미시적 잡음이 없는 surrogate 결과라고 명시. 보편적 운전 권고·seed 무관성·온도 오프셋 처방은 채택하지 않음. |
| squeezing_paper | 모델·손실식·floor·counterexample·geometry·tolerance의 집약, Pooser/Jasperse 문헌, reproducibility 제한 | 과학적 주제를 모두 위 항목에 통합. chronology·다섯 단계 구조·코드 교정 기록은 본문에서 제거. 소실된 raw grids와 원시 수치의 제한은 provenance 부록에 기록. |
| v7 전체 | SQL/주파수 규약, ideal/불균형/분포 covariance, gold power fit, 3종 optimum, 0.5 dB 예산, 7변수 GABES scan, 206/145 수렴 결과, 10개 실험·Allen tolerance·carrier/fiber 사례, 한계 | 모든 과학적 절·식·수치표·그림 5개를 보존. 기존 한계·개발 순서 절은 현재 물리 검증/실험 요구조건으로 갱신. 개발 문체를 없애고 각 결과의 데이터·보정 계약을 정의. |
| v7 science theory_notes | finite-seed exact moments, asymmetric-loss/weighted readout, joint uncertainty/Hessian, 지연 예산 | 본문 또는 부록에 독립적으로 읽히는 식을 복원. 파일을 선행학습으로 요구하지 않음. |

## 그림과 수치의 보존 원칙

- v7의 `theory_gold_comparison.pdf`, `theory_gain_absorption.pdf`, `theory_tolerance.pdf`, `coupled_maps.pdf`, `tolerance.pdf` 전부 재사용한다. 현재 소스 재실행 결과로 재명명하거나 과거 소스 해시를 변경하지 않는다.
- v1–v6/paper의 반복 지도는 해당 모델의 잡음 식·이중 흡수·heuristic gate에 묶여 있다. 그 그림을 현재 예측으로 재출판하지 않는다. 주요 수치·저각도/공차 주제는 범위가 표시된 진단 표로 보존하고, 그림 원본은 그대로 보존한다. 같은 내용을 가진 지도 수를 유지하는 것은 과학적 내용 보존의 기준이 아니다.
- 새 열적 앙상블 그림은 저장된 p2 세 seed의 실제 행렬 비교를 재사용한다. 샘플 수 증가 또는 경로 정확도가 곧 광학 squeezing 수렴이라는 오해를 막기 위해 경로 오차·ensemble 차이를 별도 표기한다.
- `v7_science` 평균장 fit(전체 residual 0.42504822, target gain 15), 별도 pole/Ultra mixing fit(0.5594938027, target 15.5), microscopic reduced-RMS fixture를 별도 계산으로 명시한다. 계수를 섞거나 하나의 no-fit 결과로 합치지 않는다.

## 갱신할 현재 증거

1. 독립 QRT/원자 source 검증과 reduced single-velocity readout의 성과·적용 범위.
2. 새로운 local null/loss probe: SQL 한계 1.11e−16, 대칭 손실 법칙 3.33e−16; 38 assumed scalar로 no-fit 감사 실패. 정량 예측 정확도가 아닌 내부 검증.
3. Thermal p2 72경로·360계산, 3/9 grids, seed 비교 여섯 방향 0/6 통과, 5% 기준 대비 metric 변화 5.70–43.61%. 최종 optical gain/SQL 미구성.
4. 고정된 gain 보정으로 Sim 15.5, Liu 6.445 대 8, McCormick 1.926 대 9. 조건 변경 시 외삽 실패와 실측 입력 부재를 결과로 제시.
5. source covariance·미시적 모델·검출 정의·독립 입력·blind holdout을 별도 판정. Sobol seed 독립성은 실험 조건 독립성과 다름.

## 명시적으로 바로잡을 내용

- 대칭 vacuum loss는 `R→ηR+(1−η)`로 SQL에 접근한다. 비대칭 loss/가중치 불일치는 큰 단일 팔 잡음을 차분에 남겨 SQL 위로 올릴 수 있으므로 “어떤 손실도 anti-squeezing 불가”라는 일반화는 사용하지 않는다.
- detuning convention은 Hz와 rad/s를 구분하고 − branch beat는 `2π(νHF−δν)`다. 동일 OPD 숫자를 다른 ground-hyperfine 기준에 그대로 비교하지 않는다.
- pump cap은 에너지 수지의 사후 제한이며 self-consistent three-field depletion이 아니다. explicit transport와 phenomenological transit reset을 이중 적용하지 않는다.
- floor/positivity/수치 수렴은 모델 충분성·실험 정합성의 증명이 아니다. No-fit은 독립 입력 조건이며 target을 fit하지 않았다는 진술만으로 성립하지 않는다.
- 광역/저각도 surrogate frontier를 실험 최적점으로 승인하지 않는다. 차분 PSD·양 출력·실제 손실·잡음의 공동 검증이 필요하다.

## 검증 계획

새 PDF 생성 전에 이 인벤토리를 작성했다. XeLaTeX 반복 빌드, undefined reference/overfull/누락 그림 검사, PDF 전 페이지 렌더·시각 확인, 본문 수치와 기존 JSON 교차검사, 입력 파일 해시와 산출물 기록을 수행한다. 새 heavy simulation은 실행하지 않는다. 전체 pytest는 주 에이전트가 담당한다.

## 개별 축약·제외의 과학적 이유

| 원본 내용 | 처리와 이유 |
|---|---|
| v1 `squeezing_map`, `squeezing_theory_max`, paper fig2 | 잘못된 `5/(8G²)` 잡음 환산과 소실된 raw grid. 해당 극값·온도별 정밀표를 유효한 계산으로 복원하지 않음. floor·plateau·비용전선의 정의는 본문에 보존. |
| v2/v3 `squeezing_frontier`, `squeezing_frontier_od`, paper fig3 | 같은 OD/frontier 주제의 반복 이미지이며 잘못된 잡음 환산에 결합. 올바른 분포 공분산의 2D 그림과 유한 optimum 표로 물리 내용을 보존. |
| v3 `squeezing_frontier_ideal_qe1`, paper의 −81.95 dB 사례 | 19,389점 탐색 조건, 실제 −81.9469 dB 좌표·두 gain을 반례로 독립 본문에 보존. 결과를 source limit라고 주장한 부분은 제외. |
| v4의 −31.35 dB 표·formula 교정 전후 지도, paper fig4 | 일반 손실 매질에서 gap 식도 covariance를 대신하지 못함. old/gap/canonical 식의 적용 조건과 실패 메커니즘은 보존. v4와 paper의 동일-pair 비교 dB는 서로 불일치하여 독립적으로 확정되지 않은 숫자를 추가 정답으로 선택하지 않음. |
| v5 강화모형의 −2.00 GHz optimum·0.28 dB lobe 간격 및 fig5 | beat 규약 변경 이전의 중복 탐색이며 현미시적 noise를 인증하지 않음. 양 팔 손실·scatter 비용과 lobe 근접성의 과학적 질문을 보존하고, 같은 beat 규약의 대표 저장 자료를 사용. |
| v6 finite/ideal/OD 및 paper fig6–7 | finite와 η=1 구분, 두 lobe·gold·상세기하·gain-gap 극단점의 좌표와 dB를 한 표/설명에 보존. 장식적인 반복 heatmap은 싣지 않되 원본은 보존. |
| v6 TPD/angle·low-angle 및 paper fig8 | 기본 angle·저각도 ridge·탈락 raw 최저점의 실제 좌표·gap·값을 보존. −20.70 dB는 검증된 source limit가 아닌 조건부 surrogate 진단. |
| v6 tolerance 및 paper fig9 | +0.25/+0.5/+1 dB의 OPD·TPD·온도·loss 전 열, TPD 재조정 폭, seed 1–64 µW 변화, 5 MHz 대 17.5 MHz grid, loss slope·QE 조건을 보존. 운전점의 보편적 권고는 제외. |
| v1–paper의 “손실로 양의 IDS 불가” | 동효율 vacuum-loss 명제로 제한하고, 비대칭 손실의 일반 covariance를 추가. |
| paper의 “gain 크고 OD=0이면 모순”, “고정 population라 광펌핑 없음”, “self-consistent depletion” | 현재 모델·양자채널의 사실과 맞지 않는 일반화. explicit Hamiltonian/state·source·전파·에너지 수지 검사로 대체. |
| 실행 시간·커밋·버전별 변경 요약, 중복 build commands | 과학적 observable이 아니므로 보고서 본문에서는 제거. 원본·본 인벤토리·새 build/QA 로그에 기록. |

최신 보고서의 부록은 이전 판을 열지 않아도 surrogate 목적함수의 의미, 대표 좌표, 검출 조건, 탈락 조건과 공차를 이해할 수 있도록 작성했다. 과거 이미지는 삭제·덮어쓰기하지 않았다.


## 편집 실행 기록

- 저장된 증거만 사용했으며 새로운 heavy simulation을 실행하지 않았다. 보고서 작성 때문에 원자 엔진의 식을 변경하지 않았다. 이 두 문장은 과학적 본문에서 제거하고 여기에만 기록한다.
- Thermal manifest의 square full width는 346.410 µm이고 half-width는 173.205 µm다. T=373 K, n=1e18 m⁻³, L=12.5 mm, pump 600 mW/530 µm, 무편극 유입 및 측벽 pump intensity 65.2–80.8%를 본문에 명시했다.
- build_report.txt는 초기 조립 스냅샷이며, 추가 공차·경계조건 교정을 포함한 최종 source는 squeezing_report_v8.tex다. 최종 PDF는 source를 직접 XeLaTeX으로 빌드한다.

- 새 build/QA 기록 79개를 `docs/grand_challenge/publication_refresh_2026_09_18/squeezing_build/`로 이동했다. 되돌릴 수 있는 정확한 원·목적지와 해시는 `squeezing_move_manifest.json`에 기록했다. 과학 보고서 디렉터리에는 최신 TeX/PDF와 실제 과학 자료만 남긴다.

- 고정 Cmix 비교는 alternative_results.json의 pole/Ultra·1 m/s·4σ·64구간 조건을 명시했다. 임시 초안의 별도 화면-scan 131→141°C 문장은 같은 artifact의 116/121/126°C 및 5.305/15.500/61.018 수치로 교체하여 서로 다른 실행 계약을 합치지 않았다. 이 수정은 이전 v1–v7 과학 내용을 삭제하지 않는다.

## 최종 검증 결과

- 23쪽 PDF를 XeLaTeX 두 pass로 생성했다. 최종 log에서 LaTeX/font warning, overfull/underfull, undefined reference, missing glyph가 모두 0이다.
- Poppler 95 dpi로 23쪽 전부 렌더하고 시각 검토했다. 누락·잘림·겹침 없이 수식·표·그림을 확인했다.
- 보존 입력 31개 SHA256이 모두 일치했다. v7 equation 23/23 및 align 3/3 블록, 모든 label과 그림 5개가 그대로 보존됐고, 열적 수렴 그림을 추가했다.
- 검증된 build PDF와 게시 PDF는 byte-identical이다. build 명령·최종 해시·검사 상세는 `squeezing_build/validation.json`에 기록했다.
