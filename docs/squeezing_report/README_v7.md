# 스퀴징 보고서 v7

최신 TeX/PDF만으로 이론, 최적점, 다변수 예측, 문헌 정합성, 허용범위와 코드 한계를 읽을 수 있는 독립 보고서다. 변경 이력과 과거 분석의 평가는 squeezing_report_history.md에 분리했다.

## 읽기

- squeezing_report_v7.pdf / squeezing_report_v7.tex: 과학 보고서.
- squeezing_report_history.md: 2026-07-22 이후 26개 커밋, 코드 변화, 과거 분석 기록, 레퍼런스 CSV의 출처 정정.
- v7_science/theory_notes.md: 분포 Gaussian 모형의 상세 수식과 가정.
- 원본 v6 및 references/fwm_squeezing_paper_parameters.csv는 수정하지 않았다.

## 새 계산과 해석

1. Gold 평균 출력 111/109 µW와 seed 8 µW만으로 분포 이득·probe 흡수 모형을 제약했다. 이 출력을 셀 출력으로 간주하고 후단 효율을 별도로 적용한다. 잡음은 적합하지 않았다. 모형의 G=16.02594, T=.821165에서 IDS는 η=.865일 때 −7.70926 dB, η=.8694일 때 −7.81890 dB다. 출력·효율의 기준면과 반올림 불확도가 작은 잔차보다 클 수 있다.
2. GABES의 residual ℓ을 gold probe gain=15 한 점에 맞춰 .42504822로 보정했다. 이후 모든 변수 스캔은 이 값을 고정한다. 실제 원자 잡음 대신 canonical bright-seed 이상식을 명시적으로 적용한 조건부 예측이다.
3. 주 206개 scan의 NF3/NF2 비교가 모두 CONVERGED였다. dv=2.5 m/s, cutoff=5σ 정밀화의 145회 평가도 모두 통과했다. 이는 일차원 모델의 수치 판정이며 실험 검증과 구별된다.
4. 고정 TPD의 +0.5 dB 경계는 refinement.json을 우선한다. TPD 연결구간은 −18.40872..+61.54507 MHz다. scan 끝에 걸린 반대쪽 경계는 물리적 한계가 미확정이다.
5. 실측 tolerance는 Allen 2026에서 별도 인용했다. GABES 조건부 구간, 분포 이론의 계수 구간, 실제 장치의 구간은 서로 다른 기준이다.

## 재현

저장소 루트의 PowerShell에서 순서대로 실행한다.

~~~powershell
python docs/squeezing_report/analyze_squeezing_v7.py
python docs/squeezing_report/refine_squeezing_v7.py
python docs/squeezing_report/distributed_squeezing_v7.py
python -m pytest -q
Set-Location docs/squeezing_report
xelatex -disable-installer -interaction=nonstopmode -halt-on-error squeezing_report_v7.tex
xelatex -disable-installer -interaction=nonstopmode -halt-on-error squeezing_report_v7.tex
~~~

의존성은 저장소의 numpy/matplotlib/scipy와 XeLaTeX, kotex, Malgun Gothic이다.
실행 당시 HEAD와 gabes 소스 hash를 results.json에 기록한다. 엔진이 바뀌면 현재 입력으로 새 계산이 되므로 역사적 재현에는 동일 소스가 필요하다.
재계산은 입력과 소스에 따른 cache를 사용한다. cache와 정밀화 자료의 기준이 다르면 후속 생성기가 알려주므로 생성 순서를 지킨다.

## 파일별 역할

| 파일 | 내용 |
|---|---|
| analyze_squeezing_v7.py | 일점 gain 보정, 7변수 scan, OPD–TPD 및 OPD–온도 지도 |
| refine_squeezing_v7.py | 고정 TPD 교차경계 이분법, TPD 연결구간, cap/floor 분리 |
| distributed_squeezing_v7.py | 분포 quantum covariance, gold 출력 적합, 유한 최적점과 잡음예산 |
| v7_science/results.json | 기본 조건, 소스 hash, 주 scan, 보간 tolerance |
| v7_science/refinement.json | 정밀 경계·탐색구간·수렴판정·원본/동적/최종 gain 구분 |
| v7_science/theory_results.json | 독립 Gaussian 이론 결과와 극한해·불확정성 검사 |
| v7_science/validation.json | 최종 산출물 hash, pytest, PDF 검수 기록 |
| v7_assets/ | 처음 수행한 감사 계산·커밋 수집 자료; 최신 본문 그림은 v7_science 사용 |

results.json의 legacy 진단 키 capped_fraction은 전체 scan에서의 |Gs−Gs_smallsignal|/|Gs_smallsignal| 최대값이다. 이름과 달리 최종 pump cap만의 감소율이 아니며, 동적 전파와 수동 gain의 floor 처리도 섞인다. 개별 영향은 refinement.json의 selected_points 항목에서 확인한다.

검증: 전체 pytest 652 passed (147.49 s). 분포 모형은 무손실 이상식, coherent 수동손실, 전체 공분산 불확정성 및 독립 구간 분할 수렴 검사를 통과했다.
