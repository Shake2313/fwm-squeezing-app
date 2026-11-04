# Rydberg-EIT 전기장 센싱 문헌 DB

사용자 요청일 2026-11-03. 실험 레퍼런스 논문 탐색 후 [FWM DB](../fwm_experiments/README.md)·[PSR DB](../psr_experiments/README.md)와 같은 구조로 구축.

**Ju et al. (2026)을 당분간 Gold reference로 채택.** 실험 구도·조건의 메인 논문으로 삼는 사용자 결정. [채택 기록](reference_policy.json)·[논문 원본 JSON](papers/ju_2026_preprint.json)·[원문](https://arxiv.org/abs/2606.04354).

**[오프라인 검색·비교 열람본](index.html)** · [논문 CSV](papers.csv) · [실험·계산 조건 CSV](experiments.csv) · [수치·출처 CSV](observations.csv) · [주요 결과 CSV](findings.csv) · [전체 JSON](database.json)

## ⟦AI⟧ 수록 범위와 열람

2007–2026년 26편. 실험 23편, 이론 2편, 이론·실험 1편. 39개 실험·계산·스캔 기록, 171개 수치·미지정 항목, 42개 주요 결과. 기록 수는 독립 실험 표본 수와 다름.

열증기 RF 전기장 센싱 21편을 `primary`, 초기 EIT·OD 최적화·감도 이론·냉각 원자 SQL 비교 5편을 `context_only`로 구분. 초기 밝은 공명·AT 보정, 벡터 편광, homodyne·FM·RF mixer, superheterodyne, repump, laser noise, flat-top·원자수, comb·대역폭·다중 채널, MEMS·비공명 저주파까지 포함.

26편 모두 원문 확보. `sources/`의 PDF마다 URL·판본·SHA-256·페이지 수 기록. 출판 서지와 실제 검수한 preprint·저자 원고의 판본은 `publication`·`sources`에서 각각 보존. 출판 24편, 독립 출판 DOI를 확인하지 않은 preprint 2편. 공개 저장소와 브라우저 기관 접근 병용. 접근 이력은 [SEARCH.md](SEARCH.md).

`index.html`은 데이터 내장형. 파일을 직접 열어 검색·주제 필터·원문 위치 확인 가능. Gold 논문은 기본 목록 맨 위. DOI·원문 웹 링크는 인터넷 필요. CSV는 Excel 호환 UTF-8 BOM. PDF·추출문은 기존 DB처럼 로컬 보관하고 Git에서는 제외.

## ⟦AI⟧ Gold 논문의 보고 조건

아래 값은 Ju 논문의 보고값 재서술. 원문 위치는 [논문 JSON](papers/ju_2026_preprint.json)의 항목별 `locator` 참조.

| 항목 | 원문 조건 |
|---|---|
| 광학 전이 | ⁸⁵Rb 5S₁/₂(F=3) → 5P₃/₂(F′=4) → 40D₅/₂ |
| MW 전이 | 40D₅/₂ → 39F₇/₂, 약 37 GHz |
| 셀·온도 | 길이 50 mm, 실제 실험 약 20 °C |
| 광학 배치 | Probe·coupling 역방향, 선형·직교 편광 |
| 광파워 | Probe 6 µW, coupling 30 mW |
| 빔 크기 | 직경 0.15 mm. 1/e²·FWHM 정의 미지정 |
| RF LO·IF | LO 37.07104 GHz, ΩLO/2π=3.7 MHz, IF 40 kHz |
| 검출 | Balanced photodetector, QE 약 95% |
| Coherence 제어 | 3축 자기장 보상, 두 레이저 ULE cavity 고정 |
| EIT 선폭 | 좁은 보상 스펙트럼 1.6 MHz는 Pp 미상. Fig. 2(b)의 6 µW에서는 약 3.0 MHz |
| 측정 감도 | 12.5(8) nV·cm⁻¹·Hz⁻¹ᐟ². RF 장 RMS 환산은 식 (1) 규약 |
| 계산 PSN | 11.2 nV·cm⁻¹·Hz⁻¹ᐟ² |
| 온도 계산 | 약 33 °C의 PSN 최적점. 실제 약 20 °C 실험과 별도 기록 |

원문의 ⁸⁵Rb 전이 선택과 동위원소 순도 구분. 순수 ⁸⁵Rb 셀, 절대 probe detuning +40 MHz, 절대 동작점 OD는 이 논문 기록에 추가 가정으로 넣지 않음. 직경 0.15 mm에 AutoOD의 1/e² 정의를 적용하려면 별도 가정 필요. 대화에서 계산한 AutoOD 값은 논문의 실측 OD와 별도.

원문 내부 표기 차이도 보존. 본문 2절은 F=3 → F′=4, 그림 1(a)는 F=2 → F′=3. 메인 전이에는 본문 표기를 채택하고 그림 표기를 별도 항목으로 기록. 현재 보관한 Gold 원문은 main article 15쪽. 별도 Supplementary Materials는 미확보.

## ⟦AI⟧ 수치 해석

- `field_sensitivity`는 전기장/√Hz. 최소 검출장, 적분 시간, RBW·ENBW, IF는 별도 quantity. 긴 적분의 pV/cm 값을 nV/cm/√Hz로 직접 비교하지 않음.
- 측정 감도·잡음 스펙트럼 기반 환산·광자 PSN·원자 SQL은 정의 보존. Cai 2022의 직접 응답 35.75와 잡음 평가 12.50 nV/cm/√Hz는 다른 기록.
- 단위의 값·배율 보존. `uV/m/sqrt(Hz)`와 `nV/cm/sqrt(Hz)`는 같은 열에서 숫자만 비교하지 않음. CSV 별칭 정규화는 배율 환산을 수행하지 않음.
- Cell input·detector·원자 위치장·외부 입사장 구분. 원문이 정의하지 않은 기준면은 `not_specified`. RF detuning과 optical probe detuning도 구분.
- 빔 반경·직경·FWHM·spot size는 원문 정의 유지. 미지정은 `null`/`not_reported`. 다른 스캔·최적값을 하나의 가상 운전점으로 합치지 않음. Gold 그림 점 추출·모델 비교는 [별도 분석 보고서](../../analysis/rydberg_gold_validation_20261103/README.md)에 보관.
- `apparatus_group`은 ⟦AI⟧ 잠정 연구실·장치 계열 태그. 독립 장치 여부 인증이 아님. `findings.claim`은 출처 재서술, `model_implication`은 향후 AI 해석을 구분하기 위한 열이며 현재 비워 둠.

[스키마](SCHEMA.md) · [검색·접근 기록](SEARCH.md) · [검수](REVIEW.md) · [검증 결과](validation.json).

## ⟦AI⟧ 재생성·검증

저장소 루트에서 실행.

```powershell
python references/rydberg_experiments/build_db.py
python references/rydberg_experiments/build_db.py --check
python references/rydberg_experiments/verify_db.py
```

원본은 `papers/*.json`. CSV·통합 JSON·HTML·manifest는 생성물. 일상 갱신은 원본 편집 후 생성기 실행. `_curate_records.py`는 최초 수작업 큐레이션 기록, `_collect_sources.py`는 선언된 공개 URL의 최초 확보 보조 도구.

검증 범위는 문헌 DB·출처·CSV·검색 HTML. 수치 엔진·GC 캠페인은 이 데이터 작업의 검증 대상에서 제외.
