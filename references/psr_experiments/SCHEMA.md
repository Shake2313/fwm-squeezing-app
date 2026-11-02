# PSR 문헌 DB 구조 · version 1

사용자 요청: PSR 관련 실험·이론·잡음·검출·매질·공간 모드 문헌을 폭넓게 확인하고, `references/fwm_experiments` 구조를 참조해 저장.

## ⟦AI⟧ 원본과 생성물

- 원본: `papers/<paper_id>.json`, UTF-8. 파일명과 ID 일치.
- 생성물: `database.json`, `papers.csv`, `experiments.csv`, `observations.csv`, `findings.csv`, `index.html`.
- CSV: Excel 호환 UTF-8 BOM. JSON·HTML: UTF-8. HTML은 데이터 내장; 네트워크 없이 검색 가능.
- `schema.json`: JSON Schema 2020-12. 선택적 추가 필드는 허용. FWM의 gain·seed·legacy coverage 요건은 PSR에 적용하지 않음.

논문 기본 필드: `schema_version`, `paper_id`, `title`, `authors`, `year`, `publication`, `scope`, `study_type`, `topics`, `sources`, `experiments`, `findings`. 제목은 원제 보존. `study_type`: `experiment`, `theory`, `mixed`, `review`. `scope.inclusion`: `primary`, `context_only`, `excluded`. 이론 논문의 `experiments`는 빈 배열 가능. 실험이 없는 논문을 측정 표본으로 세지 않음.

`experiments[]`는 `experiment_id`, `label`, `kind`, `conditions`, `results` 포함. 실제 운전점·공동 스캔·이론 계산·과거 문헌 인용을 `kind`로 구별. 서로 다른 최적값을 하나의 가상 운전점으로 조합하지 않음. `literature_summary`는 독립 실험 결과 아님.

`experiments.csv`의 `species`, `line`, `configuration`, `measurement`, `seed_mode`는 실험별 명시값을 우선하고, 없으면 논문 `scope`에서 상속. `apparatus_group`은 실험별 값 → 논문 값 순서. 열람본은 논문의 입력 모드·장치 그룹과 실험별 재정의를 함께 표시. 서로 다른 seed·source·sensor 구성을 같은 조건으로 취급하지 않음.

## ⟦AI⟧ 출처와 근거

`sources[]`: `source_id`, `url`, `access`, `version`. 원문 파일은 DB 기준 상대 경로 `path`, `sha256`, PDF의 `pages`를 함께 기록. 기존 저장소 원문은 `../` 경로로 재사용 가능. 저장소 밖의 절대 경로·경로 이탈은 금지. 새 원문은 `sources/` 보관. `accessed_on`, `access_method`, `access_note` 권장.

`access`: `full_text`, `metadata_only`, `abstract_only`, `unavailable`. 브라우저로 원문을 확인했으나 파일을 확보하지 못했으면 `full_text`, `path: null`, `access_method`, `access_note` 기록. 검증기는 로컬 snapshot 부재를 경고. 초록만 읽은 경우 `full_text`로 표시 금지.

각 수치·finding은 `source_id`와 페이지·절·그림을 지정하는 `locator` 포함. `full_text` 근거의 기본 `verification`은 `full_text_checked`. 원문 미확보 출처의 수치·finding은 `verification: metadata_only`를 **명시**해야 검증 통과. CSV·열람본에도 표시. 초록의 수치를 원문 검수값으로 승격하지 않음. 다른 출처의 원문을 확보했다면 그 출처 ID로 연결.

`findings[]`: `finding_id`, `claim`, `source_id`, `locator`, `category`, `model_implication`. `claim`은 원문 내용의 간결한 재서술. `model_implication`은 GABES에 대한 **⟦AI⟧ 해석·개선 제안**, 이미 구현·검증됐다는 뜻 아님. 문헌의 관측과 AI 제안은 별도 열 유지.

## ⟦AI⟧ 수치 규약

`conditions[]`·`results[]`의 공통 필드:

```json
{"quantity":"pump_power","value":10,"unit":"mW","value_type":"reported","basis":"total power","reference":null,"plane":"cell_input","source_id":"main","locator":"PDF p. 2, experimental setup","evidence":"Concise source paraphrase.","notes":""}
```

허용 `value_type`: `reported`, `approximate`, `range`, `lower_bound`, `upper_bound`, `derived`, `digitized`, `not_reported`. 미보고는 `value: null`; 0으로 대체 금지. `range`는 `value: null`, 수치 `lower`·`upper` 필요. `derived`는 `derivation.formula`와 입력 근거 유지. 모든 수치는 유한값. 보고하지 않은 불확도·빔 반경 정의·기준면은 추정해 채우지 않음.

`squeezing_measured`, `squeezing_electronic_noise_subtracted`, `squeezing_loss_corrected`는 **SQL 아래 감소량을 양수 dB**로 기록. 예: 3 dB below SQL → `value: 3`. 부호를 가진 정규화 잡음 자체는 별도 `noise_level_min` 등으로 기록하고 기준을 설명. 전자잡음 차감·손실 역보정·직접 측정은 다른 quantity. 서로 대체하지 않음. Anti-squeezing, 유한 분석 주파수, 위상, LO 모드·손실, 검출 기준면도 분리.

`basis`, `reference`, `plane`은 원문 정의 보존. `unit`은 원문 표기; CSV `normalized_unit`은 µW/uW, µm/um, °C/degC, %/percent 별칭만 정규화하며 값의 배율은 바꾸지 않음. 그림 축 범위는 측정점 아님. 그림 수치 추출은 `digitized`와 방법·해상도·불확도 기록 필요.

## ⟦AI⟧ 재생성·검증

```powershell
python references/psr_experiments/build_db.py
python references/psr_experiments/build_db.py --check
python references/psr_experiments/verify_db.py
```

`--check`는 원본·출처 검증 후 모든 생성물의 예상 바이트를 비교. 파일을 쓰지 않음. 원본·문서·생성기 변경 후에는 재생성 필요. `manifest.json`은 입력·생성물·원문의 SHA-256 기록. 자기 자신과 검증 영수증은 순환 해시에서 제외.

검증 범위: schema·ID·출처 연결·페이지 수·파일 해시·범위·환산식·원문 접근 상태·CSV round-trip·stale view·잘못된 데이터 거부. 과학적 의미는 별도 원문 검수 대상. 이 자료 구축은 GABES 수치 예측·실험 검증·보정 학습·미열람 holdout을 생성하지 않음.
