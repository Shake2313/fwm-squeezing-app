# Rydberg opacity·현행 체크리스트 검수, 2026-10-23

사용자 요청: 독립 에이전트의 FWM·SABES 변경 후 체크리스트 현행성 재검수. 유효한 Rydberg opacity부터 반영. 기본값·기존 완료 범위 유지, 패치·검증·잔여 P1 보고.

## ⟦AI⟧ 결과

Rydberg opacity **완료**. 기존43개 done객체·기준·공식GC0/4 보존. 새 실험용FWM 완료범위1개 별도 등록. 전체102개 중45done, 잔여P1 **20개**: ready14·in_progress2·blocked_external2·needs_decision1·parked1.

숨은0.001 → `effective_opacity_participation`(0..1, navigate-only). 원자 line strength1·atom-number participation과 분리, 정적 투과/finite-IF/설정형 electrometry에 한 번 적용. 기존 `raw.ls`는 deprecated snapshot alias로 보존. Live읽기는 새 명명값 사용.

출처는 legacy 가정, **UNCALIBRATED**. 실측 표준불확도unknown. 기본 ±50%는 가정 시나리오이며 신뢰구간 아님. 온도·파워 sweep/표본최적점 **PREDICTED**. 원자 상태/IF phasor 재사용, 시나리오 추가 atomic solve0. Doppler-on 감도 시나리오는 명시한 nominal 단일 detuning 고정. 무응답 최적점null, unavailable 시나리오와 유한 범위를 구분.

기본 EIT linewidth **1.6130585827288382MHz**, T(0)**0.9400166150823321**. AT 총ASD **6.0823872151659275nV/cm/√Hz**, optimumΔ **−2.22MHz**. 스캔·χ·투과·χphys·IFphasor·ASD 배열 모두 **array-exact**. [수치/추가solve 검사](scheme_validation.json).

가정 opacity .0005/.001/.0015에서 sampled 온도최적점55/47/45°C. 범위끝55°C 포함한 제한된 격자 모델결과, 실험 운전점 권고 아님. 기본값이 같아도 가정 범위가 최적점 해석을 바꿈.

## ⟦AI⟧ 독립 변경 후 재검수

실험FWM Jasperse·양팔 gain/peak/width 보정·측정밖/조건변경 외삽 기구현. 현재v22에서 회귀검사. 과거 boundedPCHIP/10·16Jasperse closure는 역사범위 그대로 유지. 열린9개 항목의 stale 현재상태/전면 fit금지 문구만 최신 사용자 정책에 맞춰 정리; 이전 문구는 history 보존. 실험모델의 경험계수·명시한 유효covariance 허용, GC의 독립 미시/no-fit 관문은 별도 유지.

SABES 장치/profile/3card 개선 보존. `bridge.noise_observables`가 PassiveL1을 강제하고 결과속성이 rawgain/PassiveL1을 읽는 통합공백은 `sabes-mean-field-trust-boundary`에 구체 기준으로 남김. 광범위 source/SQL/RF 검증 완료로 승격하지 않음.

OD/SAS Transmission-only 선택은 의도한 현행 UI. 계산된 A10/tau/export와10·21 완료범위 유지. Rydberg/다른스킴의 기존 P1 계획은 해당범위에서 유효.

## ⟦AI⟧ 검증

- 집중: scheme/opacity **55PASS**, CLI/adapter **20PASS**. 경계0·상한·dedup·단일적용·cache·무응답·유효band·strictJSON·axial 포함.
- 실제 앱 AppTest 단일render2개: 기본 .001과 저장값 .0005, exception0, T(0).940016615→.969544540. Native browser gesture 검증과 구분.
- **종합 coverage 2747PASS·3skip·0미해결**, 최종 collection **2750개**. 단일 최종소스 전체실행PASS 주장은 아님. 최초 C전체2727PASS·3skip·1deselect·3FAIL·8ERROR(exit1) 보존. 모든11건은 AppTest registry 관련; defaults/liveplot의 실제 test-runner missing-registry만 처리. runtime/사용자callback 오류는 재발생. 전체9개 UI 소비파일과 변경된 provenance/SABES 테스트파일 재검사, thermal source/cache 통합2건 추가. 중복실행 수 합산 없이 고유 최종case 기준. [범위·원본실패·재검사receipt](full_pytest_receipt.json).
- 검사 중 독립에이전트3파일 변경: byte-keyed source-provenance cache·그 경계2검사·실제SABES Ultra fixture 재사용. 원래 물리/노이즈식 동일, 변경파일과 관련통합 재검사. 최종 Python **435파일** 원본/C byte동일·재검사 동안 불변. baseline→최종8파일 차이(새UI guard1 포함) 명시. 원본경로에 봉인된1건은 원본1PASS 증거 유지.
- 원본 불완전2회와 최초UI재검사 중간실패는 최종PASS집계 제외. Git metadata/filesystem 지연에 C검증 전환; 초기 사본 **6631파일** byte대조. 최초 original/C ordered2742 collection 동일, 최종은 동일Python 사본에서 collection. 원본 최종recollect는 지연으로 중단. 실패/불완전log 보존.
- 기존 done43·default모든기존값·SABES/quantum41Python·GC dependency34+8·봉인37개 보존. [검증](validation.json).

Rydberg readout 수렴/receiver 계약·실측 opacity/held-out 검증은 남음. 원자 감쇠 생성자·pointwise balance 변경 없음. 이번 완료는 opacity/sweep trust 범위.

## ⟦AI⟧ 남은 P1

| 항목 | 상태 |
|---|---|
| `collisional-coefficient-provenance-and-pressure-shift` | ready |
| `species-line-rabi-and-saturation-provenance` | ready |
| `lambda-spectrum-validity` | ready |
| `rydberg-readout-validity-guardrails` | needs_decision |
| `fwm-gain-only-trust-boundary` | in_progress |
| `fwm-balanced-readout-covariance-correction` | ready |
| `fwm-option-a-and-loss-ledger` | ready |
| `fwm-coupling-normalization-audit` | ready |
| `sabes-mean-field-trust-boundary` | ready |
| `sabes-calibration-trust-schema` | ready |
| `sabes-balanced-detection-and-sql-calibration` | blocked_external |
| `fwm-squeezing-spectrum-heldout-validation` | blocked_external |
| `fwm-reference-parameter-provenance` | ready |
| `fwm-finite-seed-angular-doppler-production` | ready |
| `fwm-rf-squeezing-spectrum-and-bandwidth` | parked |
| `fwm-tpd-gain-physical-replacement` | in_progress |
| `od-sas-normalization-reference-provenance` | ready |
| `scheme-current-documentation-and-examples` | ready |
| `magneto-paraffin-exchange-dephasing-contract` | ready |
| `fwm-fast-displayed-readout-error-certificate` | ready |

다음 후보: **SABES 선택모델/보정 전달 parity → Fast 표시전점 인증 → OD/SAS normalization/reference 출처**. AI 우선순위 제안, 새 사용자 제한 아님.

총책임자: PASS(명시한 유지보수 범위). GC 미시공분산·열 ensemble·공식관문 변화0. 기존 daily 취합cutoff10·18 유지. 새 batch 없음. 다음 질문: SABES의 cheap/full/selected/export가 동일 선택모델·gain보정 결과를 전달하는가?

## ⟦AI⟧ 채팅 중간 커밋 검증

별도 미커밋 FWM/SABES 변경을 제외한 커밋 구성은 [334고유검사 PASS](commit_checkpoint.json).
기존 HEAD seeded-FWM v20에 이 채팅의 biphoton 근거·읽기 변경만 결합. 독립 v22 DEFAULTS/preset/UI는 작업 트리에 보존.
공용 컴포넌트 guard의 독립 browser transport·작은 CSV 원본은 필수 의존성으로 포함.
첫333PASS·사본 SVG누락1FAIL 후 기존 자산 export 보완·해당1PASS; 제품 source/test 변경 없음.
앞의2747 종합 검증과 v22 현행성 검수는 당시 working-tree 기록이며, 이 커밋 전체의 동일바이트 검사로 승격하지 않음.
