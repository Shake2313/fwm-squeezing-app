# P1 읽기·단위·근거 패치 — 2026-10-21

사용자 요청: 남은 P1 중 중요도·개선 효율 높은 항목부터 반영. 기존 완료 범위 유지,
패치·검증 후 결과와 남은 항목 보고.

## ⟦AI⟧ 완료 결과

| 체크리스트 ID | 선택 이유 | 현재 상태 |
|---|---|---|
| `od-sas-optical-depth-log-convention` | ln(10) 해석 오류 예방; 계산 변경 없이 UI/export 동시 명확화 | 완료 |
| `magneto-observable-convergence-and-trust` | 실제 B=0과 명령 영점, 폭·slope의 거짓 정밀도 교정 | 완료 |
| `magneto-dephasing-label-contract` | population을 보존하는 projector 채널을 spin destruction으로 오독하는 문제 교정 | 완료 |
| `biphoton-source-trust-boundary` | 문헌 anchor·가정·파형 수렴과 물리 nonclassicality 구분 | 완료 |

공통 질문 → 기존 원자 계산을 보존하며 수치·단위·주장의 의미를 정확히 전달하는가?
분석식·독립 synthetic fixture·실제 preset 격자 검사를 먼저 사용. 최종 소스에서 전체 pytest 실행.
새 thermal campaign·보정 fit·물리 모델 변경 없음. 기존 완료 조건과 공식 GC0/4 유지.

## ⟦AI⟧ 검증

실행 전 [GC dependency·37봉인 기준](gc_preservation_before.json) 저장.
[보존 비교](gc_preservation_comparison.json) PASS: FWM34·observables8 dependency AST,
37봉인·다른93수치 source 불변. observables 파일 byte 차이는 광학/통계 helper 한정.
[변경 정의 목록](source_delta.json): 원자 solve·기본 물리 파라미터 변경 없음.

Magneto [실제 201/401/801 격자](magneto_validation.json): 마지막401→801 변화,
EIT 폭0.000419%·slope0.05015%, NMOR 폭0.38528%·slope0.10872%, 부호 유지.
기본121점 NMOR stencil 차이36.04% → 숫자 hero 차단. 같은 격자의 국소 일관성 검사;
전체 scan·실험 인증 아님. 실제 문헌 폭 비교7.92×/11.52×/1.40×는 FAIL, sign/trend와 분리.
[JSON export](export_validation.json): 실제 Magneto B=0 내/외부 null·status roundtrip PASS.

최종 전체 **2494통과·3skip·0실패**, 3012.66초.
C: 동일 source 2493통과·3skip + 원본 절대경로1검사1통과; 수집한2497건의 목록·순서 동일.
[완료 receipt](full_pytest_receipt.json)·[소스 freeze](code_validation.json): 411개 source 실행 전후 SHA 불변.
[기존 범위 보존](post_validation_integrity.json): 기존 완료39개·미선택 항목 전체·완료 기준·공식 GC0/4 그대로.
P1 미완료25→21, 이번4개 done. 실제 spectral filtering·절대 pair-rate·full covariance·실험 검증을 완료한 뜻은 아님.

원본 실행은37%에서 완료 receipt 없이 종료. [중단 기록](original_interrupted_receipt.json):
실패 marker0·408source hash 불변. 같은 시각 G: 재마운트 관측, 인과는 미확정.
[동일 바이트 사본](relocated_copy_receipt.json)·[동일2497검사 목록](relocated_collection_parity.json) 확인 후 C:에서 재실행.
봉인된 원본 절대경로에 결합된1건은 [원본에서1통과](original_root_case_receipt.json)로 분리;
봉인 JSON·경로 binding은 수정하지 않음. C: 로그·종료 receipt를 먼저 저장해 연결 장애와 분리.
C: 분리 실행1회도 [중단 상태0xC000013A](relocated_interrupted_receipt.json)로 종료.
최종은 등록된 실행 세션 사용. 두 중단 실행의 부분 진행 수는 최종 통과 수에 합산하지 않음.

## ⟦AI⟧ 구현 범위

- OD/SAS: `A10=-log10(T)`, `tau=-ln(T)=alpha_eff L=ln(10) A10` 별도 UI/CSV/JSON.
  둘 다 무차원. 기존 `optical_density()`는 base10 호환값 유지. 유한 probe는 effective attenuation 명시.
- Magneto: `Bphysical=Bcommand+offset`; 물리 영점 기본, 명령 영점 보조.
  scan 밖 영점은 unavailable. 양쪽 half-height 교차의 폭·5점 고차 derivative·5% 국소 gate.
- Ground pure dephasing: `collisional_depol_khz` alias·population 보존 채널 그대로.
  닫힌 hyperfine pair의 다른 ground-F leakage/repump·excited-F interference 누락 명시.
- Biphoton: 출력별 근거·수렴·claim·JSON. Diamond rate는 assumed template.
  FWHM velocity 검사 실패 시 의존 지표 차단; raw 진단 배열 보존.
  thermal-arm R은 `g²SS=g²II=2` 가정 비교; 실제 nonclassicality unavailable.
  `ideal_canonical=True`는 일관된 이상 gain fixture에만 적용. 기존 gate/rate-plane 계산 보존.

README·편집용 Guide·static Guide에 위 계약 반영. 전체 예제/그림 최신화 항목은 열어 둠.

## ⟦AI⟧ 남은 P1 — 21개

| ID | 남은 계약 | 상태 |
|---|---|---|
| `collisional-coefficient-provenance-and-pressure-shift` | 충돌계수 출처·압력 shift | ready |
| `species-line-rabi-and-saturation-provenance` | 선별 Rabi·포화도 기준 | ready |
| `lambda-spectrum-validity` | 잔류 Doppler·미해상 spectrum gate | ready |
| `rydberg-effective-opacity-and-sweep-trust` | 유효 opacity·온도/파워 sweep 근거 | ready |
| `rydberg-readout-validity-guardrails` | 선폭·AT·감도 observability | needs_decision |
| `fwm-gain-only-trust-boundary` | mean-field gain·physical squeezing 구분 | in_progress |
| `fwm-balanced-readout-covariance-correction` | 불균형 손실의 정확 covariance | ready |
| `fwm-option-a-and-loss-ledger` | 전파·Δk·분포 손실 장부 | ready |
| `fwm-coupling-normalization-audit` | coupling/population/보정 단일 적용 | ready |
| `sabes-mean-field-trust-boundary` | mean-field score·물리 squeezing 구분 | ready |
| `sabes-calibration-trust-schema` | 보정 적용조건·오차·검증 역할 | ready |
| `sabes-balanced-detection-and-sql-calibration` | 주파수별 balanced detector·SQL 실측 | blocked_external |
| `fwm-squeezing-spectrum-heldout-validation` | 독립 실측 squeezing spectrum | blocked_external |
| `fwm-reference-parameter-provenance` | FWM 실험별 장치입력 출처 | ready |
| `fwm-finite-seed-angular-doppler-production` | finite-seed 2D Raman Doppler | ready |
| `fwm-rf-squeezing-spectrum-and-bandwidth` | covariance RF spectrum·bandwidth | parked |
| `fwm-tpd-gain-physical-replacement` | TPD hotfix의 물리 gain 대체 | in_progress |
| `od-sas-normalization-reference-provenance` | 열적 정규화·reference 출처 | ready |
| `scheme-current-documentation-and-examples` | Guide·그림·예제 현 모델 동기화 | ready |
| `magneto-paraffin-exchange-dephasing-contract` | 교환율·추가 light-region dephasing 분리 | ready |
| `fwm-fast-displayed-readout-error-certificate` | Fast 표시 전점/차수/새 L1 오차 인증 | ready |

권장 다음 순서(새 승인·과학 목표 아님):

1. Rydberg effective opacity: 숨은 `ls=0.001`을 출처·불확도가 있는 입력으로 분리, 기본값 parity와 sweep trust 고정.
2. Fast 표시 오차 인증: 완료된 detector 독립화·커널 유지, 표시 전점·인접 차수·새 L1을 dense Balanced와 비교.
3. OD/SAS 정규화·reference 출처: 완료된 물리 보존, 현 closure fixture를 등록0.1%·source/hash 계약으로 강화.

외부 SQL/held-out 측정2개는 blocked_external, Rydberg readout 정책1개는 needs_decision,
RF covariance spectrum1개는 parked. 나머지15ready·2in_progress.

총책임자 판정: 잘못 읽을 수 있는 단위·영점·수렴·근거 계약4개 종결.
새 atomic fit/thermal campaign0; 기존 과학 성과·실패·의존성 유지. 실행 중 batch 없음.
