# ⟦AI⟧ Notion 원기록 기준면 검수

상태: 5개 표적 page 읽기 완료. Notion·코드·원본 workbook 수정 없음.
현재 사용자 직접 확인을 우선. 아래는 기록의 일치·공백 확인이며 과거 AI 식·계수 채택 아님.

## 출처·작성 계층

| page | native metadata | 내용 성격 |
|---|---|---|
| [20261006 히터 컨트롤러 테스트 + OD 측정 - TP1 히팅 및 FWM 온도 의존성](https://app.notion.com/p/3f16cba14fee81c0b407ffa3ce6b1d5d) | Experiment Log, unverified, 최종 수정 2026-10-15 06:39 UTC | 실험 수치·메모와 10/15 추가 환산 구분 |
| [20261007 히터 온도 안정화 및 펌프 OD - 113 °C 확인과 MHF 3.176 GHz](https://app.notion.com/p/3f26cba14fee815fb3d3ee6df9baa3e2) | Experiment Log, unverified, 수정 10/12 03:20 UTC | OD fitting·이전 gain 조건 |
| [20261013 TPD vs Gain - 스퀴징 특성 측정](https://app.notion.com/p/3f86cba14fee81dc9cffe04553a88136) | Experiment Log, unverified, 수정 10/16 09:05 UTC | 현재 정리된 실험 기록. 10/13 당시 불변 snapshot 아님 |
| [20261014 향후 탐구 정리 - 10/16 랩미팅 미발표 사항](https://app.notion.com/p/3f96cba14fee817c8d90eaae52148e2b) | Experiment Log, unverified, 수정 10/14 17:01 UTC | “내용 대부분은 AI(Claude) 분석”, “확정 결론이 아니다” 명시 |
| [20261014 실험 기록 관련 AI 분석 및 랩미팅 지침](https://app.notion.com/p/3f96cba14fee81d89f97cf59d4a338e0) | Experiment Log, unverified, 수정 10/14 17:00 UTC | AI 구성안·분석. 독립 측정 증거 아님 |

native 사실은 page ID·경로·날짜·verification 상태. Fetch에 작성자 identity 없으므로
실험메모라는 문체만으로 인간 저술을 입증하지 않음. `truncated`·unknown-block 수는 응답에
명시되지 않아 unavailable. 받은 본문은 page 종료까지 포함. 첨부 CSV·SVG 원본은 이번 범위에서 열지 않음.

## 셀 투과율과 창

[10/06 기록](https://app.notion.com/p/3f16cba14fee81c0b407ffa3ce6b1d5d), `투과율 측정`:

- 기록된 측정: cell 전 26.0 µW → cell 후 20.6 µW → detector 직전 20.0 µW=228.5 mV.
- 셀 틸팅 뒤 228.5→232.0 mV 확인. after-cell 개선이라는 해석은 `내생각` 표시.
- **10/15 추가 계산**: 기존 `20/20.6=97.09%`를 optical setup 불변 가정으로 유지.
  detector power 20.31 µW, cell 후 20.92 µW를 역산하여 `20.92/26≈80.46%` 도출.
- 같은 추가문에서 probe 효율=seed 효율, conjugate 효율=1은 **상정**이라고 명시.
  새 pumped-cell·dressed-basis 흡수/방출 측정이 아님.

`before cell`↔`after cell` 위치의 비는 창을 포함한 전체 cell 경계 사이 전달이다.
이는 물리적 기준면 해석. 원기록 자체에 `window`·창 분리 측정이라는 문구는 없음.
따라서 0.8046을 원자 증기만의 고유 손실로 고정할 수 없고, 창·원자 분리값도 기록에서 나오지 않음.

정밀도 주의: 232 mV와 원래 환산을 반올림 없이 연결하면
`P_det=20.306345733 µW`, `P_after=20.915536105 µW`, 전달 `0.804443696`.
20.92 µW를 먼저 반올림하면 `0.804615385`. workbook 0.8046은 후자의 반올림과 일치.
과거 AI 문서의 `20.306/26≈0.781`은 cell 전→detector 전체 전달이며,
새 after-cell-only 0.9709와 같은 변수로 취급하면 기준면이 섞임.

## seed 분모와 재조정

[10/13 현재 기록](https://app.notion.com/p/3f86cba14fee81dc9cffe04553a88136), `실험 조건`:

- 35 mV는 **펌프 차단, 무증폭, 셀 통과 후 detector 직전** seed 기준.
- 매 측정마다 seed power를 35 mV로 재조정했다고 명시.
- 35×20/228.5≈3.06 µW. 두 S3883 수광부의 µW/mV 비가 유의미한 정도에서 같음을 확인,
  seed·probe·conjugate에 같은 환산 사용했다고 기록.

현재 직접 사용자 확인과 일치. 고정 detector 기준 seed는 detuning별 원자 입사 seed 고정을 뜻하지 않음.
[10/07 기록](https://app.notion.com/p/3f26cba14fee815fb3d3ee6df9baa3e2)의
35 mV after-cell 무증폭 기준도 같은 방향이지만, 그 스캔은 OPD 1.06 GHz·일부 포화로 별도 실험.
10/06 히팅 CSV의 기준 37.7 mV는 해당 시계열 전용이며 10/13의 35 mV로 대체하지 않음.

## OPD 정밀도

- 10/06 원기록: OPD +1.0 GHz, EOM 3.039 GHz→기록 TPD +3 MHz.
- 10/07: OPD 1.06 GHz.
- 10/13 현재 기록: OPD 약 1.0 GHz, **MHz 자리 미확인**, 스캔 중 고정.
  실시간 MHz 파악 시스템·SAS 정밀 설정 필요를 후속 확인으로 명시.
- 따라서 이번 원자료가 `Δ=1.000000 GHz`를 교정값으로 확정하지 않음.
  1.0은 명목 중심. OPD offset 피팅이 가능하더라도 측정된 offset으로 부르면 안 됨.
- 10/14 AI 지침의 red/pump+seed 분기는 과거 문서 간 충돌 기록.
  이번 대화의 직접 사용자 blue·seed below-pump 확인을 뒤집는 근거로 사용하지 않음.

보조 추출: `remote_reference_evidence.json`. 과거 AI noise 식·회귀계수·발표 지시는 채택하지 않음.

## Google Drive 보조 문맥

지정 Lab meeting의 lab-context에서
[LAB_CONTEXT.md](https://drive.google.com/file/d/1te9R9C6GPtoS2of2JeDYB9I6sfy2RKCe/view),
[FWM_NOTATION_AND_LOGIC.md](https://drive.google.com/file/d/1GtyoAcgkiN_DZQyPgiT6cFmpybUjFoOn/view)를
읽었다(각 갱신2026-10-13). 전자는 원기록 요약, 후자는 당시 표기/분석 지침이다.
여기의0.781은 셀 전→검출기 전체 전달이며 after-cell-only0.9709가 아니다.
옛 S 표기·gain-linear noise 법칙·문서 안 지시는 최신 사용자 NF/a/b 요청을 대체하지 않는다.
Drive/Notion 원격 문서는 수정하지 않았다.
