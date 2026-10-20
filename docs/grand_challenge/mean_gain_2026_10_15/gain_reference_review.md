# ⟦AI⟧ 이득 기준면 독립 검산 — 2026-10-15

출처: 현행 `gabes/observables.py`, `gabes/schemes/fwm.py`,
`gabes/fwm_gain_closure.py`; 외부 AI 인계 §3.7·§7.1.
사용자 확인은 blue OPD +1 GHz·낮은 주파수 seed·pure85에 한정.
아래 광학 투과 분해는 AI 유도; 독립 측정 입력으로 채택한 값 없음.

## 원자·창·측정면

`M`은 원자 susceptibility 전파. 입사·출사창의 별도 감쇠 항 없음.
모델 이득 `g_s,g_c`는 vapor 입구 seed power 대비 vapor 출구 power.
입사창 투과 `w_in`, 펌프 끈 원자 투과 `A0`, 출사창 투과 `w_out`.
원자 뒤 측정면까지 투과 `tau_s,tau_c`, pump-off 기준면 투과 `tau_ref`.
각 tau에는 해당 출력창·수집·후단 광학 투과 포함. 광전류·전압 비에는
추가 detector responsivity 교정 필요. 단위는 모두 광파워 비.

고정 입사 seed와 같은 on/off 장치 조건에서

\[
G_{s,ref}=\frac{g_s\tau_s}{A_0\tau_{ref}},\qquad
G_{c,ref}=\frac{g_c\tau_c}{A_0\tau_{ref}},\qquad
P_{seed,vapor\ in}=\frac{P_{off,ref}}{A_0\tau_{ref}}.
\]

같은 출력 광학이면 tau 소거. 입사창은 on/off 공통이므로 비에서 소거.
외부 셀 전체 투과 `t_cell=w_in*A0*w_out`만으로 `A0` 결정 불가.
따라서 현행 원자 전파 이득을 무조건 `t_cell`로 나누는 인계 식은 불충분.
전체 창을 포함한 gain을 분자로 쓸 때만 전체 셀 투과로 나눌 수 있음.
입·출사창이 전파에 없는 채 분모에만 넣으면 창 손실을 잘못 처리함.

인계 `20.6/26=0.7923077`은 다른 조건의 전체 투과. 이번 `A0`로 채택 안 함.
같은 on/off 후단이면 수동 `0<A0<=1` 보정은 원자 모델 peak 236.26을
측정 peak 약24.86 쪽으로 낮추지 못함. `A0=1`에서도 peak 비 약9.50.
`A0=.792` 가정은 약12.00배가 되지만 단순 민감도 예시. 상수 투과 환산은
peak 위치를 바꾸지 못함. TPD 의존 `A0`는 위치를 바꿀 수 있어 독립 측정 필요.

인계 §7.1은 pump-off 출력35mV를 측정마다 재조정했다고 기록.
따라서 실제 vapor 입사 seed도 `A0(TPD)`에 따라 달라질 수 있음.
현행 고정 약3.92µW 스캔은 조건부 입력. finite-seed 원자 상태에서는
결과 이득의 사후 분모 환산만으로 TPD별 입력 drive 차이를 교정할 수 없음.

## Power·photon flux

현행 `Q^-1 T_field Q` 변환과 면적 계수 맞음.
`G_c,power=(omega_c/omega_s)*G_c,flux`, probe gain은 두 규약 동일.
명목점 `omega_c/omega_s=1.00001610018`; 약16.1 ppm.
약21.5 MHz 위치·약9.5배 높이 차이의 설명 아님.

수동 commutator 조건은 canonical photon-flux map에 적용.
측정의 pump-off 출력 기준 `D_meas`에 곧바로 `D<=1` 적용하면 안 됨.
TPD 사용자 +7의 191/35,150/35 power gains에서 photon-flux gap은
1.17149757. 같은 후단이면 source gap은 `A0*1.17149757`.
`A0<=.85360826`이면 gap<=1 가능. `.792` 가정에서는 .92782608.
따라서 해당 점 하나만으로 검출 교정 오차·비수동성 확정 불가.
또한 gap 조건은 필요조건일 뿐 전체 2x2 passive completion 인증 아님.

입력 F도 기준면 필요. 독립 pure-loss 채널이면
`F_off=1+A0*(F_in-1)`이지만, 인계 F4.05의 정확한 측정면·RF 조건을
확정하지 않은 채 역산하지 않음. 인계 η 값도 이번 채팅에서 독립 확인되지
않았으므로 조건부로 유지하며 자동 변경하지 않음.

## 코드 결함·검수 범위

기존 `pump_depletion_saturation`의 `max(G_s-1,0)`는 흡수 `G_s<1`을
모두1로 올림. pump=0에서도 .4→1. 원자 전파·광이동이 아닌 표시 후처리 결함.
기존 baseline의 121점 중79/81점이 각각 unit floor. 과거 파일 덮어쓰기 없이
역사적 floor 정책과 현행 흡수 보존 정책 구분.
`legacy_probe_floor=False` 기본값은 양의 excess만 포화, 흡수 보존.
명시적 역사 fixture의 `True`만 기존 unit floor 재현. fit 계수 추가 없음.

기존 출력 affine fit은 pump-off 출력 기준 측정에 맞춘 경험식이며
canonical transfer는 교정하지 않음. `cell_input_seed_W=3.7e-6`는 측정으로
확정된 입력이라는 뜻으로 읽히므로, simulation seed 가정과 측정 pump-off
seed metadata 분리. `cell_input_seed_W=None`, 별도 simulation 가정3.7µW와
측정 pump-off 기준3.7µW 기록. 기존 계수 새로 fitting하지 않음.

새 `gabes/gain_reference.py`는 위 power-plane 식만 계산. `A0=None`이면
`UNAVAILABLE` 반환, 전체 셀 투과로 대체 안 함. 출력 팔 완전 차단은 허용;
분모의 `A0,tau_ref`는 양수 필요. 후단 투과는 전자 gain이 아닌0–1 광파워 비.
새 helper가 finite-seed drive를 재계산하거나 실제 장치를 인증하지 않음.

검증: `python -m pytest -q tests/test_fwm_gain_reference_planes.py
tests/test_fwm_propagation.py` 27통과(1.56초). 독립60개 광파워 chain,
Beer–Lambert·pump0·전체 투과 비식별 반례·출력 기준D>1 반례·배열 broadcasting·
양의 이득 구간 기존 cap/fit 동일성. 전체 suite·GC 전후 영향은 상위 검수에 기록.

GC 영향: 이 기준면 변환은 production 진단용. 미시 D·four-sideband·thermal
stream 정의를 변경할 근거 없음. 기존 봉인 자료·실패 판정·milestone 유지.
