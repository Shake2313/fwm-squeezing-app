# ⟦AI⟧ 평균장 광이동·이득 기준면 재검증 — 2026-10-15

사용자 요청: 평균장 이득 불일치 우선 검토, 광이동·이득 기준면 독립 재검증,
임의 계수 fitting 없이 확인 가능한 코드 결함 패치, Grand Challenge 영향 검수.
직접 확인된 장치 조건은 **⁸⁵Rb F=2→F′=3 기준 blue OPD+1 GHz,
seed<pump, pure85(87약0.1% 가능)**. 외부 AI 인계의 나머지 수치는 조건부 입력이다.

## 판정과 실제 변경

| 항목 | 독립 검증 결과 | 처리 |
|---|---|---|
| 광이동 부호·단위 | 기존 4준위 Hamiltonian에 포함; 독립 복소 응답 상대오차3.09×10⁻¹⁴ | 임의 detuning offset 추가 없음 |
| 흡수 후처리 | pump0에서도 `Gs<1`을1로 올리던 오류 | cap와 기존 출력 fit 모두 흡수 보존 |
| 광학 carrier 원점 | scan은 F2→F′3, 상수는 fine-structure centroid | Option-A 파수·photon energy에 누락1.921502256 GHz 반영 |
| 측정 gain 환산 | atomic map에 창 손실 없음; total-cell T만으로 atomic A0 식별 불가 | 명시적 기준면 helper·누락 입력의 `UNAVAILABLE`·raw metadata 추가 |
| 절대 gain·peak | 이번 결함 수정으로 실험 불일치 해결 안 됨 | 모델·입력 의존성 유지, 새 fitting0 |

코드: `gabes/observables.py`, `gabes/fwm_gain_closure.py`,
`gabes/gain_reference.py`, `gabes/schemes/fwm.py`.
새 계수 fitting이나 기존 affine/participation 계수 재조정은 하지 않았다.
검산 스캔은 gain closure OFF지만 inherited residual0.74와 조건부 seed·η를
사용하므로 parameter-free 장치 예측으로 승격하지 않는다.

## 광이동 독립 경로

[유도·조건](independent_light_shift.md), [독립 구현](independent_light_shift.py),
[수치](independent_light_shift.json). 직접 matrix-action Lindblad와 별도 trace-zero
basis로 production과 분리하여 검산했다. Δ=ωpump−ωF2→F′3,
δ=ωseed−(ωpump−ωHF); 실험 beat축은 −δ다.

294.8mW, pump1/e² 반경530µm, Δ/2π=+1GHz, v=0에서
Ω/2π=495.958975MHz. 정확한 dressed-state Raman 이동52.269937MHz,
소산 포함 pole52.269382MHz. 2차 섭동62.511844MHz는 이 세기에서 충분히 정확하지 않다.
이 단일 속도 pole을 ensemble gain peak에 그대로 빼는 보정은 근거가 없다.
0.1µW도 좁은 개별 Raman pole 근처에서는 infinitesimal response와 약17% 차이가
남으므로 전 구간 weak-seed 기준으로 보증할 수 없다. 더 작은 seed로의 수렴은 검증했다.

기존 lumped CG² 축약의 strong-pump 정당화는 미완료다. signed 24준위
Hamiltonian의 한정된 q=0,B=0,v=0 진단은 same-m 이동37.330–44.202MHz를
주지만, 원자상태·편광·소산·gain 전파를 닫은 대체 모델은 아니다.
따라서 scalar CG 계수 재조정으로 실험을 맞추지 않았다.

carrier 원점 보정은 hyperfine Casimir 고유에너지와 degeneracy centroid에서
`(7ωHF+5ωexcHF)/12`로 독립 확인했다.
[Steck ⁸⁵Rb 원자상수](https://steck.us/alkalidata/rubidium85numbers.pdf).
이는 optical 파수 약5.1ppm 보정이며 atomic Δ를1.9215GHz 이동하는 패치가 아니다.
현재 atomic Doppler의 centroid 단일-k 근사와 `PHASE_LEGACY`의 fixed-k 근사는 유지한다.

## 이득 기준면과 이전 해석 정정

[독립 유도·반례](gain_reference_review.md),
[power-chain tests](../../../tests/test_fwm_gain_reference_planes.py).
atomic map은 vapor 입구 seed 대비 vapor 출구 power gain `g_j`다.
펌프를 끈 원자 투과 A0와 각 출력 측정면까지 투과 τj,τref에 대해

\[
G_{j,ref}=g_j\frac{\tau_j}{A_0\tau_{ref}},\qquad
P_{seed,vapor\ in}=\frac{P_{off,ref}}{A_0\tau_{ref}}.
\]

각 τ는 출력창·수집·후단 광학을 포함한다. 같은 후단이면 소거되며 입사창도
on/off 공통 비에서 소거된다. 전체 셀 투과 `win*A0*wout`만으로 A0를 알 수 없다.
측정마다 pump-off 출력을35mV로 재조정했다는 인계에 따르면 입력 seed도
A0(TPD)에 따라 변할 수 있다. 고정 seed 스캔을 사후 분모 환산만으로
장치와 일치시킬 수 없다. 검출 전압에는 별도 responsivity 교정도 필요하다.

이전 [잡음 판독 검수](../noise_readout_2026_10_15/README.md)의 `model/t_cell≈298`
예시는 현재 원자 전파의 올바른 측정면 환산으로 채택하지 않는다.
또 인계의 measured gap1.17에 곧바로 수동 source 조건 `D<=1`을 적용할 수 없다.
TPD+7의 output power 비191/35,150/35를 photon flux로 환산하면 gap1.17149757;
같은 출력 광학에서 A0≤0.85360826이면 source gap≤1 가능하다.
이 한 점으로 검출 교정 오차나 비수동성을 확정하지 않는다.
기존 canonical M 자체의 passive-completion 판정과는 다른 문제다.

생산 raw에는 small-signal atomic gains의 기준면과 미확정 A0를 기록한다.
기존 empirical 출력 fit을 켠 표시 gain은 별도로 출력 비·물리 기준면 미확정으로
표기한다. `cell_input_seed_W=3.7µW`라는 기존 calibration metadata도
입력 미확정, simulation 가정3.7µW, 측정 pump-off 기준3.7µW로 분리했다.
기존 fit이 명시적 η knob를 사용하지 않았다는 사실은 측정면 독립성을 보증하지 않는다.
인계 η=.706817033과 F4.05의 RF·측정면 역시 이번 채팅에서 독립 확인되지 않았다.

## 패치 후 수치·호환성

[401점 비교](production_comparison.json),
[현행 배열](current_production_spectrum.npz),
[재현 코드](production_comparison.py).
T386.15K,pump294.8mW,seed약3.918µW,530/330µm,L12.5mm,angle0.32°,
NF3/2,Balanced pole,closure OFF,internalδ±60MHz 조건이다.

| 조건부 결과 | 수정 전 | 수정 후 |
|---|---:|---:|
| 표시 peak Gs | 236.261616 | 236.269121 |
| 실험 beat축 peak(MHz) | −23.7 | −23.7 |
| Gs<1인 표시 점 | 0/401 | 89/401 |
| canonical passive-completion 불가 점 | 41/401 | 41/401 |

상정한 고정 vapor-input seed약3.918µW·동일 on/off 출력 throughput 조건에서
인계 peak870/35≈24.857,beat−2.2MHz와 비교하면 높이 최소9.505배,
위치21.5MHz 차이가 남는다. 이 조건의 수동 A0≤1 환산은
이 원자 gain을 낮추지 못한다. A0의 TPD 의존성·진짜 입력 seed는 별도다.
알려지지 않은 실제 seed와 unequal throughput까지 포함한 보편적 하한은 아니다.
원시 실험 자료를 이번에 새로 취득한 것은 아니다.

`legacy_mean_field_fixture=True`는 역사적 centroid 원점·흡수 unity floor를
명시적으로 복원하고 historical badge를 붙인다. 과거 snapshot 파일을 재생성하지
않았고 기준 스캔의38배열이 수정 전과 bitwise 일치했다. 새 기본값과 섞지 않는다.
후보 F′2의 과거 finite-grid 결과도 당시 source의 역사적 결과로 유지한다.
새 carrier에서 기존81점 안정성 검사를 반복한 것은 아니다.

Manley–Rowe 후처리 cap 자체는 기존 near-degenerate power 근사·독립 팔 포화다.
흡수 보존 수정이 이를 self-consistent depletion이나 exact photon-pair 불변량으로
만들지는 않는다. `physical_squeezing_dB=None`을 유지한다.

## 검증과 Grand Challenge 영향

새 power-plane·carrier·실제 pump0 limit·역사 snapshot·L1 연동 관련100검사 통과,
fast tier29검사 통과. 극단 random drift 스트레스(최대gain약3.2×10⁴⁶)의
NumPy/compiled 차이는 [독립60/90자리 expm](scorer_numerical_check.json)으로
확인했다. NumPy오차2.54×10⁻⁸,compiled1.20×10⁻⁹여서 compiled 결함이 아니었다.
기존 스트레스는 유지하고 근거 있는1e−7 허용오차를 사용한다.
별도 canonical drift 세 조건에는 더 엄격한1e−12 parity 검사를 추가했다.

[GC 최종 대조](final_gc_impact_comparison.json):21진입점 import,
34FWM 의존 binding·8observables 의존 node·상수 보존,4native 보고서의 전체
수치 payload 동일,37봉인 artifact/thermal ZIP byte 보존,기존 strict validator
거부 동작 보존. 묶음 import의 무관 alias 추가는 원본 AST 차이도 남긴 채
binding별로 비교했다. [변조 검출5control](gc_ast_scope_review_v2.json)도 통과했다.
GC 경로는 주어진 optical frequency/atomic detuning을 사용하며 이번 생산 carrier
helper를 호출하지 않는다. off-path observables 파일 hash 변화는 보존·공개한다.

전체 `python -m pytest -q` **2313통과·4skip·0실패**.
[검증 요약](code_validation.json),[완료 receipt](full_pytest_receipt.json),
[로그](full_pytest.log). 실행 전후 252개 source hash가 동일하다.
확인된 코드 결함·기준면 진단 범위는 종결하고 절대 이득 불일치는 열린 문제로 유지한다.

이번 production 비교6회(배열 dtype 비교기 실패 후 동일3회 재실행 포함),
소수 단일 원자·limit 검산,새 thermal transport solve0.
WP-INPUT/PREDICT에 기준면 비식별성과 strong-pump 축약의 열린 의존성을 연결한다.
기존6개 미시 WP·L1의 선언 범위, frozen v2FAIL, 공식 milestone **0/4** 유지.
다음 판별 질문은 같은 vapor-input 조건의 Doppler-averaged infinitesimal response와
finite-seed response가 gain의 높이·위치 차이를 얼마나 만드는가다.
실험 대조에는 독립 A0(TPD),팔별 throughput,seed input과 편광 계측이 필요하다.
