# Quotient structure 과학 내용 보존·검토 기록

기준일: 2026-09-18. 개발·판본 비교 기록은 이 파일에만 둠. 과학 문서는 [독립 HTML](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html>).

## 산출물과 재현

- `fwm_quotient_structure_v4_builder.py`: 기존 61개 노드·98개 연결의 독립 데이터와 HTML 생성기.
- `fwm_quotient_structure_v4_science.py`: 순서론·조건부 미시 이론·열원자 transport 정의와 확장 데이터.
- `fwm_quotient_structure_v4_sources.json`: 현행 analytic reconstruction의 188개 label–section 제목 사본. 이전 판본 실행 의존 없음.
- `fwm_quotient_structure_v4_view.js`, `fwm_quotient_structure_v4_view.css`: 독립 탐색 UI 자산. 생성 HTML 안에 내장.
- `fwm_quotient_structure_v4.html`: 모든 87개 노드·142개 연결·수식·조건·근거를 정적 본문과 내장 graph data에 동시 저장. JS·네트워크·동반 파일 없이 과학 내용 열람 가능. 추가 근거 TeX/PDF/코드 링크는 선택 사항.

재현: 저장소 루트에서 `python "docs/FWM physics and analytic reconstruction/fwm_quotient_structure_v4_builder.py"`. 동반 analytic v3 TeX가 있으면 사용된 TeX key가 실제 source label에 모두 존재하는지 엄격 검사. 없으면 독립 catalog로 생성 가능. v2/v3 builder·HTML·notes import/read 없음.

## 보존 기준과 정정

기존 v3의 61개 node ID, 98개 edge record 전체 보존. 방향·종류·label·pending 값도 동일. 모든 기존 수식·입력·출력·차원·노드 제목 보존. 본문 범위 정정·추가 6개 노드와 근거 key 추가 1개 노드는 아래 표에 표시. 나머지 54개 노드는 모든 field 동일.

| 정정·추가 대상 | 근거와 처리 |
|---|---|
| `i-drive` | 생산의 temperature→density와 연구의 독립 temperature/density 입력 구분. 파워·반경·Rabi의 역할 유지. |
| `d-chi` | 생산 1/12·경험 0.74와 연구 uniform-Zeeman-RMS의 manifold별 평균 구분. 기존 수식 보존. |
| `m-drift` | 미시 drift 일반 미구현처럼 읽히는 표현 수정. 조건부 single-velocity/two-band/stationary 경로와 미완료 moving-thermal 연결 분리. |
| `m-reservoir` | atomic jump-product diffusion 구현과 QRT 대조 반영. commutator로 symmetrized noise를 정할 수 없다는 조건 유지. |
| `o-claim` | 조건부 microscopic readout 존재와 실제 장치 미검증 동시 명시. |
| `r-geometry` | 기존 본문 유지 후 2D rms 수치 추가: 해석 1.380173564 MHz, 40×40·5σ quadrature 1.380163304 MHz, 차이 −7.43×10⁻⁴%. |
| `r-poles` | 본문·식 유지. `tab:analytic-verification` 근거 key 추가. |

그룹 `micro`, `mean`, `output`의 개요·검증 범위도 위 조건과 맞춤. 기존 7행 검증 표, 유한 시드 4행 표, 수학 설명 6개, 용어집, double-Λ 도식, TMSV·독립 손실·seeded 통계 보존. Scoped analytic fixture 2행 추가: spectral resolvent 1.0×10⁻¹², 별도 Floquet/literal-SI 대조 ≤1.64×10⁻¹⁴. N_F=1→2에서 G_p 16.8%·계수 갭 부호 반전도 조건이 있는 수치 반례로 복원.

v2의 과학 개념 54개 전부 현행 node로 매핑. 미확립 주장은 사실처럼 보존하지 않음:

- 공통 규약 패널로 옮기면 자동으로 N-free/SP quotient가 된다는 주장 제거. D<X, D<F, E<F의 induced N 반례와 명시적 준순서·동치류 정리로 대체.
- 검증 차이를 자연변환·2-cell 진폭으로 부르던 표현 제거. 정의된 사상·합성 법칙 없는 범주론적 주장이었음. 종류·영역·가정·오차 기준이 있는 증거 관계로 대체.
- 생산 finite-seed Floquet와 pump-only 정상상태 참조를 별도 경로로 명시. v2의 혼합 도식 순서를 생산 실행 순서로 보존하지 않음.
- K→microscopic D와 K→vacuum D를 동등한 물리 대안으로 해석하지 않음. K는 commutator 제약, D는 reservoir correlation 입력. vacuum completion은 대수 시험.
- rate+Sylvester 제거는 선언한 scale hierarchy의 조건부 근사로 보존. Ω_pump/(2π)=0.708 GHz, excited splitting 361.58 MHz에서 약한 혼합 제거가 정당화되지 않는 반례 포함. 모든 rate model 실패로 일반화하지 않음.
- `(G_p−g_c)²/(G_p+g_c)` 점수와 `(0.2,0.1)→0.0333` 반례 복원. matched SQL의 physical squeezing으로 사용하지 않음.
- v2의 검증 사다리·프리즘·개정 내역·이전 버전 오류 서술은 과학 문서에서 제외. 실제 수치 검증과 반례는 본문에 재정의하여 유지.

## 순서론의 정확한 범위

비교 context `𝒦=(D,Γ,관측 interfaces,허용 Π)` 고정. `P_A=π∘P_B`를 만족하는 허용 회수 map이 있으면 `A≼B`. Π는 항등·합성에 닫힘. 입력별 임의 보상·target fit 불허. 반사·추이의 준순서이며 전체 순서 아님. `A∼B ⇔ A≼B ∧ B≼A`; `[A]≤[B] ⇔ A≼B`는 대표 선택에 무관한 부분순서. 표현 변환과 channel 합성을 함께 몫내려면 interface와 congruence 별도 확인.

가정의 함의는 `Mod(Γ₂)⊆Mod(Γ₁)`로 별도 정의. 수치 설정 세분은 동일 모델에서만 비교. 실험 근거는 주장별 조건부 증거. 물리 확장·수치 세분·증거 추가를 하나의 정확도 순서로 놓지 않음. 임의 meet/join·격자·오차 단조 감소·CP/QRT→실험 정확도 주장 없음. ε-근접은 일반적으로 비추이적이므로 정확한 몫 증명에 사용하지 않음.

Loewner 잡음 순서는 같은 transfer/input/readout/SQL 조건에만 적용. `X=0`, `Y_s=diag(s,s⁻¹)/2`의 CP prepare channels는 서로 비교 불가. 이 집합에 보편적 least noise 없음. least/minimal/특정 scalar 최적화 구분.

## 검증 범위

정적 검증: 기존 node·edge 전체 보존, 모든 endpoint 유효, 중복 HTML ID 없음, 내부 anchor 누락 없음, 87개 정적 node record 모두 존재. 네트워크 runtime script/link 없음. 외부 view JS 및 내장 JS 2개 모두 `node --check` 통과. 본문 수식은 원래 61개 모두 동일. 추가 수식·수치는 명시한 GC 유도·checkpoint·원시 ensemble artifact와 대조.

최종 catalog: 188개 theory label 중 106개 직접 참조, logic·GC key 20개 포함 총 126개 근거 key. 참조 label 누락 0. 현행 TeX SHA256 `20b052e8fdbdb10fd3015a4f9725ff9e6171861a36d6b35ff2cd470ea9592527`; HTML 생성 주석에도 원본 바이트 hash 내장. 고정된 관측 계약의 정보 동치와 모든 물리 관측의 동치를 구별하고, 가정 함의는 논리적 동치로 몫내기 전 준순서임을 명시.

최종 브라우저 검토: 루트 agent가 현행 HTML의 전체 지도, L/`q-quotient`, T/`t-convergence`, 검증 보기, detail→정적 본문 이동, 기본 screenshot 확인. JS 오류 없음. 정적·탐색 배지 일치와 동반 PDF 경로 수정 확인. 전체 pytest·동반 산출물 종합 검사는 루트 agent 담당. 이 문서 변경으로 새 numerical campaign 실행하지 않음.

동반 링크 검사: analytic v3 TeX/PDF, 코드, GC 유도·원시 ensemble JSON 모두 존재 확인. 검토 시점 `../squeezing_report/squeezing_report_v8.pdf`만 동시 제작 중이었음. 경로는 담당 agent 확인값. 최종 파일 생성 후 루트 agent의 종합 검증에서 재확인 필요.

아래 inventory는 생성 HTML의 내장 JSON에서 직접 산출. 보존 여부가 완료율·물리 정확도 백분율을 의미하지 않음.

## v3 node 61개 보존 inventory

| 그룹 | ID·현행 위치 | 제목 | field 보존 |
|---|---|---|---|
| input | [i-levels](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#i-levels>) | 네 준위와 double-Λ | 전체 동일 |
| input | [i-frequency](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#i-frequency>) | 광학 carrier와 RF 축 | 전체 동일 |
| input | [i-drive](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#i-drive>) | 입사 파워와 원자 환경 | 위 근거에 따라 body 보완 |
| input | [i-modes](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#i-modes>) | 수집 모드와 입력 평균장 | 전체 동일 |
| input | [i-state](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#i-state>) | 입력 양자 상태 | 전체 동일 |
| input | [i-detection](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#i-detection>) | 셀 이후 측정 조건 | 전체 동일 |
| atomic | [a-assembly](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#a-assembly>) | 구동과 산일자 조립 | 전체 동일 |
| atomic | [a-probe](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#a-probe>) | probe 기준장 해 | 전체 동일 |
| atomic | [a-conjugate](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#a-conjugate>) | conjugate 기준장 해 | 전체 동일 |
| atomic | [a-readout](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#a-readout>) | 가중 polarization 읽기 | 전체 동일 |
| atomic | [a-average](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#a-average>) | 생산 1D 속도 평균 | 전체 동일 |
| atomic | [a-order](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#a-order>) | 인접 차수 전체 경로 | 전체 동일 |
| atomic | [a-gate](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#a-gate>) | 전체 스캔 수렴 판정 | 전체 동일 |
| reference | [r-frame](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-frame>) | minus 정적 pump 프레임 | 전체 동일 |
| reference | [r-state](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-state>) | trace-one 정상상태 | 전체 동일 |
| reference | [r-gauge](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-gauge>) | 프레임·seed 극한 검증 | 전체 동일 |
| reference | [r-source](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-source>) | trace-zero 구동과 주파수 | 전체 동일 |
| reference | [r-direct](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-direct>) | 구속 직접 선형해 | 전체 동일 |
| reference | [r-poles](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-poles>) | 쌍직교 pole와 residue | 위 근거에 따라 refs 보완 |
| reference | [r-defective](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-defective>) | 결함 행렬의 대안 | 전체 동일 |
| reference | [r-response](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-response>) | Nambu 응답 조립 | 전체 동일 |
| reference | [r-geometry](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-geometry>) | 속도와 비공선 detuning | 위 근거에 따라 body 보완 |
| reference | [r-average](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-average>) | 복소 응답의 속도 구적 | 전체 동일 |
| reference | [r-voigt](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-voigt>) | 단일극 Faddeeva 극한 | 전체 동일 |
| mean | [d-chi](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#d-chi>) | 무차원 χ로 변환 | 위 근거에 따라 condition 보완 |
| mean | [d-geometry](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#d-geometry>) | bare k와 기하 위상부정합 | 전체 동일 |
| mean | [d-maxwell](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#d-maxwell>) | Maxwell 행렬과 정준 변환 | 전체 동일 |
| mean | [d-constant](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#d-constant>) | 상수 drift의 행렬지수 | 전체 동일 |
| mean | [d-segment](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#d-segment>) | z 의존 전파 | 전체 동일 |
| mean | [d-amplitude](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#d-amplitude>) | 출력 복소 평균장 | 전체 동일 |
| mean | [d-gains](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#d-gains>) | 두 이득과 정규화 | 전체 동일 |
| gain | [x-flux](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#x-flux>) | 광자속 gap | 전체 동일 |
| gain | [x-canonical](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#x-canonical>) | 전체 정준성 조건 | 전체 동일 |
| gain | [x-ideal](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#x-ideal>) | 밝은 seed 이상 비교 | 전체 동일 |
| micro | [m-pair](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#m-pair>) | Nambu pairing | 전체 동일 |
| micro | [m-drift](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#m-drift>) | 원자 요동 drift 구성 | 위 근거에 따라 body, condition 보완 |
| micro | [m-reservoir](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#m-reservoir>) | 미시 저장고 상관 | 위 근거에 따라 body, condition 보완 |
| micro | [m-k](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#m-k>) | 교환관계 채널 분류 | 전체 동일 |
| micro | [m-diffusion](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#m-diffusion>) | symmetrized diffusion | 전체 동일 |
| micro | [m-vacuum](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#m-vacuum>) | 독립진공 완성 시험 | 전체 동일 |
| covariance | [v-propagator](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#v-propagator>) | 요동 전파자 | 전체 동일 |
| covariance | [v-input](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#v-input>) | 입력 요동의 전파 | 전체 동일 |
| covariance | [v-integral](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#v-integral>) | 셀 내부 분포 잡음 | 전체 동일 |
| covariance | [v-lyapunov](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#v-lyapunov>) | 상수계수 Lyapunov 문제 | 전체 동일 |
| covariance | [v-phi](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#v-phi>) | 좌표무관 φ₁ 해 | 전체 동일 |
| covariance | [v-eigen](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#v-eigen>) | 조건부 고유기저 전개 | 전체 동일 |
| covariance | [v-vanloan](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#v-vanloan>) | Van Loan 블록지수 | 전체 동일 |
| covariance | [v-segments](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#v-segments>) | z 의존 국소 채널 합성 | 전체 동일 |
| covariance | [v-output](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#v-output>) | 두 항의 합류 | 전체 동일 |
| covariance | [v-quadrature](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#v-quadrature>) | companion과 quadrature 조립 | 전체 동일 |
| detector | [f-loss](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#f-loss>) | 공분산의 외부 순수 손실 | 전체 동일 |
| detector | [f-carrier](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#f-carrier>) | 검출 carrier와 위상 | 전체 동일 |
| detector | [f-current](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#f-current>) | 광전류 선형화 | 전체 동일 |
| detector | [f-sql](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#f-sql>) | 같은 검출 조건의 SQL | 전체 동일 |
| detector | [f-measure](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#f-measure>) | 검출 covariance 투영 | 전체 동일 |
| detector | [f-tmsv](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#f-tmsv>) | 이상 TMSV의 공유 pair 수 | 전체 동일 |
| detector | [f-thinning](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#f-thinning>) | 두 팔의 독립 순수 손실 | 전체 동일 |
| detector | [f-difference](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#f-difference>) | 차동의 열적 요동 상쇄 | 전체 동일 |
| detector | [f-states](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#f-states>) | seeded 상태와 thermal 환경 | 전체 동일 |
| output | [o-spectrum](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#o-spectrum>) | 정규화와 dB | 전체 동일 |
| output | [o-claim](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#o-claim>) | 예측의 성립 범위 | 위 근거에 따라 condition 보완 |

## v3 연결 98개 보존 inventory

각 행은 source·target뿐 아니라 label·kind·pending까지 기존 record와 동일. AND 필요 입력, 대안, 축약, 검증을 일반 순서관계 하나로 해석하지 않음.

| # | 출발 | 도착 | 종류 | 연결 의미 | 미완료 |
|---|---|---|---|---|---|
| 1 | [i-levels](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#i-levels>) | [a-assembly](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#a-assembly>) | data | 준위·전이 강도 | 아니오 |
| 2 | [i-frequency](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#i-frequency>) | [a-assembly](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#a-assembly>) | data | Δ·δ·branch·beat | 아니오 |
| 3 | [i-drive](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#i-drive>) | [a-assembly](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#a-assembly>) | data | Rabi·N·완화율 | 아니오 |
| 4 | [i-drive](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#i-drive>) | [a-average](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#a-average>) | data | 속도 가중치 | 아니오 |
| 5 | [i-frequency](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#i-frequency>) | [d-geometry](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#d-geometry>) | data | bare carrier 주파수 | 아니오 |
| 6 | [i-modes](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#i-modes>) | [d-geometry](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#d-geometry>) | data | L·θ·면적 | 아니오 |
| 7 | [i-modes](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#i-modes>) | [d-amplitude](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#d-amplitude>) | data | 입력 복소 carrier | 아니오 |
| 8 | [i-state](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#i-state>) | [v-input](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#v-input>) | data | V_in | 아니오 |
| 9 | [i-state](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#i-state>) | [v-segments](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#v-segments>) | data | 초기 V_in | 아니오 |
| 10 | [i-state](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#i-state>) | [v-quadrature](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#v-quadrature>) | data | companion 입력 상태 | 아니오 |
| 11 | [i-detection](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#i-detection>) | [f-loss](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#f-loss>) | data | η·vacuum port | 아니오 |
| 12 | [i-detection](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#i-detection>) | [f-carrier](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#f-carrier>) | data | 동일 η | 아니오 |
| 13 | [i-detection](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#i-detection>) | [f-current](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#f-current>) | data | H·w | 아니오 |
| 14 | [i-detection](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#i-detection>) | [f-sql](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#f-sql>) | data | SQL 보정 | 아니오 |
| 15 | [i-detection](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#i-detection>) | [f-measure](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#f-measure>) | data | S_el | 아니오 |
| 16 | [a-assembly](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#a-assembly>) | [a-probe](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#a-probe>) | data | probe 기준장 블록 | 아니오 |
| 17 | [a-assembly](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#a-assembly>) | [a-conjugate](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#a-conjugate>) | data | conjugate 기준장 블록 | 아니오 |
| 18 | [a-probe](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#a-probe>) | [a-readout](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#a-readout>) | data | ρ₀ᵃ·ρ₊ᵃ | 아니오 |
| 19 | [a-conjugate](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#a-conjugate>) | [a-readout](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#a-readout>) | data | ρ₀ᵇ·ρ₊ᵇ | 아니오 |
| 20 | [a-readout](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#a-readout>) | [a-average](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#a-average>) | data | 네 복소 χ̄ | 아니오 |
| 21 | [a-average](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#a-average>) | [d-chi](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#d-chi>) | data | 평균 reduced response | 아니오 |
| 22 | [i-levels](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#i-levels>) | [r-frame](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-frame>) | data | 공유 원자 조립 | 아니오 |
| 23 | [i-drive](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#i-drive>) | [r-frame](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-frame>) | data | pump·완화율 | 아니오 |
| 24 | [r-frame](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-frame>) | [r-state](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-state>) | data | ℒ_pump | 아니오 |
| 25 | [r-state](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-state>) | [r-source](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-source>) | data | ρ_ss | 아니오 |
| 26 | [i-frequency](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#i-frequency>) | [r-source](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-source>) | data | δ·독립 Ω_SA | 아니오 |
| 27 | [r-source](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-source>) | [r-direct](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-direct>) | data | trace-zero 구동 | 아니오 |
| 28 | [r-direct](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-direct>) | [r-response](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-response>) | data | 직접 polarization | 아니오 |
| 29 | [r-response](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-response>) | [r-average](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-average>) | data | 속도별 복소 도함수 | 아니오 |
| 30 | [r-geometry](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-geometry>) | [r-average](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-average>) | data | 1D / opt-in 2D 가중치 | 아니오 |
| 31 | [i-drive](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#i-drive>) | [r-geometry](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-geometry>) | data | T·σ_v | 아니오 |
| 32 | [i-modes](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#i-modes>) | [r-geometry](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-geometry>) | data | θ·k | 아니오 |
| 33 | [d-chi](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#d-chi>) | [d-maxwell](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#d-maxwell>) | data | physical χ | 아니오 |
| 34 | [d-geometry](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#d-geometry>) | [d-maxwell](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#d-maxwell>) | data | k·Δk·Q | 아니오 |
| 35 | [d-maxwell](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#d-maxwell>) | [d-constant](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#d-constant>) | data | 상수 M_cl | 아니오 |
| 36 | [d-constant](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#d-constant>) | [d-amplitude](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#d-amplitude>) | data | T_cl | 아니오 |
| 37 | [d-constant](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#d-constant>) | [d-gains](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#d-gains>) | data | 정준 T | 아니오 |
| 38 | [d-segment](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#d-segment>) | [d-amplitude](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#d-amplitude>) | data | 순서 있는 전파자 | 아니오 |
| 39 | [d-segment](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#d-segment>) | [d-gains](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#d-gains>) | data | 정규화한 작은 신호 T | 아니오 |
| 40 | [d-gains](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#d-gains>) | [x-flux](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#x-flux>) | data | G_p·G_c | 아니오 |
| 41 | [d-gains](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#d-gains>) | [x-canonical](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#x-canonical>) | data | full complex T | 아니오 |
| 42 | [d-amplitude](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#d-amplitude>) | [f-carrier](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#f-carrier>) | data | 복소 출력 carrier | 아니오 |
| 43 | [m-pair](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#m-pair>) | [m-drift](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#m-drift>) | data | 정준 paired basis | 아니오 |
| 44 | [m-drift](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#m-drift>) | [m-k](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#m-k>) | data | M_q·J | 아니오 |
| 45 | [m-reservoir](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#m-reservoir>) | [m-k](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#m-k>) | data | J_F·B 교환관계 | 아니오 |
| 46 | [m-reservoir](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#m-reservoir>) | [m-diffusion](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#m-diffusion>) | data | B·N_res | 아니오 |
| 47 | [m-k](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#m-k>) | [m-vacuum](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#m-vacuum>) | data | K와 추가 vacuum 선택 | 아니오 |
| 48 | [v-propagator](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#v-propagator>) | [v-input](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#v-input>) | data | T_q(L,0) | 아니오 |
| 49 | [v-propagator](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#v-propagator>) | [v-integral](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#v-integral>) | data | T_q(L,z) | 아니오 |
| 50 | [v-input](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#v-input>) | [v-output](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#v-output>) | data | 입력 전파 항 | 아니오 |
| 51 | [v-output](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#v-output>) | [v-quadrature](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#v-quadrature>) | data | paired output | 아니오 |
| 52 | [f-loss](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#f-loss>) | [f-measure](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#f-measure>) | data | V_R,det | 아니오 |
| 53 | [f-carrier](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#f-carrier>) | [f-current](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#f-current>) | data | 검출 진폭·위상 | 아니오 |
| 54 | [f-carrier](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#f-carrier>) | [f-sql](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#f-sql>) | data | 동일 검출 파워 | 아니오 |
| 55 | [f-current](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#f-current>) | [f-measure](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#f-measure>) | data | 측정 벡터 m | 아니오 |
| 56 | [f-sql](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#f-sql>) | [f-measure](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#f-measure>) | data | 같은 SQL | 아니오 |
| 57 | [f-measure](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#f-measure>) | [o-spectrum](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#o-spectrum>) | data | PSD·SQL | 아니오 |
| 58 | [o-spectrum](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#o-spectrum>) | [o-claim](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#o-claim>) | data | 주장 범위 | 아니오 |
| 59 | [f-tmsv](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#f-tmsv>) | [f-thinning](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#f-thinning>) | data | 공유 pair number | 아니오 |
| 60 | [f-thinning](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#f-thinning>) | [f-difference](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#f-difference>) | data | count moments | 아니오 |
| 61 | [a-average](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#a-average>) | [a-order](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#a-order>) | validation | 두 차수 응답 | 아니오 |
| 62 | [d-gains](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#d-gains>) | [a-order](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#a-order>) | validation | 두 차수 정준 전파 | 아니오 |
| 63 | [a-order](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#a-order>) | [a-gate](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#a-gate>) | validation | full-scan 비교량 | 아니오 |
| 64 | [r-state](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-state>) | [r-gauge](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-gauge>) | validation | gauge map 대상 | 아니오 |
| 65 | [a-average](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#a-average>) | [r-gauge](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-gauge>) | validation | 별도 finite-seed 극한 비교 | 아니오 |
| 66 | [r-direct](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-direct>) | [r-gauge](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-gauge>) | validation | 참조 응답 비교 | 아니오 |
| 67 | [r-source](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-source>) | [r-poles](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-poles>) | alternative | 대각화 가능시 | 아니오 |
| 68 | [r-poles](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-poles>) | [r-response](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-response>) | alternative | 같은 응답의 spectral 표현 | 아니오 |
| 69 | [r-defective](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-defective>) | [r-response](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-response>) | alternative | 안정한 응답 평가 | 아니오 |
| 70 | [v-lyapunov](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#v-lyapunov>) | [v-phi](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#v-phi>) | alternative | 좌표무관 행렬함수 | 아니오 |
| 71 | [v-lyapunov](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#v-lyapunov>) | [v-eigen](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#v-eigen>) | alternative | 조건 좋은 고유기저 | 아니오 |
| 72 | [v-lyapunov](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#v-lyapunov>) | [v-vanloan](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#v-vanloan>) | alternative | 권장 블록지수 | 아니오 |
| 73 | [v-phi](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#v-phi>) | [v-output](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#v-output>) | alternative | 같은 W | 아니오 |
| 74 | [v-eigen](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#v-eigen>) | [v-output](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#v-output>) | alternative | 같은 W · 조건수 guard | 아니오 |
| 75 | [v-vanloan](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#v-vanloan>) | [v-output](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#v-output>) | alternative | 같은 W | 아니오 |
| 76 | [v-segments](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#v-segments>) | [v-output](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#v-output>) | alternative | z 의존 전체 V | 아니오 |
| 77 | [r-source](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-source>) | [r-defective](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-defective>) | limit | 결함 행렬 조건 | 아니오 |
| 78 | [r-response](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-response>) | [r-voigt](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-voigt>) | limit | 단일 선형 pole·상수 residue | 아니오 |
| 79 | [d-maxwell](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#d-maxwell>) | [d-segment](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#d-segment>) | limit | z 의존 계수 | 아니오 |
| 80 | [x-flux](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#x-flux>) | [x-ideal](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#x-ideal>) | limit | 이상 가정 추가 | 아니오 |
| 81 | [x-canonical](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#x-canonical>) | [x-ideal](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#x-ideal>) | limit | 정준성·입력·검출 조건 추가 | 아니오 |
| 82 | [v-integral](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#v-integral>) | [v-lyapunov](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#v-lyapunov>) | limit | 상수 M_q·D | 아니오 |
| 83 | [i-levels](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#i-levels>) | [m-drift](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#m-drift>) | data | 구동 원자 모델 | 예 |
| 84 | [r-average](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-average>) | [m-drift](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#m-drift>) | data | 주파수 응답 재료·상태 대응 필요 | 예 |
| 85 | [m-drift](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#m-drift>) | [v-propagator](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#v-propagator>) | data | 미시 M_q 공급 | 예 |
| 86 | [m-diffusion](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#m-diffusion>) | [v-integral](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#v-integral>) | data | 미시 D 공급 | 예 |
| 87 | [m-drift](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#m-drift>) | [v-segments](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#v-segments>) | data | 국소 M_k | 예 |
| 88 | [m-diffusion](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#m-diffusion>) | [v-segments](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#v-segments>) | data | 국소 D_k | 예 |
| 89 | [m-reservoir](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#m-reservoir>) | [v-quadrature](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#v-quadrature>) | data | companion reservoir 상관 | 예 |
| 90 | [v-quadrature](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#v-quadrature>) | [f-loss](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#f-loss>) | data | 물리적 V_R,out | 예 |
| 91 | [d-constant](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#d-constant>) | [v-eigen](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#v-eigen>) | analogy | trace–discriminant 형태만 유사 | 아니오 |
| 92 | [f-states](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#f-states>) | [f-tmsv](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#f-tmsv>) | analogy | 입력 상태와 환경 비교 | 아니오 |
| 93 | [f-thinning](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#f-thinning>) | [f-loss](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#f-loss>) | analogy | 순수 손실의 count / covariance 표현 | 아니오 |
| 94 | [i-levels](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#i-levels>) | [i-drive](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#i-drive>) | data | 가중 dipole·원자종 | 아니오 |
| 95 | [i-frequency](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#i-frequency>) | [i-modes](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#i-modes>) | data | carrier 에너지·정규화 | 아니오 |
| 96 | [i-drive](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#i-drive>) | [m-reservoir](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#m-reservoir>) | data | 실제 구동 상태와 reservoir 모델 | 예 |
| 97 | [a-probe](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#a-probe>) | [m-drift](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#m-drift>) | data | 유한 시드 상태 주위의 요동 선형화 필요 | 예 |
| 98 | [a-probe](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#a-probe>) | [m-reservoir](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#m-reservoir>) | data | 실제 구동 원자 상관의 미시 연결 | 예 |

## 새 노드 26개

| 그룹 | ID·본문 | 역할 |
|---|---|---|
| structure | [q-objects](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-objects>) | 비교의 객체와 관측 계약 |
| structure | [q-information](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-information>) | 정확한 정보 회수의 준순서 |
| structure | [q-assumptions](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-assumptions>) | 가정의 함의 순서 |
| structure | [q-equivalence](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-equivalence>) | 표현의 정확한 동치 |
| structure | [q-quotient](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-quotient>) | 몫의 정의와 보존되는 합성 |
| structure | [q-numerics](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-numerics>) | 수치 세분과 오차 |
| structure | [q-evidence](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-evidence>) | 주장마다 다른 증거 |
| structure | [q-noise](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-noise>) | 잡음의 부분순서와 한계 |
| structure | [q-and](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-and>) | 결론은 전제들의 합류 |
| structure | [q-counterexamples](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-counterexamples>) | 축약과 단일 점수의 반례 |
| microscopic | [q-gksl](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-gksl>) | 명시적 원자 generator |
| microscopic | [q-diffusion](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-diffusion>) | Jump products의 ordered diffusion |
| microscopic | [q-qrt](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-qrt>) | 독립 QRT와 극한 검사 |
| microscopic | [q-rms](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-rms>) | Pump·weak dipole의 일관성 |
| microscopic | [q-periodic](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-periodic>) | Finite-seed periodic noise |
| microscopic | [q-spatial](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-spatial>) | Stationary spatial mode projection |
| microscopic | [q-gaussian](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-gaussian>) | Gaussian field와 photocurrent |
| microscopic | [q-uncertainty](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-uncertainty>) | 입력·오차·검출의 공동 조건 |
| transport | [t-characteristic](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#t-characteristic>) | 유한 atomic characteristic |
| transport | [t-boundary](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#t-boundary>) | 입사 상태와 경계 noise |
| transport | [t-poisson](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#t-poisson>) | Maxwell 유입과 Poisson 합 |
| transport | [t-solvers](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#t-solvers>) | 경로 해와 독립 참조 |
| transport | [t-thermal](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#t-thermal>) | 조건부 열린 Rb 기둥 |
| transport | [t-convergence](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#t-convergence>) | Ensemble 수렴의 별도 관문 |
| transport | [t-nonlocal](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#t-nonlocal>) | Moving-atom Maxwell의 조건 |
| transport | [t-prediction](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#t-prediction>) | 절대 예측과 실험 판정 |

## 새 연결 44개

| 출발 | 도착 | 종류 | 연결 의미 |
|---|---|---|---|
| [q-objects](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-objects>) | [q-information](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-information>) | logic | 공통 domain·관측·허용 회수 map |
| [q-objects](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-objects>) | [q-assumptions](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-assumptions>) | logic | 물리 가정의 해석 |
| [q-information](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-information>) | [q-quotient](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-quotient>) | logic | 양방향 정확 회수의 동치류 |
| [q-equivalence](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-equivalence>) | [q-quotient](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-quotient>) | logic | 표현 교체 시 관측과 합성 보존 |
| [q-assumptions](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-assumptions>) | [q-and](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-and>) | logic | 필요 조건의 동시 충족 |
| [q-numerics](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-numerics>) | [q-evidence](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-evidence>) | evidence | 고정 모델에서의 측정된 수치 오차 |
| [q-noise](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-noise>) | [q-and](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-and>) | logic | 같은 m·SQL에서만 잡음 비교 |
| [q-counterexamples](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-counterexamples>) | [q-information](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-information>) | evidence | gain 정보와 covariance 정보 구별 |
| [q-counterexamples](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-counterexamples>) | [q-and](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-and>) | logic | 모든 실제 입력선 보존 |
| [q-gksl](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-gksl>) | [q-diffusion](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-diffusion>) | data | 같은 state와 jump operators |
| [q-gksl](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-gksl>) | [q-qrt](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-qrt>) | data | 같은 ℒ의 독립 Liouville 경로 |
| [q-diffusion](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-diffusion>) | [q-qrt](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-qrt>) | validation | ordered spectrum 비교 대상 |
| [q-rms](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-rms>) | [q-gksl](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-gksl>) | data | reciprocal Hamiltonian·readout |
| [q-gksl](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-gksl>) | [q-periodic](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-periodic>) | data | periodic state·explicit bath |
| [q-diffusion](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-diffusion>) | [q-periodic](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-periodic>) | data | source별 harmonic noise |
| [q-periodic](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-periodic>) | [q-gaussian](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-gaussian>) | limit | 두-band field covariance와 밝은 carrier |
| [q-spatial](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-spatial>) | [q-uncertainty](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-uncertainty>) | data | 조건부 spatial prediction Jacobian |
| [q-gaussian](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-gaussian>) | [q-uncertainty](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-uncertainty>) | data | 같은 출력과 SQL의 sensitivity |
| [q-diffusion](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-diffusion>) | [m-diffusion](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#m-diffusion>) | limit | 원자 제거·mode 사상 조건 필요 |
| [q-periodic](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-periodic>) | [m-drift](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#m-drift>) | limit | two-band fixed-pump 조건부 drift |
| [q-spatial](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-spatial>) | [m-drift](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#m-drift>) | limit | stationary-center 조건부 drift |
| [q-gaussian](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-gaussian>) | [f-current](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#f-current>) | limit | bright approximation 또는 quadratic 보충 |
| [q-uncertainty](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-uncertainty>) | [o-claim](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#o-claim>) | evidence | 독립 입력의 적용 범위 |
| [q-qrt](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-qrt>) | [q-evidence](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-evidence>) | evidence | 동일 GKSL의 대수 비교 |
| [r-frame](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-frame>) | [q-equivalence](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-equivalence>) | evidence | minus pump-only gauge와 domain |
| [r-direct](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-direct>) | [r-poles](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-poles>) | equivalence | 대각화 가능한 동일 resolvent의 표현 |
| [v-phi](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#v-phi>) | [v-vanloan](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#v-vanloan>) | equivalence | 동일 상수 M,D,L에서 같은 W |
| [m-vacuum](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#m-vacuum>) | [q-noise](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-noise>) | evidence | 허용 noise의 한 선택 |
| [i-modes](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#i-modes>) | [q-objects](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-objects>) | logic | mode·normalization·observable contract |
| [q-gksl](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-gksl>) | [t-characteristic](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#t-characteristic>) | data | 위치·시간 의존 H와 같은 jumps |
| [t-characteristic](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#t-characteristic>) | [t-boundary](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#t-boundary>) | data | full finite-time evolution |
| [t-boundary](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#t-boundary>) | [t-poisson](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#t-poisson>) | data | connected packet와 mean pulse |
| [t-characteristic](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#t-characteristic>) | [t-solvers](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#t-solvers>) | validation | 동일 경로의 독립 계산 대상 |
| [t-thermal](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#t-thermal>) | [t-characteristic](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#t-characteristic>) | data | boundary·state·beam·signed k |
| [t-thermal](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#t-thermal>) | [t-poisson](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#t-poisson>) | data | Maxwell incoming measure |
| [t-poisson](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#t-poisson>) | [t-convergence](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#t-convergence>) | data | grid별 source-resolved spectra |
| [t-solvers](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#t-solvers>) | [t-convergence](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#t-convergence>) | evidence | 경로 수치 gate는 필요조건 |
| [q-numerics](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-numerics>) | [t-convergence](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#t-convergence>) | logic | consecutive grids·independent seeds |
| [t-boundary](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#t-boundary>) | [t-nonlocal](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#t-nonlocal>) | data | cross-position source/readout 정보 |
| [t-convergence](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#t-convergence>) | [t-prediction](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#t-prediction>) | evidence | 수렴한 thermal 자료 필요 |
| [t-nonlocal](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#t-nonlocal>) | [t-prediction](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#t-prediction>) | logic | field closure와 canonical covariance |
| [q-uncertainty](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-uncertainty>) | [t-prediction](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#t-prediction>) | logic | independent ledger·detector uncertainty |
| [q-and](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-and>) | [t-prediction](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#t-prediction>) | logic | 모든 전제의 동시 충족 |
| [t-convergence](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#t-convergence>) | [q-evidence](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-evidence>) | evidence | observed refinement/scramble 오류 |

## v2 내부 과학 개념 54개 대응

ID 표기 변경·노드 재분할 허용. 계산 내용·조건·반례의 현행 본문 위치 모두 명시. 앞서 적은 잘못된 구조 주장은 정정했으며 단순 복제하지 않음.

| v2 그룹·ID | 내용 | 현행 위치 |
|---|---|---|
| Bref · β̃₁ | 프레임 변환 | [r-frame](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-frame>) |
| Bref · β̃₄ | 16×16 영공간 직접 해 | [r-state](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-state>) |
| Bref · β̃ᵍ | 게이지 사상 → Floquet | [r-gauge](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-gauge>) |
| Bref · β̃⊥ | plus 브랜치 패리티 | [r-frame](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-frame>) · [r-gauge](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-gauge>) |
| Cref · γ̃₁ | 대각합-0 레졸번트 | [r-direct](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-direct>) |
| Cref · γ̃′ | Nambu 도함수 | [r-response](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-response>) |
| Cref · 𝒟̃ | 텐서 Gauss–Legendre | [r-geometry](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-geometry>) |
| Cref · σ̃ | 해석 rms 대조 | [r-geometry](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-geometry>) |
| B · β₁ | 프레임 변환 U₋ | [r-frame](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-frame>) |
| B · β₂ | 펌프 Hamiltonian | [a-assembly](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#a-assembly>) |
| B · β₃ | 산일자 + 열 재장전 | [a-assembly](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#a-assembly>) · [r-state](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-state>) |
| B · β₄ | 16×16 영공간 | [r-state](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-state>) |
| B · β₄′ | rate + Sylvester | [q-counterexamples](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-counterexamples>) |
| B · β₅ | Floquet 블록 (일반 N_F) | [a-order](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#a-order>) |
| B · β₆ | N_F 수렴 게이트 | [a-gate](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#a-gate>) |
| C · γ₁ | 대각합-0 레졸번트 | [r-direct](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-direct>) |
| C · γ₅ | 비공선 Doppler (1D 생산) | [r-geometry](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-geometry>) · [a-average](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#a-average>) |
| C · γ₂ | 쌍직교 고유기저 | [r-poles](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-poles>) |
| C · γ₃ | 단순극 합 (대각화 가능시) | [r-poles](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-poles>) |
| C · γ₃′ | Jordan · 등고선 | [r-defective](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-defective>) |
| C · γ₄ | 응답행렬 χ(Ω_SA) | [r-response](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-response>) |
| C · γ₆ | 속도 평균 ⟨·⟩_v | [r-average](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-average>) |
| C · γ₇ | Faddeeva · Voigt | [r-voigt](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-voigt>) |
| C · γ₈ | 속도평균 χ̄(Ω_SA) | [r-average](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#r-average>) |
| D · δ₀ | 정준 재정규화 | [d-maxwell](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#d-maxwell>) · [i-modes](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#i-modes>) |
| D · δ₁ | 기하 위상부정합 | [d-geometry](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#d-geometry>) |
| D · δ₂ | Maxwell 전기장 행렬 | [d-maxwell](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#d-maxwell>) |
| D · δ₃ | 대칭 회전틀 | [d-maxwell](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#d-maxwell>) |
| D · δ₄ | Cayley–Hamilton | [d-constant](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#d-constant>) |
| D · δ₅ | 두-이득 규약 | [d-gains](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#d-gains>) |
| D · δ₆ | 두 갭 · 광자속 항등식 | [x-flux](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#x-flux>) |
| E · ε₀ | Nambu 벡터 | [m-pair](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#m-pair>) |
| E · ε₁ | 양자 Langevin | [m-drift](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#m-drift>) |
| E · ε₃ | K 채널 분류기 | [m-k](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#m-k>) |
| E · ε₂ | 전파자 | [v-propagator](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#v-propagator>) |
| E · ε₄ | 미시 확산행렬 D | [m-diffusion](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#m-diffusion>) · [q-diffusion](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-diffusion>) |
| E · ε₄′ | 독립진공 완성 (테스트용) | [m-vacuum](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#m-vacuum>) |
| E · ε₅ | 미분 Lyapunov | [v-lyapunov](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#v-lyapunov>) |
| E · ε₆ | 좌표무관 닫힌 해 | [v-phi](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#v-phi>) |
| E · ε₆′ | 고유기저 trace–disc. | [v-eigen](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#v-eigen>) |
| E · ε₇ | Van Loan 블록지수 | [v-vanloan](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#v-vanloan>) |
| E · ε₈ | 공분산 전파 | [v-output](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#v-output>) |
| E · ε₈′ | 분절 누적 (z-의존) | [v-segments](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#v-segments>) |
| E · ε₉ | 4-모드 조립 | [v-quadrature](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#v-quadrature>) |
| X · ξ₀ | 컨쥬게이트 광자속 계수 | [x-flux](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#x-flux>) |
| X · ξ₁ | 정준 이상 IDS | [x-ideal](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#x-ideal>) |
| X · ξ_L | 레거시 스캔 점수 | [q-counterexamples](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#q-counterexamples>) |
| X · ξ₂ | 대칭 손실 | [x-ideal](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#x-ideal>) |
| X · ξ₃ | 비대칭 가중 차동 | [f-current](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#f-current>) · [f-sql](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#f-sql>) |
| X · ξ₄ | 유효성 게이트 | [x-canonical](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#x-canonical>) |
| F · φ₁ | 셀 이후 손실만 | [f-loss](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#f-loss>) |
| F · φ₂ | 광전류 선형화 | [f-current](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#f-current>) |
| F · φ₃ | SQL 정규화 | [f-sql](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#f-sql>) |
| F · φ₄ | 최종 관측량 | [f-measure](<../../FWM physics and analytic reconstruction/fwm_quotient_structure_v4.html#f-measure>) |

## 구조 비교 hash

검토한 이전 문서의 SHA256. 바이트 불변 검사는 루트의 `historical_artifacts_before.json`과 `verify_publications.py`가 별도 수행.

| 파일 | SHA256 |
|---|---|
| fwm_quotient_structure_v2.html | `3b66ecaec66a0177c08bba429173e5619fe9c992e93d55ec912c17c87c77484e` |
| fwm_quotient_structure_v3.html | `72fb7ba20f07fb0706c5fdf2c4aea70c2c4c397287155f55290b8c03991c2600` |
| fwm_quotient_structure_v3_builder.py | `c70a8947d1437c76d8ef9ab7b96d53a5d6297f5cdf7bbc6880ccabd81940445a` |
| fwm_quotient_structure_v3_view.js | `fe76918aa2482b22d7992458d1e4181d2e27d902aed6e79b1c30871c5713ae37` |
| fwm_quotient_structure_v3_view.css | `d475b3181913b9c1e4d091467c64bf8ded5565e306a73d12b070419fb8f96b51` |
| fwm_quotient_structure_v3_notes.md | `adcc2629a5a57d720ee7c430e583f322d117e6a0a59a14df2d43695081866ed6` |
