# 고정 열적 격자 실행·합산 계약

2026-09-16 작성, 09-17 확장. 기존 `thermal_campaign_v1` 유지. 수치 소스 108개·물리 입력·오차 기준 불변.
실행 제어만 `tools/`에 추가. 실행한 제어 파일도 원시 SHA-256별 별도 보존.

## 경로 배치

```powershell
python tools/thermal_campaign_batch.py --run docs/grand_challenge/thermal_campaign_v1 --grid 1 11 --output docs/grand_challenge/thermal_campaign_v1/batch-p1-s11.json --workers 4
```

- 하나의 process pool에서 격자 전체 처리. 긴 계산부터 제출; worker 최대 4, BLAS 각 1.
- 새 interpreter는 봉인 ZIP 사용. 제어 파일은 capsule 루트에 복사; 수치 소스 목록에 추가 안 됨.
- 기존 cache 먼저 검증. 같은 물리 경로·solver 요청만 재사용.
- 경로마다 CF4 3해상도 + 독립 adjoint 2해상도. 다섯 결과가 있어야 경로 판정 저장.
- 보고서·cache 덮어쓰기 금지. 기존 보고서는 현재 오차·경로 identity·5개 cache 참조와 대조.
- 실패 시 미제출 계산 중단. 실행 중 worker의 정상 cache 저장은 기다림.
- 전체 열적 수렴·광학 squeezing 플래그는 경로 성공으로 변경하지 않음.

중단 뒤 같은 명령 재개 가능. 이미 완료된 배치 보고서가 있으면 새 `--output` 필요.
재개 여부는 파일 존재만으로 판정하지 않음. 봉인 hash·원시 packet·실제 solver work 재검사.

## 반복 solve 생략의 수학적 근거

단일 원자의 미가중 해를

$$F = F(\mathbf r_{\rm in},\mathbf v,\tau; H,L,\rho_{\rm in},O,\Omega)$$

로 표기. Sobol 격자 번호·도착률은 이 운동방정식에 없음.
같은 물리 경로·모델·소스·solver 설정이면 같은 $F$.
중첩 격자의 동일 경로를 다시 풀 이유 없음. 수치 오차 증거와 합산 가중치는 현재 격자에서 재구성.

Face마다 $2^p$개 노드. 매핑은 `face * 2**new_power + offset`.
전역 index의 단순 배수 변환은 일반적으로 틀림.
P0→p1에서는 face 내부 offset이 0뿐이므로 `i → 2*i` 성립.

Pump-only 공통 입사 위상: $s=(1,-1,-1,1)$, $D(\phi)=\mathrm{diag}(e^{is_j\phi})$.
$D X D^\dagger$의 위상 harmonic은 $0,\pm2$뿐.
따라서 $\phi=0,\pi/2,\pi,3\pi/2$ 네 점 평균은 균일 위상 적분과 정확히 같음.
합산 코드의 charge mask와 별개로 직접 네 위상 변환을 계산해 대조.
위상별 원자 solve 추가 불필요. **Finite seed에는 적용 불가.**

Poisson number 항은 $\langle m m^\dagger\rangle_\phi$.
$\langle m\rangle_\phi\langle m\rangle_\phi^\dagger=0$으로 대체 금지.
Density는 경계 도착률에 한 번 포함. $nV$, occupancy, 가중치 합으로 재정규화 금지.
위 근거를 코드 주석에도 명시. 후속 에이전트의 불필요한 반복 solve 복원 방지.

## 격자 합산

```powershell
python tools/thermal_campaign_grid.py --run docs/grand_challenge/thermal_campaign_v1 --grid 0 11 --output docs/grand_challenge/thermal_campaign_v1/grid-p0-s11.json
python tools/thermal_campaign_grid.py --run docs/grand_challenge/thermal_campaign_v1 --grid 1 11 --output docs/grand_challenge/thermal_campaign_v1/grid-p1-s11.json
```

Cache 읽기 전용. 미계산 경로를 자동 생성하거나 누락 가중치를 보정하지 않음.
경로 실패·누락 또는 독립 가중합 불일치: **격자 전체 spectrum 보류**.
모든 RF·named source·두 ordering·복소 response·Poisson 항을 별도 비교.
독립 합산 상대 오차 기준 $10^{-10}$; 원자 경로 오차 기준과 별도.

한 격자 내부 수치 통과와 열적 앙상블 수렴은 다른 판정.
기존 ensemble gate: 최소 3격자 × 3독립 scramble, 두 refinement와 scramble 비교, 기준 5%.
P0·p1의 seed11 결과만으로 gate 통과 불가. 미검증 조합의 성공을 만들어 넣지 않음.

## 중첩 격자 비교

```powershell
python tools/thermal_campaign_compare.py --run docs/grand_challenge/thermal_campaign_v1 --coarse docs/grand_challenge/thermal_campaign_v1/grid-p0-s11.json --fine docs/grand_challenge/thermal_campaign_v1/grid-p1-s11.json --output docs/grand_challenge/thermal_campaign_v1/comparison-p0-p1-s11.json
```

입력 보고서 seal만 신뢰하지 않음. 두 격자의 실제 cache·경로 gate·합산을 다시 읽어 대조.
공통 6경로 × 5해상도의 8개 raw metric bitwise 동일성 확인.
현재 source 이름·packet digest 재결합, 도착률 반감도 검사.
비교 제어 파일과 합산 helper의 실제 byte를 함께 보존.

격자 합산의 선형성으로

$$S_{p+1}=\tfrac12 S_p + \sum_{i\in\mathrm{new}}\lambda_i^{(p+1)}F_i$$

성립. 두 ordering·source별 항·Poisson number·복소 response 모두 직접 대조.
이 식의 수치 통과는 합산·재사용의 일관성 검사. 격자 수렴 인증과 별도.
실제 p→p+1 변화는 기존 SI floor의 RF/source별 상대 Frobenius 오차로 비교.
CLI 성공은 감사 계산 성공. `refinement_diagnostic.passed`와 `convergence_gate.passed` 별도 확인.

## 여러 격자의 선언 전체 감사

```powershell
python tools/thermal_campaign_ensemble.py --run docs/grand_challenge/thermal_campaign_v1 --grids docs/grand_challenge/thermal_campaign_v1/grid-p0-s11.json docs/grand_challenge/thermal_campaign_v1/grid-p1-s11.json docs/grand_challenge/thermal_campaign_v1/grid-p2-s11.json --output docs/grand_challenge/thermal_campaign_v1/ensemble-seed11-p0-p2.json
```

명시한 입력마다 원본 cache·경로 gate·독립 가중합 재검증. 새 원자 solve 0회.
캠페인 선언의 powers·seeds 유지. 입력한 seed만으로 필수 목록을 줄이지 않음.
다른 캠페인·중복 격자·변조 보고서·누락 cache 거부. 입력 보고서·제어 파일도 종료 시 재검사.
실행 제어와 두 helper 원시 byte 보존. 실제 import한 수치 모듈의 capsule 위치·hash 기록.

- `audit_passed`: **제출한** 모든 격자의 원본 증거·합산 감사 통과.
- `coverage`: 선언한 전체 격자, 제출한 격자, 누락 격자 명시.
- `pairwise_diagnostics`: 현재 있는 연속 refinement·양방향 scramble 비교. 자료 부족이어도 실제 실패 보존.
- `convergence_gate`: 기존 전체 수렴 함수의 실제 판정. 수식·5% 기준 불변.
- `thermal_ensemble_converged`: 선언 전체 자료와 기존 gate가 모두 통과해야 true.
- `physical_optical_prediction`: 이 원자 감사에서는 false. Optical Maxwell·SQL·실험 검증 별도.

격자 상대 변화는 coarse/reference에 대한 차이. 연속 열적 적분의 참오차 상계 아님.
경로 수치 통과 후 격자 변화가 커도 즉시 누락 physics로 결론 내리지 않음.
RF별 그림은 각 주파수에서 source 최대값 유지; source 평균으로 실패를 숨기지 않음.
독립 seed는 같은 열적 적분의 quadrature 검사. 다른 Δ·δ·T·pump power 실험의 out-of-sample 검증 아님.

## 뒤이을 격자

V1 선언: 9격자, 126경로 항목, 72고유 경로·360고유 계산 요청.
Seed11 p2의 24경로 + seed211 p1의 12경로, 총36고유 경로·180계산 검증 완료.
V1의 계획상 미수행180요청과 기존 실패 기록 유지. 후속 계산은 확장 v2에서 수행 가능.

- V2 p2/seed11: 기존24경로·120계산 이전 후 새 경로 gate·합산 재검증 완료.
- V2 p2/seed211: 24경로·120계산 검증 완료. 기존60개 재사용·새60개 계산.
- V2 p2/seed811: 24경로·새120계산 검증 완료. 세 p2 격자의 여섯 방향 공동 감사 완료.
- V2 p3/seed11: 48/48경로·합산 감사 완료. 공통24경로·120계산 재사용, 새24경로·120계산.
- 남은 p3/seed211·811: 신규48경로·240계산. 이후 p4 신규144경로·720계산. 미수행 합계960계산.
- V2 전체1440요청 중180원본 이전·새300계산, 96고유 경로·480계산 확보. 복사본을 새 원자 계산으로 세지 않음.

2026-09-28 현재 [4격자 공동 감사](thermal_campaign_v2/ensemble-p2-three-seeds-p3-s11.json) 통과, 5격자 누락.
P2→p3/seed11의 5% refinement는 실패. 고정 v2 `[2,3,4]`는 이 필수 edge를 포함하므로
남은 격자를 채워도 현 선언 그대로 전체 gate 통과 불가. 기존 선언·실패 보존.
다음은 p3의 남은 독립 seed 감사와 별도 더 미세한 창(예: `[3,4,5]`)의 설계·선언.
새 창도 두 연속 refinement·독립 scramble의 실제 검증 필요; 통과 보장 없음.

Occupancy가 우연히 $nV$에 가까운 격자만 선택해 종료하지 않음.
기존 선언의 occupancy/$nV$: p0/seed11 0.617675, p1 0.945434, p2 0.830227.
이 값은 경계 기하·속도 quadrature 진단. 원자 spectrum 수렴 또는 실험 일치 증거 아님.

## 더 미세한 격자 선언·정확한 cache 이전

```powershell
python tools/thermal_campaign_extend.py --run docs/grand_challenge/thermal_campaign_v1 --target NEW_CAMPAIGN --powers 2 3 4
```

`NEW_CAMPAIGN`은 아직 없는 별도 경로. 상위 폴더는 존재해야 함.
기존 선언·실패 보고서 불변. 새 선언도 동일한 세 seed·모델·오차 기준 유지.
P2·p3·p4: 9격자, 504경로 항목, 288고유 물리 경로·1440고유 계산 요청.
더 미세한 두 refinement를 검사할 공간 확보. 통과 보장 아님.

- 기존 ZIP의 **실제 코드**를 새 interpreter에서 실행해 선언 생성. 현재 앱 코드 혼입 거부.
- 수치 소스 108개·원시 ZIP byte·모델·환경·solver 설정·예산 동일성 검사.
- 부모와 새 선언의 모든 경로·도착률·occupancy 재구성. Cache 없는 격자도 검사.
- 시작 시 존재하는 cache 목록·원시 hash 고정. 이후 추가된 정상 결과는 이번 이전 대상 밖.
- 동일 ODE key만 원본 byte 복사. 반복 solve·hash 수정·기존 path/grid 통과 판정 이전 없음.
- 새 선언에 묶인 두 번째 source capsule에서 `provider=None`으로 다시 검증.
- `transfer.json`에 부모 선언·선택 cache·복사 원본 hash·제어 파일·누락 요청 기록.
- 임시 폴더의 전체 검증 후 독점 rename으로 공개. 대상 충돌·실패 시 새 캠페인 공개 안 함.

새 계획의 `historical_numerical_results_reused`는 실제 import 여부 기록.
이전 성공은 raw 계산 동일성 검사. 경로 gate·현재 가중합·열적 수렴은 새 캠페인에서 별도 재구성.
`transfer_passed=true`만으로 `thermal_ensemble_converged`나 optical 판정 변경 금지.
Source/model/path/spec/environment가 같은 미가중 방정식이라는 재사용 근거를 코드 주석에 명시.

## 선언이 다른 중첩 격자 비교

```powershell
python tools/thermal_campaign_cross_compare.py --parent-run PARENT --run CHILD --coarse PARENT_P1.json --fine CHILD_P2.json --output NEW_JSON
```

직접 확장한 부모·자식 캠페인, 동일 seed·연속 power만 허용.
V1 p1→v2 p2 진단 가능. 부모 p1을 자식 선언에 추가하거나 과거 보고서의 선언 hash를 바꾸지 않음.

- 부모·자식 각각 자기 source ZIP·선언 marker·새 interpreter에서 원본 경로 gate·합산 검증.
- 확장 계획·transfer·부모 snapshot·검증 기록의 seal 및 실제 cache byte 대조.
- Source ZIP·108개 소스·모델·환경·solver·예산·seed 동일성 확인. 양쪽 저장소 각각 읽기.
- 같은 물리 경로의 5계산 key·record SHA·payload digest·8metric 원본 동일성 검사.
- Face 안 Sobol offset 유지, 현재 source/digest 재결합, 도착률 절반 및 분할합 검증.
- `S_fine = 0.5*S_coarse + 신규 경로 합`: 기존10⁻¹⁰ 기준. Density 한 번, Poisson mean-outer 유지.
- 동일 미가중 ODE의 수학적 동일성으로 반복 solve 생략. 근거 주석 유지; provider 전달 금지.
- 읽은 cache·보고서·선언·ZIP·소스·제어 파일 마지막 byte 재검사. 무관한 새 cache 추가는 허용.
- 원본 캠페인의 제어 파일도 수정하지 않음. 비교에 사용한 제어4파일 원본은 결과에 base64로 보존.
- 과거 절대 위치는 기록용. 파일 이동 후에도 부모 계획 raw hash·ZIP·record 내용으로 계보 검증.
- 감사 실패 결과 보존·종료 코드 실패. 5% refinement 진단 실패와 감사 실패는 별개.

비교 결과는 한 refinement edge 진단. V2의 필수 p2→p3·p3→p4 및 독립 seed gate 대체 불가.
출력은 독점 생성. 전체 thermal·optical 판정 false 유지.

그림 명령:

```powershell
python -m docs.grand_challenge.plot_thermal_grid PARENT_P1.json CHILD_P2.json NEW.png --cross-audit NEW_JSON
```

두 grid의 record SHA·raw file hash·선언·selection과 통과한 감사 연결 필수.
그림에서 재계산한 6metric 최대값도 감사 원본과 일치해야 출력. RF별 source 최대값 유지.
새 과학적 인증 기능 아님. 자체계산·검증은 위 CLI 담당.

## 적용 범위

조건부 열린 정사각 기둥·373 K·연속 Gaussian pump·reduced Rb·pump-only 원자 응답.
RF 0.1 / 1 / 4 MHz. 입력 구조의 seed power는 이 pump-only 상태 방정식의 유한 seed 구동 아님.
열적 전체 수렴·finite seed·nonlocal Maxwell·full atom·검출 SQL·절대 실험 squeezing: 별도 후속.
Fitted Fast/Balanced gain·excess-noise 계수를 이번 미시적 원자 계산에 도입하지 않음.

## 수렴 변화의 층별 분해 — 2026-09-28

```powershell
python tools/thermal_campaign_diagnose.py --run docs/grand_challenge/thermal_campaign_v2 --reference docs/grand_challenge/thermal_campaign_v2/grid-p2-s11.json --candidate docs/grand_challenge/thermal_campaign_v2/grid-p2-s211.json --output NEW_DIAGNOSTIC.json
```

같은 seed의 연속 격자 또는 같은 power의 다른 seed 두 개만 비교. 출력은 독점 생성.
봉인 source capsule에서 두 격자의 모든 원시 packet·경로 gate·네 위상 합산 재검증.
Provider 전달 없음; 새 원자 solve **0회**. 누락 자료를 계산·보간·재가중하지 않음.

입사면, 체류시간, 종속도/σ, 횡속도/σ의 네 고정 분할. 각 분할은 모든 경로를 한 번씩 포함.
서로 다른 분할은 같은 경로를 다시 보는 것이므로 분할 종류 사이의 기여를 합산하지 않음.
실제 도착률·Poisson mean outer·모든 복소 RF/source 행렬 유지. 각 집단 합이 원본 격자와
두 격자의 차이를 재현하는지 기존 10⁻¹⁰ 기준으로 검사.

한 RF/source에서 Δ=S_candidate−S_reference, 집단별 D_k의 합이 Δ라 할 때,

$$a_k=\frac{\operatorname{Re}\operatorname{Tr}(\Delta^\dagger D_k)}{\|\Delta\|_F\max(\|S_{reference}\|_F,128\epsilon\,s_{candidate})}.$$

Σa_k는 해당 방향의 상대 변화. Δ=0이면 a_k=0으로 정의. 집단 norm은 일반적으로 더해지지 않음.
음의 a_k는 상쇄 기여. 참 적분 오차의 상계나 독립적인 양의 noise 비율이 아님.
각 metric의 최악 RF/source에서 순위를 보고하되, 전체 RF/source 배열도 결과에 보존.

최초 native 검증은 유지. 검증 전에 잡은 원시 byte의 봉인 hash를 실제 반환된 native record와
연결하고 읽기 직후 같은 byte인지 검사. 후속 읽기와 최종 저장 직전에도 원시 SHA-256 대조.
동일한 byte에 대한 마지막 native 재검사를 이 확인으로 대체: 경로당 native 읽기 12→7회.
Source·실행 제어·보존한 제어 파일·계획·입력 보고서 최종 검사 유지. 실제 가속률 미측정.
파일 변조·누락·이름 변경·검증 중 A→B→A 교체 거부 검사 포함. 기존 기록의 제어 파일 교체 없음.

`audit_passed`는 증거·분해 산술 검증. `directed_diagnostic.passed`는 선택한 한 방향의 5% 판정.
전체 thermal·optical 인증은 false 유지. [실제 p2·p2→p3 감사·독립 검산·그림](progress_2026_09_28/README.md).

## 후속 수치 개선 후보 — 검토만 완료

2026-09-17 독립 코드 검토. 현재 고정 캠페인 수정 없음. 후보의 가속률 미측정.

현재 Gaussian 경로는 $H(t)=H_0+f(t)H_1$. 광학 carrier를 그대로 시간 적분하는 구조 아님.
Ballistic Doppler는 경로마다 상수. Pump frame·adjoint demodulation·RF/source 공통 계산은 이미 적용.
다른 회전 좌표계를 추가한다고 가속이 보장되지는 않음.
CF4의 시간 순서 오차와 atomic commutator 보존 검사는 별개. 보존 검사만으로 step 수렴 인증 금지.

1. **경로별 refinement 사다리 별도 선언.** 겹치는 계산 재사용, 실패 경로만 세분화.
   두 연속 edge·독립 비교·독립 refinement, 8metric·RF/source·SI floor·기존 오차 기준 유지.
   Fine trajectory를 subsampling한 값을 독립 coarse solve로 세지 않음. 새 evidence 계약 필요.
2. **상수 연산 사전 계산.** 고정 port Gaussian 모델의 $L_0,L_1$·drift·readout/drive 대수 검토.
   실제 운동 위상·Gaussian envelope 유지. Conditioning 검사·fallback 보존.
   항등식, 기존 구현의 모든 packet metric, 독립 adjoint 대조 후 채택. Density eigensystem 공유는 이미 구현됨.

소스·수치 계획 변경 시 별도 선언과 출처 검증. 기존 기록의 hash·실패 판정 교체 금지.
Production Ultra의 정상상태/Floquet 항등식을 유한 Gaussian 경로에 자동 적용하지 않음.
근거: `smooth_transport.py`, `reference/exponential_transport.py`, `reference/adjoint_transport.py`.

## 코드 검증

- 새 배치·합산·중첩 비교 검사 93개 통과. Synthetic fixture만 사용; 원자 예측 증거와 구분.
- 실제 별도 interpreter 실행, 누락/변조/재봉인한 허위 보고서, source 재결합, 모든 metric 변경 거부 포함.
- 독립 검토에서 대기 작업 취소·기존 보고서 cache 참조·제어 파일 hash/content race 보완.
- 2026-09-16 `python -m pytest -q`: **1816 passed, 3 skipped, 1 failed**, 422.09초.
- 유일 실패: 기존 삭제된 `FWM_physics.tex` 참조. Skip 3건: Windows symlink 권한. 기존 상태 유지.

2026-09-17 추가:

- 선언 전체 감사 전용 37검사 통과. 전체 3×3 synthetic 통과·자료 부족·실제 refinement 실패 분리.
- 별도 프로세스가 실제 cache 검증까지 도달했는지 stderr 조건 보강. 이후 37검사 재통과.
- 실제 p0·p1 원본 감사 통과: [새 검사기 실행 기록](thermal_campaign_v1/ensemble-seed11-p0-p1.json).
  새 solve 0회, 누락 격자 7개, 수렴 false. 고정 수치 모듈 62개 import 출처 기록.
- 실제 p0·p1·p2 원본 감사 통과: [세 격자 기록](thermal_campaign_v1/ensemble-seed11-p0-p2.json).
  새 solve 0회, 선언 9격자 중 3검증·6누락. 실제 두 refinement 진단 모두 5% 초과; 전체 수렴 false.
- 전체 `python -m pytest -q`: **1853 passed / 3 skipped / 기존 문서 누락 1 failed**, 486.90초.
- RF별 그림의 최대값은 기존 봉인 비교 원본과 동일. 0.1·1·4 MHz별 source 최대값 보존.

2026-09-17 격자 확장 도구 추가:

- 전용 **31검사 통과**, 352.31초. 독립 코드 검토에서 수정 필요 결함 없음.
- Byte 보존·현재 선언 재결합·추가 cache 허용·변조/누락 거부·복사 실패·대상 충돌 검사.
- 실제 별도 interpreter: synthetic 원자 증거 거부, 빈 캠페인 이전 성공, 잘못된 target marker 거부.
- 비어 있지 않은 unit-test 이전은 synthetic fixture. 실제 원자 cache 이전 검증과 구분.
- 전체 `python -m pytest -q`: **1884 passed / 3 skipped / 기존 문서 누락 1 failed**, 868.32초.
- 독점 폴더 공개는 Windows/Linux 지원. 다른 플랫폼에서는 지원 없이 덮어쓰기하지 않고 중단.

2026-09-17 실제 seed211·확장 이전 완료:

- [결과·원본·두 그림](thermal_campaign_seed211_v1.md). Seed211 12/12경로·새60계산 통과.
- P0/p1 직접 합산 최대 2.136374517e-16; p0→p1 여섯 metric 최대값은 모두5% 초과.
- 5/9격자 원본 감사 통과. 실제refinement3개·양방향scramble4개 모두5% 기준 초과.
- [V2 실제 이전](thermal_campaign_v2/transfer.json): 동일ZIP·180raw record byte 보존, 새solve0회.
- [새 p2/seed11](thermal_campaign_v2/grid-p2-s11.json): 24경로·120hit·0miss, 경로 gate·합산 재구성.
  기존8개 spectrum의 canonical JSON byte 동일. Source·모델·spec·환경 계약 유지.
- DriveFS의 링크 수0 지원 보완. 2이상·symlink·비정규 파일 거부와 독립복사/byte/native 검증 유지.
  관련6회귀검사 통과. 확장 도구 전용37검사 포함해 전체검사 재실행.
- 전체 `python -m pytest -q`: **1890 passed / 3 skipped / 1 failed**, 858.11초.
  유일실패: 기존 삭제 FWM_physics.tex. Skip3건: Windows symlink 권한.

2026-09-18 p2/seed211·선언 간 비교 완료:

- [결과·원본·그림](thermal_campaign_p2_seed211_v2.md). 24/24경로·새60계산 통과.
- 독립 가중합 최대3.136425967e-16. 부모p1→자식p2 공통60계산 재사용·분할합 잔차4.714798084e-16.
- V2 두 p2 원본 공동 감사 통과. 2/9격자·7누락. 전체 thermal·optical false.
- 선언 간 비교 전용52검사 통과. 실제 두 capsule·양쪽 native cache 검증도 통과.
- 독립 검토의 실패 종료 코드 보완, byte 불변인 캠페인 이동 지원 및 회귀 검사.
- 전체 `python -m pytest -q`: **1942 passed / 3 skipped / 1 failed, 1079.41초**. 기존 FWM_physics.tex 누락 실패·Windows symlink skip3건.


## 2026-09-18 — p2 세 독립 seed 완료

[결과·원본·그림](thermal_campaign_p2_seed811_v2.md). Seed811 24/24경로·새120계산 통과.
독립 합산 잔차 4.528539686e-16. 세 p2 격자 원본 공동 감사 통과, 새 solve0회.
독립 seed 여섯 방향 0/6통과, metric 최대값 범위5.70–43.61%.
3/9격자·6누락, 검증된 v2 refinement0개. 전체 thermal·optical false.
V2 72고유 경로·360계산 확보. 다음 p3/seed11 새120계산; p3 전체360·p4 전체720계산 남음.
수치 source·모델·허용오차·도착률 불변. 기존 증거·실패 판정 보존.
병렬 그림 배치 수정·독립 행렬 산술 대조 완료.
전체 `python -m pytest -q`: **1942 passed / 3 skipped / 1 failed, 946.00초**. 기존 문서 누락1실패·Windows symlink skip3건.


## 2026-09-28 — p3/seed11·4격자 공동 감사 완료

[결과·원본·그림](progress_2026_09_28/README.md),
[배치](thermal_campaign_v2/batch-p3-s11.json), [격자](thermal_campaign_v2/grid-p3-s11.json).
같은 고정 수치 소스108개·모델·환경·오차 기준. P3/seed11 **48/48경로 통과**.
기존120계산 재사용·새120계산, 배치10,730.423초. 격자 직접 합산 잔차3.69606571×10⁻¹⁶.
[경로 기록 별도 집계](progress_2026_09_28/batch-path-union-verification.json):
누적96고유 경로·480계산, 선언 잔여960계산. 저장 기록 집계와 native 재감사 구분.
96경로·8metric 최대값: primary 두 edge 5.46665188×10⁻⁶ / 2.18536636×10⁻⁷,
독립 비교8.40851049×10⁻⁸, 독립 refinement1.87343815×10⁻⁷. 기존 경로 기준 모두 통과.

[정식 중첩 비교](thermal_campaign_v2/comparison-p2-p3-s11.json):
공통24경로·120원시 계산 동일성, source/digest·반감 도착률 재결합 확인.
분할합 `S_p3 = 0.5*S_p2 + 신규 경로 합` 최대 잔차5.299764275×10⁻¹⁶.
새 solve0회. **산술 감사 통과, 5% refinement 실패**.
Greater10.514126%, lesser7.847664%, source별13.788644%/12.330948%,
Poisson number3.677599%, response36.169445%. 6metric 중5개 기준 초과.

[4격자 공동 감사](thermal_campaign_v2/ensemble-p2-three-seeds-p3-s11.json) 통과: 4/9격자·5누락.
누락은 p3/seed211·811 및 p4/seed11·211·811. 실제 평가한 refinement1개 중0통과,
p2 독립 seed 여섯 방향도0/6통과. 전체 thermal·optical false 유지.
P2→p3 실패를 포함하는 현 v2 창은 남은 격자만 채워 통과시킬 수 없음.
더 미세한 별도 창 설계와 남은 p3 seed 감사로 후속 검증 범위 확보. 기존 실패 소급 변경 금지.

[층별 진단](progress_2026_09_28/diagnostic-p2-p3-s11.json)과
[별도 검산](progress_2026_09_28/diagnostic-p2-p3-validation.json) 일치.
입사면·체류시간·속도별 모든 RF/source·복소 행렬 보존. 합산·차이 복원 최대9.009749614×10⁻¹⁶.
Signed 기여의 상쇄 유지; 분할 종류 사이 합산 금지. 참 적분 오차 상계·물리적 원인 인증 아님.
기존 p2 진단은 당시 제어 파일로 보존. 최종 byte guard 실행 증거는 새 p2→p3 진단에 귀속.

진단 전용45검사·독립 검토 통과. 삭제된 이론 파일의 고정 참조를 현행 manifest 검사로 교정.
최종 `python -m pytest -q`: **2077 passed / 3 skipped / 0 failed, 1156.52초**.
[코드 hash·검증 요약](progress_2026_09_28/test_summary.json). 동시 배치 실행 시간을 가속률로 비교하지 않음.
경계 이력은 설계 단계. 전체 thermal·nonlocal Maxwell·절대 squeezing·실험 holdout 미완료,
공식 milestone **0/4** 유지.
