# ⟦AI⟧ 새 통합문서 독립 감사

원본 데이터 검토. 원본 문서의 메모는 자료로만 취급. 계산·코드 변경 지시로 채택하지 않음.
상태: 원본 읽기 완료. 피팅·원본 수정·재계산·xlsx 출력 없음.

원본: `D:/20261013_TPD_vs_Gain 및 스퀴징 특성 추정.xlsx`.
SHA-256 `2653fe73abca75550b2b8665da8c8fff16f5c28f155e029f1e6da155a9409ed0`, 47,107 bytes.
읽기 전후 hash 일치. Artifact Tool로 2개 sheet·표·기호·수식·cached 값·Solver 이름 확인,
source 영역 2개 render 검수. 원본 저장 없음. `cells.json`은 보조 추출물.
importer가 shared-formula 자식 수식을 일부 노출하지 않음. 최종 수식 근거는 root의 OOXML 해석과 대조 필요.

## 이득·축·운전 조건

| 원본 위치 | 확인 내용 | 해석 범위 |
|---|---|---|
| Sheet1!A1:B1 | 초기 Seed 35 mV | 전압 분모. 이 workbook 자체에 pump-off·셀 앞뒤 기준면 설명 없음 |
| Sheet1!D1 | `=B1*20/228.5` → 3.0634573304 µW | 표시 전압의 파워 환산 |
| Sheet1!F5, G5 | `=C5/$B$1`, `=D5/$B$1` | detector 전압 비. 원자 입·출력 photon-flux gain을 직접 관측한 값 아님 |
| Sheet1!B2 | `=76*3.879` → 294.804 mW | 원본의 펌프 파워 환산 |
| Sheet1!C2:D2 | OPD `~1.0` GHz | reference transition·blue/red 정의는 기존 사용자 확인과 별도 연결 필요 |
| Sheet1!B5 | `=A5-3036` → TPD +8 MHz | workbook 축은 `f_EOM−3036 MHz` |
| Sheet1!A15:B16, F15:F16 | 3034/3033 MHz, workbook TPD −2/−3 MHz에서 probe gain 모두 24.85714286 | 1 MHz 간격 plateau. sub-MHz peak를 workbook만으로 특정 불가 |
| Sheet1!E1:F2 | in-cell seed 0.8046, after-cell 0.9709, 모두 `10/06` 표기 | 새 스캔과 동시 측정인지 미기재 |
| Sheet1!H1 | `=0.58*1239.8/795` → QE 0.9045081761 | responsivity·파장 환산 |
| Sheet2!B7:X9 | seed 효율 0.8046, after-cell 0.9709, QE 0.9047 | Sheet1 QE와 0.000191824 차이. 외부 η는 각각 0.8781869882 / 0.8783732300 |
| Sheet2!B15:B16 | `B7*B4/B6`, `B7*B5/B6` | 측정 gain에 과거 in-cell 효율을 곱한 계산값 |

현행 GABES의 `δ_GABES=ν_HF−f_EOM`를 사용하면
`δ_GABES [MHz]=−TPD_workbook−0.267561` (`ν_HF=3035.732439 MHz`).
원본 축을 그대로 solver 축으로 쓰면 부호와 원점이 함께 어긋남.
sheet 전역의 숫자·라벨에서 온도·waist·cell length·각도 기록 확인되지 않음.
과거 handoff 조건을 이번 workbook의 명시 입력으로 승격하지 않음.

## 스퀴징의 측정·추산 구분

- Sheet1!B30: EOM 3.038 GHz에서 FWM와 SQL 차이 **3 dB 직접 측정**.
  SQL이 optical power에 비례한다는 성질로 나머지 degree 추산했다고 명시.
- Sheet1!B31: probe·conjugate 전압, SA에서 읽은 최소 FWM noise, seed·pump·OPD는 측정.
  나머지는 환산이라고 명시. 매 detuning 동일 RF·RBW·VBW인지 workbook에 없음.
- Sheet1!H11=−79 dBm, I11=−76 dBm. Sheet1!I5는
  `=$I$11+10*LOG10(E5/$E$11)`, K5는 `=H5-I5`.
  K5:K23은 공통 SQL anchor에서 추산한 차이. 독립 SQL 실측 19개로 취급 불가.
- Sheet1!B34: seed drift 기록. B35: noise fluctuation에 의한 약 ±0.5 dB 불확도.
  표준편차·포함인자·반복 수·파워 교정 불확도·공통 anchor 공분산 미기재.
  모든 점의 독립 Gaussian σ=0.5 dB로 확정할 근거 없음.
- Sheet1!J5 `=LOG10((F5-G5)^2/(F5+G5))`.
  J열 이름 `Squeezing Factor`는 선형 NF도 `10log₁₀ NF`도 아님.

## Sheet2 계산·피팅 문제

- B10:X11의 S,T는 수치 입력. Solver `solver_adj=Sheet2!$C$10:$X$11`,
  `solver_opt=Sheet2!$Y$17`. Y17=`SUM(C17:X17)`.
  S,T는 gain 역산의 fitting 변수. 독립 원자 계산·측정으로 취급 불가.
- B17=`(B13^2/B15-1)^2+(B14^2/B16-1)^2`는 이득 상대제곱오차.
  V10:X10은 S=10⁻⁵ 경계, V17:X17은 각 약 1. Y17≈3.
  마지막 3점의 두 arm 동시 재현 실패가 남아 있음.
- B18은 `Jasperse Eq. (16) NF`라는 계산식. `η=B8*B9`,
  `ξ=sqrt(S²+(ln T/4)²)`, `χ=atanh(ln T/(4ξ))` 사용.
  마지막 항에 `ln(T)² sinh(ξ)⁴` 포함. 논문 원식과 독립 대조 전 검증된 NF로 채택 불가.
- B23:T23은 Sheet1에서 추산한 dB 수치, B24는 `10^(B23/10)`.
  B25는 이 값에서 모델 NF를 뺀 ΔF. 새 측정 데이터가 추가된 영역 아님.
- B27과 W27은 `ΔF/(G−1)` 또는 그 회귀계수. 현행 seed RIN a `[Hz⁻¹]`와 기호·단위 다름.
  W27≈0.200740522를 새 a 입력에 넣을 근거 없음.
- **B28=`B17+$W$27*B26`**. 행 이름은 fitted F_experiment지만,
  NF baseline B18 대신 이득 objective B17을 더함. B28:T29의 fitted NF/dB에 잘못된 기준항 사용.
- B30=`SLOPE(B25:X25,B26:X26)`≈−0.0113549,
  D30=`INTERCEPT(...)`≈3.624689. 이 회귀의 a,b는 현행 PSD a,b와 다름.
- B33은 `η(α₁²−α₂²)²/(α₁²+α₂²)`, B34는 ΔF/B33.
  `F_seed−1`은 모델·전달계수 가정에 따른 역산. seed noise 직접 관측 아님.

독립 감사 권고: gain 교정에는 원본 전압·EOM 주파수·분모와 이전 사용자 기준면 확인을 연결.
스퀴징 교정은 직접 SQL anchor, 추산 degree, 모델·역산·회귀 결과를 구분.
이 문서는 새 fit·production 보정계수 채택을 승인하지 않음.
