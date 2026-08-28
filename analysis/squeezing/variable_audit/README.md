# FWM 동작 변수의 gain-referred 진단 감사

실험자가 현장에서 체감한 변수별 감각을 세 소스와 대조한다. 엔진에서 읽는
`gain_referred_noise_dB`는 mean-field gain에 ideal twin-beam law를 결합한
**gain-referred compatibility diagnostic**이다.

> **Claim gate — `MEAN_FIELD_DIAGNOSTIC`:** GABES에는 주파수 의존 microscopic
> Langevin diffusion/covariance와 동일 조건의 측정 SQL이 없다. 따라서 이 감사는
> 물리적 squeezing, squeezing bandwidth, 또는 above/below-SQL을 예측하지 않는다.
> `S_dB`와 `xi_*`는 기존 소비자를 위한 deprecated compatibility aliases다.

- 엔진 — GABES hardened Ultra (η = 0.8694, δ 격자 5 MHz, gap gate [0.5, 1.5])
- 이론 — `docs/squeezing_report/squeezing_report_v6.tex` (특히 §tolerance)
- 문헌 — `references/fwm_squeezing_paper_parameters.csv` (10편)

```bash
python analysis/squeezing/variable_audit/scan_variable_audit.py
```

약 20 초. `generated/variable_audit.md` (한글 보고서) + `.json` (원시 수치).

## 대상 변수와 판정

| 변수 | 판정 | 핵심 수치 |
|---|---|---|
| TPD | 진단 추세 지지 | δ* ∝ Ω^1.5, 파워 ×2 → δ* +24 MHz vs gain-ref. 진단 창 폭 13.8 MHz |
| 펌프 파워 | 맞음 | 겹침이 파워보다 8.8배 강한 지렛대 (= 1 / (dln q/dln P)) |
| 시드 파워 | 진단상 둔감 | 1–200 µW에서 gain-ref. 진단 변화 0.11 dB |
| 시드 누설 | 조건부 fixture | ideal covariance 가정에서 DC 재균형은 진단값을 거의 바꾸지 않음 |
| 온도 | gain/진단 추세만 | δ 고정 진단 최저점 120 °C; physical squeezing·SQL 판정 불가 |
| 셀 길이 | 보충 | qL·OD 양쪽에 선형 → 밀도와 N·L로 묶임 |
| 교차각 | 절반 | 겹침이 아니라 Δk_z가 지배. 모델은 구조상 0°를 선호 |
| Loss | ideal-law 진단 | 0.20 dB/%p; `10log10(1−η)`는 물리적 floor가 아닌 대수적 점근선 |
| 외부 빔 유입 | 조건부 항등식 | assumed ideal covariance에서만 S = ηS_ideal + (1−η) ≤ 1; 측정 SQL 판정 아님 |

## 이 감사에서 나온 독립 결과

**문헌 동작점 산포와 v6 gain-referred 진단 창은 정성적으로 일치한다.** 8편이
독립적으로 고른 TPD가 전부 13.7 MHz 폭 안에 있는데, v6가 1D 스캔으로 계산한
+0.5 dB 진단 창
(−9.6/+4.2 MHz)의 폭이 13.8 MHz다 (비 0.99). OPD·온도는 창의 약 2배.
반대로 시드 파워는 40배, 펌프 파워는 4.4배(세기로는 5.1배) 흩어져 있다 —
v6의 조건부 난이도 순서 `TPD ≈ OPD ≫ T ≫ probe power`와 같은 그림이다. 이는
물리적 squeezing tolerance의 검증이 아니다.

## 판독 규약 · 주의

- **엔진의 절대 이득은 예측값이 아니다** (v2 「Model point」: gold 동작점에서
  +7500%). 읽을 것은 민감도와 순서. 파워/겹침 지렛대 비는 연쇄법칙으로
  절대 이득이 약분되도록 구성했다.
- 이득 기준이 필요한 곳(§7 loss)은 엔진 대신 **측정 이득 G = 111/8**을 쓴다.
  이때 G_c는 반드시 G−1로 둘 것 — 원자료 109/8을 넣으면 gap이 0.25로
  twin-beam 유효 범위를 벗어나 −26 dB짜리 가짜 값이 나온다
  (v2 §Gain convention).
- JSON의 primary 결과 키는 `gain_referred_*`, `ideal_law_*`,
  `above_normalized_unity`다. `S_dB`, `xi_*`, `above_sql`은 호환용 alias이며
  물리적 squeezing·SQL claim을 뜻하지 않는다.
- §6 저각 수치는 설계 근거로 쓰지 말 것. Option A 규약에서 Δk_geom이 θ=0에서
  항등적으로 0이라 엔진이 구조적으로 0°를 선호한다 (v6 「해석」 4항).

## 관련

- `analysis/squeezing/low_pump_power/` — 저출력 펌프 mean-field/ideal-law 분석
- `.claude/skills/fwm-squeezing-frontier/` — (Δ, T) 평면의 효율 전선
