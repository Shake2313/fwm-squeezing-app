# 저출력 펌프 FWM gain/ideal-law 진단

Jain, Choi, Hull, Marino, *Opt. Lett.* **50**, 5165 (2025) — 135 mW 펌프로
셀 직후 −7.2 dB (파이버 출력 −4.4 dB) 세기차 스퀴징 — 의 저출력 동작점을 GABES와
`docs/FWM physics and analytic reconstruction/squeezing_analytic_reconstruction_v2.tex`의
평균장 민감도 및 조건부 ideal twin-beam 법칙과 대조하는 계산이다. 미시적
quantum-Langevin 잡음이 없으므로 GABES 자체의 물리적 squeezing 예측이나
독립 검증으로 읽지 않는다.

```bash
python analysis/squeezing/low_pump_power/scan_low_pump_power.py
```

약 110 초. 결과는 `generated/`에 떨어진다.

| 파일 | 내용 |
|---|---|
| `generated/low_pump_power.md` | 한글 보고서 (모든 표) |
| `generated/low_pump_power.json` | 원시 수치 |
| `generated/low_pump_power.png` | 4-패널 그림 |

## 주장의 세 다리

1. **P/w² 불변성** — 펌프는 Ω² ∝ P/w²로만 이론에 들어간다. GABES legacy에서
   (600 mW, 530 µm)와 (192 mW, 300 µm)의 `G_s_smallsignal` 배열이 **비트 단위로
   동일**함을 확인한다. 엔진의 유일한 절대 파워 의존성은 Manley–Rowe 광자
   예산(`observables.pump_depletion_saturation`)이고, 실제 변환율은 ~0.3%다.
2. **발표된 동작점은 모두 같은 세기** — Dowran 140, Sim 136, Jain 95–177 W/cm².
   Jain의 135–250 mW 창이 550–600 mW 실험들을 라비 주파수에서 감싼다.
3. **필요한 이득이 작고, 이득은 파워에 대해 포화** — −7 dB에 G ≈ 8이면 되고
   (η = 0.8694), 그 영역에서 dln(G−1)/dlnP ≈ 0.1–0.16.

## 층 분리 (중요)

- **평균장 이득**: GABES legacy fidelity `G_s_smallsignal` (δ 최대점).
  legacy를 쓰는 이유는 P/w² 불변성이 그 층에서 정확하기 때문이다. Ultra의
  겹침·동적 고갈 보정은 `invariance` 절에서 따로 잰다.
- **잡음**: 폐형식 twin-beam 법칙 `observables.intensity_difference_squeezing_dB`.
  이는 측정 이득을 입력으로 받는 조건부 ideal-law 진단이다. v2 문서의 layer
  표대로 미시적 quantum-Langevin 층은 아직 없다.

**엔진의 절대 이득은 예측값이 아니다** (v2 「Model point and experimental
comparison」: Sim 동작점에서 +7500%). 보고서의 헤드라인 수치는 모두 (a) 구조적
불변성, (b) 로그 민감도, (c) 측정 이득을 입력으로 받는 잡음 법칙 중 하나이며
절대 이득 눈금을 필요로 하지 않는다.

## 관련

- `.claude/skills/fwm-squeezing-frontier/` — (Δ, T) 평면의 효율 전선
- `analysis/squeezing/analytic_reconstruction/` — v2 문서의 수치 감사
