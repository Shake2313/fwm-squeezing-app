# FWM Squeezing Fast/Balanced 리마스터 — 완료 기준 (Gate v1)

동결: 2026-09-11, 구현 착수 전. 이후 기준을 바꾸려면 이 문서에 사유·날짜를 추가하고
이전 기준에 대한 판정도 함께 남긴다. 기준 완화로 PASS를 만들지 않는다.

## 공통 조건

- **Reference**: 같은 작업트리의 Ultra (`FWM_FIDELITY[Ultra]`, 1 m/s·4σ 격자, 401점, N_F=3/2 전체 감사).
- **Case matrix (12)**: default, v6_opt, frontier, near_res, blue_res, blue_far, red_edge,
  cold_weak, hot, low_pump, bright_seed, long_cell (`gate.py`의 `CASES`, 앱 기본값에 덮어쓰기).
- **측정 경로**: 앱과 같은 `FWMScheme.compute(params)`.
- **시간**: numba 16 threads, JIT/캐시 로드 제외 warm 상태, 다른 벤치/테스트 동시 실행 없음, 5회 중앙값.
- **개형 창**: |δ| ≤ 500 MHz (앱 그림의 x 범위), Ultra와 같은 401점.
- tube 오차 = 가로 ±1 MHz 여유 안에서의 최소 세로 오차 (가파른 측면에서 한 격자 어긋남을 과대평가하지 않기 위함).

## G1 속도

| ID | 기준 |
|---|---|
| G1.1 | Fast: default ≤ **0.12 s**, 12 cases 모두 ≤ 0.25 s |
| G1.2 | Balanced: default ≤ **0.25 s**, 12 cases 모두 ≤ 0.40 s |
| G1.3 | 변경 전 동결 측정 대비 default 가속: Fast ≥ **4.5×** (0.58 s 기준), Balanced ≥ **7×** (1.87 s 기준) |
| G1.4 | 새 프로세스(디스크 numba 캐시 존재)의 첫 Fast 호출 ≤ 3 s |

## G2 개형 (각 케이스 모두 만족)

| 지표 | Fast | Balanced |
|---|---:|---:|
| S_dB tube max | ≤ 0.10 dB | ≤ 0.05 dB |
| S_dB p95 (점별) | ≤ 0.03 dB | ≤ 0.015 dB |
| log10 G_s tube max | ≤ 0.010 | ≤ 0.005 |
| 동작점 δ=tpd의 \|ΔS\| | ≤ 0.05 dB | ≤ 0.02 dB |
| 동작점 \|ΔG_s\|/G_s | ≤ 1 % | ≤ 0.5 % |
| S_dB 최솟값 차, Ultra 최솟점에서의 차 | ≤ 0.05 dB | ≤ 0.03 dB |
| G_s 최댓값 상대차 | ≤ 2 % | ≤ 1 % |
| Pearson(S_dB) (Ultra 범위 ≥ 0.01 dB일 때) | ≥ 0.998 | ≥ 0.9995 |
| Pearson(log10 G_s) (Ultra 범위 ≥ 1e-3 dex일 때) | ≥ 0.998 | ≥ 0.9995 |

## G3 무결성

| ID | 기준 |
|---|---|
| G3.1 | Fast/Balanced `probe_axis_GHz`가 Ultra와 `array_equal` |
| G3.2 | G_s, G_c, S_dB 전 점 finite |
| G3.3 | N_F=3 vs 2 감사가 **풀이한 모든 δ 노드**에서 수행되고 12 cases 모두 CONVERGED (Balanced는 401 전 점) |
| G3.4 | 극점 guard: 풀이한 모든 행에서 pole form(Δ_eff=Δ) vs 기존 compiled kernel 상대오차 ≤ 1e-6 |
| G3.5 | 정확성: Faddeeva 평균 vs 0.5 m/s·8σ 격자 ≤ 1e-8 (단위 테스트, 양 branch, 양 seed, 양 차수) |
| G3.6 | provenance: raw에 방법·풀이 δ 수·Doppler 방식·감사 범위, Fast claim gate에 보간 사유 명시 |
| G3.7 | Ultra 불변: default Ultra의 G_s/G_c/S_dB가 변경 전 캡처와 bit-identical, Ultra tier dict 불변 |
| G3.8 | `GABES_DISABLE_NUMBA=1` 경로의 Fast/Balanced가 numba 경로와 상대 1e-9 이내 |

## G4 테스트

| ID | 기준 |
|---|---|
| G4.1 | `python -m pytest -q`: 변경 전 기준선 대비 신규 실패 0 |
| G4.2 | 신규 테스트 파일 통과 (numba on / off) |
| G4.3 | `tests/baseline_focused.npz` 재캡처 없음 |

## G5 기록

| ID | 기준 |
|---|---|
| G5.1 | `analysis/fwm_lite/DEVLOG.md`: 시도·기각·수치·다른 scheme 재사용 가이드 |
| G5.2 | `python -m analysis.fwm_lite.gate` → JSON 보고서, 전 항목 PASS일 때만 exit 0 |
