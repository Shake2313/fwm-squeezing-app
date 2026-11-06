# PSR S0–S4 코드 리뷰

사용자 요청: 여기까지 코드 리뷰 수행. 문제 발견 시 수정.
검토일: 2026-11-06.

## ⟦AI⟧ 판정·범위

추가 수정이 필요한 P1/P2 오류 미발견. 이번 리뷰의 코드 수정 없음.
S0–S4 코드·호출 관계 정적 검토와 아래 재검증 기준. 모든 운전점의 물리적 정확성이나 실험 일치 보증은 아님.

| 검토 경로 | 확인 내용 |
|---|---|
| `psr_conventions`, `psr_full_obe` | Jones·Stokes 부호, detuning·속도·reset 단위, 구동 정상상태, 방출·reset의 같은 generator 사용, 약한 여기 비교 |
| `psr_spectrum`, `psr_homodyne` | ±RF 반사와 두 ordered noise, SQL 정규화, 출력 Jones 직교 모드, 공통 위상, 검출 손실 1회 적용 |
| `psr_thermal` | drift·noise를 속도 평균한 뒤 전파, 적분 차수·cutoff·출력 covariance 검사, 유한 평가 예산 |
| `psr_propagation` | 중간점 원자·잡음 재계산, 공유 구간 합성, 균일·얇은 셀 극한, 길이·속도 수렴, 미수렴 진단 경로 |
| PSR scheme·비교·JSON | 숨김 입력 배제, 실제 적용 입력과 보관 상태, 선택점·위상 곡선·covariance 일치, 검출 cache 재사용, 실패 수치 보류 |

S4 benchmark의 수치 코드·공유 의존성 SHA-256 모두 현재 파일과 일치.
측정 후 설명만 바꾼 benchmark script는 기존 `metadata_correction`의 SHA-256과 일치.
기존 측정 기록·봉인 자료 변경 없음.
S0–S4 문서와 PSR README의 로컬 링크 누락 0건.

## ⟦AI⟧ 재검증

```text
python -m pytest -q tests/test_psr_propagation.py tests/test_psr_propagated_scheme.py tests/test_psr_thermal.py tests/test_psr_thermal_scheme.py tests/test_psr_spectrum.py tests/test_psr_homodyne.py tests/test_psr_rf_scheme.py tests/test_psr.py tests/test_psr_noise.py tests/test_psr_scheme.py tests/test_psr_conventions.py tests/test_psr_full_obe.py tests/test_squeezing_mode_selection.py tests/test_defaults_store.py tests/test_polarization_comparison.py tests/test_magneto_doppler.py tests/test_magneto_multilevel.py tests/quantum/test_traveling.py tests/quantum/test_atomic_diffusion.py tests/quantum/test_sideband_readout.py tests/quantum/test_contracts_and_channels.py --durations=15
```

**329 passed**, 21 warnings, 135.14 s.
경고는 기존 JSON export의 `datetime.utcnow()` deprecation.
PSR S0–S4와 연결된 모드 선택·설정·편광 비교, 영향을 받는 magneto 계산, 공유 quantum 전파·확산·검출 계약 포함.
코드 수정 없는 검토이므로 전체 pytest와 별도 GC 연구 캠페인은 실행하지 않음.
브라우저는 재실행하지 않음. 이번에는 그림 생성·JSON 일치 회귀를 재실행했고, 실제 다운로드 기록은 [S4 검증](s4.md)의 2026-11-05 기록 참조.

추가 경계 운전점: Δ=400 MHz, RF=±400 kHz, n=10¹⁶ m⁻³, L=25 mm, ψ=0.3 rad.

| 입력 ε [rad] | I [mW/cm²] | reset [s⁻¹] | 수렴 구간 | 출력 I [mW/cm²] |
|---|---:|---:|---:|---:|
| +π/4 | 0.1 | 1000 | 4 | 0.0953723875530 |
| −π/4 | 12.5 | 100000 | 4 | 12.0353927431 |
| 0 | 100000 | 10 | 4 | 99998.9142518 |
| 0.17 | 0.001 | 10 | 4 | 0.000951709172858 |

네 점 모두 원자·전파 검사와 길이 수렴 통과. 원편광의 방위각 미정 처리는 유지.

## ⟦AI⟧ 판정의 한계

모델 범위는 [S4 구현·근사](s4.md) 참조. 국소 정상상태·독립 속도군·1차원 전파 가정 유지.
열평균 대체 kernel 검사는 후보 평균 펌프 경로를 고정한 경험적 오차 추정이며, 완전히 결합된 속도·길이 적분의 엄밀한 오차 상한은 아님.
기존 S5/S6 계획은 `docs/checklist.json`의 `_psr_staged_model_2026_11_05`에 유지.

## ⟦AI⟧ 커밋 범위 분리 검증

2026-11-06 사용자 요청: Claude 작업과 충돌 없이 PSR S0–S4 구현·테스트·문서만 커밋·푸시.

34개 파일로 분리. PSR 정상상태·잡음·열평균·길이 전파와 관련 문서·benchmark·테스트 포함. 공유 변경은 기존 grid를 보존하는 `velocity_quadrature`의 선택적 cutoff/RF knot 인수와, JSON의 실제 적용 입력·보관 상태 구분만 포함.

Claude의 model/spectral comparison adapter·화면 레이아웃·공통 제어 변경과 기존 스테이징은 제외. PSR 파일의 비교 hook과 두 테스트의 비교 adapter 연동 구간도 별도 작업으로 남김. 펌프 그림 안내는 기존 figure selector에서도 성립하도록 표기.

선택된 Git tree와 기존 커밋의 의존성만 내보낸 독립 스냅샷에서 위 명령 중 `tests/test_polarization_comparison.py`를 제외해 재검증. **313 passed**, 21 warnings, 133.85 s. 나머지 16개 비교 기능 테스트는 앞선 결합 작업 트리의 329개 검증에 해당함. 공유 독립 QRT reference는 기존 커밋의 파일을 사용했으며 별도 GC 캠페인은 실행하지 않음. 전체 pytest는 미실행.

커밋 직전 선택 파일·문서 링크·diff 형식과 기존 스테이징·작업 파일 보존을 확인.
