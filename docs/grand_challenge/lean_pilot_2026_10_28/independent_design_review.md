# ⟦AI⟧ 독립 수학 설계 검토

검토자: 별도 에이전트 lean_pilot_review. 이 기록은 수학적 범위·반례 설계 검토이며 컴파일 실행 증거는 run/receipt.json이 소유한다. 별도 에이전트의 소스 읽기는 완료되지 않았으므로 독립 소스 감사 완료로 주장하지 않는다. 소스 검토·실행은 주 에이전트가 수행했다.

- C(g,y)(v)=gv+y와 합성(g2g1,g2y1+y2) 확인. g는 covariance/power 배율이며 field amplitude가 아님.
- 독립 입력/저장고·구간 저장고 가정 필요.
- 두 eta=1/2 손실: 올바른 y=3/8, vout=1/2. 잘못된 단순합은 vout=5/8, 첫 noise 누락은3/8.
- 일반 합성·결합·항등·loss closure와 유리수 반례는 작은 core/Std 시범에 적합.
- 진공 분산만 확인해 CP를 주장하면 안 됨. 예: g=2,y=0은 진공 출력v=1이지만 phase-insensitive CP 잡음 하한을 위반.
- 오류 대조의 비정상 종료만으로 성공 판정 금지. 최종 runner는 정확한 반례와 grind의 False propositions 진단을 함께 확인함.
- 일반 real-matrix CP·production 대응·GC acceptance는 본 시범 범위 밖.
