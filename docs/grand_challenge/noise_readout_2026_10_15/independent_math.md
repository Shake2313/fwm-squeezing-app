# ⟦AI⟧ L1 판독식 독립 유도·반례 검산

2026-10-15. 출처: 독립 수학 검토 에이전트. 인계 문서 열기 전 입력·출력
연산자에서 재유도한 뒤 대조. 아래는 AI 검토 결과이며 새 사용자 목표 아님.
대상: production 평균장 전송의 조건부 판독 진단. GC 미시 확산 대체 아님.

## 가정·유도

광자 플럭스 기저 \(v=(a,b^\dagger)^T\), \(J=\operatorname{diag}(1,-1)\).
입력 평균 \(\langle a\rangle=\alpha\), \(\langle b\rangle=0\).
켤레 입력 진공, 시드 진폭 quadrature 분산 \(F\), 진공 분산 1.
밝은 시드 선형화로 spontaneous \(O(1)\) 세기 항 생략. \(M\)은 같은
평균·잡음 경로의 전송이어야 함. 두 팔의 사후 효율 \(\eta\) 동일 가정.

\(a=M_{00}a_{in}+M_{01}b_{in}^\dagger+n_0\),
\(b^\dagger=M_{10}a_{in}+M_{11}b_{in}^\dagger+n_1\).
\(G_s=|M_{00}|^2\), \(G_c=|M_{10}|^2\), \(D=G_s-G_c\), \(T=G_s+G_c\),
\(q=M_{00}^*M_{01}-M_{10}^*M_{11}\), \(\ell=(M_{00}^*,-M_{10}^*)\).
시드 위상을 흡수하면

\[
\delta I_- /|\alpha|=D\,\delta X_a+q b_{in}^\dagger+q^*b_{in}
 +\ell n+(\ell n)^\dagger.
\]

\(D\) 실수이므로 시드 위상·진폭/위상 공분산의 기여 0.
이 결론은 zero-mean 켤레 입력에 의존. 두 입력 모두 coherent mean이 있으면
동일한 위상 소거를 일반화할 수 없음.
또한 같은 유효 M이 carrier와 요동을 전파한다는 가정 필요.
유한 RF에서 평균은 \(M(0)\), sideband는 \(M(\pm\Omega)\)를 쓰므로
carrier 방향으로 투영한 시드 계수가 일반적으로 실수 아님.
단일 DC 2×2 M으로 주파수별 시드 위상 잡음 소거를 인증할 수 없음.

소멸형 진공 bath만 \(n=Cf\)라면
\(CC^\dagger=K=J-MJM^\dagger\succeq0\),
\(\ell K\ell^\dagger=D-D^2+|q|^2\). 따라서

\[
S_{src}=\frac{D+(F-1)D^2+2|q|^2}{T},\qquad
S_{det}=1-\eta+\eta S_{src},
\]
\[
F_p=(F-1)G_s+2|M_{01}|^2+1,\quad
F_c=(F-1)G_c+2|M_{11}|^2-1,
\]
\[
\Lambda=S_{src}T-D(F_p-F_c-1)=2|\det M|^2.
\]

독립 위상 무관 creation bath \(Ef^\dagger\) 추가 시 분자 증가량
\(2\|\ell E\|^2\ge0\). 이 조건에서만 signed L1식을 하한으로 사용 가능.
입력·bath의 상관 또는 위상 민감 bath에는 보장 없음.
\(G_c=0\)이면 반환 \(F_c\)는 형식값이며 밝은 켤레 팔의 측정 Fano 아님.

## 범위·반례

- \(K\succeq0\)는 endpoint 소멸형 진공 완성 가능성. 국소 프로브 손실만의
  전파를 인증하지 않음. \(M=\operatorname{diag}(0.9,10)\): \(K\succeq0\),
  \(|\det M|^2=81\), \(\Lambda=162\). 따라서 \(\Lambda\le2\)는 이 검사만으로 주장 불가.
- pair Hamiltonian+국소 프로브 손실만의 명시 모델에서는
  \(|\det M|^2=e^{-\int\alpha_p dz}\le1\), \(D\le1\),
  \(r=|M_{01}|^2-|M_{11}|^2\le-1\). \(0\le D\le1\)에서
  \(S_{src}\le[2-D+(F-1)D^2]/T\). 인계의 정리는 이 모델 범위에서 일치.
- 위상 민감 bath 반례: 프로브 투과율 \(\tau=0.4\), 진폭 압축 bath \(r=1\).
  같은 \(M\)에서 L1=1이나 실제 \(S=\tau+(1-\tau)e^{-2r}=0.48120117\).
  Gaussian CP 행렬·출력 uncertainty 모두 양정/준양정. 보편 하한 아님.
- \(M=\begin{psmallmatrix}0&0\\1&0\end{psmallmatrix}\): 형식 \(S_{src}=-1\),
  \(\lambda_{min}(K)=-2\). 음의 값 보존하되 dB 판독·추정치로 표시 금지.
  비양정 K에서 형식값이 양수여도 수동 완성 예측으로 표시 금지.
- 고이득 오판 반례: \(M=1.0001 B(r=8)\).
  \(\lambda_{min}(K)=-2.00010\times10^{-4}\).
  기존 후보 \(10^{-10}\sum|M|^2\) 허용오차는 \(8.88789\times10^{-4}\)여서
  실패를 통과시킴. \(64\epsilon_{mach}\max(1,\sum|M|^2)\)로 좁혀 반례 거부.
  이는 부동소수점 후방오차 규모 검사이며 물리 오차 허용치 아님.

## 독립 수치 검산·종결

참조는 물리 좌표 \((x_a,p_a,x_b,p_b)\)에서 pair Hamiltonian·손실률로
drift와 양의 diffusion을 직접 구성. Van Loan의
\(\exp L\begin{psmallmatrix}A&Q\\0&-A^T\end{psmallmatrix}\) 사용.
별도 DOP853 공분산 ODE와도 일치. \(M\)에서 diffusion 역산하지 않음.

160개 복소 결합·위상 부정합·1–5구간 프로브 손실 프로파일에서
S/Fp/Fc 상대오차 최대 \(2.78\times10^{-15}/1.08\times10^{-15}/1.82\times10^{-15}\).
위상 분산 1→50 변화 최대 \(8.0\times10^{-15}\),
\(\Lambda\) 절대오차 최대 \(1.78\times10^{-15}\).
켤레 손실 추가·진폭/위상 공분산·D=0 sub-SQL·상한·CP·API 오류·broadcast도 검사.

검증: `python -m pytest -q tests/test_fwm_passive_l1_reference.py` → **20 passed**.
후속 소비자 검토: masked NaN을 JSON null로 내보내는 strict JSON 검사 통과.
현재 F=4.05 판독과 cached F=1, invalid flag·선택점을 구분해 보존.
P1/P3 식·P6 flag 필요성 통과. 양의 formal 값과 음의/비양정 값을 구분해
표시해야 함. 다음 검토: production 소비자에서 같은 M·기준면·효율을 사용하는지 확인.
