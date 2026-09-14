# Passive finite-velocity bosonic transport bridge

This independent reference is a **two-oscillator passive testbed**. It verifies
canonical inflow, transport memory, reservoir noise, and optical complete
positivity (CP). It does not implement a finite-dimensional Rb atom, FWM
squeezing, or a general moving-vapor closure. It imports neither the production
transport module nor the Lindblad characteristic solver. Its role is to provide
an exact bosonic control for the separate `gabes/quantum/transport.py` work.

## Ports, units, and exact channel

The forward, flux-normalized annihilation ports are \(u=(b,a)^T\), with
\([u_i(0),u_j^\dagger(0)]=\delta_{ij}\). A common frequency delta function is
suppressed; equivalently use normalized frequency modes. The optical and atomic
inflows and local reservoir are independent. At a fixed real \(\Omega\),

\[
\partial_z u=K(z)u+Bf(z),\qquad
K=\begin{pmatrix}i\Omega/c&-i\kappa(z)\\
-i\kappa(z)^*&(i\Omega-\gamma)/v\end{pmatrix},\quad
B=\binom{0}{\sqrt{2\gamma/v}},\quad Q=BB^\dagger.
\]

Here \(c,v>0\), \(\gamma\ge0\), and
\([f(z),f^\dagger(z')]=\delta(z-z')\). With length unit \(\ell\) and time
unit \(t\), velocities have units \(\ell/t\), \(\Omega,\gamma\) have units
\(1/t\), and \(\kappa\) has units \(1/\ell\). The example uses arbitrary
declared units; its velocities and coupling are not Rb parameters. Flux
normalization and the conjugate coupling make the exact identity
\(K+K^\dagger+Q=0\) possible without velocity weights in the commutator.

Writing \(U(z,s)\) for the ordered fundamental solution gives

\[
u(L)=T u(0)+\int_0^L U(L,s)Bf(s)\,ds,\quad T=U(L,0),\qquad
W=\int_0^L U(L,s)Q U(L,s)^\dagger\,ds.
\]

The exact covariance equations are \(T'=KT\) and
\(W'=KW+WK^\dagger+Q\), with \(T(0)=I,W(0)=0\). Consequently
\(M=TT^\dagger+W\) obeys \(M'=KM+MK^\dagger+Q\), whose unique solution
is \(M=I\). This proves preservation of both output commutators, including
their vanishing cross commutator.

For a constant segment of length \(h\), the implementation computes
\(E=e^{Kh}\) and the Gramian separately. It exponentiates the augmented affine
Lyapunov generator with upper blocks
\([I\otimes K+K^*\otimes I,\operatorname{vec}_F Q]\) and a zero last row.
This is the exact integral, works with a singular Lyapunov map, and does not
subtract two nearly equal matrices. The identity \(W=I-EE^\dagger\) is used
only as a check. Noise is never fitted, eigenvalue-clipped, or completed from a
commutator residual.

For segments numbered from input to output,

\[
T_{j+1}=E_jT_j,\qquad W_{j+1}=E_jW_jE_j^\dagger+W_j^{\rm local}.
\]

Later segments multiply on the left. Earlier noise is propagated through every
later segment; it cannot simply be added at the output without propagation.

## Optical reduction, vacuum, and thermal inflow

Trace the **atomic output** to obtain the optical channel. Its input operator
remains

\[
b_{\rm out}=T_{bb}b_{\rm in}+T_{ba}a_{\rm in}+b_{\rm res}.
\]

Define \(\tau=|T_{bb}|^2\), \(\beta=|T_{ba}|^2\), and
\(r=W_{bb}\). All three are nonnegative, and

\[
\boxed{\tau+\beta+r=1.}
\]

The reported reservoir Gramian is a **commutator** contribution, independent of
thermal occupation. For quadratures \([x,p]=i\) and \(V_{\rm vac}=I/2\),
independent thermal atomic and reservoir occupations \(n_a,n_r\ge0\) give

\[
X=\begin{pmatrix}\Re T_{bb}&-\Im T_{bb}\\\Im T_{bb}&\Re T_{bb}\end{pmatrix},
\quad Y=\big[(n_a+\tfrac12)\beta+(n_r+\tfrac12)r\big]I.
\]

The CP matrix is \(Y+i(J-XJX^T)/2\), where
\(J=\left(\begin{smallmatrix}0&1\\-1&0\end{smallmatrix}\right)\).
Since \(XJX^T=\tau J\), its smallest eigenvalue is
\(n_a\beta+n_rr\ge0\). For vacuum optical input the output occupation is
the same expression, and its variance is \(1/2+n_a\beta+n_rr\).
With all inputs vacuum, the output is exactly vacuum. Passivity introduces no
creation operators or anomalous moments, so this control predicts no squeezing.

Dropping atomic inflow changes the commutator to \(\tau+r=1-\beta\).
Even though the remaining vacuum noise \(rI/2\) is PSD, the CP matrix then has
minimum eigenvalue \(-\beta/2\). Its apparent variance
\((1-\beta)/2\) is an invalid output for a declared canonical optical port.
This is an analysis-only negative control; no repair noise is introduced.

## Elimination creates spatial memory and correlated noise

Let \(\lambda=(i\Omega-\gamma)/v\) and \(g(z)=e^{\lambda z}\) for
\(z\ge0\). Eliminating the atomic equation as an exact integral identity gives

\[
a(z)=F(z)-i\int_0^z g(z-s)\kappa(s)^*b(s)\,ds,\qquad
F(z)=g(z)a_{\rm in}+\sqrt{2\gamma/v}\int_0^z g(z-s)f(s)\,ds.
\]

The resulting optical equation retains a causal memory term:

\[
b'(z)=\frac{i\Omega}{c}b(z)-i\kappa(z)F(z)
-\kappa(z)\int_0^z g(z-s)\kappa(s)^*b(s)\,ds.
\]

For \(m=\min(z,z')\), the two source commutator kernels are

\[
C_a(z,z')=g(z)g(z')^*,\qquad
C_r(z,z')=e^{\lambda(z-m)+\lambda^*(z'-m)}
                 (1-e^{-2\gamma m/v}).
\]

Both are PSD kernels and are kept separately as full complex matrices. Their
sum is exactly

\[
C_a(z,z')+C_r(z,z')=
e^{-\gamma|z-z'|/v+i\Omega(z-z')/v}.
\]

The symmetrized kernel is \((n_a+1/2)C_a+(n_r+1/2)C_r\). Only equal
occupations allow a common thermal multiplier. At zero damping the reservoir
vanishes but the shared inflow still gives a rank-one, spatially correlated
kernel. \(F\) is the freely propagated atomic source entering the memory
equation; it is not the full atom \(a\), which also includes optical feedback.

The exact optical response to forcing at \(z\) is \(U_{bb}(L,z)\).
Thus, with \(h(z)=-i\kappa(z)U_{bb}(L,z)\),

\[
b_{\rm out}=T_{bb}b_{\rm in}+\int_0^L h(z)F(z)\,dz,
\quad \beta=\iint h(z)C_a(z,z')h(z')^*\,dz\,dz',
\quad r=\iint h(z)C_r(z,z')h(z')^*\,dz\,dz'.
\]

These identities follow by applying the same causal resolvent to the Volterra
equation; they do not posit an independent noise at each position. The constant
example checks them with 48, 96, 192, and 384 Gauss-Legendre nodes. The reservoir
kernel has a derivative kink on \(z=z'\), so the unsplit double quadrature
converges at approximately second order. Its finite quadrature error is reported
separately from the exact segment identity.

An explicitly labelled negative control discards only off-diagonal entries of
the sampled reservoir kernel, keeping atomic inflow intact. The remaining
kernel is PSD but loses the shared reservoir history and fails the optical
commutator. Its diagonal entries are **not** delta-function densities; multiplying
them by an inverse grid spacing would invent a different reservoir. No such
rescaling or local noise fit is performed.

## Reproducible controls and precision

`build_control()` returns a strict JSON-serializable dictionary and writes no
files. It exposes declared parameters, separate boundary/reservoir contributions,
complex matrices as `{"real": ..., "imag": ...}`, numerical errors, CP
eigenvalues, refinement records, negative controls, thresholds, and `passed`.
The same result is available as `expected_controls_passed`; top-level `source`,
`frequency`, and `units` identify provenance and the `exp(-i Omega t)` convention.
The main audit can import it directly and hash the three sidecar sources.

```python
from analysis.grand_challenge.reference.ballistic_linear_channel import build_control
control = build_control()
assert control["passed"]
```

The other callable interfaces are `Parameters`, `Segment`, `segment_channel`,
`propagate_segments`, `direct_profile_ode`, `optical_diagnostics`, and
`atomic_source_covariance`. `PassiveChannel.transfer` includes both inflows;
`PassiveChannel.reservoir` includes only the distributed reservoir commutator.
No public propagation API has a switch to remove inflow or diagonalize noise.

Default values are \(L=1.3,c=4,v=0.65,\Omega=0.8,\gamma=0.35\), with
constant \(\kappa=0.9+0.35i\). The vacuum control uses \(n_a=n_r=0\);
the separate thermal control uses \(n_a=0.7,n_r=0.3\).

The nonconstant control uses \(u=z/L\) and
\(\kappa(z)=[0.85+0.20\cos(2\pi u)]
\exp(i[0.35+1.1u+0.25\sin(2\pi u)])\). Midpoint constant segments with
8, 16, 32, 64, and 128 cells are compared to an independent direct adaptive
amplitude/covariance ODE using DOP853 (`rtol=2e-12`, `atol=2e-14`). The ODE
uses the continuous profile and imposes no commutator identity. Matrix errors
are maximum absolute entry errors.

| Check | Absolute gate / observation |
|---|---|
| Segment commutator, Hermiticity, vacuum variance; PSD/CP tolerance | `5e-12` |
| Constant exponential/Gramian versus direct ODE | `2e-10` |
| Finest nonconstant transfer and reservoir errors | each below `2e-5` |
| Successive midpoint refinement order | greater than `1.9` |
| Finest double-kernel quadrature commutator error | below `5e-6` |
| Negative-control commutator deficit | greater than `0.05` |

A representative float64 run gives:

| Contribution / check | Value |
|---|---:|
| Optical input \(\tau\) | 0.294558587349037 |
| Atomic inflow \(\beta\) | 0.397925485039272 |
| Reservoir \(r\) | 0.307515927611691 |
| Full segment commutator error | \(2.3\times10^{-16}\) or less |
| Constant ODE transfer / reservoir error | \(5.2\times10^{-14}\) / \(6.9\times10^{-14}\) |
| 128-segment transfer / reservoir error | \(1.32\times10^{-5}\) / \(4.42\times10^{-6}\) |
| Final observed midpoint order | 2.00055 |
| 384-node correlated-kernel commutator error | \(1.42\times10^{-6}\) or less |
| Omitted-inflow commutator deficit | 0.397925485039272 |
| Omitted-inflow CP minimum eigenvalue | -0.198962742519636 |
| Diagonalized-reservoir commutator deficit, inflow retained | 0.305682832027507 |

Targeted verification, leaving the full suite to the main worker:

```powershell
python -B -m pytest -q -p no:cacheprovider tests/quantum/test_ballistic_linear_channel.py
```

Tests also cover the analytic lossless beam splitter, zero length, zero
coupling, free atomic attenuation, complex couplings and both frequency signs,
independent Gramian quadrature, ordered discontinuous segments against a split
ODE, thermal covariance against a full covariance ODE, spatial kernel quadrature,
zero-damping shared inflow, and strict JSON serialization.
