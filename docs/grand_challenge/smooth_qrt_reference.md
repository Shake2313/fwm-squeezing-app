# Independent smooth full-density QRT reference

`analysis/grand_challenge/reference/smooth_qrt.py` exports:

```python
smooth_qrt(h0, h1, envelope, reservoirs, boundary_state, duration_s,
           readouts, frequencies_rad_s, *, rtol=2e-11, atol=2e-14,
           max_step_s=None)
```

The prescribed Hamiltonian is `H(t) = h0 + envelope(t) * h1`. Both matrices
have shape `(n,n)` and units rad/s; each must be Hermitian. `envelope(t)` is
a real scalar callback evaluated in **physical seconds**, including when the
solver internally uses dimensionless time. Reservoirs are one fixed
`ExplicitReservoirs` instance. The boundary state is a physical trace-one
density `(n,n)` and may be nonstationary. The duration is positive in seconds.

Readout operators are fixed complex matrices `(p,n,n)`. Angular frequencies
are a real array `(nf,p)` in rad/s; ports in a row may have unequal frequencies
of either sign. For each row, the finite pulse is

\[
Y_j=\int_0^T e^{i\omega_j t}O_j(t)\,dt,\qquad \mu_j=\langle Y_j\rangle.
\]

Here the time dependence of the operator is its physical evolution; the
supplied readout matrix is fixed. The reservoir channels and readout matrices
cannot vary with the envelope through this API. No cross-frequency-row
moments are returned.

| Result key | Shape or type | Meaning |
|---|---|---|
| `greater` | `(nf,p,p)` | Connected \(\langle\delta Y_j\delta Y_k^\dagger\rangle\) |
| `lesser` | `(nf,p,p)` | Connected \(\langle\delta Y_k^\dagger\delta Y_j\rangle\) |
| `raw_greater`, `raw_lesser` | `(nf,p,p)` | Moments before global mean subtraction |
| `mean_pulse` | `(nf,p)` | \(\mu\), in readout units times seconds |
| `exit_state` | `(n,n)` | Full-density state at \(T\) |
| `residence_time_s` | scalar | \(T\) |
| `evaluations` | integer | DOP853 RHS evaluation count |
| `qrt_block_dimension` | integer | Total dimension, sharing one density across all rows |
| `affine_generator_nnz` | pair of integers | Nonzero entries in the compiled CSR matrices \(G_0,G_1\) |
| `stored_time_points` | integer | One: only the endpoint is retained |
| `solver_method` | string | `DOP853` |
| `rtol`, `atol`, `max_step_s` | scalars | Requested tolerances and effective maximum step |

Second moments have readout-product units times seconds squared. The result
does not include retarded response; the main microscopic implementation owns
that calculation. No drift `A`, diffusion `D`, stationary-noise result, or
covariance propagator enters this reference.

## Raw QRT equations

Let \(U(t,s)\) be the density propagator of the supplied Hamiltonian and
reservoirs. For greater ordering, define

\[
X_k(t)=\frac1T\int_0^t U(t,s)[O_k^\dagger\rho(s)]e^{-i\omega_k s}\,ds,
\qquad x_k(t)=e^{i\omega_k t}X_k(t).
\]

For lesser ordering use \(\rho(s)O_k^\dagger\) as the source. These are raw
sources: the operator means are not removed here. Define a normalized
positive-time triangle and its rotating coordinate by

\[
B_{jk}(t)=\frac1T\int_0^t e^{i\omega_j v}\operatorname{Tr}[O_j X_k(v)]\,dv,
\qquad Z_{jk}(t)=e^{-i(\omega_j-\omega_k)t}B_{jk}(t).
\]

In physical time,

\[
\dot\rho=L(t)\rho,\qquad
\dot x_k=(L(t)+i\omega_k)x_k+
\begin{cases} O_k^\dagger\rho/T,&>,\\ \rho O_k^\dagger/T,&<,\end{cases}
\]
\[
\dot Z_{jk}=i(\omega_k-\omega_j)Z_{jk}+\operatorname{Tr}(O_jx_k)/T,
\qquad
\dot m_j=-i\omega_jm_j+\operatorname{Tr}(O_j\rho)/T.
\]

Initially `x = Z = m = 0`. At exit,

\[
Q_{jk}=T^2e^{i(\omega_j-\omega_k)T}Z_{jk}(T),\qquad
\mu_j=Te^{i\omega_jT}m_j(T).
\]

For each ordering, the raw moment is \(Q+Q^\dagger\), and the connected moment
is this raw moment minus \(\mu\mu^\dagger\). The Hermitian partner is the other
time integration half-plane. Equal times have zero integration measure for
these ordinary pulses. This construction applies to unequal port frequencies
without dropping their difference phases.

## Sparse affine propagation

With row-major density vectorization,

\[
L_0=\mathcal D-i(h_0\otimes I-I\otimes h_0^T),\qquad
L_1=-i(h_1\otimes I-I\otimes h_1^T),\qquad
L(t)=L_0+f(t)L_1.
\]

The dissipator comes from the explicit reservoir channels. `L1` is assembled
directly from `h1`, avoiding cancellation from subtracting two dissipative
generators. All rotating sources, triangles and means are linear in the full
state. Using `u=t/T`, the lifted equation is

\[
\frac{dy}{du}=[G_0+f(uT)G_1]y.
\]

The code compiles both matrices to CSR once from small matrix blocks. It does
not allocate a dense lifted matrix. Each RHS evaluates `envelope(u*T)` and
`G0 @ y + envelope(u*T) * (G1 @ y)`; no generator assembly or tensor contraction
occurs during integration. All frequency rows share the same density, giving
the total dimension

\[
n^2+nf\,(2pn^2+2p^2+p).
\]

This is 82 for two levels, two ports and three frequency rows, or 508 for four
levels, four ports and three rows. The reported dimension differs from the
segmented reference's dimension per frequency row. Unlike the segmented
reference's separate density exponential, this density shares DOP853 error
control with the auxiliary variables; changing the readouts or frequency grid
can change its numerical trajectory within integration error.

`solve_ivp(method="DOP853", t_eval=[1.], dense_output=False)` retains only the
endpoint and fixed solver workspace. The rotating variables are an exact
coordinate transformation of these QRT equations. The integrator still
resolves their fast modes: no frequency term or Hamiltonian coupling is
discarded, and no additional rotating-wave or envelope-segmentation
approximation is made.

`max_step_s=None` means `T/64`. For an envelope with narrower features, supply
a sufficiently small maximum step and check refinement. The solver cannot
infer unsampled features of an arbitrary callback. The real scalar and finite
value checks apply at every RHS evaluation. Nonfinite input arrays, invalid
shapes, non-Hermitian Hamiltonians, unphysical entry states and invalid scalar
duration/tolerance/step values are rejected.

The tolerances apply to the normalized density, source, triangle and mean
coordinates. They do **not** guarantee the same relative accuracy for small
physical pulse moments. In particular, restoring seconds squared and
subtracting the global mean outer product can expose absolute-error limits
and cancellation. Outputs are not clipped, renormalized or repaired. For
Poisson arrivals, `poisson_beam_spectrum` accepts this dictionary and restores
the raw mean outer product responsible for arrival number noise.

## Focused verification and scope

The owned tests are run with:

```powershell
python -B -m pytest -q -p no:cacheprovider tests/quantum/test_smooth_qrt.py
```

Final focused run on 2026-09-14: **40 passed in 24.07 seconds**. The source
and test files are frozen for the main audit with these SHA-256 digests:

| File | SHA-256 |
|---|---|
| `analysis/grand_challenge/reference/smooth_qrt.py` | `385d5097bb9013ce84f79e0471d3834c9085a583d0f583ad0bdc83c8c45ec8b7` |
| `tests/quantum/test_smooth_qrt.py` | `ac2b1d2d14be7ca9ff5bf9c802c2a5ccbb7a2e94d82952b1ff0254fad27e4bff` |

The older `characteristic_qrt` supplies an independent smooth comparison. It
integrates centered sources and explicit Fourier phases in physical seconds,
assembling the full density generator at every RHS evaluation. To exercise
unequal port frequencies through its common-frequency API, its readout is
`exp(1j*(omega_j-common_omega)*t)*O_j`. Constant-envelope cases are compared
with `segmented_qrt` exponentials. Other checks cover an analytic damped
finite Lorentzian window, SI time rescaling, batching, complex identity
offsets, Poisson number noise and invalid inputs.

Measurements on 2026-09-14, Python 3.14.2 / NumPy 2.4.1 / SciPy 1.17.0, using
one BLAS thread:

| Control | Observed error or cost |
|---|---|
| Smooth noncommuting two-level drive, complex readouts, three unequal signed frequency rows | Maximum relative error `6.60e-16` across both connected moments, mean and exit density against the older centered-source QRT |
| Driven identity readouts | Maximum absolute connected residual `1.33e-15` in the toy time units |
| Smooth reduced 85Rb D1, 40 ns, opposite hyperfine-frequency ports | Relative errors: greater `2.01e-10`, lesser `9.52e-11`, mean `8.14e-9`, exit density `2.99e-11` against the older QRT |
| Same 40 ns case, tightening `(rtol, atol)` from `(2e-11, 2e-14)` to `(2e-12, 2e-15)` | Changes: greater `1.45e-10`, lesser `8.05e-11`, mean `8.78e-9`, exit density `3.91e-11` |
| 40 ns sparse solve at default tolerances | About `0.45 s`; `30,449` RHS evaluations; dimension `90`; CSR nonzeros `(156, 320)` |
| 40 ns sparse solve at tighter tolerances | `40,601` RHS evaluations |

The short Rb fixture uses 60 MHz peak pump Rabi frequency, 0.9 GHz detuning,
the supplied radiative reservoirs and entry population `(5/12, 7/12, 0, 0)`.
Its readout frequencies are `-OMEGA_HF + 2*pi*0.4e6` and
`OMEGA_HF - 2*pi*0.7e6`. The real envelope is
`exp(-((t/T - 0.45)/0.3)**2)`. Norms are approximately `8.65e-19 s^2` for the
greater moment, `5.46e-19 s^2` for the lesser moment and `2.28e-13 s` for the
mean. Errors are checked per frequency row so a DC row cannot mask a small
carrier-frequency result. Runtime is a local measurement, not a guarantee.

This sidecar validates the reference only on the stated controls. It does not
certify the main 2 microsecond moving-Rb fixture, a thermal ensemble, a
self-consistent optical-field model or experimental squeezing. The main task
owns the full 2 microsecond audit and repository-wide pytest run. Only this
new reference, its new tests and this document are changed by the sidecar.
