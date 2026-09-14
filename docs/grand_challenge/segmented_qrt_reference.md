# Independent segmented full-density QRT reference

`analysis/grand_challenge/reference/segmented_qrt.py` exports:

```python
segmented_qrt(hamiltonians, reservoirs, boundary_state, durations_s,
              readouts, frequencies_rad_s)
```

Inputs are Hamiltonians `(S,n,n)` in rad/s, one shared `ExplicitReservoirs`, a
physical trace-one boundary density `(n,n)`, strictly positive segment durations
`(S,)` in seconds, complex operators `(S,p,n,n)`, and real angular frequencies
`(nf,p)`. Segments are ordered from entry to exit; `T = sum(durations_s)`.
Zero-duration segments, an empty protocol, invalid shapes, nonfinite data,
non-Hermitian Hamiltonians and unphysical boundary states are rejected.
Hamiltonians and readouts may change at every boundary. The entry state may
be nonstationary, including in a closed system without a unique steady state.

For each frequency row independently, define

\[
Y_j = \int_0^T e^{i\omega_j t}O_j(t)\,dt,\qquad
\mu_j=\langle Y_j\rangle.
\]

The result is a dictionary:

| Key | Shape | Meaning |
|---|---|---|
| `greater` | `(nf,p,p)` | Connected \(\langle\delta Y_j\delta Y_k^\dagger\rangle\) |
| `lesser` | `(nf,p,p)` | Connected \(\langle\delta Y_k^\dagger\delta Y_j\rangle\) |
| `raw_greater`, `raw_lesser` | `(nf,p,p)` | Corresponding moments before subtracting \(\mu\mu^\dagger\) |
| `mean_pulse` | `(nf,p)` | \(\mu\), in operator units times seconds |
| `exit_state` | `(n,n)` | Full-density state at exit |
| `residence_time_s` | scalar | \(T\) |
| `qrt_block_dimension` | integer | \(n^2+2pn^2+2p^2+p\) |
| `qrt_exponentials`, `density_exponentials` | integers | \(S\,nf\) augmented and \(S\) density exponentials |

Second moments have operator-product units times seconds squared. No frequency
rows are mixed with other rows. No drift `A`, diffusion `D`, stationary noise
object, or propagated covariance is consumed. The shared density generator is
assembled from the supplied Hamiltonian and explicit reservoir channels.

## Rotating source and triangle equations

Let \(U(t,s)\) be the full-density propagator. For greater ordering, define

\[
X_k(t)=\frac1T\int_0^t U(t,s)[O_k(s)^\dagger\rho(s)]e^{-i\omega_k s}\,ds,
\qquad x_k(t)=e^{i\omega_k t}X_k(t).
\]

For lesser ordering the source is \(\rho(s)O_k(s)^\dagger\). Define the normalized
positive-time triangle and its rotating coordinate by

\[
B_{jk}(t)=\frac1T\int_0^t e^{i\omega_j u}\operatorname{Tr}[O_j(u)X_k(u)]\,du,
\qquad Z_{jk}(t)=e^{-i(\omega_j-\omega_k)t}B_{jk}(t).
\]

With raw, uncentered source operators, the equations are linear:

\[
\dot\rho=L\rho,\qquad
\dot x_k=(L+i\omega_k)x_k+\begin{cases}
O_k^\dagger\rho/T & >,\\
\rho O_k^\dagger/T & <,
\end{cases}
\]
\[
\dot Z_{jk}=i(\omega_k-\omega_j)Z_{jk}+\operatorname{Tr}(O_jx_k)/T,\qquad
\dot m_j=-i\omega_jm_j+\operatorname{Tr}(O_j\rho)/T.
\]

Initially `x = Z = m = 0`. All coordinates carry continuously across segment
boundaries: their rotating phases refer to global time zero, so a readout jump
changes a source or derivative without introducing a boundary phase factor.
At exit the physical triangle is
\(Q_{jk}=T^2e^{i(\omega_j-\omega_k)T}Z_{jk}(T)\), and
\(\mu_j=Te^{i\omega_jT}m_j(T)\). Each raw moment is \(Q+Q^\dagger\).
Subtract \(\mu\mu^\dagger\) once, after integrating the whole protocol.
The Hermitian partner supplies the other time half-plane; it is not a PSD
repair. Equal times have zero integration measure for these finite pulses.

Using `u=t/T` makes this a constant homogeneous block system within each
segment. Its semigroup is exactly `expm(G * duration_s/T)` up to floating-point
error. No sampling or adaptive stepping of the GHz Fourier phases occurs.
The block dimension is 30 for two levels/two ports, and 90 for four levels/two
ports. A separate small density exponential makes the reported exit state
independent of the requested readouts and frequency grid.

No moments or eigenvalues are clipped, fitted, renormalized or repaired.
The ordinary floating-point cancellation in raw-minus-mean subtraction is
retained. For independent identical Poisson arrivals of rate `r`, the existing
`poisson_beam_spectrum` accepts this dictionary and gives `r * raw_greater`
and `r * raw_lesser`. In particular, `O=I` has zero internal connected noise,
but retains the raw mean outer product responsible for number noise.

## Verification and scope

Run only the owned tests:

```powershell
python -B -m pytest -q -p no:cacheprovider tests/quantum/test_segmented_qrt.py
```

On 2026-09-14, Python 3.14 / NumPy 2.4.1 / SciPy 1.17.0: **31 passed**.
The original `characteristic_qrt` supplies an independent physical-time
`solve_ivp` comparison, using the callback
`exp(1j*(omega_j-common_omega)*t)*O_j(segment)` for unequal port frequencies.
It integrates centered sources; the new reference integrates raw sources.

| Control | Observed error |
|---|---|
| Two-level, changed H/readouts, nonzero mean, signed unequal frequencies | Maximum relative Frobenius error `8.59e-14` across both connected moments, mean and exit state |
| Reduced 85Rb D1, two segments totaling 0.36 ns, opposite hyperfine-frequency ports | Maximum relative error `4.27e-11` across the same quantities |
| Closed two-level system, 2 microseconds and up to 3.0357 GHz | Maximum per-frequency moment relative error `1.36e-10` against direct Heisenberg operator integrals; mean error `1.65e-12` |
| Driven identity readout | Maximum absolute connected residual `4.45e-16` in the one-second test units |

Additional controls cover segment subdivision, disjoint-position cross
correlations, piecewise complex identity shifts, the analytic damped finite
Lorentzian window, time-unit rescaling and invalid inputs. The Rb comparison
shares only the reduced pump Hamiltonians, radiative reservoirs and readout
operators; it calls no stationary-noise or covariance solver. This is a finite
atomic wavepacket reference for a prescribed piecewise-constant protocol, not
a self-consistent optical-field model or experimental squeezing validation.
The full repository suite is left to the main task as requested.
