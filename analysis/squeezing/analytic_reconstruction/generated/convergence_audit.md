# Reduced-model Floquet and velocity convergence audit

This is an analysis-only audit of `analysis/squeezing/analytic_reconstruction/ref_solver.py`; it does not modify the LaTeX report or production solver.

**Classification warning:** the initial Floquet/velocity gain tables inherit the archived `ref_solver.py` propagation, which uses dressed optical wave numbers and a refractive phase mismatch together. They isolate Floquet/velocity numerics inside that shared implementation; they are **not** predictions of the corrected no-double-count Option-A propagation. A separate corrected Option-A literature-point diagnostic appears later.

## Solver self-checks

| check | max absolute difference |
|---|---:|
| `n_f_1_rho0_max_abs_difference` | 0.000e+00 |
| `n_f_1_rho1_max_abs_difference` | 0.000e+00 |
| `n_f_2_direct_block_rho0_max_abs_difference` | 1.978e-13 |
| `n_f_2_direct_block_rho1_max_abs_difference` | 4.798e-14 |
| `n_f_3_direct_block_rho0_max_abs_difference` | 5.372e-13 |
| `n_f_3_direct_block_rho1_max_abs_difference` | 1.290e-13 |

## Floquet truncation at the common operating point

Common point: $\Delta/2\pi=-1.50$ GHz, $T=110$ C, $\delta/2\pi=-280$ MHz. The velocity grid is 5 m/s to 3 sigma. Phases are modulo 360 degrees.

| N_F | G_s | G_c | G_s-G_c | arg chi_sc (deg) | arg chi_cs (deg) |
|---:|---:|---:|---:|---:|---:|
| 1 | 22.5913 | 21.9687 | 0.622609 | 178.142 | -179.333 |
| 2 | 19.3415 | 19.4472 | -0.105684 | 178.173 | -179.366 |
| 3 | 19.3415 | 19.4472 | -0.105684 | 178.173 | -179.366 |

| change | rel. G_s | rel. G_c | rel. gap | phase sc | phase cs | delta-star shift |
|---|---:|---:|---:|---:|---:|---:|
| 1 to 2 | 16.8023% | 12.966% | 689.122% | 0.0310223 deg | 0.0331477 deg | 5 MHz |
| 2 to 3 | 3.792e-07% | 3.53297e-07% | 4.3872e-06% | 5.02297e-10 deg | 1.23822e-09 deg | 0 MHz |

## Legacy reduced-objective minimizer

`delta_star` below is only the minimizer of the reconstruction's legacy, gain-only squeezing objective (including its gap gate); it is not a commutator-preserving quantum prediction. Scan spacing is 5 MHz.

| N_F | delta_star (MHz) | G_s | G_c | gap | arg chi_sc | arg chi_cs | legacy xi finite |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | -280 | 22.5913 | 21.9687 | 0.622609 | 178.142 | -179.333 | -8.10236 dB |
| 2 | -275 | 8.97871 | 8.4684 | 0.510314 | 178.402 | -179.582 | -7.94704 dB |
| 3 | -275 | 8.97871 | 8.4684 | 0.510314 | 178.402 | -179.582 | -7.94704 dB |

## One-dimensional velocity-step refinement

All rows use N_F=3, cutoff 5 sigma, and the fixed $\delta/2\pi=-280$ MHz point. Errors are relative to the last row.

| dv (m/s) | points | G_s | G_c | gap | G_s err. | gap err. | phase-sc err. |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 10 | 195 | 18.9366 | 19.0243 | -0.087738 | 0.00156072% | 0.0132176% | 7.08312e-09 deg |
| 5 | 389 | 18.9363 | 19.024 | -0.0877264 | 2.65898e-07% | 4.08418e-05% | 3.98035e-09 deg |
| 2.5 | 777 | 18.9363 | 19.024 | -0.0877264 | 9.57873e-08% | 2.26003e-05% | 2.40198e-09 deg |
| 1.25 | 1551 | 18.9363 | 19.024 | -0.0877264 | 0% | 0% | 0 deg |

## One-dimensional velocity-cutoff refinement

All rows use N_F=3, dv=2.5 m/s, and the fixed $\delta/2\pi=-280$ MHz point. Errors are relative to the last row.

| cutoff (sigma) | points | G_s | G_c | gap | G_s err. | gap err. | phase-sc err. |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 2 | 311 | 22.3584 | 22.6171 | -0.258671 | 18.0716% | 194.861% | 0.00780145 deg |
| 2.5 | 389 | 20.0906 | 20.2336 | -0.143031 | 6.09568% | 63.042% | 0.00147851 deg |
| 3 | 467 | 19.3619 | 19.4686 | -0.106622 | 2.24781% | 21.5389% | 0.000198897 deg |
| 3.5 | 545 | 19.1153 | 19.2104 | -0.0951388 | 0.945165% | 8.44943% | 2.02091e-05 deg |
| 4 | 621 | 18.9821 | 19.0717 | -0.0895745 | 0.242127% | 2.10666% | 1.3015e-06 deg |
| 4.5 | 699 | 18.9364 | 19.0241 | -0.0877267 | 0.000643727% | 0.000266004% | 3.07159e-06 deg |
| 5 | 777 | 18.9363 | 19.024 | -0.0877264 | 0% | 0% | 0 deg |

## Corrected Option-A literature-point diagnostic

This section is separate from the archived tables above. It uses bare, frequency-specific optical wave numbers in the susceptibility terms and only vacuum/geometric phase mismatch. No refractive-index contribution is inserted into the mismatch.

Operating point: $\Delta/2\pi=+0.900$ GHz, $\delta/2\pi=-8.000$ MHz, $T=121.0$ C, pump=600 mW, seed=8 uW, $L=12.5$ mm, $\theta=0.32$ deg. The one-dimensional velocity grid has 1967 points ($dv=1.0$ m/s, cutoff 5 sigma). Angular two-photon Doppler broadening is not included.

All rows below are evaluated at the same fixed $\delta/2\pi=-8$ MHz literature point; no detuning optimization is mixed into this table.

| N_F | G_s | G_c | G_s-G_c | arg chi_sc (deg) | arg chi_cs (deg) |
|---:|---:|---:|---:|---:|---:|
| 1 | 5.60949515 | 4.66124517 | 0.948249981 | -179.880518 | 179.797297 |
| 2 | 5.61434133 | 4.66815963 | 0.946181694 | -179.879968 | 179.796859 |
| 3 | 5.61434133 | 4.66815963 | 0.946181694 | -179.879968 | 179.796859 |

### Weak-field reference-amplitude check

The atomic response is still the inherited approximate four-level model, and its finite seed/reference field enters the steady solve. The following N_F=3 check changes only that reference from 2 to 8 uW; it tests numerical weak-field linearity, not microscopic or experimental validity.

| seed (uW) | Omega_s/2pi (MHz) | G_s | G_c | gap | dG_s vs 2uW | dgap vs 2uW |
|---:|---:|---:|---:|---:|---:|---:|
| 2 | 2.074717 | 5.614368372 | 4.668742917 | 0.945625455 | +0.000000% | +0.000000% |
| 4 | 2.934093 | 5.614359421 | 4.668548520 | 0.945810900 | -0.000159% | +0.019611% |
| 8 | 4.149435 | 5.614341328 | 4.668159634 | 0.946181694 | -0.000482% | +0.058822% |

Vacuum/geometric mismatch: $\Delta k_{\rm vac}=246.535124464\,\mathrm{m^{-1}}$.

**M** [m^-1]

| row | column 1 | column 2 |
|---:|---:|---:|
| 1 | `-1.33032779+242.349727i` | `-0.410000095+195.707754i` |
| 2 | `0.693783566-195.680004i` | `0.341352189-171.01377i` |

**T = exp(M L)** [dimensionless]

| row | column 1 | column 2 |
|---:|---:|---:|
| 1 | `-0.376875777+2.33929604i` | `-0.929013535+1.95099446i` |
| 2 | `0.931714213-1.94937638i` | `1.59337708-1.7776828i` |

## Minimum-vacuum mathematical commutator completion

**This is a mathematical dilation only.** It is not microscopic atomic diffusion, is not frequency dependent, and must not be presented as a Langevin-corrected squeezing-spectrum prediction. It additionally assumes the classical two-mode amplitudes have canonical photon-flux normalization.

For $J=\mathrm{diag}(1,-1)$, the eigenvalues of $K=-(MJ+JM^\dagger)$ are 0.189390818 and 3.153969133 m^-1. Both are positive, so $J_f=I_2$ for the displayed eigenfactor.

**K** [m^-1]

| row | column 1 | column 2 |
|---:|---:|---:|
| 1 | `2.66065557` | `-1.10378366+0.0277496366i` |
| 2 | `-1.10378366-0.0277496366i` | `0.682704377` |

**B** [m^-1/2]

| row | column 1 | column 2 |
|---:|---:|---:|
| 1 | `0.177525001` | `-1.62146244` |
| 2 | `0.397210236+0.00998605079i` | `0.724221775+0.0182072736i` |

**D_min = B B^dagger / 2** [m^-1]

| row | column 1 | column 2 |
|---:|---:|---:|
| 1 | `1.33032779` | `-0.55189183+0.0138748183i` |
| 2 | `-0.55189183-0.0138748183i` | `0.341352189` |

**V_out for vacuum input** [dimensionless]

| row | column 1 | column 2 |
|---:|---:|---:|
| 1 | `5.16944552` | `-4.94851817+1.45717846i` |
| 2 | `-4.94851817-1.45717846i` | `5.19900667` |

| commutator/diffusion check | max residual |
|---|---:|
| Bare transfer $TJT^\dagger-J$ | 5.510e-02 |
| Factorization $BJ_fB^\dagger-K$ | 4.441e-16 m^-1 |
| Completed output commutator | 1.332e-15 |
| 200-to-400 point covariance-integral change | 1.672e-16 |

Here `max residual` means the entrywise norm $\max_{ij}|R_{ij}|$.

| bright-seed diagnostic | S_- | dB | classification |
|---|---:|---:|---|
| Bare T only | 0.0872051 | -10.595 | invalid; commutator not restored |
| Minimum-vacuum dilation, unweighted | 0.0922962 | -10.348 | mathematical diagnostic |
| Minimum-vacuum dilation, DC-balanced | 0.0977365 | -10.099 | mathematical diagnostic |
| Ideal Bogoliubov matched to G_s | 0.0977643 | -10.098 | counterfactual benchmark |
| Minimum dilation after external eta=0.8694 | 0.2108423 | -6.760 | external-loss diagnostic |

Against the repository literature benchmark (`README.md`, Sim et al. 85Rb optimum), the corrected Option-A mean-field result is $G_s=5.614$ versus approximately 15.5 (-63.8%). The mathematical dilation after external loss is -6.760 dB versus the reported scale near -7.8 dB. No bandwidth comparison is available because the static model has no spectrum-analyzer frequency.

## Interpretation and API limitations

- The production APIs `gabes.core.floquet_solve(...)` and `gabes.kernels.floquet_chi_grid(...)` are fixed to N_F=1. `gabes.schemes.fwm.chi_matrix_table(...)` selects that fused kernel when Numba is available and exposes no truncation-order argument.
- `gabes.schemes.fwm.compute_spectrum(...)` exposes `velocity_step` and `velocity_cutoff`, but not N_F. Its one-dimensional path calls `gabes.doppler.velocity_grid(...)`, `build_Delta_eff_axis(...)`, and `doppler_average(...)`.
- The N_F=1 result is not Floquet-converged: N_F=2 changes the common-point gains substantially. N_F=2 and N_F=3 agree to the precision reported here.
- Velocity refinement here converges only the existing collinear integral $\Delta_{eff}=\Delta-kv$. The current susceptibility API keeps delta independent of velocity and therefore cannot represent the crossing-angle two-photon Doppler distribution.
- At theta=0.32 deg, sigma_v=193.694 m/s, and lambda=794.979 nm, the omitted one-sigma angular width is 1.361 MHz.
- A correct geometry extension must separate the lab beat frequency from the velocity-shifted atomic two-photon detuning. Reusing the current `floquet_chi_grid` with `delta_eff` would incorrectly shift both because that kernel computes `omega_beat = omega_hf + branch*delta` internally.
- The continued-fraction extension assumes the same periodic Hamiltonian with only +/-1 Fourier couplings. It is exact for that finite truncation, but it does not repair the pump steady-state, quantum-Langevin, or four-level-model limitations.
- The initial/common-point convergence tables retain `ref_solver.py`'s dressed-k plus refractive-mismatch convention. Use them only to diagnose convergence of the archived calculation. The separate Option-A section above supplies the corrected bare-k/vacuum-mismatch literature-point calculation.
