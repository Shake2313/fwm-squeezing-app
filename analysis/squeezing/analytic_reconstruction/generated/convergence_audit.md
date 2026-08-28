# Reduced-model Floquet and velocity convergence audit

This is an analysis-only audit of `analysis/squeezing/analytic_reconstruction/ref_solver.py`; it does not modify the LaTeX report or production solver.

**Classification warning:** the initial Floquet/velocity gain tables inherit the archived `ref_solver.py` propagation, which uses dressed optical wave numbers and a refractive phase mismatch together. They isolate Floquet/velocity numerics inside that shared implementation. Its dissipator also applies the inherited ground-coherence rate as coherence-only damping, without the production thermal transit reload channel. They are **not** predictions of the corrected no-double-count Option-A propagation. A separate corrected Option-A literature-point diagnostic appears later.

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

This section is separate from the archived tables above. It uses bare, frequency-specific optical wave numbers in the susceptibility terms and only vacuum/geometric phase mismatch. No refractive-index contribution is inserted into the mismatch. The trace-normalized rho_ss supplies the manifold population once; the external structural factor is 1/[2(2I+1)].

Operating point: $\Delta/2\pi=+0.900$ GHz, $\delta/2\pi=-8.000$ MHz, $T=121.0$ C, pump=600 mW, seed=8 uW, $L=12.5$ mm, $\theta=0.32$ deg. The one-dimensional velocity grid has 1967 points ($dv=1.0$ m/s, cutoff 5 sigma). This base finite-seed table is one-dimensional; the separate slow two-dimensional reference below includes angular two-photon Doppler broadening.

All rows below are evaluated at the same fixed $\delta/2\pi=-8$ MHz literature point; no detuning optimization is mixed into this table.

| N_F | G_s | G_c power | G_c flux | flux gap | arg chi_sc (deg) | arg chi_cs (deg) |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 401.708332 | 403.862425 | 403.855906 | -2.14757357 | -179.886152 | 179.814904 |
| 2 | 403.535497 | 405.853591 | 405.84704 | -2.31154273 | -179.885623 | 179.814496 |
| 3 | 403.535497 | 405.853591 | 405.84704 | -2.31154273 | -179.885623 | 179.814496 |

### Weak-field reference-amplitude check

The production atomic response is the current four-level finite-Floquet density-matrix model, and its finite seed/reference field enters the steady solve. The following N_F=3 check changes only that reference from 2 to 8 uW; it tests numerical weak-field linearity, not microscopic or experimental validity.

| seed (uW) | Omega_s/2pi (MHz) | G_s | G_c power | flux gap | dG_s vs 2uW | dflux gap vs 2uW |
|---:|---:|---:|---:|---:|---:|---:|
| 2 | 2.074717 | 403.752755516 | 406.102206349 | -2.342895387 | +0.000000% | -0.000000% |
| 4 | 2.934093 | 403.680318836 | 406.019312600 | -2.332439657 | -0.017941% | -0.446274% |
| 8 | 4.149435 | 403.535497293 | 405.853591460 | -2.311542734 | -0.053810% | -1.338201% |

### Pump-only stationary state and infinitesimal response reference

This is a separate slow reference, not the production finite-seed gain path. For the standard minus Raman branch, the explicit rotating frame $U=\exp(+i\Omega_{\rm beat}t|g_2\rangle\langle g_2|)$ makes both pump couplings static and gives $\delta+\Omega_{\rm beat}=\omega_{\rm HF}$. The inherited plus-branch frame is not gauge-equivalent and is therefore reported as unsupported rather than silently mapped.

| pump/reference check | value |
|---|---:|
| max normalized pump residual | 4.002e-17 |
| max pump trace error | 2.887e-15 |
| minimum pump-state eigenvalue | 2.370583480e-03 |
| max trace-zero response residual | 7.312e-17 |
| max response trace error | 6.478e-25 |
| static-frame vs N_F=3 pump Floquet max difference | 8.366e-15 |
| projected DC relative frequency | 0.000e+00 rad/s |
| projected DC residual | 1.794e-15 |

The finite-reference susceptibility approaches the analytic complex-field derivative as the reference amplitude is reduced. The table reports the empirical full-scan error; it does not impose the nominal quadratic law:

| Rabi fraction | seed power (uW) | worst normalized chi error |
|---:|---:|---:|
| 1.000 | 8.000000 | 2.502685924e-03 |
| 0.300 | 0.720000 | 5.093108808e-04 |
| 0.100 | 0.080000 | 9.253710408e-05 |
| 0.030 | 0.007200 | 9.137899826e-06 |

Fitted Rabi-amplitude error slope: 1.748955; pairwise slopes 1.552366, 1.922950. The asymptotic perturbative expectation is 2, but the moving full-scan maximum is reported without relabelling it as exact quadratic convergence.

| N_F at smallest reference amplitude | worst normalized chi error |
|---:|---:|
| 1 | 4.484861781e-03 |
| 2 | 9.137899962e-06 |
| 3 | 9.137899826e-06 |

Direct resolvent and diagonalizable Liouvillian pole/residue expansion:

| analysis frequency (MHz) | max direct/pole difference (s) | eigenvector condition | stationary residue max |
|---:|---:|---:|---:|
| 0.000 | 2.306e-24 | 3.389562 | 1.544e-17 |
| 0.100 | 2.312e-24 | 3.389562 | 1.544e-17 |
| 1.000 | 2.329e-24 | 3.389562 | 1.544e-17 |
| 4.000 | 2.447e-24 | 3.389562 | 1.544e-17 |

### Two-dimensional non-collinear Raman-Doppler reference

This is another separate slow reference, not the production finite-seed scan. It tensor-averages the pump-only trace-zero response over independent Maxwellian $v_z$ and $v_x$. The laboratory optical beat and spectrum-analyzer frequency remain fixed while only the atomic detunings are velocity shifted:

$\Delta_{\rm eff}=\Delta-k_pv_z$,  $\delta_{\rm eff}=\delta+(k_p-k_s\cos\theta)v_z$-k_s\sin\theta\,v_x$.

| Raman-Doppler width budget | value |
|---|---:|
| angle | 0.320 deg |
| analytic total rms | 1.380173564 MHz |
| analytic transverse rms | 1.380161171 MHz |
| analytic axial rms at angle | 5.848714 kHz |
| collinear residual rms | 1.994567 kHz |
| order-40 quadrature rms | 1.380163304 MHz |
| quadrature error | -0.000743% |

Grid refinement uses the same 0.05 MHz refined feature scan and a five-sigma velocity cutoff:

| order/axis | velocity pairs | feature delta (MHz) | feature G_s | feature G_c | G_s err. vs order 40 | shift vs order 40 | runtime (s) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 24 | 576 | 36.150000 | 4517.402502479 | 4585.295223347 | +0.014954% | +0.200000 MHz | 0.508 |
| 32 | 1024 | 36.000000 | 4517.048080097 | 4584.672226523 | +0.007107% | +0.050000 MHz | 0.932 |
| 40 | 1600 | 35.950000 | 4516.727074904 | 4584.259386638 | +0.000000% | +0.000000 MHz | 1.526 |

Cutoff refinement uses order 40 per velocity axis:

| cutoff (sigma) | feature delta (MHz) | feature G_s | G_s err. vs 5 sigma | shift vs 5 sigma |
|---:|---:|---:|---:|---:|
| 4.5 | 35.950000 | 4519.082750675 | +0.052154% | +0.000000 MHz |
| 5.0 | 35.950000 | 4516.727074904 | +0.000000% | +0.000000 MHz |

| nominal -8 MHz comparison | value |
|---|---:|
| one-dimensional G_s | 403.825149966 |
| two-dimensional G_s | 404.304420146 |
| G_s change | +0.118683% |
| max response normalized residual | 2.181e-16 |
| max response trace error | 5.825e-26 |

All grid/cutoff acceptance flags pass. Remaining exclusions are the finite-seed two-dimensional production path, a beam angular distribution, segmentwise pump-state recomputation, and microscopic Langevin diffusion.

Vacuum/geometric mismatch: $\Delta k_{\rm vac}=246.535124464\,\mathrm{m^{-1}}$.
Equal declared Gaussian collected-mode areas: $A_p=A_c=1.710597200e-07\,\mathrm{m^2}$.

**Q (E = Q a)** [field per sqrt(photon/s)]

| row | column 1 | column 2 |
|---:|---:|---:|
| 1 | `3.31753644e-05` | `0` |
| 2 | `0` | `3.31756322e-05` |

**M_field** [m^-1]

| row | column 1 | column 2 |
|---:|---:|---:|
| 1 | `-1.79544422+8.518638i` | `-0.593341782+297.227397i` |
| 2 | `0.962246255-297.203293i` | `0.508739931+60.4654553i` |

**T_field** [dimensionless]

| row | column 1 | column 2 |
|---:|---:|---:|
| 1 | `18.8396938+6.97147291i` | `-8.60352176+18.2180186i` |
| 2 | `8.62545956-18.2059067i` | `17.4837706+10.2248631i` |

**M_canonical = Q^-1 M_field Q** [m^-1]

| row | column 1 | column 2 |
|---:|---:|---:|
| 1 | `-1.79544422+8.518638i` | `-0.593346571+297.229796i` |
| 2 | `0.962238488-297.200894i` | `0.508739931+60.4654553i` |

**T_canonical = Q^-1 T_field Q** [dimensionless]

| row | column 1 | column 2 |
|---:|---:|---:|
| 1 | `18.8396938+6.97147291i` | `-8.6035912+18.2181656i` |
| 2 | `8.62538994-18.2057598i` | `17.4837706+10.2248631i` |

The production Maxwell/Q path is independently reconstructed from the stored reduced susceptibilities using literal SI prefactors and `scipy.linalg.expm`:

| independent reference comparison | max absolute difference |
|---|---:|
| `M_field_max_abs_difference` | 0.000e+00 |
| `T_field_max_abs_difference` | 1.638e-14 |
| `Q_max_abs_difference` | 0.000e+00 |
| `M_canonical_max_abs_difference` | 0.000e+00 |
| `T_canonical_max_abs_difference` | 1.638e-14 |

## Parameter provenance and illustrative sensitivity

These sweeps are machine-readable under `parameter_provenance` and `sensitivity` in the JSON artifact. They are illustrative local tests, not parameter uncertainties or fits.

Current atomic construction: `gabes.schemes.fwm.collisional_atom -> gabes.schemes.fwm.thermal_transit_reset_superoperator + gabes.atoms.double_lambda_rb85`.

| parameter | current value | source | status |
|---|---:|---|---|
| `ell_s` | 0.74 dimensionless | `LITERATURE_ARGS.line_strength_residual` | inherited unrefitted residual; not first-principles |
| `kappa` | 0.1 dimensionless | `gabes.schemes.fwm.HARDENED_PUMP_SCATTER_KAPPA` | uncalibrated technical-noise diagnostic coefficient |
| `gamma_gg_floor_over_2pi_kHz` | 100 kHz | `gabes.constants.GAMMA_GG_2PI` | thermal transit-reset rate; density collision dephasing is added separately in the current atomic path |

| table/block | atomic or numerical solver provenance |
|---|---|
| `floquet_convergence` | atomic_response: analysis.squeezing.analytic_reconstruction.convergence_audit.floquet_solve_truncated; atomic_model: gabes.schemes.fwm.collisional_atom -> gabes.schemes.fwm.thermal_transit_reset_superoperator + gabes.atoms.double_lambda_rb85; propagation: gabes.observables Option-A path, checked against literal SI scipy.linalg.expm reference |
| `seed_reference_linearity` | atomic_response: recomputed N_F=3 for each finite seed reference; atomic_model: gabes.schemes.fwm.collisional_atom -> gabes.schemes.fwm.thermal_transit_reset_superoperator + gabes.atoms.double_lambda_rb85; propagation: same independently checked Option-A path |
| `pump_only_weak_response_reference` | atomic_response: gabes.schemes.fwm.pump_only_weak_response_reference; state: static physical pump frame, trace-one 16x16 null solve; response: two-column complex Nambu derivative; projected trace-zero DC and ordinary finite-frequency resolvent; production_default: False |
| `noncollinear_doppler_reference` | atomic_response: gabes.schemes.fwm.pump_only_weak_response_noncollinear_reference; velocity_average: tensor Gauss-Legendre quadrature of the truncated Maxwell distribution over independent (v_z, v_x); frequency_contract: Omega_beat,lab and Omega_SA fixed; only Delta_eff and delta_eff depend on velocity; production_default: False |
| `ell_s_propagation_only` | atomic_response: stored current N_F=3 reduced susceptibility; atomic_model: gabes.schemes.fwm.collisional_atom -> gabes.schemes.fwm.thermal_transit_reset_superoperator + gabes.atoms.double_lambda_rb85; propagation: Option-A rerun only; atomic response held fixed |
| `gamma_gg_floor_mean_field` | atomic_response: N_F=3 recomputed at 90/100/110 kHz floor; atomic_model: production collisional_atom thermal-reset path with sensitivity-only transit-rate override; propagation: same independently checked Option-A path |
| `kappa_pump_scatter_diagnostic` | solver: gabes.schemes.fwm._pump_scatter_noise -> gabes.schemes.absorption._hyperfine_alpha; atomic_covariance: not computed |
| `algebraic_dilation_fixture` | solver: static eigenfactorization plus Gauss-Legendre integral; atomic_diffusion: not computed; physical_bound: False |

### Propagation-only ell_s sweep

Illustrative +/-10% sweep; fixed stored/current n_f=3 reduced chi; not a fit or uncertainty interval.

| ell_s | effective line strength | G_s | G_c power | G_c flux | flux gap |
|---:|---:|---:|---:|---:|---:|
| 0.666 | 0.055500000 | 188.447070742 | 189.053808218 | 189.050756445 | -0.603685703 |
| 0.740 | 0.061666667 | 403.535497293 | 405.853591460 | 405.847040027 | -2.311542734 |
| 0.814 | 0.067833333 | 857.835745355 | 863.690860021 | 863.676918017 | -5.841172662 |

### Ground-coherence-floor mean-field sensitivity

Actual n_f=3 atomic-response recomputation through the current thermal-transit-reset path with all non-gamma inputs fixed; central row reuses the identical current solve.

| gamma_gg transit floor / 2pi (kHz) | collision dephasing / 2pi (kHz) | total coherence rate / 2pi (kHz) | G_s | G_c power | G_c flux | flux gap | source |
|---:|---:|---:|---:|---:|---:|---:|---|
| 90.0 | 2.881854 | 92.881854 | 440.249924433 | 442.733154620 | 442.726007865 | -2.476083432 | recomputed N_F=3 atomic solve |
| 100.0 | 2.881854 | 102.881854 | 403.535497293 | 405.853591460 | 405.847040027 | -2.311542734 | reused current N_F=3 atomic solve |
| 110.0 | 2.881854 | 112.881854 | 370.625582907 | 372.789955981 | 372.783938273 | -2.158355366 | recomputed N_F=3 atomic solve |

### Pump-scatter kappa diagnostic

`N_ps = kappa * (1 - exp(-OD_pump))`. SQL-normalized technical pump-scatter diagnostic only; not physical squeezing or microscopic diffusion. At kappa=0.1 and OD_pump=0.149609275, dN_ps/dkappa=0.138955658 and N_ps=0.0138955658.

## Algebraic commutator/diffusion fixture

**This is an algebraic dilation fixture with one chosen vacuum-reservoir covariance.** It is not a Caves bound, microscopic atomic diffusion, or frequency dependent, and must not be presented as a Langevin-corrected squeezing-spectrum prediction. The canonical basis is constructed explicitly above from the declared frequencies and mode areas.

For $J=\mathrm{diag}(1,-1)$, the eigenvalues of $K=-(MJ+JM^\dagger)$ are 0.285202198 and 4.323166101 m^-1. Both are positive, so $J_f=I_2$ for the displayed eigenfactor.

**K** [m^-1]

| row | column 1 | column 2 |
|---:|---:|---:|
| 1 | `3.59088844` | `-1.55558506+0.0289018754i` |
| 2 | `-1.55558506-0.0289018754i` | `1.01747986` |

**B** [m^-1/2]

| row | column 1 | column 2 |
|---:|---:|---:|
| 1 | `0.227422334` | `-1.88126753` |
| 2 | `0.48311563+0.00897601043i` | `0.885284159+0.0164480703i` |

**D_vacuum_fixture = B B^dagger / 2** [m^-1]

| row | column 1 | column 2 |
|---:|---:|---:|
| 1 | `1.79544422` | `-0.77779253+0.0144509377i` |
| 2 | `-0.77779253-0.0144509377i` | `0.508739931` |

**V_out for vacuum input** [dimensionless]

| row | column 1 | column 2 |
|---:|---:|---:|
| 1 | `406.423341+1.17624413e-14i` | `35.8550356+406.492771i` |
| 2 | `35.8550356-406.492771i` | `409.73006+5.65169026e-15i` |

| commutator/diffusion check | max residual |
|---|---:|
| Bare transfer $TJT^\dagger-J$ | 3.388e+00 |
| Factorization $BJ_fB^\dagger-K$ | 4.441e-16 m^-1 |
| Completed output commutator | 3.011e-13 |
| 200-to-400 point covariance-integral change | 9.948e-14 |

Here `max residual` means the entrywise norm $\max_{ij}|R_{ij}|$.

| bright-seed diagnostic | S_- | dB | classification |
|---|---:|---:|---|
| Bare T only | 0.0200977 | -16.969 | invalid; commutator not restored |
| Algebraic vacuum fixture, unweighted | 0.0241363 | -16.173 | algebraic diagnostic |
| Algebraic vacuum fixture, DC-balanced | 0.0052331 | -22.812 | algebraic diagnostic |
| Ideal Bogoliubov matched to G_s | 0.0012406 | -29.064 | counterfactual benchmark |
| Algebraic fixture after external eta=0.8694 | 0.1515841 | -8.193 | external-loss diagnostic |

Against the repository literature benchmark (`README.md`, Sim et al. 85Rb optimum), the corrected Option-A mean-field result is $G_s=403.535$ versus approximately 15.5 (2503.5%). The algebraic fixture diagnostic after external loss is -8.193 dB versus the reported scale near -7.8 dB. No bandwidth comparison is available because the static model has no spectrum-analyzer frequency.

## Interpretation and API limitations

- Production exposes `gabes.core.floquet_solve_truncated(...)` and an independently assembled dense-block reference for arbitrary finite N_F. `gabes.kernels.floquet_chi_grid(...)` and `gabes.schemes.fwm.chi_matrix_table(...)` carry the same order argument; tests pin compiled/reference parity at N_F=1,2,3.
- `gabes.schemes.fwm.compute_spectrum(...)` defaults to N_F=3 and compares the entire reported scan with N_F=2 using complex-response/transfer, gain, wrapped-phase, and optimum-shift criteria. It also exposes `velocity_step` and `velocity_cutoff`. Its one-dimensional path calls `gabes.doppler.velocity_grid(...)`, `build_Delta_eff_axis(...)`, and `doppler_average(...)`.
- The archived N_F=1 result is not Floquet-converged: N_F=2 changes the common-point gains substantially. N_F=2 and N_F=3 agree to the precision reported here; production no longer promotes a single-point comparison.
- The production velocity refinement still converges only the collinear integral $\Delta_{eff}=\Delta-kv$. The opt-in two-dimensional reference above represents the crossing-angle Raman-Doppler distribution without changing the laboratory beat frequency.
- At theta=0.32 deg, sigma_v=193.694 m/s, and lambda=794.979 nm, the one-sigma angular width is 1.361 MHz; the separate reference resolves this width explicitly.
- The two-dimensional reference separates the lab beat frequency from the velocity-shifted atomic two-photon detuning. Reusing the production `floquet_chi_grid` with `delta_eff` remains invalid because that kernel also derives `omega_beat = omega_hf + branch*delta` from the same value.
- The continued-fraction extension assumes the same periodic Hamiltonian with only +/-1 Fourier couplings. It is exact for that finite truncation, but it does not repair the pump steady-state, quantum-Langevin, or four-level-model limitations.
- The initial/common-point convergence tables retain `ref_solver.py`'s dressed-k plus refractive-mismatch convention. Use them only to diagnose convergence of the archived calculation. The separate Option-A section above supplies the corrected bare-k/vacuum-mismatch literature-point calculation.
