# Change log: squeezing_analytic_reconstruction_v1.tex

This log maps the requested corrections to the pre-revision report. Original
line numbers refer to the version in git before this revision.

| Request | Original location | Revised location and disposition |
|---|---|---|
| 1. Bosonic commutators | Noise section, old equations io and S, lines 597–628 | Part V replaces the gain-only two-mode map and cov_0 formula with quantum-Langevin propagation, local/integrated commutator identities, and a covariance diffusion integral. The minimum-vacuum dilation passes numerically; microscopic atomic diffusion remains pending. |
| 2. Analysis frequency | Beat-frequency subsection and scalar noise result, lines 207–214 and 597–664 | Parts I, III, V, and VI distinguish Omega_beat from Omega_SA, Fourier transform fluctuations, and define the spectrum S_-(Omega_SA). Static-model bandwidth limitations are explicit. |
| 3. Pump steady state | Old carrier rate equations and Sylvester block, lines 255–305 | Part II uses the trace-constrained pump Liouvillian null space. The former split is retained only as a labeled controlled approximation with assumptions. |
| 4. Floquet truncation | Old Floquet expansion, lines 216–247 | Parts III and VII report N_F = 1, 2, 3 for both gains, the gain gap, cross-response phases, and delta-star. N_F = 1 is classified qualitative and unconverged. |
| 5. Non-collinear Doppler | Old one-dimensional Doppler section, lines 455–494 | Part III adds vector one- and two-photon detunings, a two-dimensional velocity average, the zero-angle limit, and the 1.36 MHz angular Raman width. Existing numerical data are explicitly only a one-dimensional convergence test. |
| 6. Operating points | Old input and summary sections, lines 113–188 and 667–742 | Part I classifies the negative-detuning point as an archived reduced-model optimum. Parts IV and VII evaluate the documented +0.9 GHz, -8 MHz, 121 C experimental point and report discrepancies. |
| 7. Validation claims | Abstract, summary, and error budget, lines 55–67 and 667–718 | Part VII separates algebraic, implementation, truncation, microscopic, and experimental validation. Shared-model agreement is no longer physical validation. |
| 8. Dimensions and notation | Old macroscopic coupling, Voigt, phase matching, and mixed-frequency notation | Parts I, III, and VII correct the time dimension of reduced susceptibility, remove the extra 1/k in the Voigt identity, use omega_c = 2 omega_p - omega_s, and enforce angular-frequency notation. |
| 9. Dispersion double counting | Old propagation equations Mentries and phasematch, lines 498–524 | Part IV selects Option A: diagonal complex susceptibilities carry medium dispersion; Delta k_z contains only vacuum/geometric mismatch. |
| 10. Absorption double counting | Old arm-OD and detection-noise sections, lines 614–645 | Parts V and VI keep distributed atomic absorption with Langevin reservoirs and apply only post-cell path/detector loss afterward. |
| 11. Empirical parameters | Old ground decoherence, line factor, pump noise, and parameter summary | Appendix B gives definition, source, class, range/fitting status, and gain/noise sensitivity for ell_s, kappa, gamma_gg, and external losses. |
| 12. Four-level limitations | Old four-level and Clebsch-Gordan definitions, lines 88–106 | Appendices C and D add the normalization ledger, derive the stated averaged factors by definition, flag possible repeated weights, and list omitted Zeeman/polarization/noise physics. |
| 13. Valid retained results | Old mean-field, Faddeeva, propagation, and gain-gap sections | Parts III and IV retain the corrected base Faddeeva integral, classical response chain, constant-matrix exponential, and classical gain identity with explicit assumptions. |
| 14. Organization | Whole report | Reorganized into the seven requested parts plus audit appendices. |
| 15. Limiting cases | Previously absent | Part VII adds no-coupling, ideal-amplifier, bright-seed, external-loss, zero-angle, and zero-dissipation checks. |
| 16. Title and conclusion | Original title and summary, lines 55–75 and 667–726 | Title is restricted to mean-field reconstruction. The conclusion separates analytic, matrix-solved, inherited, phenomenological, experimental, and unvalidated content and invalidates the physical status of the old -8.10 dB and -15.62 dB values. |
| 17. Deliverables | Previously absent as an auditable set | The report now contains the assumptions, provenance, dimensional, frequency, convergence, commutator, comparison, limitations, and acceptance tables. Reproducible numerical artifacts are under analysis/squeezing/analytic_reconstruction/generated. |

The revision is intentionally provisional with respect to physical squeezing:
the current repository does not contain a microscopic frequency-dependent
atomic diffusion matrix, a full non-collinear velocity implementation, or
measured detector transfer functions. The report marks those requirements
pending instead of substituting a gain-only squeezing formula.
