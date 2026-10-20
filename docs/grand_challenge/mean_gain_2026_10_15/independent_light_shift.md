# ⟦AI⟧ Independent pump-light-shift and weak-response control

2026-10-15. Diagnostic evidence, not an adopted apparatus correction. Source:
[reference script](independent_light_shift.py), [numbers](independent_light_shift.json).
Reproduce from the repository root:
`python -m docs.grand_challenge.mean_gain_2026_10_15.independent_light_shift`.
No production source edits or ensemble/transport batch was performed by this audit.

**Sign and magnitude within the existing four-level model.** In angular units,
the physical pump frame is

\[
H/\hbar=\operatorname{diag}(0,\omega_g,-\Delta-\omega_e,-\Delta)+V,
\quad V_{ge}=V_{eg}=\Omega_g c_{ge}/2,
\quad c_{ge}^2=3C_{FF'}^2.
\]

Here \(\Delta=\omega_p-\omega_{F2\to F'3}\),
\(g_1=F2\), \(g_2=F3\), and
\(\delta=\omega_s-(\omega_p-\omega_g)\) for branch −1.
At second order in the pump,

\[
s_1=\frac{\Omega_A^2}{4}\left[
\frac{c_{12}^2}{\Delta+\omega_e}+\frac{c_{13}^2}{\Delta}\right],\qquad
s_2=\frac{\Omega_B^2}{4}\left[
\frac{c_{22}^2}{\Delta+\omega_e+\omega_g}
+\frac{c_{23}^2}{\Delta+\omega_g}\right].
\]

The shifted ground gap is \(\omega_g+s_2-s_1\), hence the Raman probe
coordinate is **\(\delta_R=s_1-s_2>0\)** at the stated blue detuning.
Laboratory beat detuning \((\omega_p-\omega_s)-\omega_g\) is \(-\delta\).
The sign does not establish that a measured gain maximum equals this pole.

At pump 294.8 mW, Gaussian 1/e² radius 530 µm, Δ/2π=1 GHz, v=0:
Ω/2π=495.958975 MHz. The second-order ground shifts are
96.441471 and 33.929626 MHz; their difference is 62.511844 MHz.
Exact diagonalization gives 86.171065 and 33.901128 MHz and a differential
shift **52.269937 MHz**. The dissipative Raman pole is **52.269382 MHz**,
HWHM 0.357507 MHz at 386.15 K with the production reset/collision model.
The second-order result is already 10.242 MHz above exact diagonalization;
it is inadequate for precision peak inference at this pump strength.

The independently constructed matrix uses direct matrix-action commutators,
explicit radiative jumps, dephasing and thermal reset. Its stationary state
and complex response use an independent trace constraint and a 15-dimensional
trace-zero QR basis. Against production, generator relative error is
8.22×10⁻¹⁷, state maximum error 9.11×10⁻¹⁵ and response relative error
3.09×10⁻¹⁴. The four-level response contains the pump light shift; no missing
light-shift term, extra factor two or MHz/rad/s conversion error was found here.

**Velocity and finite-seed controls.** At v=−400,−200,0,+200,+400 m/s the
same Raman poles are 30.963,39.699,52.269,71.271,101.732 MHz. A Doppler-averaged,
spatially propagated power-gain optimum depends on the complex residues,
absorption, phase mismatch, finite seed and depletion. Therefore the reported
ensemble maximum δ=23.7 MHz cannot be corrected by subtracting the v=0
52.269 MHz pole; the experimental −2.2 MHz versus model −23.7 MHz beat-axis
gap remains unresolved by this local control.

For six declared detunings and Δ_eff/2π=750,1000,1250 MHz, finite N_F=3
responses converge to the independent infinitesimal reference as seed power
tends to zero. Relative Frobenius errors at the sampled v=0 Raman-pole
detuning are 86.5%,17.0%,2.03%,0.207%,0.0207% for reference seed powers
3.918,0.1,0.01,0.001,0.0001 µW. The final two errors decrease proportionally
to power, as expected from a cubic field correction divided by seed amplitude.
At δ=23.7 MHz, the same three-velocity error is only 0.197% at 3.918 µW.
Thus **0.1 µW is not a globally certified weak-seed setting**, and the local
error near an individual pole does not quantify the ensemble gain error.
Zero seed cannot be passed to the finite-ratio production extractor, which
divides by reference Rabi frequency; the pump-only trace-zero response gives
the well-defined zero-seed limit directly.

**CG reduction and patch scope.** The existing signed dipole reference
[uncoupled_d1.py](../../../analysis/grand_challenge/reference/uncoupled_d1.py)
gives uniform-manifold per-atom strengths
\(d_J^2\,3C_{FF'}^2/(2F+1)\). The production saturation-intensity Rabi
has implicit base dipole approximately \(d_J/\sqrt3\), so its effective
four-level squared drive dipole is approximately \(d_J^2 C_{FF'}^2\).
The ratio to the uniform-manifold mean is consequently 5/3 for F2 and 7/3
for F3. The pre-existing
[normalization audit](../normalization_derivation.md) already records that
the strong-pump four-level reduction is unproved. A scalar correction does
not restore optical pumping, signed loop phases or polarization interference.

As a diagnostic only, exact diagonalization of the 24-level Hamiltonian with
the signed I⊗J dipoles, linear q=0 pump, B=0 and v=0 gives same-m Raman
shifts 37.330,42.493,44.202,42.493,37.330 MHz for m=−2,…,+2.
This has no probe-polarization/dissipator/gain solve and is not substituted
into production. It demonstrates that agreement of a lumped CF² sum does
not certify its light-shift spectrum.

A separate provable small bookkeeping defect is the production optical
carrier origin: D_GHz is referenced to F2→F′3 while `OMEGA_D1` is the centroid.
The missing offset is \((7\omega_g+5\omega_e)/12\), or 1.921502256 GHz.
This changes optical wave numbers/photon energies by about 5.1×10⁻⁶; it
does not move the already transition-referenced atomic detuning by 1.9215 GHz
or explain a 21.5 MHz Raman mismatch. Production patch ownership remains with
the parent review.

**Head-manager checkpoint.** Local Raman-sign/unit/reference parity passes;
strong-pump model validity and measured gain agreement remain open. The next
discriminating comparison is an ensemble weak-reference/finite-reference gain
comparison at the same input and output reference planes, not an added fitted
light-shift offset. The sealed GC microscopic/thermal results are untouched.
