"""
Phase-4 checks for Zeeman magneto-optics.

  CG / manifold     : angular-momentum plumbing stays normalized.
  Polarized Hanle   : QWP angle switches a paraffin-cell dip into a peak.
  Two-region Ramsey : wall coherence narrows the central feature.
  Buffer cell       : single-region buffer mode keeps a broad Hanle feature.
  NMOR              : rotation remains antisymmetric around B = 0.

    python tests/test_magneto.py    # or: pytest tests/test_magneto.py
"""
import sys
from pathlib import Path
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from gabes import atoms, constants, core, observables, schemes, zeeman  # noqa: E402
from gabes.constants import GAMMA  # noqa: E402
from gabes.schemes.magneto import MagnetoScheme  # noqa: E402

cg = zeeman.clebsch_gordan


def test_clebsch_gordan_known_values():
    assert abs(cg(1, 0, 1, 0, 2, 0) - np.sqrt(2 / 3)) < 1e-9
    assert abs(cg(1, 1, 1, -1, 2, 0) - np.sqrt(1 / 6)) < 1e-9
    assert abs(cg(1, 1, 1, -1, 0, 0) - 1 / np.sqrt(3)) < 1e-9
    assert abs(cg(0.5, 0.5, 0.5, -0.5, 1, 0) - 1 / np.sqrt(2)) < 1e-9
    s = sum(cg(2, m1, 1, 1 - m1, 2, 1) ** 2 for m1 in range(-2, 3))
    assert abs(s - 1.0) < 1e-9


def test_manifold_emission_normalised():
    # Spontaneous emission now lives in polarization-grouped jump operators Σ_q.
    # Σ_q^†Σ_q summed over q must give Γ on every excited level (total decay = Γ),
    # and the per-channel `decay` list is reserved for incoherent transit reload.
    atom = zeeman.zeeman_manifold(2, 3)
    assert atom.n_levels == 12
    assert atom.decay == ()                       # no transit_rate => no rate channels
    T = sum(op.conj().T @ op for op in atom.emission_ops)
    for e in atom.excited:
        assert abs(T[e, e].real / GAMMA - 1.0) < 1e-9
        assert abs(T[e, e].imag) < 1e-12


def test_manifold_emission_transfers_ground_coherence():
    # The grouped Σ_q (unlike one jump per channel) must source a ground Zeeman
    # coherence |g⟩⟨g'| from an excited coherence ρ_{ee'}. Check the dissipator has
    # a nonzero matrix element coupling an excited coherence into a ground one.
    atom = zeeman.zeeman_manifold(2, 1)
    n = atom.n_levels
    ng = len(atom.ground)
    # pick two ground levels and the two excited levels reached from them by the
    # same emitted polarization; their ground coherence must be fed by TOC.
    g0, g1 = 0, 2                                  # m_g = -2, 0
    e0, e1 = ng + 0, ng + 2                        # m_e = -1, +1  (q=+1 and q=+1)
    src = atom.rho_index(e0, e1)                   # excited coherence ρ_{e0 e1}
    dst = atom.rho_index(g0, g1)                   # ground coherence ρ_{g0 g1}
    assert abs(atom.lindblad[dst, src]) > 1e-6


def test_elastic_optical_dephasing_preserves_lifetime_branching_and_toc():
    gamma_natural = 2 * np.pi * 5.746e6
    gamma_elastic = 2 * np.pi * 39.1e6
    base = zeeman.zeeman_manifold(
        2, 1, gamma=gamma_natural, gamma_gg=0.0,
        transit_rate=0.0, optical_dephasing=0.0)
    broadened = zeeman.zeeman_manifold(
        2, 1, gamma=gamma_natural, gamma_gg=0.0,
        transit_rate=0.0, optical_dephasing=gamma_elastic)

    # Elastic broadening cannot alter natural lifetime or polarization branching.
    assert base.collapse_ops == ()
    assert len(broadened.collapse_ops) == 1
    for op_base, op_broad in zip(base.emission_ops, broadened.emission_ops):
        assert np.array_equal(op_base, op_broad)
    total_decay = sum(op.conj().T @ op for op in broadened.emission_ops)
    for excited in broadened.excited:
        assert np.isclose(total_decay[excited, excited], gamma_natural,
                          rtol=1e-12, atol=1e-12)

    delta = broadened.lindblad - base.lindblad
    n = base.n_levels
    # Pure optical dephasing leaves every population equation unchanged.
    for level in range(n):
        population = base.rho_index(level, level)
        assert np.allclose(delta[:, population], 0.0, rtol=0.0, atol=0.0)
    # It adds exactly half the collisional FWHM to each optical coherence.
    optical = base.rho_index(base.ground[0], base.excited[0])
    assert np.isclose(delta[optical, optical], -gamma_elastic,
                      rtol=1e-14, atol=1e-9)

    # The pressure-independent spontaneous-emission TOC source is unchanged.
    ng = len(base.ground)
    src = base.rho_index(ng + 0, ng + 2)
    dst = base.rho_index(0, 2)
    assert base.lindblad[dst, src] != 0
    assert broadened.lindblad[dst, src] == base.lindblad[dst, src]


def test_transit_reload_and_collisional_dephasing_are_distinct_channels():
    gamma_natural = 2 * np.pi * 5.746e6
    gamma_transit = 2 * np.pi * 20e3
    gamma_depol = 2 * np.pi * 2e3
    base = zeeman.zeeman_manifold(
        1, 1, gamma=gamma_natural, gamma_gg=0.0, transit_rate=0.0)
    transit = zeeman.zeeman_manifold(
        1, 1, gamma=gamma_natural, gamma_gg=0.0,
        transit_rate=gamma_transit)
    depol = zeeman.zeeman_manifold(
        1, 1, gamma=gamma_natural, gamma_gg=gamma_depol,
        transit_rate=0.0)
    combined = zeeman.zeeman_manifold(
        1, 1, gamma=gamma_natural, gamma_gg=gamma_depol,
        transit_rate=gamma_transit)

    assert base.decay == () and base.dephasing == () and base.collapse_ops == ()
    assert transit.dephasing == () and transit.collapse_ops == ()
    assert depol.decay == () and depol.dephasing == ()
    assert len(depol.collapse_ops) == len(depol.ground)

    n = base.n_levels
    rho = np.zeros((n, n), dtype=complex)
    rho[base.ground[0], base.ground[0]] = 0.55
    rho[base.ground[1], base.ground[1]] = 0.25
    rho[base.excited[0], base.excited[0]] = 0.20
    rho[base.ground[0], base.ground[1]] = 0.03j
    rho[base.ground[1], base.ground[0]] = -0.03j

    delta_transit = transit.lindblad - base.lindblad
    actual = (delta_transit @ rho.reshape(-1)).reshape(n, n)
    reservoir = np.zeros((n, n), dtype=complex)
    for ground in base.ground:
        reservoir[ground, ground] = 1.0 / len(base.ground)
    expected = gamma_transit * (np.trace(rho) * reservoir - rho)
    assert np.allclose(actual, expected, rtol=1e-12, atol=1e-9)

    delta_depol = depol.lindblad - base.lindblad
    for level in range(n):
        population = base.rho_index(level, level)
        assert np.allclose(delta_depol[:, population], 0.0,
                           rtol=0.0, atol=0.0)
    coherence = base.rho_index(base.ground[0], base.ground[1])
    assert np.isclose(delta_depol[coherence, coherence], -gamma_depol,
                      rtol=1e-14, atol=1e-9)
    optical = base.rho_index(base.ground[0], base.excited[0])
    assert np.isclose(delta_depol[optical, optical], -gamma_depol / 2.0,
                      rtol=1e-14, atol=1e-9)

    # Transit already damps a departing coherence once; depolarization adds once.
    delta_combined = combined.lindblad - base.lindblad
    assert np.isclose(delta_combined[coherence, coherence],
                      -(gamma_transit + gamma_depol),
                      rtol=1e-14, atol=1e-9)


def test_angular_momentum_commutator():
    Fx, Fy, Fz = zeeman.angular_momentum_matrices(2)
    assert np.allclose(Fx @ Fy - Fy @ Fx, 1j * Fz)


def _fast_defaults(sc):
    p = sc.defaults()
    p.update(scan_points=81, velocity_classes=1, doppler="off")
    return p


def _alpha_rot(raw, p):
    xprobe = observables.chi_phys(raw["chi_probe"], raw["N_eff"],
                                  dipole=raw["dipole"], line_strength=p["line_strength"])
    xp = observables.chi_phys(raw["chi_p"], raw["N_eff"],
                              dipole=raw["dipole"], line_strength=p["line_strength"])
    xm = observables.chi_phys(raw["chi_m"], raw["N_eff"],
                              dipole=raw["dipole"], line_strength=p["line_strength"])
    x = raw["b_ut"]
    alpha = raw["k_vec"] * np.imag(xprobe)
    rot = 0.25 * raw["k_vec"] * p["cell_mm"] * 1e-3 * np.real(xp - xm)
    return x, alpha, rot


def _feature_amp(x, alpha):
    ic = int(np.argmin(np.abs(x)))
    bg = 0.5 * (alpha[0] + alpha[-1])
    return alpha[ic] - bg


def _central_halfwidth(x, alpha):
    ic = int(np.argmin(np.abs(x)))
    bg = 0.5 * (alpha[0] + alpha[-1])
    amp = alpha[ic] - bg
    target = bg + 0.5 * amp
    sign = np.sign(amp)
    right = np.arange(ic, len(x))
    vals = (alpha[right] - target) * sign
    below = np.where(vals <= 0)[0]
    if below.size == 0 or below[0] == 0:
        return np.nan
    j = right[below[0]]
    i = j - 1
    t = (target - alpha[i]) / (alpha[j] - alpha[i])
    return abs((x[i] + t * (x[j] - x[i])) - x[ic])


def test_default_is_87rb_d1_paraffin_polarized_hanle():
    sc = schemes.get("magneto")
    p = _fast_defaults(sc)
    raw = sc.compute(p)
    assert raw["isotope"] == "87Rb"
    assert raw["line"] == "D1"
    assert raw["cell_type"] == "Paraffin coated cell"
    assert raw["Fg"] == 2 and raw["Fe"] == 1
    assert raw["valid"] is True
    assert abs(raw["gFg"] - 0.5) < 0.01
    assert isinstance(raw["light_region_occupation"], float)
    assert np.isclose(raw["light_region_occupation"], 1 / 81, rtol=1e-9)
    assert np.isclose(raw["light_state_trace"], 1.0, rtol=0.0, atol=1e-12)


def _two_level_region_fixture():
    atom = atoms.two_level(gamma=2 * np.pi * 5.75e6)
    omega = 2 * np.pi * 0.9e6
    h_light = np.array([[0.0, omega / 2], [omega / 2, 0.0]], dtype=complex)
    h_dark = np.zeros((2, 2), dtype=complex)
    light = core.build_liouvillian(h_light, atom)[None, :, :]
    dark = core.build_liouvillian(h_dark, atom)[None, :, :]
    return atom, light, dark


def test_two_region_returns_conditional_state_and_separate_occupation():
    atom, light, dark = _two_level_region_fixture()
    gamma_out = 2 * np.pi * 80e3
    gamma_in = 2 * np.pi * 1e3
    rho, occupation = MagnetoScheme._steady_state_two_region(
        light, dark, np.array([0.0, 2 * np.pi * 0.2e6]), atom.S_v,
        atom.n_levels, gamma_out, gamma_in)

    expected = gamma_in / (gamma_in + gamma_out)
    assert np.allclose(occupation, expected, rtol=1e-9, atol=0.0)
    assert np.allclose(np.trace(rho, axis1=-2, axis2=-1), 1.0,
                       rtol=0.0, atol=1e-12)


def test_zero_light_exit_recovers_one_region_light_state():
    atom, light, dark = _two_level_region_fixture()
    deff = np.array([0.0, 2 * np.pi * 0.2e6])
    rho_two, occupation = MagnetoScheme._steady_state_two_region(
        light, dark, deff, atom.S_v, atom.n_levels,
        gamma_out=0.0, gamma_in=2 * np.pi * 1e3)
    rho_one = MagnetoScheme._steady_state_buffer(
        light, deff, atom.S_v, atom.n_levels)

    scale = max(float(np.max(np.abs(rho_one))), 1e-300)
    assert np.max(np.abs(rho_two - rho_one)) / scale <= 1e-9
    assert np.allclose(occupation, 1.0, rtol=1e-9, atol=0.0)


def test_local_density_enters_susceptibility_exactly_once():
    sc = schemes.get("magneto")
    p = _fast_defaults(sc)
    p.update(scan_points=3)
    density = 2.5e15
    with patch("gabes.schemes.magneto.species.number_density", return_value=density):
        raw_one = sc.compute(p)
    with patch("gabes.schemes.magneto.species.number_density", return_value=2 * density):
        raw_two = sc.compute(p)

    assert np.allclose(raw_one["chi_probe"], raw_two["chi_probe"],
                       rtol=0.0, atol=0.0)
    x_one = observables.chi_phys(raw_one["chi_probe"], raw_one["N_eff"],
                                 dipole=raw_one["dipole"],
                                 line_strength=p["line_strength"])
    x_two = observables.chi_phys(raw_two["chi_probe"], raw_two["N_eff"],
                                 dipole=raw_two["dipole"],
                                 line_strength=p["line_strength"])
    assert np.allclose(x_two, 2 * x_one, rtol=1e-12, atol=0.0)


def test_normalization_convention_and_absolute_scale_status_are_reported():
    sc = schemes.get("magneto")
    p = _fast_defaults(sc)
    raw = sc.compute(p)
    view = sc.headless_observables(raw, p)
    metrics = {metric["label"]: metric for metric in view["metrics"]}
    derived = next(table["markdown"] for table in view["tables"]
                   if table["title"] == "Derived quantities")

    assert metrics["Absolute contrast status"]["kind"] == "status"
    assert metrics["Absolute contrast status"]["value"] == "external validation required"
    assert "trace-one conditional light-region state" in derived
    assert "diagnostic only" in derived
    assert "local vapor density N, applied once" in derived

    nmor_params = dict(p, signal_type="NMOR rotation")
    nmor_view = sc.headless_observables(raw, nmor_params)
    nmor_metrics = {metric["label"]: metric for metric in nmor_view["metrics"]}
    assert nmor_metrics["Absolute rotation/slope status"]["kind"] == "status"
    assert (nmor_metrics["Absolute rotation/slope status"]["value"]
            == "external validation required")


def test_paraffin_linear_qwp_gives_zero_field_dip():
    sc = schemes.get("magneto")
    p = _fast_defaults(sc)
    p.update(cell_type="Paraffin coated cell", qwp_deg=0.0,
             residual_transverse_b_ut=0.05)
    x, alpha, _ = _alpha_rot(sc.compute(p), p)
    assert _feature_amp(x, alpha) < 0


def test_paraffin_circular_qwp_gives_zero_field_peak():
    sc = schemes.get("magneto")
    p = _fast_defaults(sc)
    p.update(cell_type="Paraffin coated cell", qwp_deg=45.0,
             residual_transverse_b_ut=0.08)
    x, alpha, _ = _alpha_rot(sc.compute(p), p)
    assert _feature_amp(x, alpha) > 0


def test_paraffin_wall_coherence_narrows_central_feature():
    sc = schemes.get("magneto")
    base = _fast_defaults(sc)
    base.update(cell_type="Paraffin coated cell", qwp_deg=0.0,
                residual_transverse_b_ut=0.05, b_max_ut=0.8, scan_points=201)
    short = dict(base, wall_coherence_ms=0.05)
    long = dict(base, wall_coherence_ms=10.0)
    x_s, a_s, _ = _alpha_rot(sc.compute(short), short)
    x_l, a_l, _ = _alpha_rot(sc.compute(long), long)
    assert _central_halfwidth(x_l, a_l) < _central_halfwidth(x_s, a_s)
    assert abs(_feature_amp(x_l, a_l)) > abs(_feature_amp(x_s, a_s))


def test_buffer_mode_has_single_broad_hanle_and_ground_relaxation_broadens():
    sc = schemes.get("magneto")
    base = _fast_defaults(sc)
    base.update(cell_type="Buffer gas cell", qwp_deg=0.0, b_max_ut=120.0,
                ne_pressure_torr=20.0, collisional_depol_khz=0.0)
    low = dict(base, buffer_ground_relax_khz=5.0)
    high = dict(base, buffer_ground_relax_khz=80.0)
    raw_low = sc.compute(low)
    raw_high = sc.compute(high)
    x_l, a_l, _ = _alpha_rot(raw_low, low)
    x_h, a_h, _ = _alpha_rot(raw_high, high)
    assert raw_low["buffer_gamma"] == constants.neon_buffer_broadening(20.0)
    assert _feature_amp(x_l, a_l) < 0
    assert _central_halfwidth(x_h, a_h) > _central_halfwidth(x_l, a_l)


def test_buffer_rate_swap_changes_susceptibility_and_reports_split_rates():
    sc = schemes.get("magneto")
    base = _fast_defaults(sc)
    base.update(cell_type="Buffer gas cell", qwp_deg=0.0, b_max_ut=80.0,
                scan_points=41, ne_pressure_torr=20.0,
                residual_transverse_b_ut=0.0)
    transit_dominated = dict(
        base, buffer_ground_relax_khz=20.0, collisional_depol_khz=2.0)
    depol_dominated = dict(
        base, buffer_ground_relax_khz=2.0, collisional_depol_khz=20.0)
    raw_transit = sc.compute(transit_dominated)
    raw_depol = sc.compute(depol_dominated)

    for key in ("chi_probe", "chi_p", "chi_m"):
        scale = max(float(np.max(np.abs(raw_transit[key]))), 1e-300)
        relative_difference = float(np.max(np.abs(
            raw_transit[key] - raw_depol[key]))) / scale
        assert relative_difference > 1e-3

    buffer_gamma = constants.neon_buffer_broadening(20.0)
    assert np.isclose(raw_transit["gamma_natural"],
                      raw_transit["gamma"] - buffer_gamma, rtol=1e-15)
    assert raw_transit["gamma_elastic"] == buffer_gamma / 2.0
    assert np.isclose(raw_transit["gamma_transit"] / (2 * np.pi), 20e3)
    assert np.isclose(raw_transit["gamma_depol"] / (2 * np.pi), 2e3)
    assert np.isclose(raw_depol["gamma_transit"] / (2 * np.pi), 2e3)
    assert np.isclose(raw_depol["gamma_depol"] / (2 * np.pi), 20e3)


def test_buffer_steady_states_remain_trace_one_hermitian_and_positive():
    gamma_natural = 2 * np.pi * 5.746e6
    gamma_elastic = 2 * np.pi * 39.1e6
    scheme = MagnetoScheme()
    drive = {+1: 1.0, -1: 1.0}

    for transit_khz, depol_khz in ((20.0, 2.0), (2.0, 20.0)):
        atom = zeeman.zeeman_manifold(
            2, 1, gamma=gamma_natural,
            gamma_gg=2 * np.pi * depol_khz * 1e3,
            transit_rate=2 * np.pi * transit_khz * 1e3,
            optical_dephasing=gamma_elastic)
        hamiltonian = scheme._hamiltonian(
            atom, (0.03e-6, 0.0, 0.02e-6), 0.5, -1.0 / 6.0,
            2 * np.pi * 1e6, drive)
        liouvillian = core.build_liouvillian(hamiltonian, atom)[None, :, :]
        rho = scheme._steady_state_buffer(
            liouvillian, np.array([0.0, 2 * np.pi * 3e6]),
            atom.S_v, atom.n_levels)

        trace_residual = np.max(np.abs(
            np.trace(rho, axis1=-2, axis2=-1) - 1.0))
        hermiticity_residual = np.max(np.abs(
            rho - np.swapaxes(rho.conj(), -1, -2)))
        hermitian_rho = 0.5 * (rho + np.swapaxes(rho.conj(), -1, -2))
        minimum_eigenvalue = float(np.linalg.eigvalsh(hermitian_rho).min())
        assert trace_residual <= 1e-10
        assert hermiticity_residual <= 1e-10
        assert minimum_eigenvalue >= -1e-10


def test_extreme_ui_dephasing_counterexample_remains_positive():
    sc = schemes.get("magneto")
    params = _fast_defaults(sc)
    params.update(
        cell_type="Buffer gas cell", Fg=2.0, Fe=1.0,
        ne_pressure_torr=0.0, intensity_mw_cm2=5.0, qwp_deg=45.0,
        residual_transverse_b_ut=5.0, transverse_field_angle_deg=0.0,
        laser_detuning_mhz=10.0, buffer_ground_relax_khz=0.1,
        collisional_depol_khz=200.0, b_max_ut=1.0, b_offset_ut=0.0,
        scan_points=3, velocity_classes=1, doppler="off")

    captured = {}
    steady_state = MagnetoScheme._steady_state_buffer

    def capture_steady_state(*args, **kwargs):
        rho = steady_state(*args, **kwargs)
        captured["rho"] = rho
        return rho

    with patch.object(MagnetoScheme, "_steady_state_buffer",
                      new=staticmethod(capture_steady_state)):
        raw = sc.compute(params)

    # This is the exact UI-valid point that produced min eigenvalue -3.3774e-6
    # with direct coherence-only damping. The middle scan point has physical Bz=0.
    rho = captured["rho"][1, 0]
    trace_residual = abs(np.trace(rho) - 1.0)
    hermiticity_residual = float(np.max(np.abs(rho - rho.conj().T)))
    hermitian_rho = 0.5 * (rho + rho.conj().T)
    minimum_eigenvalue = float(np.linalg.eigvalsh(hermitian_rho).min())
    assert raw["buffer_gamma"] == 0.0
    assert raw["gamma_elastic"] == 0.0
    assert trace_residual <= 1e-10
    assert hermiticity_residual <= 1e-10
    assert minimum_eigenvalue >= -1e-10


def test_longitudinal_b_offset_shifts_physical_field_axis():
    sc = schemes.get("magneto")
    p = _fast_defaults(sc)
    p.update(b_offset_ut=0.25)
    raw = sc.compute(p)
    assert raw["b_offset_ut"] == 0.25
    assert np.allclose(raw["b_physical_ut"] - raw["b_ut"], 0.25)


def test_intrinsic_eia_on_cycling_transition():
    # With transfer of coherence (grouped Σ_q emission), the open Fg=1->Fe=2
    # (Fe=Fg+1) transition is EIA at linear pol / zero residual field, while the
    # Fe<=Fg transitions stay EIT. (Lezama; arXiv physics/0512199.)
    sc = schemes.get("magneto")
    base = _fast_defaults(sc)
    base.update(cell_type="Paraffin coated cell", qwp_deg=0.0,
                residual_transverse_b_ut=0.0, b_max_ut=0.5, scan_points=201)
    eia = dict(base, Fg=1.0, Fe=2.0)
    eit = dict(base, Fg=2.0, Fe=1.0)
    x_a, a_a, _ = _alpha_rot(sc.compute(eia), eia)
    x_t, a_t, _ = _alpha_rot(sc.compute(eit), eit)
    assert _feature_amp(x_a, a_a) > 0      # Fe=Fg+1 -> absorption peak (EIA)
    assert _feature_amp(x_t, a_t) < 0      # Fe=Fg-1 -> transparency dip (EIT)


def test_buffer_circular_lca_needs_transverse_field():
    # Circular light orients the ground state along the beam, an eigenstate of the
    # longitudinal B scan -> flat (no feature) without a transverse field. A small
    # transverse residual field makes the orientation precess and gives a B=0
    # level-crossing ABSORPTION peak (Yu, PRA 81, 023416).
    sc = schemes.get("magneto")
    base = _fast_defaults(sc)
    base.update(cell_type="Buffer gas cell", Fg=2.0, Fe=2.0, qwp_deg=45.0,
                b_max_ut=1.0, scan_points=201, ne_pressure_torr=20.0,
                buffer_ground_relax_khz=5.0, collisional_depol_khz=0.5,
                transverse_field_angle_deg=90.0)
    flat = dict(base, residual_transverse_b_ut=0.0)
    lca = dict(base, residual_transverse_b_ut=0.03)
    x_f, a_f, _ = _alpha_rot(sc.compute(flat), flat)
    x_l, a_l, _ = _alpha_rot(sc.compute(lca), lca)
    assert abs(_feature_amp(x_f, a_f)) < 1e-6 * max(np.abs(a_f).max(), 1e-30)
    assert _feature_amp(x_l, a_l) > 0


def test_nmor_zero_crossing():
    sc = schemes.get("nmor")
    p = _fast_defaults(sc)
    p.update(signal_type="NMOR rotation", qwp_deg=0.0)
    x, _, rot = _alpha_rot(sc.compute(p), p)
    ic = int(np.argmin(np.abs(x)))
    iL = int(np.argmin(np.abs(x + 0.5 * np.max(np.abs(x)))))
    iR = int(np.argmin(np.abs(x - 0.5 * np.max(np.abs(x)))))
    assert abs(rot[ic]) < 1e-6 * max(np.abs(rot).max(), 1e-30)
    assert rot[iL] * rot[iR] < 0


def test_cell_type_controls_visibility_metadata():
    specs = schemes.get("magneto").param_schema()
    by_name = {s.name: s for s in specs}
    assert by_name["cell_type"].control == "segmented"
    assert by_name["transit_relax_khz"].visible_if == {"cell_type": "Paraffin coated cell"}
    assert by_name["ne_pressure_torr"].visible_if == {"cell_type": "Buffer gas cell"}


def test_invalid_transition_handled():
    sc = schemes.get("magneto")
    p = _fast_defaults(sc)
    p.update(Fg=1.0, Fe=3.0)
    raw = sc.compute(p)
    assert raw["valid"] is False
    view = sc.observables(raw, p)
    assert view.get("figure") is not None
    heroes = [m for m in view["metrics"] if m.get("tier") == "hero"]
    assert [m["label"] for m in heroes] == ["Status", "Transition"]
    assert heroes[0].get("kind") == "status"


def test_unresolved_central_width_uses_status_hero():
    sc = schemes.get("magneto")
    p = _fast_defaults(sc)
    p.update(b_offset_ut=50.0)
    view = sc.headless_observables(sc.compute(p), p)
    heroes = [m for m in view["metrics"] if m.get("tier") == "hero"]
    assert [m["label"] for m in heroes] == [
        "Transmission at B=0", "Width status"]
    assert heroes[1].get("kind") == "status"
    assert all(str(m["value"]).lower() != "n/a" for m in heroes)


if __name__ == "__main__":
    test_clebsch_gordan_known_values()
    test_manifold_emission_normalised()
    test_manifold_emission_transfers_ground_coherence()
    test_angular_momentum_commutator()
    test_default_is_87rb_d1_paraffin_polarized_hanle()
    test_paraffin_linear_qwp_gives_zero_field_dip()
    test_paraffin_circular_qwp_gives_zero_field_peak()
    test_paraffin_wall_coherence_narrows_central_feature()
    test_buffer_mode_has_single_broad_hanle_and_ground_relaxation_broadens()
    test_longitudinal_b_offset_shifts_physical_field_axis()
    test_nmor_zero_crossing()
    test_cell_type_controls_visibility_metadata()
    test_invalid_transition_handled()
    print("Phase-4 magneto OK (polarized Hanle, two-region paraffin, buffer cell, NMOR).")
