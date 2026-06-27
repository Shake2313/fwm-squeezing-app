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

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from gabes import constants, observables, schemes, zeeman  # noqa: E402
from gabes.constants import GAMMA  # noqa: E402

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
