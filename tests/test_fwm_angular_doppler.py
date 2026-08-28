"""Two-dimensional non-collinear Raman-Doppler reference tests."""

import numpy as np
import pytest

from gabes import constants, doppler, observables
from gabes.schemes import fwm


T_LITERATURE = 394.15
D_GHZ = 0.9
DELTA_MHZ = -8.0
THETA_RAD = np.deg2rad(0.32)


def _geometry(delta_mhz=DELTA_MHZ):
    probe_GHz = np.atleast_1d(
        D_GHZ - constants.NU_HF / 1e9 + np.asarray(delta_mhz) * 1e-3)
    k_pump, k_probe, k_conj = fwm.seeded_option_a_wavenumbers(
        D_GHZ, probe_GHz)
    return probe_GHz, k_pump, k_probe, k_conj


def _reference(*, order=18, cutoff=5.0, theta=THETA_RAD,
               k_probe_override=None):
    probe_GHz, k_pump, k_probe, _k_conj = _geometry()
    if k_probe_override is not None:
        k_probe = np.asarray([k_probe_override], dtype=float)
    pump = fwm.rabi_freq(0.6, fwm.W_PUMP)
    delta = 2.0 * np.pi * np.asarray([DELTA_MHZ * 1e6])
    return fwm.pump_only_weak_response_noncollinear_reference(
        pump,
        pump,
        delta,
        2.0 * np.pi * D_GHZ * 1e9,
        T=T_LITERATURE,
        pump_k_rad_m=k_pump,
        probe_k_axis_rad_m=k_probe,
        crossing_angle_rad=theta,
        quadrature_order=order,
        cutoff_sigma=cutoff,
    )


def _option_a_gain(reference):
    probe_GHz, _k_pump, k_probe, k_conj = _geometry()
    number_density = fwm.hyperfine.number_density(T_LITERATURE)
    line_strength = (
        fwm.SEEDED_REFERENCE_RESIDUAL * fwm.physical_coupling_norm(-1))
    mismatch = fwm.seeded_phase_mismatch_z(
        D_GHZ, probe_GHz, angle_deg=np.rad2deg(reference.crossing_angle_rad))
    gain_probe, gain_conj, _transfer = observables.gain_from_chi(
        reference.chi_ss[0],
        reference.chi_sc[0],
        reference.chi_cs[0],
        reference.chi_cc[0],
        k_probe,
        k_conj,
        fwm.L_CELL,
        number_density,
        line_strength=line_strength,
        delta_k_z=mismatch,
    )
    return float(gain_probe[0]), float(gain_conj[0])


def test_maxwell_legendre_grid_and_analytic_raman_width_budget():
    velocity, weights = doppler.maxwell_legendre_grid(
        T_LITERATURE, order=24, cutoff_sigma=5.0)
    sigma = np.sqrt(constants.KB * T_LITERATURE / constants.MASS_85RB)
    assert np.sum(weights) == pytest.approx(1.0, abs=2e-15)
    assert np.sum(weights * velocity) == pytest.approx(0.0, abs=1e-12)
    assert np.sqrt(np.sum(weights * velocity**2)) == pytest.approx(
        sigma, rel=2e-5)

    _probe, k_pump, k_probe, _conj = _geometry()
    angular = doppler.noncollinear_raman_rms_budget(
        T_LITERATURE, k_pump, k_probe[0], THETA_RAD)
    collinear = doppler.noncollinear_raman_rms_budget(
        T_LITERATURE, k_pump, k_probe[0], 0.0)
    assert angular["total_rms_hz"] == pytest.approx(1.38e6, rel=0.01)
    assert 1e3 < collinear["axial_rms_hz"] < 3e3
    assert collinear["transverse_rms_hz"] == pytest.approx(0.0, abs=0.0)


def test_atomic_detunings_have_the_declared_transverse_sign():
    one_photon = 2.0 * np.pi * 0.9e9
    two_photon = 2.0 * np.pi * -8e6
    _probe, k_pump, k_probe, _conj = _geometry()
    Delta_eff, delta_eff = doppler.noncollinear_atomic_detunings_rad_s(
        one_photon,
        two_photon,
        vx_m_s=10.0,
        vz_m_s=20.0,
        pump_k_rad_m=k_pump,
        probe_k_rad_m=k_probe[0],
        crossing_angle_rad=THETA_RAD,
    )
    expected = (two_photon
                + (k_pump - k_probe[0] * np.cos(THETA_RAD)) * 20.0
                - k_probe[0] * np.sin(THETA_RAD) * 10.0)
    assert Delta_eff == pytest.approx(one_photon - k_pump * 20.0)
    assert delta_eff == pytest.approx(expected)


def test_collinear_equal_k_limit_matches_doppler_free_raman_reference():
    order = 10
    cutoff = 4.5
    _probe, k_pump, _k_probe, _conj = _geometry()
    two_dimensional = _reference(
        order=order, cutoff=cutoff, theta=0.0, k_probe_override=k_pump)

    velocity, weights = doppler.maxwell_legendre_grid(
        T_LITERATURE, order=order, cutoff_sigma=cutoff)
    pump = fwm.rabi_freq(0.6, fwm.W_PUMP)
    delta = 2.0 * np.pi * np.asarray([DELTA_MHZ * 1e6])
    one_dimensional = fwm.pump_only_weak_response_reference(
        pump,
        pump,
        delta,
        2.0 * np.pi * D_GHZ * 1e9 - k_pump * velocity,
        atom=fwm.collisional_atom(T_LITERATURE),
    )
    expected = np.einsum(
        "v,voi->oi", weights, one_dimensional.chi_matrix[0, 0])
    np.testing.assert_allclose(
        two_dimensional.chi_matrix[0, 0], expected, rtol=3e-13, atol=3e-24)
    np.testing.assert_allclose(
        two_dimensional.raman_shift_grid_rad_s, 0.0, atol=0.0)


def test_lab_beat_is_not_velocity_shifted_and_plus_branch_fails_closed():
    reference = _reference(order=14)
    expected_beat = fwm.seeded_sideband_beat(
        reference.lab_delta_axis_rad_s, -1)
    np.testing.assert_array_equal(
        reference.lab_optical_beat_axis_rad_s, expected_beat)
    assert reference.diagnostics["lab_beat_velocity_invariant"]
    assert np.ptp(reference.raman_shift_grid_rad_s[0]) > 0.0
    assert reference.provenance["production_default"] is False

    with pytest.raises(NotImplementedError, match="minus"):
        fwm.pump_only_weak_response_noncollinear_reference(
            1.0,
            1.0,
            [0.0],
            0.0,
            T=T_LITERATURE,
            pump_k_rad_m=1.0,
            probe_k_axis_rad_m=1.0,
            crossing_angle_rad=0.0,
            branch=+1,
        )


def test_grid_and_cutoff_refinement_stabilize_reference_gain():
    medium = _reference(order=14, cutoff=5.0)
    fine = _reference(order=18, cutoff=5.0)
    tail = _reference(order=18, cutoff=4.5)
    medium_gain = _option_a_gain(medium)
    fine_gain = _option_a_gain(fine)
    tail_gain = _option_a_gain(tail)
    for candidate, reference in ((medium_gain, fine_gain), (tail_gain, fine_gain)):
        relative = np.max(np.abs(
            (np.asarray(candidate) - np.asarray(reference))
            / np.asarray(reference)))
        assert relative < 0.01
    rms_numeric = fine.diagnostics["quadrature_raman_rms_rad_s"][0]
    rms_analytic = fine.diagnostics["analytic_raman_rms_rad_s"][0]
    assert rms_numeric == pytest.approx(rms_analytic, rel=0.01)
    assert fine.diagnostics["max_response_normalized_residual"] < 1e-12
    assert fine.diagnostics["max_response_trace_error"] < 1e-12
