"""Focused checks for the weak-SIG Rydberg electrometry primitives."""

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from gabes import atoms, constants, core  # noqa: E402
from gabes.rydberg_electrometry import (  # noqa: E402
    BalancedDetector,
    PhotodiodeChannel,
    RFDipoleCoupling,
    asd_to_rms,
    balanced_detector_noise,
    coherent_weighted_average,
    current_responsivity_from_atomic_phasor,
    electrometry_sensitivity,
    photodiode_responsivity_a_per_w,
    weak_signal_response,
)


def _driven_two_level(omega_lo=1.3):
    atom = atoms.two_level(gamma=1.0)
    h = np.zeros((2, 2), dtype=complex)
    h[0, 1] = h[1, 0] = omega_lo / 2.0
    return atom, h


def _steady_excited_population(atom, hamiltonian):
    L = core.build_liouvillian(hamiltonian, atom)
    return core.steady_state_from_liouvillian(L, atom.n_levels)[1, 1].real


def test_weak_signal_sidebands_satisfy_linearized_equations_and_hermiticity():
    atom, h = _driven_two_level()
    omega_if = 0.37
    phase = 0.41
    response = weak_signal_response(
        atom, h, omega_if, signal_transition=(0, 1), signal_phase_rad=phase)

    L0 = core.build_liouvillian(h, atom)
    hm = np.zeros((2, 2), dtype=complex)
    hm[1, 0] = 0.5 * np.exp(-1j * phase)
    hp = hm.conj().T
    cm = core.comm_super(hm)
    cp = core.comm_super(hp)
    r0 = response.steady_state.reshape(-1)
    rm = response.rho_minus_per_angular_rabi.reshape(-1)
    rp = response.rho_plus_per_angular_rabi.reshape(-1)
    eye = np.eye(4)

    assert np.linalg.norm((L0 + 1j * omega_if * eye) @ rm + cm @ r0) < 1e-11
    assert np.linalg.norm((L0 - 1j * omega_if * eye) @ rp + cp @ r0) < 1e-11
    assert abs(np.trace(response.steady_state) - 1.0) < 1e-12
    assert abs(np.trace(response.rho_minus_per_angular_rabi)) < 1e-11
    assert abs(np.trace(response.rho_plus_per_angular_rabi)) < 1e-11
    assert np.allclose(
        response.rho_plus_per_angular_rabi,
        response.rho_minus_per_angular_rabi.conj().T,
        rtol=1e-10, atol=1e-11)


def test_low_if_phasor_matches_static_rabi_derivative_and_high_if_rolls_off():
    omega_lo = 1.3
    atom, h = _driven_two_level(omega_lo)
    population = np.diag([0.0, 1.0])

    low = weak_signal_response(
        atom, h, 1.0e-4, signal_transition=(0, 1))
    q_low = low.real_observable_phasor_per_angular_rabi(population)
    high = weak_signal_response(
        atom, h, 100.0, signal_transition=(0, 1))
    q_high = high.real_observable_phasor_per_angular_rabi(population)

    eps = 1.0e-5
    hp = h.copy()
    hm = h.copy()
    hp[0, 1] = hp[1, 0] = (omega_lo + eps) / 2.0
    hm[0, 1] = hm[1, 0] = (omega_lo - eps) / 2.0
    finite_difference = (
        _steady_excited_population(atom, hp)
        - _steady_excited_population(atom, hm)) / (2.0 * eps)

    assert np.isclose(q_low.real, finite_difference, rtol=2e-4, atol=2e-6)
    assert abs(q_low.imag) < 1e-3
    assert abs(q_high) < abs(q_low) / 20.0


def test_complex_phasors_are_coherently_averaged():
    phasors = np.array([1.0 + 1.0j, 1.0 - 1.0j])
    assert coherent_weighted_average(phasors, [1.0, 1.0]) == 1.0 + 0.0j


def test_rf_dipole_conversion_round_trip_and_rms_convention():
    dipole = 3.2e-27
    peak = RFDipoleCoupling(dipole, angular_factor=0.4,
                            field_amplitude_convention="peak")
    rms = RFDipoleCoupling(dipole, angular_factor=0.4,
                           field_amplitude_convention="rms")
    field = 2.5e-6
    omega = peak.angular_rabi_from_field(field)

    assert np.isclose(peak.field_from_angular_rabi(omega), field)
    assert np.isclose(
        peak.field_from_cyclic_rabi_hz(peak.cyclic_rabi_from_field_hz(field)),
        field)
    assert np.isclose(
        rms.angular_rabi_from_field(field), np.sqrt(2.0) * omega)
    assert np.isclose(
        peak.angular_rabi_per_field_rad_s_per_v_m,
        dipole * 0.4 / constants.HBAR)

    generalized_split_hz = 5.0e6
    detuning_hz = 3.0e6
    assert np.isclose(
        peak.field_from_at_splitting_hz(
            generalized_split_hz, detuning_hz=detuning_hz),
        peak.field_from_cyclic_rabi_hz(4.0e6),
    )
    try:
        peak.field_from_at_splitting_hz(2.0e6, detuning_hz=3.0e6)
    except ValueError as exc:
        assert "at least" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("unphysical AT split smaller than detuning was accepted")


def test_photodiode_responsivity_and_balanced_noise_budget():
    wavelength = 780e-9
    eta = 0.82
    responsivity = photodiode_responsivity_a_per_w(eta, wavelength)
    expected = (eta * constants.ELEMENTARY_CHARGE * wavelength
                / (2.0 * np.pi * constants.HBAR * constants.C_LIGHT))
    assert np.isclose(responsivity, expected)

    signal = PhotodiodeChannel(1.0e-3, 0.5, dark_current_a=2.0e-9)
    reference = PhotodiodeChannel(1.0e-3, 0.5, dark_current_a=3.0e-9)
    detector = BalancedDetector(
        signal, reference, reference_weight=1.0,
        electronic_noise_current_asd_a_per_sqrt_hz=4.0e-12,
        relative_intensity_noise_per_sqrt_hz=2.0e-6,
        rin_correlation=1.0)
    budget = balanced_detector_noise(detector)

    expected_shot = np.sqrt(2.0 * constants.ELEMENTARY_CHARGE * (
        0.5e-3 + 2.0e-9 + 0.5e-3 + 3.0e-9))
    assert np.isclose(budget.shot_noise_current_asd_a_per_sqrt_hz, expected_shot)
    assert budget.rin_noise_current_asd_a_per_sqrt_hz < 1e-18
    assert budget.technical_noise_current_asd_a_per_sqrt_hz == 4.0e-12
    assert np.isclose(
        budget.total_noise_current_asd_a_per_sqrt_hz,
        np.hypot(expected_shot, 4.0e-12))


def test_atomic_to_detector_chain_produces_absolute_sensitivity_without_anchor():
    coupling = RFDipoleCoupling(2.0e-27)
    coupling_rms = RFDipoleCoupling(
        2.0e-27, field_amplitude_convention="rms")
    atomic_phasor_per_rabi = 3.0e-7 + 4.0e-7j
    current_per_atomic = 2.0e-3
    current_phasor = current_responsivity_from_atomic_phasor(
        atomic_phasor_per_rabi, current_per_atomic, coupling)
    current_phasor_rms = current_responsivity_from_atomic_phasor(
        atomic_phasor_per_rabi, current_per_atomic, coupling_rms)
    # A linear responsivity is invariant when both field and photocurrent are
    # switched together from peak to RMS amplitude conventions.
    assert np.isclose(current_phasor_rms, current_phasor)

    detector = BalancedDetector(
        PhotodiodeChannel(20.0e-6, 0.6),
        electronic_noise_current_asd_a_per_sqrt_hz=1.0e-12)
    noise = balanced_detector_noise(detector)
    sensitivity = electrometry_sensitivity(noise, abs(current_phasor))

    assert sensitivity.total_field_asd_v_m_per_sqrt_hz >= (
        sensitivity.psn_field_asd_v_m_per_sqrt_hz)
    assert np.isclose(
        sensitivity.total_field_asd_nv_cm_per_sqrt_hz,
        sensitivity.total_field_asd_v_m_per_sqrt_hz * 1.0e7)
    assert np.isclose(
        asd_to_rms(sensitivity.total_field_asd_v_m_per_sqrt_hz, 25.0),
        5.0 * sensitivity.total_field_asd_v_m_per_sqrt_hz)


def test_invalid_dc_and_missing_calibrations_are_rejected():
    atom, h = _driven_two_level()
    try:
        weak_signal_response(atom, h, 0.0, signal_transition=(0, 1))
    except ValueError as exc:
        assert "if_angular_frequency" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("finite-IF solver accepted DC")

    try:
        RFDipoleCoupling(0.0)
    except ValueError as exc:
        assert "transition_dipole" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("zero transition dipole was accepted")


if __name__ == "__main__":
    test_weak_signal_sidebands_satisfy_linearized_equations_and_hermiticity()
    test_low_if_phasor_matches_static_rabi_derivative_and_high_if_rolls_off()
    test_complex_phasors_are_coherently_averaged()
    test_rf_dipole_conversion_round_trip_and_rms_convention()
    test_photodiode_responsivity_and_balanced_noise_budget()
    test_atomic_to_detector_chain_produces_absolute_sensitivity_without_anchor()
    test_invalid_dc_and_missing_calibrations_are_rejected()
    print("Rydberg electrometry linear-response checks OK.")
