"""Physics-limit tests for seeded-FWM propagation and normalization."""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm

from gabes import constants, hyperfine, observables
from gabes.schemes import fwm


def _dimensionless_chi(chi_physical, density, line_strength=1.0):
    coupling = (
        -2.0 * density * line_strength * constants.DIPOLE_D1**2
        / (constants.EPS_0 * constants.HBAR)
    )
    return np.asarray(chi_physical, dtype=complex) / coupling


def test_positive_imaginary_susceptibility_obeys_beer_lambert():
    density = 2.0e15
    k = 7.9e6
    length = 0.0125
    chi_physical = np.array([1.3e-6 + 2.1e-6j])
    chi_bar = _dimensionless_chi(chi_physical, density)
    zero = np.zeros(1, dtype=complex)

    probe_gain, conjugate_gain, transfer = observables.gain_from_chi(
        chi_bar, zero, zero, chi_bar,
        k, k, length, density, line_strength=1.0,
    )

    expected_intensity = np.exp(-k * np.imag(chi_physical) * length)
    expected_probe_amplitude = np.exp(0.5j * k * chi_physical * length)
    expected_conjugate_amplitude = np.exp(
        -0.5j * k * np.conj(chi_physical) * length)
    np.testing.assert_allclose(probe_gain, expected_intensity, rtol=1e-13)
    np.testing.assert_allclose(transfer[:, 0, 0], expected_probe_amplitude,
                               rtol=1e-13, atol=1e-15)
    np.testing.assert_allclose(transfer[:, 1, 1], expected_conjugate_amplitude,
                               rtol=1e-13, atol=1e-15)
    np.testing.assert_allclose(conjugate_gain, 0.0, atol=1e-15)
    assert probe_gain[0] < 1.0


def test_mismatch_rotation_matches_unrotated_fwm_phases():
    """Pin the sign implied by exp[i(kz-wt)] and delta_k=2kp-ks-kc."""
    density = 3.0e17
    k_probe, k_conj = 7.90e6, 7.91e6
    delta_k = 73.0
    length = 0.007
    chi_physical = (
        0.7e-6 + 0.2e-6j,
        -0.4e-6 + 0.3e-6j,
        0.2e-6 - 0.5e-6j,
        -0.1e-6 + 0.4e-6j,
    )
    reduced = tuple(
        _dimensionless_chi([value], density) for value in chi_physical)
    M0 = observables._gain_matrix_from_chi(
        *reduced, k_probe, k_conj, density, constants.DIPOLE_D1, 1.0)[0]
    M_rot = observables._gain_matrix_from_chi(
        *reduced, k_probe, k_conj, density, constants.DIPOLE_D1, 1.0,
        delta_k_z=delta_k)[0]
    expected_rot = M0 + np.diag([-0.5j * delta_k, 0.5j * delta_k])
    np.testing.assert_allclose(M_rot, expected_rot, rtol=0.0, atol=1e-13)

    def unrotated_rhs(z, field):
        phase = np.exp(1j * delta_k * z)
        matrix = np.array([
            [M0[0, 0], M0[0, 1] * phase],
            [M0[1, 0] / phase, M0[1, 1]],
        ])
        return matrix @ field

    initial = np.array([1.0 + 0.2j, -0.3 + 0.1j])
    solved = solve_ivp(
        unrotated_rhs, (0.0, length), initial, rtol=2e-12, atol=2e-13)
    rotate_at_exit = np.diag([
        np.exp(-0.5j * delta_k * length),
        np.exp(0.5j * delta_k * length),
    ])
    np.testing.assert_allclose(
        rotate_at_exit @ solved.y[:, -1], expm(M_rot * length) @ initial,
        rtol=2e-10, atol=2e-11)


def test_zero_coupling_is_identity_for_homogeneous_and_segmented_paths():
    chi = np.array([1.0 + 2.0j, -0.4 + 0.7j])
    zero = np.zeros_like(chi)
    eye = np.broadcast_to(np.eye(2, dtype=complex), (chi.size, 2, 2))

    for kwargs in (
        {},
        {"propagation_segments": 8, "segment_profile": np.ones(8)},
    ):
        probe_gain, conjugate_gain, transfer = observables.gain_from_chi(
            chi, zero, zero, chi, constants.K_VEC, constants.K_VEC,
            0.01, 0.0, **kwargs,
        )
        np.testing.assert_allclose(transfer, eye, atol=1e-15)
        np.testing.assert_allclose(probe_gain, 1.0, atol=1e-15)
        np.testing.assert_allclose(conjugate_gain, 0.0, atol=1e-15)


def test_constant_medium_one_segment_matches_many_segments():
    rng = np.random.default_rng(17)
    chi = [rng.normal(size=3) * 1e-8 + 1j * rng.normal(size=3) * 1e-8
           for _ in range(4)]
    common = dict(
        k_probe=7.90e6,
        k_conj=7.91e6,
        L=0.0125,
        N_atoms=4.0e17,
        line_strength=0.08,
        delta_k_z=np.array([20.0, 30.0, 40.0]),
    )
    one = observables.gain_from_chi(*chi, **common)
    many = observables.gain_from_chi(
        *chi, **common, propagation_segments=64, segment_profile=np.ones(64))
    for lhs, rhs in zip(one, many):
        np.testing.assert_allclose(lhs, rhs, rtol=2e-12, atol=2e-13)


def test_segment_probe_od_is_reported_without_silent_clipping():
    density = 2.0e15
    target_od = 7.0
    chi_physical = 1j * target_od / (constants.K_VEC * fwm.L_CELL)
    chi_bar = _dimensionless_chi([chi_physical], density)
    profile, reported_od = fwm._uniform_segment_profile_and_probe_od(
        chi_bar, density, 1.0, nseg=8)
    np.testing.assert_allclose(profile, 1.0, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(reported_od, target_od, rtol=2e-15, atol=0.0)


def test_coupling_ledger_matches_declared_auto_od_structure():
    """Gate the coefficient ledger; this is not a pump-state rho_ss test."""
    for branch, ground in ((-1, fwm.G2), (+1, fwm.G1)):
        ledger = fwm.physical_coupling_ledger(branch)
        ground_F = fwm.GROUND_F[ground]
        rho_ground = hyperfine.GROUND_POP[ground_F]

        assert ledger["manifold_population_source"] == "trace-normalized rho_ss"
        assert ledger["external_manifold_population_factor"] == 1.0
        assert ledger["macroscopic_coupling_norm"] == (
            1.0 / hyperfine.N_GROUND_SUBLEVELS)
        assert fwm.physical_coupling_norm(branch) == (
            1.0 / hyperfine.N_GROUND_SUBLEVELS)

        for excited in fwm.EXCITED_STATES:
            transition_factor = fwm.TRANSITION_DIPOLE_SCALE[ground, excited]**2
            reduced_strength = (
                transition_factor * rho_ground
                * ledger["macroscopic_coupling_norm"])
            auto_od_strength = (
                3.0 * hyperfine.CF2[(ground_F, fwm.EXCITED_F[excited])]
                * hyperfine.GROUND_POP[ground_F]
                / hyperfine.N_GROUND_SUBLEVELS)
            np.testing.assert_allclose(reduced_strength, auto_od_strength,
                                       rtol=1e-15, atol=1e-15)


def test_unequal_mode_q_recovers_canonical_bogoliubov_map():
    squeeze = 0.7
    phase = 0.31
    mu = np.cosh(squeeze)
    nu = np.exp(1j * phase) * np.sinh(squeeze)
    transfer_canonical = np.array(
        [[mu, nu], [np.conj(nu), mu]], dtype=complex)
    omega_probe = constants.OMEGA_D1 + 2 * np.pi * 0.9e9
    omega_conj = omega_probe - 2 * constants.OMEGA_HF
    area_probe = 0.5 * np.pi * (330e-6)**2
    area_conj = 0.5 * np.pi * (470e-6)**2
    q_probe = np.sqrt(
        2.0 * constants.HBAR * omega_probe
        / (constants.EPS_0 * constants.C_LIGHT * area_probe))
    q_conj = np.sqrt(
        2.0 * constants.HBAR * omega_conj
        / (constants.EPS_0 * constants.C_LIGHT * area_conj))
    Q_literal = np.diag([q_probe, q_conj])
    Q_api = observables.photon_flux_mode_matrix(
        omega_probe, omega_conj, area_probe, area_conj)
    np.testing.assert_allclose(
        observables.gaussian_mode_area(330e-6), area_probe, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(Q_api, Q_literal, rtol=2e-15, atol=0.0)
    transfer_field = Q_literal @ transfer_canonical @ np.linalg.inv(Q_literal)

    result = observables.canonical_transfer_diagnostics(
        transfer_field, omega_probe, omega_conj, area_probe, area_conj)
    np.testing.assert_allclose(result["transfer_canonical"], transfer_canonical,
                               rtol=1e-15, atol=1e-15)
    np.testing.assert_allclose(result["photon_flux_gap"], 1.0,
                               rtol=1e-15, atol=1e-15)
    np.testing.assert_allclose(result["commutator_defect"], 0.0, atol=1e-15)
    expected_power_gain = (
        omega_conj / omega_probe * np.abs(transfer_canonical[1, 0])**2)
    np.testing.assert_allclose(result["conjugate_power_gain"],
                               expected_power_gain, rtol=1e-15, atol=1e-15)


def test_stored_literature_chi_matches_independent_option_a_golden():
    """Independent SI/expm anchor for sign, normalization, Q, and full T."""
    reduced = np.array([
        -2.825216423278075e-11 - 3.156193807767438e-13j,
        -4.6431534133265106e-11 - 9.727224915826937e-14j,
        -4.6424201142453104e-11 + 1.6459703154104033e-13j,
        -1.132757317884393e-11 + 8.098427191674804e-14j,
    ])
    density = 2.149563741450344e19
    line_strength = 0.74 / 12.0
    k_probe = 7903541.441086388
    k_conj = 7903669.024915995
    delta_k = 246.53512446396053
    length = 0.0125
    prefactor = (
        -2.0 * density * constants.DIPOLE_D1**2 * line_strength
        / (constants.EPS_0 * constants.HBAR))
    physical = prefactor * reduced
    M_literal = np.array([
        [0.5j * k_probe * physical[0] - 0.5j * delta_k,
         0.5j * k_probe * physical[1]],
        [-0.5j * k_conj * np.conj(physical[2]),
         -0.5j * k_conj * np.conj(physical[3]) + 0.5j * delta_k],
    ])
    T_literal = expm(M_literal * length)
    area = 0.5 * np.pi * (330e-6)**2
    omega_probe = constants.C_LIGHT * k_probe
    omega_conj = constants.C_LIGHT * k_conj
    Q_literal = np.diag([
        np.sqrt(2 * constants.HBAR * omega_probe
                / (constants.EPS_0 * constants.C_LIGHT * area)),
        np.sqrt(2 * constants.HBAR * omega_conj
                / (constants.EPS_0 * constants.C_LIGHT * area)),
    ])
    T_canonical_literal = np.linalg.inv(Q_literal) @ T_literal @ Q_literal

    chi_args = tuple(np.asarray([value]) for value in reduced)
    _, _, T_api = observables.gain_from_chi(
        *chi_args, k_probe, k_conj, length, density,
        line_strength=line_strength, delta_k_z=delta_k)
    diagnostics = observables.canonical_transfer_diagnostics(
        T_api, omega_probe, omega_conj, area, area)
    np.testing.assert_allclose(T_api[0], T_literal, rtol=2e-13, atol=2e-13)
    np.testing.assert_allclose(
        diagnostics["transfer_canonical"][0], T_canonical_literal,
        rtol=2e-13, atol=2e-13)
    np.testing.assert_allclose(
        [diagnostics["probe_power_gain"][0],
         diagnostics["conjugate_power_gain"][0],
         diagnostics["conjugate_photon_flux_gain"][0],
         diagnostics["photon_flux_gap"][0]],
        [1054.3059596473247, 1061.885343068735,
         1061.8682017387723, -7.562242091447615],
        rtol=2e-13, atol=2e-12)


def test_manley_rowe_cap_is_energy_bounded_and_not_a_linear_map():
    pump_power = 0.6
    seed_power = 8e-6
    generated = np.array([1e2, 1e5, 1e9])
    probe_gain, conjugate_gain = observables.pump_depletion_saturation(
        1.0 + generated, generated, pump_power, seed_power)
    cap = 0.5 * pump_power
    assert np.all((probe_gain - 1.0) * seed_power <= cap)
    assert np.all(conjugate_gain * seed_power <= cap)
    np.testing.assert_allclose(probe_gain - 1.0, conjugate_gain)


def test_option_a_wavenumbers_and_mismatch_are_bare_vacuum_values():
    detuning = 0.9
    probe = np.array([-2.15, -2.14, -2.13])
    angle = 0.32
    k_pump, k_probe, k_conj = fwm.seeded_option_a_wavenumbers(detuning, probe)
    expected_probe = (
        constants.OMEGA_D1 + 2 * np.pi * probe * 1e9) / constants.C_LIGHT
    expected_conj = (
        constants.OMEGA_D1 + 2 * np.pi * (2 * detuning - probe) * 1e9
    ) / constants.C_LIGHT
    expected_mismatch = (
        2 * k_pump - (expected_probe + expected_conj)
        * np.cos(np.radians(angle)))

    np.testing.assert_allclose(k_probe, expected_probe, rtol=0.0, atol=1e-12)
    np.testing.assert_allclose(k_conj, expected_conj, rtol=0.0, atol=1e-12)
    np.testing.assert_allclose(
        fwm.seeded_phase_mismatch_z(detuning, probe, angle_deg=angle),
        expected_mismatch, rtol=0.0, atol=1e-12)
