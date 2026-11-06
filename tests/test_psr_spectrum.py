"""Independent PSR atom-to-field checks in the fixed and outgoing-carrier frames."""
import numpy as np
import pytest
from scipy.integrate import solve_ivp

from analysis.grand_challenge.reference.atomic_qrt import qrt_ordered_spectrum
from gabes import constants, core, species
from gabes.psr_full_obe import solve_psr_pump
from gabes.psr_homodyne import psr_homodyne_readout
from gabes.psr_spectrum import psr_local_generator, psr_uniform_transfer
from gabes.quantum.channels import covariance_uncertainty_minimum
from gabes.quantum.traveling import constant_segment


DETUNING_MHZ = 400.
INTENSITY_MW_CM2 = 12.5
DENSITY_M3 = 1.e16
RESET_PER_S = 1.e5
RF_HZ = np.array([-4.e5, 0., 4.e5])


@pytest.fixture(scope="module")
def driven():
    pump = solve_psr_pump(DETUNING_MHZ, intensity_mw_cm2=INTENSITY_MW_CM2,
                          density_m3=DENSITY_M3, ground_relaxation_hz=RESET_PER_S)
    return pump, psr_local_generator(pump, RF_HZ)


def test_static_nambu_jacobian_matches_independent_maxwell_field_derivatives(driven):
    pump, local = driven
    np.testing.assert_array_equal(local.signs, [1, 1, -1, -1])
    assert local.mode_labels == ("x", "y", "x_dagger", "y_dagger")
    k = 2*np.pi*pump["conditions"]["optical_frequency_hz"]/constants.C_LIGHT
    carrier = pump["conditions"]["peak_field_v_m"]
    linear_x = np.array([1., 0.])
    step = 1.e-5

    def emitted_field(relative_field):
        # Perturb the physical envelope. A real parallel perturbation changes
        # intensity as well as the unit Jones vector in the S1 API.
        norm = np.linalg.norm(relative_field)
        solution = solve_psr_pump(
            DETUNING_MHZ, intensity_mw_cm2=INTENSITY_MW_CM2*norm**2,
            jones_xy=relative_field/norm, density_m3=DENSITY_M3,
            ground_relaxation_hz=RESET_PER_S)
        return .5j*k*solution["polarization_xy_c_m2"]/(constants.EPS_0*carrier)

    # Four independent real perturbations constrain every complex entry.
    for unit in np.eye(2):
        for direction in (unit, 1j*unit):
            finite = (emitted_field(linear_x + step*direction)
                      - emitted_field(linear_x - step*direction))/(2*step)
            vector = np.r_[direction, direction.conj()]
            expected = np.r_[finite, finite.conj()]
            np.testing.assert_allclose(local.drift[1] @ vector, expected,
                                       rtol=3e-6, atol=3e-7)


def test_two_time_qrt_matches_both_projected_langevin_orderings(driven):
    # This reference propagates delta(operator)*rho through the original
    # Liouvillian. It never reads production adjoint rows or jump diffusion.
    pump, local = driven
    atom = species.level_structure("Rb87", "D1")
    lowered = np.zeros((2, 16, 16), complex)
    lowered[0, :8, 8:] = (atom.dipoles[1] + atom.dipoles[-1])/np.sqrt(2)
    lowered[1, :8, 8:] = 1j*(atom.dipoles[1] - atom.dipoles[-1])/np.sqrt(2)
    operators = []
    for dipole in lowered:
        operators.extend(((dipole + dipole.conj().T)/2,
                          (dipole - dipole.conj().T)/(2j)))
    with core.blas_single_thread():
        ordered = qrt_ordered_spectrum(pump["generator_s_inverse"], pump["rho"],
                                       operators, 2*np.pi*RF_HZ)

    # Maxwell emission into a photon-flux envelope is proportional to i*d,
    # into the conjugate envelope to -i*d†. For plane waves n A |g|² =
    # n*k*d0²/(2 epsilon0 hbar), independent of the auxiliary quantization area.
    transform = np.array([[1j, -1., 0., 0.], [0., 0., 1j, -1.],
                          [-1j, -1., 0., 0.], [0., 0., -1j, -1.]])
    k = 2*np.pi*atom.frequency_hz/constants.C_LIGHT
    prefactor = DENSITY_M3*k*atom.dipole_sq/(2*constants.EPS_0*constants.HBAR)
    expected_greater = prefactor*(transform @ ordered @ transform.conj().T)
    # Reverse time and operator order on the Hermitian QRT observables.
    expected_lesser = prefactor*(transform @ ordered[::-1].transpose(0, 2, 1) @ transform.conj().T)
    for expected, actual in ((expected_greater, local.noise_greater),
                             (expected_lesser, local.noise_lesser)):
        for reference, calculated in zip(expected, actual):
            scale = np.linalg.norm(reference)
            assert np.linalg.norm(calculated - reference) < 1.e-7*scale
        assert np.linalg.norm(expected[:, :2, 2:]) > .001*np.linalg.norm(expected)


def test_passive_independent_sidebands_follow_shifted_susceptibility_and_preserve_vacuum():
    # No drive gives a thermal ground mixture and no optical emission.
    def passive(delta):
        return solve_psr_pump(delta, intensity_mw_cm2=0., density_m3=DENSITY_M3,
                              ground_relaxation_hz=RESET_PER_S)

    pump = passive(DETUNING_MHZ)
    local = psr_local_generator(pump, RF_HZ)
    k = 2*np.pi*pump["conditions"]["optical_frequency_hz"]/constants.C_LIGHT
    fixed_transfer = constant_segment(local, .025)
    greater, lesser = fixed_transfer.vacuum_output()
    for i, f in enumerate(RF_HZ):
        positive = .5j*k*passive(DETUNING_MHZ + f/1.e6)["chi_probe"]
        mirrored = .5j*k*passive(DETUNING_MHZ - f/1.e6)["chi_probe"]
        np.testing.assert_allclose(local.drift[i], np.diag([positive, positive,
                                                           mirrored.conjugate(), mirrored.conjugate()]),
                                   rtol=3.e-10, atol=1.e-10)
        np.testing.assert_allclose(greater[i], np.diag([1, 1, 0, 0]), atol=3.e-10)
        np.testing.assert_allclose(lesser[i], np.diag([0, 0, 1, 1]), atol=3.e-10)


def test_outgoing_carrier_gauge_matches_rotating_laboratory_slice_coefficients(driven):
    pump, local = driven
    k = 2*np.pi*pump["conditions"]["optical_frequency_hz"]/constants.C_LIGHT
    beta = .5*k*pump["chi_probe"].real
    length = .10/abs(beta)  # Resolve an appreciable mean phase increment.
    transfer = psr_uniform_transfer(local, pump, length)
    index = 2
    drift = local.drift[index]
    noises = (local.noise_greater[index], local.noise_lesser[index])

    def phase(z):
        return np.exp(1j*beta*z*local.signs)

    def propagate_lab(z, packed):
        u = phase(z)
        def rotate(matrix):
            return u[:, None]*matrix*u.conj()[None, :]
        lab = rotate(drift)
        t, ng, nl = packed.reshape(3, 4, 4)
        return np.array([lab @ t, lab @ ng + ng @ lab.conj().T + rotate(noises[0]),
                         lab @ nl + nl @ lab.conj().T + rotate(noises[1])]).ravel()

    initial = np.zeros((3, 4, 4), complex)
    initial[0] = np.eye(4)
    integration = solve_ivp(propagate_lab, (0., length), initial.ravel(), rtol=1.e-10, atol=1.e-12)
    assert integration.success
    t_lab, ng_lab, nl_lab = integration.y[:, -1].reshape(3, 4, 4)
    u = phase(length)
    np.testing.assert_allclose(transfer.transfer[index], u.conj()[:, None]*t_lab, atol=2.e-10)
    for expected, actual in ((ng_lab, transfer.noise_greater[index]),
                             (nl_lab, transfer.noise_lesser[index])):
        np.testing.assert_allclose(actual, u.conj()[:, None]*expected*u[None, :], atol=2.e-10)


def test_zero_density_removes_all_atomic_drift_and_noise():
    pump = solve_psr_pump(DETUNING_MHZ, intensity_mw_cm2=INTENSITY_MW_CM2, density_m3=0.)
    local = psr_local_generator(pump, [-1e6, 1e6])
    for matrix in (local.drift, local.noise_greater, local.noise_lesser):
        np.testing.assert_array_equal(matrix, 0.)
    transfer = psr_uniform_transfer(local, pump, .075)
    np.testing.assert_allclose(transfer.transfer, [np.eye(4)]*2, atol=1e-15)


@pytest.mark.parametrize("jones", [(0., 1.), (1j, 0.)])
def test_phase_referenced_path_checks_its_x_carrier_axes(jones):
    pump = solve_psr_pump(DETUNING_MHZ, intensity_mw_cm2=INTENSITY_MW_CM2,
                          jones_xy=jones, density_m3=DENSITY_M3, ground_relaxation_hz=RESET_PER_S)
    local = psr_local_generator(pump, [-4.e5, 4.e5])
    with pytest.raises(ValueError):
        psr_uniform_transfer(local, pump, .025)


def test_driven_homodyne_spectrum_matches_frequency_domain_photocurrent(driven):
    pump, local = driven
    transfer = psr_uniform_transfer(local, pump, .025)
    efficiency = .9025
    phases = np.array([0., 30., 90., 120.])
    readout = psr_homodyne_readout(transfer, phase_degrees=phases, detection_efficiency=efficiency)
    t = transfer.transfer[2]
    positive_ordering = t @ np.diag([1, 1, 0, 0]) @ t.conj().T + transfer.noise_greater[2]
    for shift, trace in ((0., "noise_a"), (90., "noise_b")):
        reference = []
        for theta in np.deg2rad(phases + shift):
            # i(+Ω) combines the annihilator at +Ω and creator at -Ω.
            row = np.array([0., np.exp(-1j*theta), 0., np.exp(1j*theta)])
            power = row @ positive_ordering @ row.conj()
            assert abs(power.imag) < 1.e-12
            reference.append(efficiency*power.real + 1 - efficiency)
        np.testing.assert_allclose(readout[trace][0], reference, atol=3.e-12)
    assert covariance_uncertainty_minimum(readout["source_sideband_covariance"][0]) >= -1.e-12


@pytest.mark.parametrize("detuning,intensity,reset,jones", [
    (-900., .1, 1.e3, (1., 0.)),
    (350., 100., 1.e3, (1., 0.)),
    (400., 0., 1.e5, (1., 0.)),
    (-120., 1.e3, 1.e5, (np.sqrt(.7), 1j*np.sqrt(.3))),
])
def test_optimized_noise_matches_direct_with_dc_paired_and_unpaired_rf(
        detuning, intensity, reset, jones):
    pump = solve_psr_pump(detuning, intensity_mw_cm2=intensity,
                          ground_relaxation_hz=reset, density_m3=DENSITY_M3,
                          jones_xy=jones)
    frequency = [-2.e6, -4.e5, 0., 4.e5, 1.7e6, 2.e6]
    direct = psr_local_generator(pump, frequency)
    optimized = psr_local_generator(pump, frequency, method="optimized")
    assert optimized.audit()["passed"]
    assert optimized.reservoir_names == direct.reservoir_names
    for name in ("drift", "noise_greater_by_reservoir", "noise_lesser_by_reservoir"):
        expected, actual = getattr(direct, name), getattr(optimized, name)
        # Source-by-source Frobenius errors remain meaningful near exact zero
        # entries, unlike pointwise relative comparisons of numerical zeros.
        np.testing.assert_allclose(actual, expected, rtol=1.e-10, atol=1.e-12)
        assert np.linalg.norm(actual - expected) <= 1.e-11*np.linalg.norm(expected)


def test_replacement_gram_identity_uses_declared_jumps_and_has_general_fallback():
    from gabes.psr_spectrum import _jump_gram, _replacement_distribution, _replacement_gram
    from gabes.quantum.reservoirs import thermal_reset_channels

    random = np.random.default_rng(19)
    raw = random.normal(size=(5, 5)) + 1j*random.normal(size=(5, 5))
    rho = raw @ raw.conj().T
    rho /= np.trace(rho)
    observables = random.normal(size=(4, 5, 5)) + 1j*random.normal(size=(4, 5, 5))
    jumps = np.stack([channel.operator for channel in thermal_reset_channels(
        7.3, [.5, .3, .2, 0., 0.], name="test", source="test")])
    replacement = _replacement_distribution(jumps)
    assert replacement is not None
    np.testing.assert_allclose(_replacement_gram(replacement, observables, rho),
                               _jump_gram(jumps, observables, rho), rtol=2.e-15, atol=3.e-14)
    changed = jumps.copy()
    changed[0, 0, 1] = .01
    assert _replacement_distribution(changed) is None
    assert _replacement_distribution(jumps[:-1]) is None


def test_local_method_validation_and_zero_density_optimized():
    pump = solve_psr_pump(0., density_m3=0.)
    optimized = psr_local_generator(pump, [-1.e6, 0., 1.e6], method="optimized")
    for name in ("drift", "noise_greater_by_reservoir", "noise_lesser_by_reservoir"):
        np.testing.assert_array_equal(getattr(optimized, name), 0.)
    with pytest.raises(ValueError, match="method"):
        psr_local_generator(pump, [-1.e6, 1.e6], method="unknown")
