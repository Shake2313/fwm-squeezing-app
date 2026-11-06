"""Independent physical checks on the Rb87 D1 single-velocity steady state."""

import numpy as np
import pytest

from gabes import constants, magneto_multilevel, psr, species
from gabes.psr_full_obe import solve_psr_pump


DENSITY_M3 = 1e15
RESET_PER_S = 1000.


def test_zero_drive_is_thermal_and_has_the_linear_probe_susceptibility():
    detuning = 400.
    state = solve_psr_pump(detuning, intensity_mw_cm2=0., density_m3=DENSITY_M3,
                           ground_relaxation_hz=RESET_PER_S)
    np.testing.assert_allclose(state["rho"], np.diag([1/8]*8 + [0.]*8), atol=1e-13)
    np.testing.assert_array_equal(state["polarization_xy_c_m2"], 0.)
    assert state["excited_fraction"] == pytest.approx(0., abs=1e-14)

    # Sum unpumped line strengths across the two hyperfine components. A full
    # thermal replacement dephases rho_eg at r, giving HWHM Gamma/2 + r.
    atom = magneto_multilevel.build_multilevel_atom("Rb87", "D1")
    ground_energies = np.repeat(atom.ground_energy_mhz, [3, 5])
    transition_offsets = atom.excited_energy_mhz[np.newaxis, :] - ground_energies[:, None]
    laser_mhz = detuning + magneto_multilevel.component_centroid_mhz(atom)
    widths = (transition_offsets - laser_mhz)*species.MHZ - 1j*(atom.gamma/2 + RESET_PER_S)
    x_dipole = (atom.dipoles[-1] + atom.dipoles[1])/np.sqrt(2)
    susceptibility = (DENSITY_M3*atom.dipole_sq/(constants.EPS_0*constants.HBAR)
                      * np.sum(abs(x_dipole)**2/widths)/8)
    for key in ("chi_plus", "chi_minus", "chi_probe"):
        np.testing.assert_allclose(state[key], susceptibility, rtol=2e-13)
    small_drive = solve_psr_pump(detuning, intensity_mw_cm2=1e-11, density_m3=DENSITY_M3,
                                 ground_relaxation_hz=RESET_PER_S)
    np.testing.assert_allclose(small_drive["chi_probe"], susceptibility, rtol=2e-8)


def test_stationary_lindblad_balance_and_absorbed_power():
    intensity = 12.5
    reset = 1e5
    state = solve_psr_pump(400., intensity_mw_cm2=intensity, ground_relaxation_hz=reset,
                           jones_xy=(3/5, 4j/5), density_m3=DENSITY_M3)
    rho = state["rho"]
    np.testing.assert_allclose(rho, rho.conj().T, atol=2e-12)
    assert rho.trace() == pytest.approx(1., abs=2e-12)
    assert np.linalg.eigvalsh(rho)[0] >= -1e-12
    excited_fraction = rho[8:, 8:].trace().real
    assert excited_fraction == pytest.approx(state["excited_fraction"], rel=1e-12)

    # Check the unvectorized master equation directly from the named jumps.
    # Includes their coherence decay as well as level-to-level repopulation.
    hamiltonian = state["hamiltonian_rad_s"]
    residual = -1j*(hamiltonian @ rho - rho @ hamiltonian)
    excited_projector = np.diag([0.]*8 + [1.]*8)
    emission_metric = np.zeros_like(rho)
    for channel in state["reservoirs"].channels:
        jump = channel.operator
        metric = jump.conj().T @ jump
        residual += jump @ rho @ jump.conj().T - .5*(metric @ rho + rho @ metric)
        if channel.kind == "spontaneous_emission":
            emission_metric += metric
    atom = species.level_structure("Rb87", "D1")
    np.testing.assert_allclose(emission_metric/atom.gamma, excited_projector, atol=1e-14)
    assert np.linalg.norm(residual) < 1e-9*reset

    # The steady pump replaces all excited-state losses. Both axially
    # absorbed power and independent atomic decay power have units W/m^3.
    absorbed_power = state["alpha_per_m"] * (10*intensity)
    photon_energy = 2*np.pi*constants.HBAR*atom.frequency_hz
    decay_power = DENSITY_M3*photon_energy*(atom.gamma + reset)*excited_fraction
    assert absorbed_power == pytest.approx(decay_power, rel=1e-9)

    # A coherent superposition of excited F=1 and F=2 paths can interfere
    # in each q-resolved decay into one ground F. Resolving jumps separately
    # for every excited F would still pass the total decay closure above.
    superposition = np.zeros(16, complex)
    for excited_f in (1, 2):
        sublevel = atom.excited.index((excited_f, 1))
        superposition[8 + sublevel] = 1/np.sqrt(2)
    emitted = np.zeros((8, 8), complex)
    for channel in state["reservoirs"].channels:
        if channel.kind == "spontaneous_emission":
            after = (channel.operator @ superposition)[:8]
            emitted += np.outer(after, after.conj())/atom.gamma
    expected = np.zeros_like(emitted)
    for q in (-1, 0, 1):
        after = atom.dipoles[q] @ superposition[8:]
        for ground_f in (1, 2):
            mask = np.array([f == ground_f for f, _ in atom.ground])
            resolved = np.where(mask, after, 0.)
            expected += np.outer(resolved, resolved.conj())
    np.testing.assert_allclose(emitted, expected, atol=2e-15)


@pytest.mark.parametrize("detuning", [-600., 350.])
@pytest.mark.parametrize("intensity", [.0001, .01])
def test_weak_overlap_matches_susceptibilities_absorption_and_rotation(detuning, intensity):
    common = dict(density_m3=DENSITY_M3, ground_relaxation_hz=RESET_PER_S,
                  intensity_mw_cm2=intensity)
    step_deg = .05
    weak = psr.compute_psr_response([detuning], doppler=False, cell_length_m=.01,
                                    ellipticity_step_deg=step_deg, **common)
    assert weak["conditional_model_applicable"][0]
    linear = solve_psr_pump(detuning, **common)
    np.testing.assert_allclose(linear["alpha_per_m"], weak["alpha_per_m"][0], rtol=2e-4)

    # Perturb pure ellipticity with fixed major-axis orientation. The weak
    # QWP stencil also rotates the ellipse, which is equivalent only at B=0
    # in this isotropic model. QWP beta=±h/2 produces epsilon=∓h/2.
    half_step = np.deg2rad(step_deg/2)
    rotation_rates = []
    for sign, weak_index in ((-1, 3), (1, 4)):
        jones = (np.cos(half_step), sign*1j*np.sin(half_step))
        full = solve_psr_pump(detuning, jones_xy=jones, **common)
        for key in ("chi_plus", "chi_minus"):
            expected = weak["stencil"][key][0, weak_index]
            # Separate comparisons keep a small absorptive part observable.
            np.testing.assert_allclose(full[key].real, expected.real, rtol=2e-4)
            np.testing.assert_allclose(full[key].imag, expected.imag, rtol=2e-4)
        rotation_rates.append(full["azimuth_rate_per_m"])
    slope = (rotation_rates[1] - rotation_rates[0])/(2*half_step)
    # Full OBE retains inter-ground-F coherences and optical reset dephasing.
    # At these four points the largest relative slope mismatch is ~2.1e-4.
    np.testing.assert_allclose(slope, weak["g_per_m"][0], rtol=4e-4)


def test_strong_drive_saturates_the_isolated_d1_manifold():
    # Numerical limiting case of a 16-state D1 model. Linear polarization
    # couples both entire J=1/2 manifolds. A drive dominating hyperfine
    # splitting equalizes total excited/ground populations.
    mid = solve_psr_pump(407.24, intensity_mw_cm2=1e7, density_m3=DENSITY_M3)
    high = solve_psr_pump(407.24, intensity_mw_cm2=1e9, density_m3=DENSITY_M3)
    assert .20 < mid["excited_fraction"] < .35
    assert .49 < high["excited_fraction"] < .5
    assert high["alpha_per_m"] < mid["alpha_per_m"]/40


def test_moving_atom_matches_its_rest_frame_with_separate_density_scaling():
    atom = species.level_structure("Rb87", "D1")
    velocity = 52.
    detuning = -390.
    shifted = detuning - atom.frequency_hz*velocity/constants.C_LIGHT/1e6
    moving = solve_psr_pump(detuning, velocity_m_s=velocity, intensity_mw_cm2=1.,
                            density_m3=DENSITY_M3)
    stationary = solve_psr_pump(shifted, intensity_mw_cm2=1., density_m3=2*DENSITY_M3)
    np.testing.assert_allclose(moving["rho"], stationary["rho"], atol=1e-12)
    assert stationary["alpha_per_m"] == pytest.approx(2*moving["alpha_per_m"], rel=1e-10)
    np.testing.assert_allclose(stationary["polarization_xy_c_m2"],
                               2*moving["polarization_xy_c_m2"], atol=1e-26, rtol=1e-10)
