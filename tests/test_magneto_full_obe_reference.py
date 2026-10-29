"""Selected weak-drive checks against an explicit ground+excited-state OBE.

The reference shares atomic angular data, but constructs its own full density
matrix Hamiltonian, emission dissipator, optical coherences and field readout.
It does not use effective-operator or dark-reservoir implementation helpers.
These isolated Doppler-free points do not certify coated-cell spectra or
high-intensity optical elimination.
"""
import numpy as np
import pytest
from scipy.linalg import solve

from gabes import constants, core, species
from gabes.magneto_multilevel import (
    build_multilevel_atom, component_centroid_mhz, multilevel_response,
)


def _hyperfine_g(I, J, F, gj, gi):
    ff = F * (F + 1)
    if not ff:
        return 0.0
    return (gj * (ff + J * (J + 1) - I * (I + 1))
            + gi * (ff + I * (I + 1) - J * (J + 1))) / (2 * ff)


def _full_obe_rotation(atom, delta_mhz, intensity_mw_cm2, *, field_ut,
                       reset_rate, density, length):
    """Return full rotation and its exact zero-intensity linear background."""
    iso = species.RB87
    ng, ne = len(atom.ground), len(atom.excited)
    n = ng + ne
    eye = np.eye(n)
    ground_energy = np.concatenate([
        np.full(block.stop - block.start, energy)
        for block, energy in zip(atom.ground_slices, atom.ground_energy_mhz)
    ])
    laser_mhz = component_centroid_mhz(atom) + delta_mhz
    energies = np.r_[ground_energy, atom.excited_energy_mhz - laser_mhz] * species.MHZ
    # Include optical-transition Zeeman shifts in this full reference. Their
    # intensity-independent Faraday background is subtracted below.
    gi = -0.0009951414
    for index, (F, m) in enumerate(atom.ground):
        energies[index] += (constants.MU_B / constants.HBAR * field_ut * 1e-6 * m
                            * _hyperfine_g(iso.I, 0.5, F, 2.002331070, gi))
    for index, (F, m) in enumerate(atom.excited):
        energies[ng + index] += (constants.MU_B / constants.HBAR * field_ut * 1e-6 * m
                                 * _hyperfine_g(iso.I, 0.5, F, 0.666, gi))
    hamiltonian = np.diag(energies).astype(complex)
    field = np.sqrt(2 * intensity_mw_cm2 * 10 / (constants.EPS_0 * constants.C_LIGHT))
    rabi_scale = np.sqrt(atom.dipole_sq) * field / constants.HBAR
    polarization = 1 / np.sqrt(2)
    excitation = polarization * (atom.dipoles[1].T + atom.dipoles[-1].T)
    hamiltonian[ng:, :ng] = -0.5 * rabi_scale * excitation
    hamiltonian[:ng, ng:] = hamiltonian[ng:, :ng].conj().T
    liouvillian = -1j * (np.kron(hamiltonian, eye)
                         - np.kron(eye, hamiltonian.T))
    # Photons ending in distinct ground hyperfine manifolds have distinguishable
    # frequencies. Split those channels, retaining coherent virtual Fe paths
    # within each fixed final-F and emitted-polarization channel.
    for target in atom.ground_slices:
        for q in (-1, 0, 1):
            jump = np.zeros((n, n), complex)
            jump[target, ng:] = np.sqrt(atom.gamma) * atom.dipoles[q][target]
            loss = jump.conj().T @ jump
            liouvillian += (np.kron(jump, jump.conj())
                            - 0.5 * np.kron(loss, eye)
                            - 0.5 * np.kron(eye, loss.T))
    thermal = np.zeros((n, n))
    thermal[:ng, :ng] = np.eye(ng) / ng
    liouvillian += reset_rate * (
        np.outer(thermal.ravel(), eye.ravel()) - np.eye(n * n))
    liouvillian[0] = eye.ravel()
    rhs = np.zeros(n * n)
    rhs[0] = 1
    with core.blas_single_thread():
        rho = solve(liouvillian, rhs).reshape(n, n)
    np.testing.assert_allclose(rho, rho.conj().T, atol=2e-12)
    assert np.linalg.eigvalsh(rho).min() >= -1e-12
    assert np.trace(rho) == pytest.approx(1, abs=2e-12)
    prefactor = 2 * density * np.sqrt(atom.dipole_sq) / (constants.EPS_0 * field)
    chi = {q: prefactor * np.sum(atom.dipoles[q] * rho[ng:, :ng].T) / polarization
           for q in (-1, 1)}
    wave_number = 2 * np.pi * atom.frequency_hz / constants.C_LIGHT
    # At zero drive the ground state is exactly isotropic. Each first-order
    # optical coherence has this independent scalar solution to the full OBE.
    # Include the reset-induced optical decoherence, as in the full generator.
    optical_energy = energies[ng:][None, :] - energies[:ng, None]
    denominator = optical_energy - 1j * (atom.gamma / 2 + reset_rate)
    linear_prefactor = density * atom.dipole_sq / (constants.EPS_0 * constants.HBAR * ng)
    chi_linear = {q: linear_prefactor * np.sum(atom.dipoles[q] ** 2 / denominator)
                  for q in (-1, 1)}
    return (0.25 * wave_number * length * np.real(chi[1] - chi[-1]),
            0.25 * wave_number * length * np.real(chi_linear[1] - chi_linear[-1]))


@pytest.mark.parametrize("detuning_mhz", [100.0, -250.0])
def test_effective_nonlinear_rotation_matches_selected_full_obe_points(detuning_mhz):
    atom = build_multilevel_atom("Rb87", "D1")
    intensity = 1e-4
    settings = dict(field_ut=1e-3, reset_rate=100.0, density=1e15, length=0.1)
    full, linear_background = _full_obe_rotation(atom, detuning_mhz, intensity, **settings)
    nonlinear = full - linear_background
    effective = multilevel_response(
        [detuning_mhz], isotope="Rb87", line="D1", intensity_mw_cm2=intensity,
        b_ut=settings["field_ut"], gamma_out_hz=0, doppler=False,
        ground_relaxation_hz=settings["reset_rate"],
        density_m3=settings["density"], cell_length_m=settings["length"],
    )["rotation_rad"][0, 0]
    # Subtract the exact linear response to isolate nonlinear rotation. No
    # amplitude or detuning calibration is applied.
    assert abs(nonlinear) > 1e-13
    assert effective == pytest.approx(nonlinear, rel=0.01)
