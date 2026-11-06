"""Saturated, single-velocity Rb87 D1 mean field for polarization self-rotation.

The 16-state rotating-wave Hamiltonian retains ground- and excited-hyperfine
coherences. Emission distinguishes final ground F and spherical polarization
but coherently sums excited-F and magnetic-sublevel paths, as in the weak PSR
pumping model. Thermal replacement acts on the entire density operator. The
state, Hamiltonian and named jumps are available for subsequent linearization.
"""
from functools import lru_cache

import numpy as np

from . import constants, core, magneto_multilevel, species
from .psr_conventions import atom_frame_detuning_mhz, circular_to_jones, jones_to_circular, stokes_rates
from .quantum.reservoirs import CollapseChannel, ExplicitReservoirs, thermal_reset_channels


MODEL_ID = "psr_rb87_d1_single_velocity_full_obe_v1"


@lru_cache(maxsize=4)
def _prepared_reservoirs(ground_relaxation_hz):
    """Keep the dissipator fixed while scanning the Hamiltonian or polarization."""
    atom = magneto_multilevel.build_multilevel_atom("Rb87", "D1")
    n_ground = len(atom.ground)
    n_levels = n_ground + len(atom.excited)
    channels = []
    for ground_slice in atom.ground_slices:
        ground_f = atom.ground[ground_slice.start][0]
        for q in (-1, 0, 1):
            jump = np.zeros((n_levels, n_levels), complex)
            jump[ground_slice, n_ground:] = np.sqrt(atom.gamma) * atom.dipoles[q][ground_slice]
            channels.append(CollapseChannel(
                f"emission:F={ground_f}:q={q}", jump, "spontaneous_emission",
                "Rb87 D1 signed dipoles and natural width from species.level_structure",
            ))
    thermal_populations = np.zeros(n_levels)
    thermal_populations[:n_ground] = 1.0 / n_ground
    channels.extend(thermal_reset_channels(
        ground_relaxation_hz, thermal_populations, name="reset",
        source="Isotropic transit replacement into the eight Rb87 ground sublevels",
    ))
    reservoirs = ExplicitReservoirs(n_levels, tuple(channels))
    dissipator = reservoirs.dissipator()
    dissipator.flags.writeable = False
    return reservoirs, dissipator


def _hamiltonian(atom, laser_from_centroid_mhz, rabi_rad_s, excitation):
    n_ground = len(atom.ground)
    ground_energy = np.empty(n_ground)
    for ground_slice, energy in zip(atom.ground_slices, atom.ground_energy_mhz):
        ground_energy[ground_slice] = energy
    excited_energy = atom.excited_energy_mhz - laser_from_centroid_mhz
    # Subtracting a scalar identity places the upper ground F at zero energy.
    energies = np.concatenate((ground_energy, excited_energy)) - ground_energy[-1]
    hamiltonian = np.diag(energies * species.MHZ).astype(complex)
    coupling = -0.5 * rabi_rad_s * excitation
    hamiltonian[n_ground:, :n_ground] = coupling
    hamiltonian[:n_ground, n_ground:] = coupling.conj().T
    return hamiltonian


def _state_diagnostics(generator, rho, reset_rate):
    vector = rho.ravel()
    residual = float(np.linalg.norm(generator @ vector))
    # Scale by terms actually present in L rho. In particular, a large bare
    # hyperfine splitting acting on a zero coherence cannot mask a bad solve.
    scale = max(float(np.linalg.norm(np.abs(generator) @ np.abs(vector))),
                reset_rate * float(np.linalg.norm(vector)))
    diagnostics = {
        "trace_error": float(abs(np.trace(rho) - 1)),
        "hermiticity_error": float(np.linalg.norm(rho - rho.conj().T)),
        "min_eigenvalue": float(np.linalg.eigvalsh((rho + rho.conj().T) / 2)[0]),
        "steady_state_residual_s_inverse": residual,
        "steady_state_relative_residual": residual / scale,
    }
    if (not all(np.isfinite(list(diagnostics.values())))
            or diagnostics["trace_error"] > 1e-9
            or diagnostics["hermiticity_error"] > 1e-9
            or diagnostics["min_eigenvalue"] < -1e-9
            or diagnostics["steady_state_relative_residual"] > 1e-9):
        raise np.linalg.LinAlgError(f"Rb87 D1 PSR full-OBE steady-state check failed: {diagnostics}")
    diagnostics["solver_status"] = "physical steady state"
    return diagnostics


def _linear_coherence_per_half_rabi(atom, excitation, laser_mhz, reset_rate):
    """Exact first-order rho_eg / (Omega/2) at zero drive and positive reset."""
    n_ground = len(atom.ground)
    coherence = np.empty_like(excitation)
    for ground_slice, energy in zip(atom.ground_slices, atom.ground_energy_mhz):
        detuning = (laser_mhz + energy - atom.excited_energy_mhz) * species.MHZ
        denominator = -detuning - 1j * (atom.gamma / 2 + reset_rate)
        coherence[:, ground_slice] = excitation[:, ground_slice] / (n_ground * denominator[:, None])
    return coherence


def solve_psr_pump(detuning_mhz, *, intensity_mw_cm2=0.1, jones_xy=(1., 0.),
                   velocity_m_s=0., temp_c=70., ground_relaxation_hz=1000., density_m3=None):
    """Solve a local stationary pump point for pure Rb87 D1 at zero B.

    Positive detuning is blue of the upper-ground-F line-strength centroid.
    Positive velocity is along the incident beam. Intensity is the total
    local cycle-averaged power per area and ``jones_xy`` is a complex unit
    vector. Temperature supplies the vapor density when ``density_m3`` is
    omitted. Each call describes one velocity, without a Maxwell average.

    Ground states precede excited states, each in ``species.level_structure``
    order. Hamiltonian entries are angular frequencies in rad/s, the reset
    rate is in s^-1. ``polarization_xy_c_m2`` is the macroscopic peak phasor
    with physical polarization Re[P exp(-i omega t)]. At zero intensity the
    polarization vanishes and chi is its analytic infinitesimal-probe limit.
    Circular component ratios are None when the corresponding field is absent.
    Azimuth/ellipticity angle rates are undefined at exactly circular input.
    """
    detuning, intensity, velocity, temperature, reset = map(float, (
        detuning_mhz, intensity_mw_cm2, velocity_m_s, temp_c, ground_relaxation_hz))
    if not np.all(np.isfinite([detuning, intensity, velocity, temperature, reset])):
        raise ValueError("finite scalar full-OBE drive, detuning, velocity, temperature and reset required")
    if intensity < 0 or temperature <= -273.15 or reset <= 0:
        raise ValueError("nonnegative intensity, positive temperature and positive thermal reset required")
    density = (species.number_density(species.ISOTOPES["Rb87"], temperature + 273.15)
               if density_m3 is None else float(density_m3))
    if not np.isfinite(density) or density < 0:
        raise ValueError("density_m3 must be finite and nonnegative")

    jones = np.asarray(jones_xy, dtype=complex)
    circular = jones_to_circular(jones)
    atom = magneto_multilevel.build_multilevel_atom("Rb87", "D1")
    centroid = magneto_multilevel.component_centroid_mhz(atom)
    effective_detuning = atom_frame_detuning_mhz(
        detuning, velocity_m_s=velocity, optical_frequency_hz=atom.frequency_hz)
    laser_mhz = centroid + effective_detuning
    # E is the peak field. 1 mW/cm² = 10 W/m², Omega = d0 E / hbar.
    peak_field = np.sqrt(2 * intensity * 10 / (constants.EPS_0 * constants.C_LIGHT))
    rabi = np.sqrt(atom.dipole_sq) * peak_field / constants.HBAR
    excitation = sum(circular[q] * atom.dipoles[q].T for q in (-1, 1))
    hamiltonian = _hamiltonian(atom, laser_mhz, rabi, excitation)
    reservoirs, dissipator = _prepared_reservoirs(reset)
    generator = core.comm_super(hamiltonian) + dissipator
    if not np.isfinite(generator).all():
        raise ValueError("settings exceed the finite full-OBE generator range")
    with core.blas_single_thread():
        rho = core.steady_state_from_liouvillian(generator, reservoirs.n_levels)
        diagnostics = _state_diagnostics(generator, rho, reset)

    n_ground = len(atom.ground)
    optical_eg = rho[n_ground:, :n_ground]
    coherence = (optical_eg / (0.5 * rabi) if intensity > 0 else
                 _linear_coherence_per_half_rabi(atom, excitation, laser_mhz, reset))
    prefactor = density * atom.dipole_sq / (constants.EPS_0 * constants.HBAR)
    # P_q / (epsilon0 * |E|), defined even when the incident E_q is absent.
    induced = {q: prefactor * np.einsum("ge,eg->", atom.dipoles[q], coherence) for q in (-1, 1)}
    chi = {q: complex(induced[q] / circular[q]) if abs(circular[q]) > 1e-12 else None for q in (-1, 1)}
    chi_probe = sum(np.conj(circular[q]) * induced[q] for q in (-1, 1))
    induced_xy = circular_to_jones(induced[1], induced[-1])
    wave_number = 2 * np.pi * atom.frequency_hz / constants.C_LIGHT
    epsilon, rotation, ellipticity_rate, alpha = stokes_rates(jones, 0.5j * wave_number * induced_xy)
    return {
        "model_id": MODEL_ID,
        "rho": rho, "hamiltonian_rad_s": hamiltonian, "generator_s_inverse": generator,
        "reservoirs": reservoirs, "optical_coherence_eg": optical_eg.copy(),
        "ground_levels": atom.ground, "excited_levels": atom.excited,
        "ground_populations": rho.diagonal()[:n_ground].real.copy(),
        "excited_populations": rho.diagonal()[n_ground:].real.copy(),
        "excited_fraction": float(np.trace(rho[n_ground:, n_ground:]).real),
        "chi_plus": chi[1], "chi_minus": chi[-1], "chi_probe": complex(chi_probe),
        "polarization_xy_c_m2": constants.EPS_0 * peak_field * induced_xy,
        "alpha_per_m": alpha, "ellipticity_rad": epsilon,
        "azimuth_rate_per_m": rotation, "ellipticity_rate_per_m": ellipticity_rate,
        "diagnostics": diagnostics,
        "conditions": {
            "isotope": "Rb87", "line": "D1", "b_ut": 0., "velocity_m_s": velocity,
            "detuning_mhz": detuning, "effective_detuning_mhz": effective_detuning,
            "reference_centroid_mhz": centroid, "optical_frequency_hz": atom.frequency_hz,
            "intensity_mw_cm2": intensity, "jones_xy": jones.copy(),
            "peak_field_v_m": float(peak_field), "rabi_rad_s": float(rabi),
            "temp_c": temperature, "ground_relaxation_hz": reset, "density_m3": float(density),
            "density_source": "pure-isotope vapor-pressure density" if density_m3 is None else "specified isotope density",
            "ground_f_coherences": "retained", "emission_secularization": "final Fg and q resolved, excited F coherent",
            "susceptibility_basis": "finite-drive induced polarization" if intensity > 0 else "zero-drive linear-probe limit",
        },
    }
