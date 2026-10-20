"""Ordered-noise absorption/emission for the declared reduced FWM atom.

This adapter reconstructs exactly the production four-level collisional
Lindblad model as explicit jumps, then uses the existing periodic quantum
regression engine.  Positive absorption and emission spectra are kept separate;
net diagonal gain is not clipped into a fictitious passive transmission.
"""

import copy
import functools
import math
import time

import numpy as np
from scipy.integrate import quad_vec
from scipy.optimize import brentq

from . import atoms, constants, core, hyperfine, observables
from .quantum.contracts import GeneratorFrequencyAxis
from .quantum.diffusion import traceless_hermitian_basis
from .quantum.periodic import periodic_atomic_noise
from .quantum.periodic_field import PeriodicFieldPorts, eliminate_periodic_atom
from .quantum.reservoirs import (
    CollapseChannel, ExplicitReservoirs, thermal_reset_channels)


def production_collisional_reservoirs(T, *, density=None, transit_rate=None):
    """Return explicit jumps exactly matching ``fwm.collisional_atom``.

    Ground and optical coherence dephasing rates obey gamma_opt>=gamma_gg/4
    in this construction.  Outside that domain the function rejects the model;
    it does not repair a nonpositive or differently parameterized dissipator.
    The thermal reset remains the prescribed Markov replacement approximation.
    """
    from .schemes import fwm

    temperature = float(T)
    if not np.isfinite(temperature) or temperature <= 0.0:
        raise ValueError("T must be finite and positive in kelvin")
    number_density = float(hyperfine.number_density(temperature) if density is None else density)
    reset_rate = float(constants.GAMMA_GG if transit_rate is None else transit_rate)
    if not np.isfinite([number_density, reset_rate]).all() or number_density < 0 or reset_rate < 0:
        raise ValueError("density and transit_rate must be finite and nonnegative")
    gamma_ground = float(constants.ground_coherence_dephasing(
        temperature, number_density, floor=0.0))
    gamma_optical = max(float(0.5 * (
        hyperfine.self_broadened_gamma(number_density)-constants.GAMMA)), 0.0)
    margin = gamma_optical-gamma_ground/4.0
    if margin < 0.0:
        raise ValueError("declared coherence dephasing has no positive two-jump construction")
    radiative = ExplicitReservoirs.from_atom(
        atoms.double_lambda_rb85(gamma_gg=0.0),
        source="same reduced D1 spontaneous-emission channels as production")
    channels = list(radiative.channels)
    if gamma_ground:
        channels.append(CollapseChannel(
            "collision:ground", np.sqrt(gamma_ground/2)*np.diag([1., -1., 0., 0.]),
            "pure_dephasing", "exact declared ground/optical coherence dephasing reconstruction"))
    if margin:
        channels.append(CollapseChannel(
            "collision:optical", np.sqrt(2*margin)*np.diag([0., 0., 1., 1.]),
            "pure_dephasing", "exact declared ground/optical coherence dephasing reconstruction"))
    channels.extend(thermal_reset_channels(
        reset_rate, [5/12, 7/12, 0., 0.],
        source="same prescribed Markov thermal replacement as production"))
    reservoirs = ExplicitReservoirs(4, tuple(channels))
    production = fwm.collisional_atom(temperature, number_density, transit_rate=reset_rate)
    scale = max(np.linalg.norm(production.lindblad), np.finfo(float).tiny)
    residual = float(np.linalg.norm(reservoirs.dissipator()-production.lindblad)/scale)
    if residual > 1e-12:
        raise ValueError("explicit jumps do not reproduce the production dissipator")
    return reservoirs, {
        "dissipator_relative_error": residual,
        "ground_collision_dephasing_s_inverse": gamma_ground,
        "optical_collision_dephasing_s_inverse": gamma_optical,
        "positive_dephasing_construction_margin_s_inverse": margin,
        "transit_reset_s_inverse": reset_rate,
        "reservoir_scope": "declared four-level radiative, collisional and Markov-reset model",
    }


def ordered_absorption_emission(local_generator, *, coordinate_indices=(0, 1)):
    """Extract positive ordered rates and independently check net-drift identity.

    For an annihilator alpha_abs=D_greater[j,j], alpha_em=D_lesser[j,j];
    for a creator the two orderings reverse.  Rates are per metre.  Correlated
    cross-port and anomalous noise remain in the generator; two scalar rates
    alone do not reproduce its complete quantum channel.
    """
    audit = local_generator.audit()
    if not audit["passed"]:
        raise ValueError("validated ordered microscopic field generator required")
    greater, lesser = local_generator.noise_greater, local_generator.noise_lesser
    result = []
    for coordinate in coordinate_indices:
        if (isinstance(coordinate, bool) or int(coordinate) != coordinate
                or not 0 <= coordinate < len(local_generator.signs)):
            raise ValueError("coordinate_indices must identify physical field coordinates")
        coordinate = int(coordinate)
        if local_generator.signs[coordinate] == 1:
            absorption, emission = greater[:, coordinate, coordinate], lesser[:, coordinate, coordinate]
        else:
            absorption, emission = lesser[:, coordinate, coordinate], greater[:, coordinate, coordinate]
        absorption, emission = absorption.real, emission.real
        net = -2*np.real(local_generator.drift[:, coordinate, coordinate])
        scale = np.maximum.reduce((abs(absorption), abs(emission), abs(net),
                                   np.full_like(net, np.finfo(float).tiny)))
        residual = abs(absorption-emission-net)/scale
        if (np.any(absorption < -1e-9*scale) or np.any(emission < -1e-9*scale)
                or np.any(residual > 1e-8)):
            raise ValueError("ordered absorption/emission positivity or Kubo identity failed")
        result.append({
            "alpha_abs_per_m": absorption,
            "alpha_em_per_m": emission,
            "alpha_net_per_m": net,
            "kubo_relative_residual": residual,
        })
    return tuple(result)


def _block_tridiagonal_action(diagonal, lower, upper, rhs):
    """Block Thomas solve for connected QRT initial covariance sources."""
    solved_rhs, solved_upper = [], []
    for index, center in enumerate(diagonal):
        center, target = center.copy(), rhs[index].copy()
        if index:
            center -= lower@solved_upper[-1]
            target -= lower@solved_rhs[-1]
        solution = np.linalg.solve(center, np.concatenate((target, upper), axis=-1))
        solved_rhs.append(solution[:, :target.shape[-1]])
        solved_upper.append(solution[:, target.shape[-1]:])
    result = [None]*len(diagonal)
    result[-1] = solved_rhs[-1]
    for index in range(len(diagonal)-2, -1, -1):
        result[index] = solved_rhs[index]-solved_upper[index]@result[index+1]
    return np.stack(result)


def _direct_periodic_spectra(h0, harmonic, reservoirs, beat, states,
                             lowering_operators, response_order):
    """Connected QRT spectra by two covariance sources, not a covariance fit.

    The initial sources C_h and C_h^T retain the convolution of periodic
    means through twice the mean harmonic order.  A physical lowering port
    reads harmonic zero at its physical pump-frame frequency +/-beat.
    This is equivalent to the independently audited R D R-dagger spectrum
    after harmonic convergence, while solving only two source columns.
    """
    operators = traceless_hermitian_basis(4)
    basis = operators.reshape(15, 16).T
    generators = (core.comm_super(harmonic.conj().T), reservoirs.generator(h0),
                  core.comm_super(harmonic))
    drift_minus, drift_zero, drift_plus = (basis.conj().T@g@basis for g in generators)
    means = np.einsum("iab,qba->qi", operators, states)
    moments = np.einsum("ijab,qba->qij", operators[:, None]@operators[None, :], states)
    mean_order = len(states)//2
    covariance = {}
    for harmonic_index in range(-2*mean_order, 2*mean_order+1):
        item = np.zeros((15, 15), complex)
        if abs(harmonic_index) <= mean_order:
            item += moments[harmonic_index+mean_order]
        for left in range(-mean_order, mean_order+1):
            right = harmonic_index-left
            if abs(right) <= mean_order:
                item -= np.outer(means[left+mean_order], means[right+mean_order])
        covariance[harmonic_index] = item
    ladder = np.arange(-response_order, response_order+1)
    zero, eye = np.zeros((15, 15), complex), np.eye(15)
    result = []
    for lowering, frequency in zip(lowering_operators, (beat, -beat)):
        readout = np.einsum("ab,iba->i", lowering, operators)
        sources = []
        for h in ladder:
            # Build the retarded commutator source independently of the two
            # ordered covariance sources used for absorption/emission.
            rho = states[h+mean_order] if abs(h) <= mean_order else np.zeros((4, 4), complex)
            drive = -.5j*(lowering.conj().T@rho-rho@lowering.conj().T)
            response_source = np.einsum("iab,ba->i", operators, drive)
            sources.append(np.column_stack((covariance.get(int(h), zero)@readout.conj(),
                                           covariance.get(int(h), zero).T@readout.conj(),
                                           response_source)))
        source = np.asarray(sources)
        diagonal = drift_zero[None]+1j*(frequency+ladder*beat)[:, None, None]*eye
        response = _block_tridiagonal_action(diagonal, drift_plus, drift_minus, -source)
        projected = readout@response[response_order]
        spectra = 2*np.real(projected[:2])
        # Independent commutator response per classical Rabi, rather than a
        # net rate inferred by subtracting the two ordered spectra.
        chi = projected[2]
        result.append((spectra[0], spectra[1], -4*chi.imag))
    return np.asarray(result)


def local_dressed_absorption(*, T, pump_rabi_rad_s, probe_rabi_rad_s,
                            one_photon_detuning_rad_s, two_photon_detuning_rad_s,
                            optical_omega_rad_s, uniform_area_m2,
                            weak_line_strength, number_density_m3=None,
                            transit_rate=None, conjugate_rabi_rad_s=0.0,
                            mean_order=3, response_order=3, phase_samples=96,
                            method="ordered-diffusion"):
    """Calculate both-arm absorption/emission in one prescribed finite-seed state.

    Pump/probe/conjugate Rabi frequencies are explicit local classical inputs,
    allowing exact reproduction of the current legacy drive normalization.
    ``weak_line_strength`` is the same macroscopic factor multiplying d_D1^2
    in production susceptibility propagation, including its current angular
    normalization.  It is not fitted or inferred from a measured transmission.

    The minus-branch physical pump frame and two retained optical bands are
    used.  One local homogeneous state, no Doppler/spatial averaging, optical
    propagation or pump depletion is implied.  The finite-seed response is a
    differential QRT response, whereas production finite-drive chi ratios are
    secants from two different seeded states.
    """
    from .schemes import fwm

    rabis = np.asarray((pump_rabi_rad_s, probe_rabi_rad_s, conjugate_rabi_rad_s), float)
    detunings = np.asarray((one_photon_detuning_rad_s, two_photon_detuning_rad_s), float)
    omega = np.asarray(optical_omega_rad_s, float)
    area, strength = float(uniform_area_m2), float(weak_line_strength)
    density = float(hyperfine.number_density(T) if number_density_m3 is None else number_density_m3)
    if (not np.isfinite(rabis).all() or np.any(rabis < 0.0)
            or not np.isfinite(detunings).all() or omega.shape != (2,)
            or not np.isfinite(omega).all() or np.any(omega <= 0.0)
            or not np.isfinite([area, strength, density]).all()
            or area <= 0.0 or strength <= 0.0 or density < 0.0):
        raise ValueError("finite nonnegative drives/density, positive frequencies/area/weak strength required")
    if method not in ("ordered-diffusion", "connected-qrt"):
        raise ValueError("method must be ordered-diffusion or connected-qrt")
    for value, name in ((mean_order, "mean_order"), (response_order, "response_order")):
        if isinstance(value, bool) or int(value) != value or value < 1:
            raise ValueError(f"{name} must be a positive integer")
    mean_order, response_order = int(mean_order), int(response_order)
    if (isinstance(phase_samples, bool) or int(phase_samples) != phase_samples
            or phase_samples < 4*mean_order+3):
        raise ValueError("phase_samples must resolve the mean-state harmonics")
    reservoirs, reservoir_info = production_collisional_reservoirs(
        T, density=density, transit_rate=transit_rate)
    h0 = fwm.pump_hamiltonian_at_Deff_zero(rabis[0], rabis[0])
    h0 -= detunings[0]*np.diag([0., 0., 1., 1.])
    beat = -fwm.OMEGA_HF+detunings[1]
    if beat == 0.0:
        raise ValueError("nonzero minus-branch pump/seed beat required")
    lowering = np.zeros((2, 4, 4), complex)
    for excited in fwm.EXCITED_STATES:
        lowering[0, fwm.G2, excited] = fwm.TRANSITION_DIPOLE_SCALE[fwm.G2, excited]
        lowering[1, fwm.G1, excited] = fwm.TRANSITION_DIPOLE_SCALE[fwm.G1, excited]
    harmonic = (rabis[1]*lowering[0].conj().T+rabis[2]*lowering[1])/2
    q = np.diag(observables.photon_flux_mode_matrix(*omega, area, area))
    coupling = constants.DIPOLE_D1*np.sqrt(strength)*q/(2*constants.HBAR)
    if method == "connected-qrt":
        # A positive Markov replacement rate contracts every traceless mode;
        # no repeated numerical period-map integration is needed in this path.
        if reservoir_info["transit_reset_s_inverse"] <= 0:
            raise ValueError("connected-qrt requires positive thermal reset for its stability bound")
        generator = reservoirs.generator(h0)
        plus, minus = core.comm_super(harmonic), core.comm_super(harmonic.conj().T)
        states = core.floquet_solve_truncated(generator, plus, minus, beat,
            [0.], np.zeros_like(generator), 4, n_f=mean_order, return_harmonics=True)[0]
        phases = 2*np.pi*np.arange(int(phase_samples))/phase_samples
        reconstructed = np.einsum("tq,qab->tab", np.exp(-1j*np.outer(
            phases, np.arange(-mean_order, mean_order+1))), states)
        trace_error = float(np.max(abs(np.trace(reconstructed, axis1=-2, axis2=-1)-1)))
        hermiticity = float(np.max(np.linalg.norm(
            reconstructed-reconstructed.conj().swapaxes(-1, -2), axis=(-2, -1))))
        minimum = float(np.linalg.eigvalsh(
            (reconstructed+reconstructed.conj().swapaxes(-1, -2))/2).min())
        if trace_error > 1e-9 or hermiticity > 1e-9 or minimum < -1e-10:
            raise ValueError("finite periodic mean failed trace, Hermiticity or positivity")
        spectra = _direct_periodic_spectra(h0, harmonic, reservoirs, beat, states,
                                           lowering, response_order)
        values = density*area*coupling[:, None]**2*spectra
        scales = np.maximum(np.max(abs(values), axis=1), np.finfo(float).tiny)
        residual = abs(values[:, 0]-values[:, 1]-values[:, 2])/scales
        if np.any(values[:, :2] < -1e-10*scales[:, None]) or np.any(residual > 1e-8):
            raise ValueError("connected QRT ordered positivity or Kubo identity failed")
        output = {}
        for arm, row, error in zip(("probe", "conjugate"), values, residual):
            output[arm] = dict(zip(
                ("alpha_abs_per_m", "alpha_em_per_m", "alpha_net_per_m"), map(float, row)))
            output[arm]["kubo_relative_residual"] = float(error)
        output.update({
            "reservoirs": reservoir_info,
            "atomic_diagnostics": {"maximum_state_trace_error": trace_error,
                "maximum_state_hermiticity_error": hermiticity,
                "minimum_phase_sampled_state_eigenvalue": minimum,
                "analytic_trace_zero_contraction_bound_s_inverse": reservoir_info["transit_reset_s_inverse"],
                "source_covariance_harmonic_order": 2*mean_order,
                "passed": True},
            "mean_order": mean_order, "response_order": response_order,
            "method": method,
            "mean_excited_population": float(np.trace(states[mean_order, 2:, 2:]).real),
            "scope": "same conditional reduced Lindblad; connected ordered QRT; finite harmonic convergence required",
        })
        return output
    atom = periodic_atomic_noise(h0, harmonic, reservoirs, beat,
        mean_order=mean_order, phase_samples=phase_samples)
    ports = PeriodicFieldPorts(("probe", "conjugate"), [1, -1], lowering, coupling)
    local = eliminate_periodic_atom(atom, ports,
        GeneratorFrequencyAxis([0.0], "reduced-minus-branch-Floquet-carrier"),
        response_order=response_order, linear_density_m_inverse=density*area)
    rates = ordered_absorption_emission(local)
    return {
        "probe": {key: float(value[0]) for key, value in rates[0].items()},
        "conjugate": {key: float(value[0]) for key, value in rates[1].items()},
        "reservoirs": reservoir_info,
        "atomic_diagnostics": atom.diagnostics,
        "field_audit": local.audit(),
        "mean_order": int(mean_order),
        "response_order": int(response_order),
        "method": method,
        "mean_excited_population": float(np.trace(atom.state_harmonics[mean_order, 2:, 2:]).real),
        "scope": (
            "local prescribed finite-seed four-level Lindblad QRT absorption/emission; "
            "no Doppler/spatial propagation or microscopic Zeeman reduction claim"),
    }


def selected_dressed_cell_absorption(*, D_GHz, T, P_pump, P_probe,
                                    w_pump, w_probe, delta_mhz, L,
                                    weak_line_strength, transit_rate=None,
                                    velocity_cutoff=4.0, relative_tolerance=0.005,
                                    branch=-1):
    """Cached Doppler-averaged positive absorption at one selected detuning.

    ``P_probe`` is the prescribed local input seed, not the pump-off measured
    output.  Uniform pump/input-seed means and zero coherent conjugate are
    retained at every point in the atomic medium.  This is a local homogeneous
    approximation, not integration along the generated/depleted mean fields.

    A normalized truncated Maxwell average is integrated adaptively, with
    dressed-state resonance roots supplied as breakpoints.  Separate reported
    quadrature bounds must meet ``relative_tolerance`` for both positive
    absorption coefficients.  Mean/response harmonic orders are compared at
    zero velocity and each dressed resonance; they are not silently lowered.
    The absorption-only intensity transmission is exp(-mean(alpha_abs)*L).
    Emission is independently retained: this scalar transmission corresponds
    to Jasperse's pure-loss coefficient only when the extra emission/gain
    channels can be neglected in that model.  It must not be applied again
    to a transfer/noise model already containing these atomic reservoirs.
    """
    args = (D_GHz, T, P_pump, P_probe, w_pump, w_probe, delta_mhz, L,
            weak_line_strength, constants.GAMMA_GG if transit_rate is None else transit_rate,
            velocity_cutoff, relative_tolerance, branch)
    if not np.isfinite(np.asarray(args, float)).all():
        raise ValueError("selected-point absorption inputs must be finite")
    values = tuple(float(item) for item in args[:-1])+(int(branch),)
    if branch != -1:
        raise NotImplementedError("ordered dressed absorption supports the certified minus branch")
    if (T <= 0 or P_pump < 0 or P_probe < 0 or w_pump <= 0 or w_probe <= 0
            or L < 0 or weak_line_strength <= 0 or values[9] <= 0
            or velocity_cutoff <= 0 or not 0 < relative_tolerance < 0.1):
        raise ValueError("invalid selected-point absorption powers, lengths, rate or tolerance")
    # Return a copy: downstream efficiency controls must not mutate cache facts.
    with core.blas_single_thread():
        return copy.deepcopy(_cached_selected_absorption(*values))


@functools.lru_cache(maxsize=32)
def _cached_selected_absorption(D_GHz, T, P_pump, P_probe, w_pump, w_probe,
                                delta_mhz, L, weak_line_strength, transit_rate,
                                velocity_cutoff, relative_tolerance, branch):
    from .schemes import fwm

    start = time.perf_counter()
    delta = delta_mhz*2*np.pi*1e6
    beat = -fwm.OMEGA_HF+delta
    if beat == 0.0:
        raise ValueError("nonzero minus-branch pump/seed beat required")
    density = float(hyperfine.number_density(T))
    pump_rabi = float(fwm.rabi_freq(P_pump, w_pump))
    seed_rabi = float(fwm.rabi_freq(P_probe, w_probe))
    probe_GHz = fwm.branch_center_GHz(D_GHz, branch)+delta_mhz*1e-3
    _, ks, kc = fwm.seeded_option_a_wavenumbers(D_GHz, np.array([probe_GHz]))
    optical = constants.C_LIGHT*np.array([float(ks[0]), float(kc[0])])
    area = float(observables.gaussian_mode_area(w_probe))
    center = 2*np.pi*D_GHz*1e9
    sigma = float(np.sqrt(constants.KB*T/constants.MASS_85RB))
    bound = velocity_cutoff*sigma
    excited = np.diag([0., 0., 1., 1.])
    base = fwm.pump_hamiltonian_at_Deff_zero(pump_rabi, pump_rabi)
    # Use the same inherited collinear velocity shift as the mean-field path.
    # Noncollinear atomic/Raman-Doppler corrections are not added implicitly.
    def eigenvalues(velocity):
        return np.linalg.eigvalsh(base-(center-fwm.K_VEC*velocity)*excited)

    grid = np.linspace(-bound, bound, 129)
    eigenvalue_grid = np.array([eigenvalues(v) for v in grid])
    breakpoints = []
    for left in range(4):
        for right in range(left+1, 4):
            equation = eigenvalue_grid[:, right]-eigenvalue_grid[:, left]-abs(beat)
            for index in np.flatnonzero(equation[:-1]*equation[1:] < 0):
                breakpoints.append(float(brentq(
                    lambda v: eigenvalues(v)[right]-eigenvalues(v)[left]-abs(beat),
                    grid[index], grid[index+1])))
    breakpoints = sorted(set(breakpoints))
    kwargs = dict(T=T, pump_rabi_rad_s=pump_rabi, probe_rabi_rad_s=seed_rabi,
        two_photon_detuning_rad_s=delta, optical_omega_rad_s=optical,
        uniform_area_m2=area, weak_line_strength=weak_line_strength,
        number_density_m3=density, transit_rate=transit_rate,
        phase_samples=32, method="connected-qrt")
    cache, maximum_kubo = {}, 0.0
    def evaluate(velocity, mean_order=3, response_order=3):
        nonlocal maximum_kubo
        key = (float(velocity), mean_order, response_order)
        if key not in cache:
            output = local_dressed_absorption(
                one_photon_detuning_rad_s=center-fwm.K_VEC*key[0],
                mean_order=mean_order, response_order=response_order, **kwargs)
            maximum_kubo = max(maximum_kubo, *(output[arm]["kubo_relative_residual"]
                for arm in ("probe", "conjugate")))
            cache[key] = np.array([output[arm][name] for arm in ("probe", "conjugate")
                for name in ("alpha_abs_per_m", "alpha_em_per_m", "alpha_net_per_m")])
        return cache[key]

    controls = sorted(set([0.]+breakpoints))
    floquet_checks = []
    for velocity in controls:
        coarse, fine = evaluate(velocity, 3, 3), evaluate(velocity, 4, 4)
        relative = np.max(abs(coarse-fine)/np.maximum(abs(fine), 1e-12))
        floquet_checks.append({"velocity_m_s": velocity,
                               "maximum_relative_rate_difference": float(relative)})
    maximum_floquet = max(item["maximum_relative_rate_difference"] for item in floquet_checks)
    if maximum_floquet > relative_tolerance/10:
        raise ValueError("finite-seed absorption harmonic convergence failed at a control velocity")

    normalization = sigma*np.sqrt(2*np.pi)*math.erf(velocity_cutoff/np.sqrt(2))
    reference_scale = np.maximum(abs(evaluate(0.)), 1.)
    def integrand(velocity):
        weight = np.exp(-.5*(velocity/sigma)**2)/normalization
        return evaluate(velocity)*weight/reference_scale

    tolerance = relative_tolerance/4
    for attempt in range(4):
        integral, error, information = quad_vec(integrand, -bound, bound,
            points=breakpoints, epsrel=tolerance, epsabs=1e-11, norm="max",
            full_output=True, quadrature="gk15", limit=256)
        means = integral*reference_scale
        error_bounds = float(error)*reference_scale
        absorption_relative_error = error_bounds[[0, 3]]/np.maximum(abs(means[[0, 3]]), 1e-30)
        if information.success and np.all(absorption_relative_error <= relative_tolerance):
            break
        tolerance /= 5
    else:
        raise ValueError("adaptive Doppler absorption did not reach the declared tolerance")
    if np.any(means[[0, 1, 3, 4]] < 0):
        raise ValueError("Doppler-averaged ordered absorption/emission must be nonnegative")
    result = {}
    for arm, offset in (("probe", 0), ("conjugate", 3)):
        result[arm] = {
            "alpha_abs_per_m": float(means[offset]),
            "alpha_em_per_m": float(means[offset+1]),
            "alpha_net_per_m": float(means[offset+2]),
            "in_cell_absorption_transmission": float(np.exp(-means[offset]*L)),
            "net_diagonal_transmission": float(np.exp(-means[offset+2]*L)),
            "absorption_quadrature_absolute_error_per_m": float(error_bounds[offset]),
            "absorption_quadrature_relative_error": float(absorption_relative_error[offset//3]),
        }
    result.update({
        "status": "CONDITIONAL_REDUCED_LINDBLAD_ABSORPTION",
        "method": "connected periodic QRT plus adaptive normalized truncated Maxwell quadrature",
        "velocity_cutoff_sigma": velocity_cutoff,
        "maxwell_tail_convergence_claimed": False,
        "velocity_classes_solved": len({key[0] for key in cache}),
        "local_qrts_solved": len(cache),
        "dressed_resonance_breakpoints_m_s": breakpoints,
        "quadrature_reported_scaled_error": float(error),
        "quadrature_evaluations": int(information.neval),
        "quadrature_success": bool(information.success),
        "relative_tolerance": relative_tolerance,
        "mean_order": 3, "response_order": 3,
        "floquet_control_orders": [[3, 3], [4, 4]],
        "floquet_controls": floquet_checks,
        "maximum_floquet_relative_rate_difference": maximum_floquet,
        "maximum_local_kubo_relative_residual": maximum_kubo,
        "runtime_seconds": time.perf_counter()-start,
        "cell_length_m": L,
        "delta_internal_MHz": delta_mhz,
        "source_planes": "prescribed vapor-entry pump/seed; zero coherent conjugate; atomic medium before windows",
        "scope": (
            "absorption-only loss coefficient for the declared pump-dressed four-level "
            "Lindblad; finite input seed and uniform pump prescribed locally; no mean-field "
            "path, angular Doppler, Zeeman, physical squeezing or empirical fit claim"),
        "emission_retained_separately": True,
        "already_in_atomic_transfer_do_not_reapply": True,
    })
    return result
