"""Longitudinal saturated PSR pump and microscopic field-noise propagation.

The physical pump is exp(i*phi)*q, with q in sqrt(mW/cm2). A common-phase
connection removes rapid scalar optical phase while retaining both Cartesian
Jones components. Local steady states and Langevin kernels are recomputed at
midpoints. No additional cell-loss vacuum is added to the atomic reservoirs.
"""
from collections import OrderedDict
from copy import deepcopy
from time import perf_counter

import numpy as np

from . import beam, constants, core
from .psr_conventions import jones_to_circular
from .psr_full_obe import solve_psr_pump
from .psr_homodyne import psr_output_homodyne_readout
from .psr_spectrum import PSR_NAMBU_MODES, PSR_RF_FRAME, psr_local_generator, psr_uniform_transfer
from .psr_thermal import solve_psr_thermal
from .quantum.contracts import AnalysisFrequencyAxis, GeneratorFrequencyAxis
from .quantum.traveling import NambuTransfer, compose_segments

MODEL_ID = "psr_rb87_d1_longitudinal_langevin_rf_v1"
_SIGNS = np.array([1., 1., -1., -1.])


class _RefineStep(Exception):
    pass


class _LocalFailure(Exception):
    pass


def _identity(rf):
    axis = GeneratorFrequencyAxis(2*np.pi*rf, PSR_RF_FRAME)
    shape = (len(rf), 4, 4)
    return NambuTransfer(axis, _SIGNS, PSR_NAMBU_MODES,
                         np.broadcast_to(np.eye(4), shape), np.zeros(shape), np.zeros(shape))


def _pump_flow(q, state):
    """Mean Maxwell derivative and prescribed common-phase connection."""
    intensity = float(np.vdot(q, q).real)
    if not np.isfinite(intensity) or intensity <= np.finfo(float).tiny:
        raise _LocalFailure("pump intensity is unresolved; a bright output carrier is required")
    field_scale = beam.field_amplitude(10*intensity)
    k = 2*np.pi*state["conditions"]["optical_frequency_hz"]/constants.C_LIGHT
    derivative = np.sqrt(intensity)*1j*k*state["polarization_xy_c_m2"]/(2*constants.EPS_0*field_scale)
    beta = float(np.vdot(q, derivative).imag/intensity)
    return derivative-1j*beta*q, beta


def _compact_state(state):
    return {key: deepcopy(state[key]) for key in (
        "model_id", "conditions", "chi_probe", "polarization_xy_c_m2", "alpha_per_m",
        "excited_fraction", "diagnostics")}


def _scaled_difference(before, after, rtol, atol):
    a, b = np.asarray(before), np.asarray(after)
    difference = np.abs(a-b)
    return dict(max_absolute=float(np.max(difference)),
                max_scaled=float(np.max(difference/(atol+rtol*np.maximum(np.abs(a), np.abs(b))))))


def _lab_transfer(path):
    phase = path["profile"]["carrier_phase_rad"][-1]
    gauge = np.diag(np.exp(1j*phase*_SIGNS))
    t = path["transfer"]
    return gauge@t.transfer, gauge@t.noise_greater@gauge.conj().T, gauge@t.noise_lesser@gauge.conj().T


def _path_difference(before, after, rtol, input_intensity):
    a, b = before["profile"], after["profile"]
    qa = np.sqrt(a["intensity_mw_cm2"][-1])*a["jones_xy"][-1]*np.exp(1j*a["carrier_phase_rad"][-1])
    qb = np.sqrt(b["intensity_mw_cm2"][-1])*b["jones_xy"][-1]*np.exp(1j*b["carrier_phase_rad"][-1])
    checks = {"output_field": _scaled_difference(qa, qb, rtol, 1e-8*np.sqrt(input_intensity)),
              "output_intensity": _scaled_difference(a["intensity_mw_cm2"][-1], b["intensity_mw_cm2"][-1],
                                                       rtol, 1e-8*input_intensity),
              "carrier_phase": _scaled_difference(a["carrier_phase_rad"][-1], b["carrier_phase_rad"][-1],
                                                    0., max(.001, rtol))}
    pa, pb = (np.outer(p["jones_xy"][-1], p["jones_xy"][-1].conj()) for p in (a, b))
    checks["output_polarization_projector"] = _scaled_difference(pa, pb, rtol, 1e-6)
    for name, old, new in zip(("transfer", "noise_greater", "noise_lesser"), _lab_transfer(before), _lab_transfer(after)):
        checks[name] = _scaled_difference(old, new, rtol, 1e-6)
    outputs = [psr_output_homodyne_readout(p["transfer"], p["profile"]["jones_xy"][-1]) for p in (before, after)]
    # Compare covariance in one common physical analyzer chart. Each grid's
    # own output polarization is separately checked above and in its extrema.
    common_jones = b["jones_xy"][-1]*np.exp(1j*b["carrier_phase_rad"][-1])
    common_covariances = []
    for path in (before, after):
        base = path["transfer"]
        lab = NambuTransfer(base.frequency_axis, base.signs, base.mode_labels, *_lab_transfer(path))
        common_covariances.append(psr_output_homodyne_readout(lab, common_jones)["covariance"])
    checks["covariance"] = _scaled_difference(*common_covariances, rtol, 1e-5)
    for name in ("minimum_variance", "maximum_variance"):
        checks[name] = _scaled_difference(outputs[0][name], outputs[1][name], rtol, 1e-5)
    for name in ("minimum_db", "maximum_db"):
        checks[name] = _scaled_difference(outputs[0][name], outputs[1][name], 0., max(.001, rtol))
    return dict(passed=all(c["max_scaled"] <= 1 for c in checks.values()), checks=checks)


def solve_psr_propagation(detuning_mhz, rf_hz, length_m, *, intensity_mw_cm2=.1,
                          jones_xy=(1., 0.), velocity_model="single", velocity_m_s=0.,
                          temp_c=70., ground_relaxation_hz=1000., density_m3=None,
                          segments=2, max_segments=32, rtol=.005,
                          max_atomic_evaluations=20000, thermal_options=None, freeze_pump=False):
    """Propagate a saturated pump and its four-port noise, refining cell slices.

    The mean field uses explicit midpoint integration. Every accepted slice's
    midpoint supplies its own atom, response and both ordered noise kernels.
    ``velocity_model`` is ``single`` or ``thermal``. Thermal integration uses
    S3 kernel/cutoff convergence independently at every local evaluation.
    ``freeze_pump`` is an explicit prescribed-uniform reference, not depletion.
    Refinement failures return passed=False and never claim an accurate PSD.
    """
    started = perf_counter()
    rf = AnalysisFrequencyAxis.from_hz(rf_hz).omega_rad_s/(2*np.pi)
    if not np.allclose(rf, -rf[::-1], rtol=1e-12, atol=1e-9) or not np.any(rf > 0):
        raise ValueError("a symmetric signed RF axis with a positive frequency is required")
    intensity, length, tolerance = map(float, (intensity_mw_cm2, length_m, rtol))
    if not np.isfinite(intensity) or intensity <= 0:
        raise ValueError("positive finite incident pump intensity required")
    if not np.isfinite(length) or length < 0 or not np.isfinite(tolerance) or not 0 < tolerance < 1:
        raise ValueError("nonnegative finite length and relative tolerance in (0,1) required")
    for value, label in ((segments, "segments"), (max_segments, "max_segments"),
                         (max_atomic_evaluations, "max_atomic_evaluations")):
        if isinstance(value, bool) or not np.isfinite(value) or value != int(value) or value < 1:
            raise ValueError(f"{label} must be a positive integer")
    if max_segments < segments:
        raise ValueError("max_segments must be at least segments")
    if velocity_model not in ("single", "thermal"):
        raise ValueError("velocity_model must be single or thermal")
    jones = np.asarray(jones_xy, complex)
    jones_to_circular(jones)
    q0 = np.sqrt(intensity)*jones
    inputs = dict(temp_c=temp_c, ground_relaxation_hz=ground_relaxation_hz, density_m3=density_m3)
    # Validate and determine physical density without an unnecessary thermal solve.
    initial = solve_psr_pump(detuning_mhz, intensity_mw_cm2=intensity, jones_xy=jones,
                             velocity_m_s=velocity_m_s if velocity_model == "single" else 0., **inputs)
    conditions = deepcopy(initial["conditions"])
    conditions.update(velocity_model=velocity_model, incident_intensity_mw_cm2=intensity,
                      cell_length_m=length, pump_evolution="prescribed uniform" if freeze_pump else "longitudinal midpoint",
                      common_phase_gauge="beta=Im(q_dagger F)/norm(q)^2; physical pump=exp(i phi)*q")
    if velocity_model == "thermal":
        conditions.pop("velocity_m_s", None)
        conditions.pop("effective_detuning_mhz", None)
        conditions.update(doppler=True, velocity_distribution="independent one-dimensional Maxwell classes")
    else:
        conditions["doppler"] = False
    thermal_settings = dict(thermal_options or {})
    forbidden = set(thermal_settings) - {"velocity_order", "max_order", "rtol", "atol_per_m", "cutoff_sigma",
                                        "tail_sigma", "max_classes", "method"}
    if forbidden:
        raise ValueError(f"unsupported thermal_options: {sorted(forbidden)}")
    evaluated = 1
    local_calls = 0
    cache = OrderedDict()
    health = []

    def evaluate(q):
        nonlocal evaluated, local_calls
        key = tuple(np.r_[q.real, q.imag])
        if key in cache:
            cache.move_to_end(key)
            return cache[key]
        remaining = int(max_atomic_evaluations)-evaluated
        if remaining <= 0 and not (velocity_model == "single" and np.array_equal(q, q0)):
            raise _LocalFailure("atomic evaluation budget reached")
        local_intensity = float(np.vdot(q, q).real)
        if not np.isfinite(local_intensity) or local_intensity <= np.finfo(float).tiny:
            raise _LocalFailure("pump intensity is unresolved; a bright output carrier is required")
        local_jones = q/np.sqrt(local_intensity)
        if velocity_model == "thermal":
            opts = dict(thermal_settings)
            if remaining <= 1:
                raise _LocalFailure("atomic evaluation budget reached before thermal validation")
            opts["max_classes"] = min(remaining-1, int(opts.get("max_classes", 3000)))
            result = solve_psr_thermal(detuning_mhz, rf, intensity_mw_cm2=local_intensity,
                                       jones_xy=local_jones, **inputs, **opts)
            evaluated += 1+result["performance"]["class_evaluations"]
            if not result["convergence"]["converged"]:
                raise _LocalFailure("local Maxwell integral is not converged: "+result["convergence"]["status"])
            state, local = result["pump_state"], result["local"]
            audit = deepcopy(result["convergence"])
            references = [(_compact_state(entry["pump_state"]), entry["local"])
                          for entry in result["history"][-3:-1]]
            if len(references) != 2:
                raise _LocalFailure("thermal response lacks order and cutoff reference integrals")
        else:
            state = (initial if np.array_equal(q, q0) else
                     solve_psr_pump(detuning_mhz, intensity_mw_cm2=local_intensity,
                                    jones_xy=local_jones, velocity_m_s=velocity_m_s, **inputs))
            evaluated += int(state is not initial)
            local = psr_local_generator(state, rf, method="optimized")
            audit, references = None, []
        local_calls += 1
        compact = _compact_state(state)
        health.append(deepcopy(compact["diagnostics"]))
        value = (compact, local, audit, references)
        cache[key] = value
        while len(cache) > 8:
            cache.popitem(last=False)
        return value

    def profile(points):
        q = np.asarray([p[1] for p in points])
        intensities = np.einsum("ij,ij->i", q.conj(), q).real
        return dict(z_m=np.array([p[0] for p in points]), intensity_mw_cm2=intensities,
                    jones_xy=q/np.sqrt(intensities)[:, None],
                    carrier_phase_rad=np.array([p[2] for p in points]))

    def one_grid(count):
        q, phase, total = q0.copy(), 0., _identity(rf)
        shadows = [_identity(rf), _identity(rf)] if velocity_model == "thermal" else []
        points, diagnostics = [(0., q.copy(), phase)], []
        dz = length/count
        for index in range(count):
            state, local, velocity_audit, references = evaluate(q)
            slope, beta = _pump_flow(q, state)
            if freeze_pump:
                midpoint = q
            else:
                if dz*np.linalg.norm(slope) > .25*np.linalg.norm(q):
                    raise _RefineStep("pump changes too far within a slice")
                midpoint = q+.5*dz*slope
                state, local, velocity_audit, references = evaluate(midpoint)
                slope, beta = _pump_flow(midpoint, state)
                if dz*np.linalg.norm(slope) > .5*np.linalg.norm(q):
                    raise _RefineStep("midpoint pump changes too far within a slice")
            segment = psr_uniform_transfer(local, state, dz, carrier_phase_rate_rad_per_m=beta)
            total = compose_segments(total, segment)
            for source, (shadow_state, shadow_local) in enumerate(references):
                shadow_segment = psr_uniform_transfer(shadow_local, shadow_state, dz,
                                                      carrier_phase_rate_rad_per_m=beta)
                shadows[source] = compose_segments(shadows[source], shadow_segment)
            # Unit-commutator accuracy can be lost cumulatively, even when
            # each local segment passes its own propagation precision guard.
            gain = float(np.max(np.linalg.norm(total.transfer, ord=2, axis=(-2, -1)))**2)
            added = max(float(np.max(np.abs(total.noise_greater))), float(np.max(np.abs(total.noise_lesser))))
            if np.finfo(float).eps*max(gain, added) >= 1:
                raise _LocalFailure("accumulated field gain or noise obscures one shot-noise unit")
            if not freeze_pump:
                updated = q+dz*slope
                if not np.isfinite(updated).all() or np.vdot(updated, updated).real <= np.finfo(float).tiny:
                    raise _RefineStep("midpoint pump update is unresolved")
                q = updated
            phase += dz*beta
            points.append(((index+1)*dz, q.copy(), phase))
            diagnostics.append(dict(z_mid_m=(index+.5)*dz, intensity_mw_cm2=float(np.vdot(midpoint, midpoint).real),
                                    jones_xy=midpoint/np.linalg.norm(midpoint), alpha_per_m=state["alpha_per_m"],
                                    carrier_phase_rate_rad_per_m=beta, excited_fraction=state["excited_fraction"],
                                    atomic_diagnostics=deepcopy(state["diagnostics"]), velocity_convergence=velocity_audit))
        result = dict(transfer=total, profile=profile(points), segment_diagnostics=diagnostics, segments=count)
        velocity_checks = []
        # Shadow integrals follow the accepted fine mean pump path, so their
        # output differences isolate accumulated velocity-kernel uncertainty.
        if shadows:
            variants = [dict(result, transfer=t) for t in (*shadows, total)]
            for label, first, second in zip(("order", "cutoff"), variants[:-1], variants[1:]):
                velocity_checks.append(dict(comparison=label,
                                             **_path_difference(first, second, tolerance, intensity)))
        result["velocity_convergence"] = dict(
            passed=all(c["passed"] for c in velocity_checks), checks=velocity_checks,
            scope="accumulated kernel references on the same candidate mean-pump path")
        return result

    paths, checks, final = [], [], None
    status, passed = "maximum longitudinal segment count reached", False
    with core.blas_single_thread(), np.errstate(over="raise", invalid="raise"):
        if length == 0:
            final = dict(transfer=_identity(rf), profile=profile([(0., q0, 0.)]), segment_diagnostics=[], segments=0)
            passed, status = True, "exact zero-length identity"
        else:
            counts = [1] if freeze_pump else [int(segments)]
            while not freeze_pump and counts[-1] < max_segments:
                counts.append(min(2*counts[-1], int(max_segments)))
            for count in counts:
                try:
                    candidate = one_grid(count)
                except _RefineStep as error:
                    checks.append(dict(segments=count, passed=False, status=str(error)))
                    continue
                except (_LocalFailure, ValueError, np.linalg.LinAlgError, FloatingPointError) as error:
                    status = str(error)
                    break
                if not candidate["velocity_convergence"]["passed"]:
                    final, status = candidate, "accumulated velocity-kernel output did not converge; increase thermal accuracy"
                    break
                if freeze_pump:
                    final, passed, status = candidate, True, "exact prescribed uniform-pump reference"
                    break
                paths.append(count)
                if final is not None:
                    try:
                        comparison = _path_difference(final, candidate, tolerance, intensity)
                    except (ValueError, np.linalg.LinAlgError, FloatingPointError) as error:
                        final, status = candidate, str(error)
                        break
                    checks.append(dict(segments=count, **comparison))
                    if comparison["passed"]:
                        final, passed, status = candidate, True, "pump, phase, transfer and output noise converged"
                        break
                final = candidate
    if final is None:
        final = dict(transfer=None, profile=profile([(0., q0, 0.)]), segment_diagnostics=[], segments=0)
    final.setdefault("velocity_convergence", dict(passed=bool(length == 0 or velocity_model == "single"),
        checks=[], scope="not needed at zero length or for a single velocity" if length == 0 or velocity_model == "single"
        else "not assessed without a complete propagation grid"))
    final.update(model_id=MODEL_ID, conditions=conditions,
                 output_intensity_mw_cm2=float(final["profile"]["intensity_mw_cm2"][-1]),
                 output_jones_xy=final["profile"]["jones_xy"][-1].copy(),
                 carrier_phase_advance_rad=float(final["profile"]["carrier_phase_rad"][-1]),
                 convergence=dict(passed=passed, status=status, rtol=tolerance, history=checks,
                                  completed_grids=paths, max_segments=int(max_segments),
                                  estimate_kind="successive midpoint grids; empirical convergence, not a rigorous error bound"),
                 performance=dict(elapsed_seconds=perf_counter()-started, local_evaluations=local_calls,
                                  atomic_evaluations=evaluated, max_atomic_evaluations=int(max_atomic_evaluations),
                                  retained_atomic_states=0, cached_local_entries=len(cache)),
                 pump_diagnostics=dict(local_state_checks=health))
    return final
