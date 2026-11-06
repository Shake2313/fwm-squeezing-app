"""Independent Maxwell velocity classes for the saturated PSR RF model.

Average drift and both ordered Langevin kernels before field propagation.
Only integrated 4-port arrays survive each class, keeping atomic memory bounded.
Refinement differences are empirical convergence estimates, not rigorous bounds.
"""
from copy import deepcopy
from math import erfc, sqrt
from time import perf_counter

import numpy as np

from . import constants, magneto_multilevel, species
from .psr_conventions import _unit_jones
from .psr_full_obe import solve_psr_pump
from .psr_homodyne import psr_homodyne_readout
from .psr_spectrum import psr_local_generator, psr_uniform_transfer
from .quantum.contracts import AnalysisFrequencyAxis
from .quantum.traveling import LocalNambuGenerator

PUMP_MODEL_ID = "psr_rb87_d1_maxwell_full_obe_v1"
MODEL_ID = "psr_rb87_d1_maxwell_langevin_rf_v1"


def _comparison(old, new, rtol, atol):
    a, b = np.asarray(old), np.asarray(new)
    delta = np.abs(b-a)
    return dict(max_absolute=float(np.max(delta)),
                max_scaled=float(np.max(delta/(atol+rtol*np.maximum(np.abs(a), np.abs(b))))))


def _kernel_difference(old, new, rtol, atol):
    fields = {"drift": "drift", "greater": "noise_greater_by_reservoir",
              "lesser": "noise_lesser_by_reservoir"}
    checks = {name: _comparison(getattr(old["local"], attr), getattr(new["local"], attr), rtol, atol)
              for name, attr in fields.items()}
    # Express mean susceptibility error in the same inverse-metre units.
    k = 2*np.pi*new["pump_state"]["conditions"]["optical_frequency_hz"]/constants.C_LIGHT
    checks["mean_polarization"] = _comparison(
        k*old["pump_state"]["chi_probe"]/2, k*new["pump_state"]["chi_probe"]/2, rtol, atol)
    conditions = new["pump_state"]["conditions"]
    # Preserve the established S3 x-linear acceptance rule exactly. For a
    # general Jones state, the projected susceptibility alone cannot detect
    # an error in the orthogonal polarization which drives self-rotation.
    peak_field = conditions["peak_field_v_m"]
    if peak_field > 0 and not np.array_equal(conditions["jones_xy"], [1., 0.]):
        factor = k/(2*constants.EPS_0*peak_field)
        checks["polarization_vector"] = _comparison(
            factor*old["pump_state"]["polarization_xy_c_m2"],
            factor*new["pump_state"]["polarization_xy_c_m2"], rtol, atol)
    return dict(passed=all(c["max_scaled"] <= 1 for c in checks.values()), checks=checks)


def _integrate(detuning_mhz, rf, inputs, order, cutoff, extra_shifts, method):
    atom = magneto_multilevel.build_multilevel_atom("Rb87", "D1")
    centroid = magneto_multilevel.component_centroid_mhz(atom)
    shifts, weights, sigma = magneto_multilevel.velocity_quadrature(
        atom, centroid+detuning_mhz, inputs["temp_c"]+273.15,
        order=order, cutoff_sigma=cutoff, extra_shifts_mhz=extra_shifts)
    totals = None
    mean = None
    health = dict(trace_error=0., hermiticity_error=0., min_eigenvalue=1.,
                  steady_state_relative_residual=0.)
    start = perf_counter()
    for shift, weight in zip(shifts, weights):
        velocity = shift*1e6*constants.C_LIGHT/atom.frequency_hz
        pump = solve_psr_pump(detuning_mhz, velocity_m_s=velocity, **inputs)
        local = psr_local_generator(pump, rf, method=method)
        arrays = (local.drift, local.noise_greater_by_reservoir, local.noise_lesser_by_reservoir)
        if totals is None:
            totals = [np.zeros_like(a) for a in arrays]
            mean = {key: np.zeros_like(pump[key]) for key in
                    ("chi_probe", "alpha_per_m", "excited_fraction", "ground_populations", "excited_populations",
                     "polarization_xy_c_m2")}
            conditions = deepcopy(pump["conditions"])
        for total, value in zip(totals, arrays):
            total += weight*value
        for key in mean:
            mean[key] += weight*pump[key]
        for key in health:
            reducer = min if key == "min_eigenvalue" else max
            health[key] = reducer(health[key], pump["diagnostics"][key])
    averaged = LocalNambuGenerator(local.frequency_axis, local.signs, local.mode_labels,
                                    *totals, local.reservoir_names)
    if not averaged.audit()["passed"]:
        raise ValueError("Maxwell-averaged response failed the local commutator or diffusion audit")
    for key in ("velocity_m_s", "effective_detuning_mhz"):
        conditions.pop(key, None)
    conditions.update(velocity_distribution="independent one-dimensional Maxwell classes",
                      velocity_mean_m_s=0., velocity_sigma_m_s=float(sigma*1e6*constants.C_LIGHT/atom.frequency_hz),
                      velocity_mixing="none", doppler=True)
    mean.update(model_id=PUMP_MODEL_ID, conditions=conditions, diagnostics=health)
    mean["diagnostics"]["solver_status"] = "all sampled velocity steady states physical"
    for key in ("chi_probe", "alpha_per_m", "excited_fraction"):
        mean[key] = mean[key].item()
    return dict(local=averaged, pump_state=mean, interval_order=order, cutoff_sigma=cutoff,
                class_count=len(weights), weight_sum=float(weights.sum()),
                elapsed_seconds=perf_counter()-start,
                omitted_maxwell_probability=erfc(cutoff/sqrt(2.)))


def solve_psr_thermal(detuning_mhz, rf_hz, *, intensity_mw_cm2=.1, temp_c=70.,
                      ground_relaxation_hz=1000., density_m3=None, velocity_order=4,
                      max_order=10, rtol=.005, atol_per_m=1e-6, cutoff_sigma=6.,
                      tail_sigma=7., max_classes=3000, method="optimized", reference=False,
                      jones_xy=(1., 0.)):
    """Bounded thermal integration with kernel refinement and a cutoff check.

    ``reference=True`` evaluates exactly ``velocity_order`` once. It is a direct
    quadrature reference, with convergence explicitly unassessed. ``method``
    selects the independently verified atomic implementation, not a physical
    approximation. The adaptive path increases interval order by two until
    every complex kernel component and mean polarization meet the tolerances,
    then evaluates the extended cutoff at the accepted order. The common
    finite-cutoff distribution is normalized as in velocity_quadrature.
    ``jones_xy`` is the unit complex Cartesian carrier polarization, including
    its phase in the same frame as the returned field generator. The mean
    polarization vector is retained for longitudinal pump propagation.
    Call thermal_convergence at the actual cell length for the original
    x-linear uniform-pump readout. General Jones propagation belongs to S4.
    """
    started = perf_counter()
    rf = AnalysisFrequencyAxis.from_hz(rf_hz).omega_rad_s/(2*np.pi)
    if not len(rf) or not np.any(rf > 0) or not np.allclose(rf, -rf[::-1], rtol=1e-12, atol=1e-9):
        raise ValueError("a symmetric signed RF axis is required")
    if method not in ("direct", "optimized"):
        raise ValueError("method must be direct or optimized")
    for value, label, minimum in ((velocity_order, "velocity_order", 2),
                                  (max_order, "max_order", 2), (max_classes, "max_classes", 1)):
        if isinstance(value, bool) or not np.isfinite(value) or int(value) != value or value < minimum:
            raise ValueError(f"{label} must be an integer at least {minimum}")
    if max_order < velocity_order:
        raise ValueError("max_order must be at least velocity_order")
    if not all(np.isfinite(v) for v in (rtol, atol_per_m, cutoff_sigma, tail_sigma)):
        raise ValueError("finite convergence tolerances and cutoffs required")
    if rtol <= 0 or atol_per_m <= 0 or cutoff_sigma <= 0 or tail_sigma <= cutoff_sigma:
        raise ValueError("positive tolerances and tail_sigma > cutoff_sigma > 0 required")
    inputs = dict(intensity_mw_cm2=intensity_mw_cm2, temp_c=temp_c,
                  ground_relaxation_hz=ground_relaxation_hz, density_m3=density_m3, jones_xy=_unit_jones(jones_xy).copy())
    # Validate physical inputs before potentially many class solves.
    solve_psr_pump(detuning_mhz, **inputs)
    atom = magneto_multilevel.build_multilevel_atom("Rb87", "D1")
    centroid = magneto_multilevel.component_centroid_mhz(atom)
    # The extreme sidebands bracket the RF-resonant region. Every RF component
    # participates in refinement, including internal sideband resonances.
    extra = (-float(np.max(np.abs(rf)))/1e6, float(np.max(np.abs(rf)))/1e6)
    history, used = [], 0
    comparison, tail_comparison = None, None
    status = "unassessed direct quadrature" if reference else "maximum interval order reached"
    def evaluate(order, cutoff):
        nonlocal used
        shifts, _, _ = magneto_multilevel.velocity_quadrature(
            atom, centroid+detuning_mhz, float(temp_c)+273.15, order=order,
            cutoff_sigma=cutoff, extra_shifts_mhz=extra)
        if used+len(shifts) > max_classes:
            return None
        result = _integrate(detuning_mhz, rf, inputs, order, cutoff, extra, method)
        used += result["class_count"]
        history.append(result)
        return result
    orders = [int(velocity_order)] if reference else list(range(int(velocity_order), int(max_order)+1, 2))
    if not reference and orders[-1] != max_order:
        orders.append(int(max_order))
    for order in orders:
        entry = evaluate(order, float(cutoff_sigma))
        if entry is None:
            status = "velocity-class evaluation budget reached"
            break
        if len(history) > 1:
            comparison = _kernel_difference(history[-2], entry, rtol, atol_per_m)
            if comparison["passed"]:
                status = "kernel order converged; cutoff check pending"
                break
    if not history:
        raise ValueError("velocity-class budget is smaller than the initial quadrature")
    if not reference and comparison and comparison["passed"]:
        base = history[-1]
        extended = evaluate(base["interval_order"], float(tail_sigma))
        if extended is None:
            status = "velocity-class evaluation budget reached before cutoff check"
        else:
            tail_comparison = _kernel_difference(base, extended, rtol, atol_per_m)
            status = "kernel and cutoff converged" if tail_comparison["passed"] else "extended cutoff comparison failed"
    converged = bool(not reference and comparison and comparison["passed"]
                     and tail_comparison and tail_comparison["passed"])
    final = history[-1]
    arrays_bytes = sum(a.nbytes for entry in history for a in
                      (entry["local"].drift, entry["local"].noise_greater_by_reservoir,
                       entry["local"].noise_lesser_by_reservoir))
    return dict(model_id=MODEL_ID, local=final["local"], pump_state=final["pump_state"], history=history,
                convergence=dict(converged=converged, status=status, order=comparison, cutoff=tail_comparison,
                                 rtol=rtol, atol_per_m=atol_per_m, propagation_assessed=False,
                                 interval_order=final["interval_order"], cutoff_sigma=final["cutoff_sigma"],
                                 omitted_maxwell_probability=final["omitted_maxwell_probability"],
                                 requested_cutoff_sigma=cutoff_sigma, extended_cutoff_sigma=tail_sigma,
                                 estimate_kind="successive quadrature differences, not a rigorous error bound"),
                performance=dict(method=method, elapsed_seconds=perf_counter()-started,
                                 class_evaluations=used, max_class_evaluations=max_classes,
                                 retained_integral_bytes=arrays_bytes,
                                 atomic_storage="one velocity class; full atomic arrays are not retained"))


def thermal_convergence(result, length_m, *, detection_efficiency=1.):
    """Add actual-length covariance and min/max convergence to the kernel audit.

    Source (eta=1) is always checked as well as the requested detector plane.
    No atom is solved here. This lets cell length and detection reuse the
    integrated kernel history. A failed propagation check is explicit.
    """
    if not np.isfinite(length_m) or length_m < 0 or not np.isfinite(detection_efficiency) or not 0 <= detection_efficiency <= 1:
        raise ValueError("nonnegative length and efficiency in [0,1] required")
    audit = deepcopy(result["convergence"])
    history = result["history"]
    audit.update(propagation_assessed=True, length_m=float(length_m),
                 detection_efficiency=float(detection_efficiency), output_checks=[])
    if len(history) < 2:
        audit.update(converged=False, passed=False, status="propagation convergence requires at least two quadratures")
        return audit
    rtol = audit["rtol"]
    # Absolute floors have SQL-variance units and dB respectively.
    variance_atol, db_atol = 1e-5, max(.001, rtol)
    first_pair = max(0, len(history)-3) if audit["cutoff"] is not None else len(history)-2
    previous = None
    passed = True
    try:
        for entry in history[first_pair:]:
            transfer = psr_uniform_transfer(entry["local"], entry["pump_state"], length_m)
            planes = sorted(set((1., float(detection_efficiency))))
            outputs = [psr_homodyne_readout(transfer, detection_efficiency=eta) for eta in planes]
            if previous is not None:
                for eta, before, after in zip(planes, previous, outputs):
                    checks = {key: _comparison(before[key], after[key], rtol, variance_atol)
                              for key in ("covariance", "minimum_variance", "maximum_variance")}
                    for key in ("minimum_db", "maximum_db"):
                        # dB errors are absolute, without a large-value relative allowance.
                        checks[key] = _comparison(before[key], after[key], 0., db_atol)
                    item = dict(detection_efficiency=eta,
                                interval_order=entry["interval_order"], cutoff_sigma=entry["cutoff_sigma"],
                                checks=checks, passed=all(c["max_scaled"] <= 1 for c in checks.values()))
                    passed &= item["passed"]
                    audit["output_checks"].append(item)
            previous = outputs
    except (ValueError, np.linalg.LinAlgError, FloatingPointError) as error:
        passed = False
        audit["propagation_error"] = str(error)
    audit.update(converged=bool(audit["converged"] and passed), variance_atol=variance_atol, db_atol=db_atol)
    if not passed:
        audit["status"] = "propagated covariance or noise extrema did not converge"
    elif audit["converged"]:
        audit["status"] = "kernels, cutoff, covariance and noise extrema converged"
    audit["passed"] = audit["converged"]
    return audit
