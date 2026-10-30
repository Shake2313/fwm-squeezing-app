"""Classical polarization self-rotation from the existing atomic response.

The returned slope can parameterize an explicitly conditional Gaussian shear.
It is not a field-noise covariance or a microscopic squeezing prediction.
Optical excitation is eliminated in the weak-excitation approximation; ground
optical pumping is retained. The local pump is uniform along the cell. No
finite-RF atomic response, optical saturation, or Langevin diffusion is added.
"""
from __future__ import annotations

import math

import numpy as np

from . import constants, magneto_multilevel, species


DEFAULT_INTENSITY_MW_CM2 = 0.1
DEFAULT_TEMP_C = 70.0
DEFAULT_CELL_LENGTH_M = 0.075
DEFAULT_GROUND_RELAXATION_HZ = 1e3
WEAK_EXCITATION_LIMIT = 0.05
UNIFORM_PUMP_OD_LIMIT = 0.3
ELLIPTICITY_DRIFT_LIMIT = 0.1
DERIVATIVE_RTOL = 0.01
DERIVATIVE_ATOL_PER_M = 1e-6

MODEL_SCOPE = (
    "Weak-excitation multihyperfine atomic PSR; zero magnetic field, independent "
    "velocity classes and isotropic transit reset; uniform local pump. "
    "Full optical saturation, longitudinal pump evolution, finite-RF response "
    "and atomic Langevin diffusion are omitted. A Gaussian shear based on this "
    "classical DC slope is a separate conditional model."
)


def _stokes_rates(qwp_deg, chi_plus, chi_minus, wave_number_per_m):
    """Local azimuth/ellipticity rates from dE_q/dz = i k chi_q E_q / 2.

    E_+ = (E_x-i E_y)/sqrt(2), S3 = 2 Im(E_x* E_y),
    psi = atan2(S2,S1)/2 and epsilon = asin(S3/S0)/2.
    With the shared QWP convention epsilon is -qwp_deg near linear input.
    All rates have units 1/m; angles are radians.
    """
    polarization = magneto_multilevel._polarization(qwp_deg)
    ep, em = polarization[1], polarization[-1]
    dep = 0.5j * wave_number_per_m * chi_plus * ep
    dem = 0.5j * wave_number_per_m * chi_minus * em
    x, y = (ep + em) / math.sqrt(2), 1j * (ep - em) / math.sqrt(2)
    dx, dy = (dep + dem) / math.sqrt(2), 1j * (dep - dem) / math.sqrt(2)
    s0 = abs(x)**2 + abs(y)**2
    s1, s2, s3 = abs(x)**2 - abs(y)**2, 2 * (x.conjugate()*y).real, 2 * (x.conjugate()*y).imag
    ds0 = 2 * (x.conjugate()*dx + y.conjugate()*dy).real
    ds1 = 2 * (x.conjugate()*dx - y.conjugate()*dy).real
    cross_derivative = dx.conjugate()*y + x.conjugate()*dy
    ds2, ds3 = 2 * cross_derivative.real, 2 * cross_derivative.imag
    helicity = s3 / s0
    epsilon = 0.5 * math.asin(float(np.clip(helicity, -1., 1.)))
    azimuth_rate = 0.5 * (s1*ds2 - s2*ds1) / (s1*s1 + s2*s2)
    ellipticity_rate = 0.5 * (ds3*s0 - s3*ds0) / (s0*s0*math.sqrt(1-helicity*helicity))
    return epsilon, float(azimuth_rate), float(ellipticity_rate), float(-ds0/s0)


def compute_psr_response(
        detunings_mhz, *, isotope="Rb87", line="D1",
        intensity_mw_cm2=DEFAULT_INTENSITY_MW_CM2, temp_c=DEFAULT_TEMP_C,
        cell_length_m=DEFAULT_CELL_LENGTH_M,
        ground_relaxation_hz=DEFAULT_GROUND_RELAXATION_HZ,
        doppler=True, velocity_order=4, ellipticity_step_deg=0.05,
        density_m3=None):
    """Return classical local PSR slopes and explicit approximation diagnostics.

    The detuning zero is the upper ground-F optical line-strength centroid.
    Density is per isotope in m^-3; None uses the elemental vapor-pressure
    density as a pure-isotope cell. Intensity is mW/cm^2, temperature Celsius,
    length meters and the reset rate s^-1 (no 2*pi conversion).

    Five QWP settings (0, +/-h, +/-h/2) determine d(dpsi/dz)/d epsilon
    and the neglected ellipticity drift. Each detuning is evaluated separately
    because the shared solver otherwise aggregates atomic health across its
    entire scan. Estimates outside applicability screens remain unchanged.
    Screening thresholds are engineering choices, not experimental validation.
    """
    delta = np.atleast_1d(np.asarray(detunings_mhz, dtype=float))
    scalars = np.asarray([intensity_mw_cm2, temp_c, cell_length_m,
                          ground_relaxation_hz, ellipticity_step_deg], dtype=float)
    if (delta.ndim != 1 or not delta.size or not np.all(np.isfinite(delta))
            or not np.all(np.isfinite(scalars))):
        raise ValueError("finite one-dimensional detunings and finite settings required")
    if isotope not in ("Rb85", "Rb87") or line not in ("D1", "D2"):
        raise ValueError("PSR supports Rb85/Rb87 on D1/D2")
    if intensity_mw_cm2 < 0 or cell_length_m < 0 or temp_c <= -273.15 or ground_relaxation_hz <= 0:
        raise ValueError("nonnegative intensity/length, positive temperature/reset rate required")
    if not 0 < ellipticity_step_deg <= 5.:
        raise ValueError("ellipticity_step_deg must be in (0, 5] degrees")
    if isinstance(velocity_order, bool) or int(velocity_order) != velocity_order or not 2 <= velocity_order <= 20:
        raise ValueError("velocity_order must be an integer from 2 to 20")
    if not isinstance(doppler, (bool, np.bool_)):
        raise ValueError("doppler must be boolean")
    density = (species.number_density(species.ISOTOPES[isotope], temp_c+273.15)
               if density_m3 is None else float(density_m3))
    if not np.isfinite(density) or density < 0:
        raise ValueError("density_m3 must be finite and nonnegative")

    atom = magneto_multilevel.build_multilevel_atom(isotope, line)
    wave_number = 2*np.pi*atom.frequency_hz/constants.C_LIGHT
    h = float(ellipticity_step_deg)
    angles = np.array([0., h, -h, h/2, -h/2])
    epsilons = np.array([_stokes_rates(angle, 0j, 0j, wave_number)[0] for angle in angles])
    shape = (delta.size, len(angles))
    chi_plus, chi_minus, chi_probe = [np.full(shape, complex(np.nan, np.nan)) for _ in range(3)]
    azimuth_rates, ellipticity_rates, absorption_rates = [np.full(shape, np.nan) for _ in range(3)]
    health_keys = ("virtual_excitation_bound", "excited_fraction_mean", "excited_fraction_max",
                   "trace_error", "hermiticity_error", "steady_state_relative_residual")
    health = {key: np.full(delta.size, np.nan) for key in health_keys}
    health["min_ground_eigenvalue"] = np.full(delta.size, np.nan)
    atomic_valid = np.zeros(delta.size, bool)
    numerical_reasons = [[] for _ in delta]
    class_counts = np.zeros(delta.size, int)
    settings = dict(
        isotope=isotope, line=line, intensity_mw_cm2=float(intensity_mw_cm2), temp_c=float(temp_c),
        b_ut=0., b_x_ut=0., b_y_ut=0., magnetic_derivative=False,
        cell_length_m=float(cell_length_m), gamma_out_hz=0., gamma_in_hz=0.,
        ground_relaxation_hz=float(ground_relaxation_hz), wall_relaxation_hz=1.,
        velocity_mixing="none", velocity_points=int(velocity_order), doppler=bool(doppler),
        density_m3=density,
    )
    for i, detuning in enumerate(delta):
        point_health = []
        try:
            for j, angle in enumerate(angles):
                response = magneto_multilevel.multilevel_response([detuning], qwp_deg=float(angle), **settings)
                chi_plus[i, j], chi_minus[i, j], chi_probe[i, j] = (
                    complex(response[key][0, 0]) for key in ("chi_plus", "chi_minus", "chi_probe"))
                _, azimuth_rates[i, j], ellipticity_rates[i, j], absorption_rates[i, j] = _stokes_rates(
                    angle, chi_plus[i, j], chi_minus[i, j], wave_number)
                point_health.append(response["diagnostics"])
                class_counts[i] = max(class_counts[i], max(response["velocity_quadrature"]["class_counts"]))
        except (np.linalg.LinAlgError, FloatingPointError) as error:
            numerical_reasons[i].append(f"atomic solve failed: {type(error).__name__}")
            continue
        for key in health_keys:
            health[key][i] = max(item[key] for item in point_health)
        health["min_ground_eigenvalue"][i] = min(item["min_ground_eigenvalue"] for item in point_health)
        atomic_valid[i] = all(item["solver_status"] == "physical steady state" for item in point_health)
        if not atomic_valid[i]:
            numerical_reasons[i].append("atomic steady-state consistency failed")

    def derivative(values):
        coarse = (values[:, 1]-values[:, 2])/(epsilons[1]-epsilons[2])
        fine = (values[:, 3]-values[:, 4])/(epsilons[3]-epsilons[4])
        error = np.abs(fine-coarse)
        converged = np.isfinite(error) & (error <= DERIVATIVE_ATOL_PER_M + DERIVATIVE_RTOL*np.abs(fine))
        return fine, error, converged

    g, derivative_error, g_converged = derivative(azimuth_rates)
    drift, drift_error, drift_converged = derivative(ellipticity_rates)
    alpha = absorption_rates[:, 0]
    optical_depth = alpha*cell_length_m
    with np.errstate(over="ignore", invalid="ignore"):
        transmission = np.exp(-optical_depth)
    shear, drift_length = g*cell_length_m, drift*cell_length_m
    finite = (np.isfinite(chi_plus).all(axis=1) & np.isfinite(chi_minus).all(axis=1)
              & np.isfinite(chi_probe).all(axis=1) & np.isfinite(azimuth_rates).all(axis=1)
              & np.isfinite(ellipticity_rates).all(axis=1) & np.isfinite(absorption_rates).all(axis=1)
              & np.isfinite(shear) & np.isfinite(transmission)
              & np.all(np.isfinite(np.stack(list(health.values()))), axis=0))
    numerical_valid = atomic_valid & finite
    derivative_converged = g_converged & drift_converged
    weak_excitation = health["virtual_excitation_bound"] <= WEAK_EXCITATION_LIMIT
    uniform_pump = np.abs(optical_depth) <= UNIFORM_PUMP_OD_LIMIT
    small_drift = np.abs(drift_length) <= ELLIPTICITY_DRIFT_LIMIT
    passive_noise_valid = numerical_valid & (alpha >= -1e-10*np.maximum(1., np.abs(alpha)))
    approximation_valid = numerical_valid & derivative_converged & weak_excitation & uniform_pump & small_drift
    approximation_reasons, passive_reasons, statuses = [], [], []
    for i in range(delta.size):
        if not finite[i]:
            numerical_reasons[i].append("nonfinite atomic response or propagation estimate")
        reasons = []
        if numerical_valid[i] and not weak_excitation[i]:
            reasons.append("weak-excitation bound exceeds 0.05")
        if numerical_valid[i] and not uniform_pump[i]:
            reasons.append("uniform-pump screen |optical depth| <= 0.3 exceeded")
        if numerical_valid[i] and not small_drift[i]:
            reasons.append("neglected ellipticity drift screen |rate L| <= 0.1 exceeded")
        if numerical_valid[i] and not derivative_converged[i]:
            reasons.append("ellipticity half-step derivative check failed")
        approximation_reasons.append(reasons)
        passive_reasons.append([] if passive_noise_valid[i] else [
            "numerical response unavailable" if not numerical_valid[i]
            else "negative absorption unsupported by passive-loss noise model"])
        statuses.append("numerical failure" if not numerical_valid[i] else
                        "outside declared approximation screens" if reasons else
                        "within declared approximation screens")

    return dict(
        detunings_mhz=delta, shear=shear, optical_depth=optical_depth, transmission=transmission,
        g_per_m=g, alpha_per_m=alpha, derivative_error=derivative_error,
        derivative_relative_error=derivative_error/np.maximum(np.abs(g), DERIVATIVE_ATOL_PER_M),
        derivative_converged=derivative_converged,
        ellipticity_drift=drift, ellipticity_drift_per_m=drift,
        ellipticity_drift_length=drift_length, ellipticity_derivative_error=drift_error,
        linear_azimuth_rate_per_m=azimuth_rates[:, 0],
        numerical_valid=numerical_valid, approximation_valid=approximation_valid,
        passive_noise_valid=passive_noise_valid,
        conditional_model_applicable=approximation_valid & passive_noise_valid,
        numerical_reasons=numerical_reasons, approximation_reasons=approximation_reasons,
        passive_noise_reasons=passive_reasons, approximation_status=statuses,
        atomic_health=health,
        stencil=dict(qwp_deg=angles, ellipticity_rad=epsilons, chi_plus=chi_plus,
                     chi_minus=chi_minus, chi_probe=chi_probe,
                     azimuth_rate_per_m=azimuth_rates, ellipticity_rate_per_m=ellipticity_rates),
        conditions=dict(isotope=isotope, line=line, intensity_mw_cm2=float(intensity_mw_cm2),
                        temp_c=float(temp_c), cell_length_m=float(cell_length_m), density_m3=density,
                        density_source="pure-isotope vapor-pressure density" if density_m3 is None else "specified isotope density",
                        ground_relaxation_hz=float(ground_relaxation_hz), doppler=bool(doppler),
                        velocity_order=int(velocity_order), ellipticity_step_deg=h,
                        b_ut=0., gamma_out_hz=0., gamma_in_hz=0.),
        frequency_hz=atom.frequency_hz, reference_centroid_mhz=magneto_multilevel.component_centroid_mhz(atom),
        detuning_reference="line-strength centroid of upper ground-F optical component",
        level_counts=dict(ground=len(atom.ground), excited=len(atom.excited)),
        velocity_quadrature=dict(method="resonance-centered composite Gauss-Legendre Maxwell" if doppler else "stationary atoms",
                                 interval_order=int(velocity_order), class_counts=class_counts,
                                 convergence_status="not assessed; compare interval orders"),
        approximation_thresholds=dict(weak_excitation_bound=WEAK_EXCITATION_LIMIT,
            uniform_pump_optical_depth=UNIFORM_PUMP_OD_LIMIT, ellipticity_drift_length=ELLIPTICITY_DRIFT_LIMIT,
            derivative_rtol=DERIVATIVE_RTOL, derivative_atol_per_m=DERIVATIVE_ATOL_PER_M,
            interpretation="conservative engineering screens, not experimental validation"),
        derivative_method="symmetric ellipticity differences; half-step result with full-step comparison",
        phase_convention="E+ = (Ex-i Ey)/sqrt(2); S3 = 2 Im(Ex* Ey); dE_q/dz = +i k chi_q E_q/2",
        rate_convention="ground_relaxation_hz is a reset rate in s^-1; no 2pi conversion",
        model_scope=MODEL_SCOPE,
    )
