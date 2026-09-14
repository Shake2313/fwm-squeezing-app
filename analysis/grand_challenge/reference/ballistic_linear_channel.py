"""Independent passive, finite-velocity bosonic transport control.

Two canonical forward flux modes (b, a), with an independent canonical atomic
inflow and a local atomic reservoir. This is an oscillator testbed, not a finite
Rb atom, FWM squeezing calculation, or general moving-vapor closure. No imports
from gabes.quantum.transport; no fitted, clipped, or commutator-completed noise.
"""

from dataclasses import asdict, dataclass
import json

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm


@dataclass(frozen=True)
class Parameters:
    """Declared units: length ell, time t; c,v in ell/t, omega,gamma in 1/t."""

    omega: float = 0.8
    c: float = 4.0
    v: float = 0.65
    gamma: float = 0.35

    def __post_init__(self):
        for name, value in asdict(self).items():
            if not np.isrealobj(value) or not np.isfinite(value):
                raise ValueError(f"{name} must be finite and real")
        if self.c <= 0 or self.v <= 0 or self.gamma < 0:
            raise ValueError("require c > 0, v > 0, gamma >= 0")


@dataclass(frozen=True)
class Segment:
    """One constant segment; kappa has units 1/ell."""

    length: float
    kappa: complex

    def __post_init__(self):
        if not np.isrealobj(self.length) or not np.isfinite(self.length) or self.length < 0:
            raise ValueError("length must be finite, real and nonnegative")
        if not np.isfinite(self.kappa):
            raise ValueError("kappa must be finite")


@dataclass(frozen=True)
class PassiveChannel:
    """Transfer of BOTH input ports and separately integrated reservoir commutator."""

    transfer: np.ndarray
    reservoir: np.ndarray


def drift(kappa, parameters):
    """K + K^dagger = -Q exactly: conjugate coupling is purely Hamiltonian."""
    p = parameters
    return np.array([[1j*p.omega/p.c, -1j*kappa],
                     [-1j*np.conj(kappa), (1j*p.omega-p.gamma)/p.v]], complex)


def _reservoir_density(parameters):
    return np.diag([0.0, 2*parameters.gamma/parameters.v]).astype(complex)


def segment_channel(segment, parameters=Parameters()):
    """Exact constant-segment exponential and independently evaluated Gramian.

    Column vectorization gives vec(KW+WK^dagger) =
    (I tensor K + K.conj() tensor I) vec(W). Augmenting by a constant 1
    integrates W'=KW+WK^dagger+Q without inverting a singular Lyapunov map.
    Unlike the growing block of a Van Loan construction, this block is passive.
    """
    k = drift(segment.kappa, parameters)
    generator = np.zeros((5, 5), complex)
    generator[:4, :4] = np.kron(np.eye(2), k) + np.kron(k.conj(), np.eye(2))
    generator[:4, 4] = _reservoir_density(parameters).reshape(4, order="F")
    reservoir = expm(segment.length*generator)[:4, 4].reshape((2, 2), order="F")
    # W is computed above, never assigned I - T T^dagger; that is a check only.
    return PassiveChannel(expm(segment.length*k), reservoir)


def propagate_segments(segments, parameters=Parameters()):
    """Segments are ordered from inflow to outflow; later maps multiply on LEFT."""
    transfer = np.eye(2, dtype=complex)
    reservoir = np.zeros((2, 2), complex)
    for segment in segments:
        local = segment_channel(segment, parameters)
        # Independent increments: W_new = E W_old E^dagger + W_local.
        reservoir = local.transfer @ reservoir @ local.transfer.conj().T + local.reservoir
        transfer = local.transfer @ transfer
    return PassiveChannel(transfer, reservoir)


def direct_profile_ode(kappa, length, parameters=Parameters(), *, rtol=2e-12, atol=2e-14):
    """Independent adaptive amplitude/covariance ODE for a continuous profile.

    Shares only the stated physical K,Q, not segment exponentials, midpoint
    samples, ordered products or the commutator complement. No identity is
    imposed on the independently integrated reservoir covariance.
    """
    Segment(length, 0.0)  # The same physical length domain applies to both solvers.
    if length == 0:
        return propagate_segments([], parameters)
    initial = np.concatenate((np.eye(2, dtype=complex).ravel(), np.zeros(4, complex)))
    q = _reservoir_density(parameters)

    def rhs(z, state):
        k = drift(kappa(z), parameters)
        transfer, reservoir = state[:4].reshape(2, 2), state[4:].reshape(2, 2)
        return np.concatenate(((k @ transfer).ravel(),
                               (k @ reservoir + reservoir @ k.conj().T + q).ravel()))

    solution = solve_ivp(rhs, (0.0, length), initial, method="DOP853", rtol=rtol, atol=atol)
    if not solution.success:
        raise RuntimeError(solution.message)
    final = solution.y[:, -1]
    return PassiveChannel(final[:4].reshape(2, 2), final[4:].reshape(2, 2))


def _occupation(value):
    if not np.isrealobj(value) or not np.isfinite(value) or value < 0:
        raise ValueError("thermal occupation must be finite, real and nonnegative")
    return float(value)


def optical_diagnostics(channel, *, atomic_occupation=0.0, reservoir_occupation=0.0):
    """Trace the atomic OUTPUT, retaining its independent INPUT port.

    Quadratures obey [x,p]=i and V_vac=I/2. With no anomalous moments,
    Y=[(n_a+1/2)|T_ba|^2+(n_r+1/2)W_bb] I. The CP matrix is
    Y+i(J-X J X^T)/2, not Y alone. Optical input is vacuum for output diagnostics.
    """
    na, nr = _occupation(atomic_occupation), _occupation(reservoir_occupation)
    t = channel.transfer[0, 0]
    tau = float(abs(t)**2)
    boundary = float(abs(channel.transfer[0, 1])**2)
    bath = float(channel.reservoir[0, 0].real)
    x = np.array([[t.real, -t.imag], [t.imag, t.real]])
    y_boundary, y_reservoir = (na+0.5)*boundary, (nr+0.5)*bath
    y = (y_boundary+y_reservoir)*np.eye(2)
    j = np.array([[0.0, 1.0], [-1.0, 0.0]])
    cp = y + 0.5j*(j-x @ j @ x.T)
    output_variance = tau/2 + y_boundary + y_reservoir
    return {
        "optical_input_commutator": tau,
        "atomic_boundary_commutator": boundary,
        "reservoir_commutator": bath,
        "commutator_sum": tau+boundary+bath,
        "commutator_error": abs(tau+boundary+bath-1),
        "atomic_occupation": na, "reservoir_occupation": nr,
        "boundary_quadrature_noise": y_boundary,
        "reservoir_quadrature_noise": y_reservoir,
        "X": x.tolist(), "Y": y.tolist(),
        "cp_eigenvalues": np.linalg.eigvalsh(cp).tolist(),
        "optical_vacuum_input_output_variance": output_variance,
        "optical_vacuum_input_output_occupation": na*boundary+nr*bath,
        "variance_occupation_identity_error": abs(output_variance-0.5-na*boundary-nr*bath),
    }


def atomic_source_covariance(positions, parameters=Parameters(), *,
                             atomic_occupation=0.0, reservoir_occupation=0.0):
    """Spatial kernel of the freely propagated atomic source F, BEFORE feedback.

    a(z)=F(z)-i integral_0^z g(z-s) kappa(s)* b(s) ds,
    F(z)=g(z)a_in+sqrt(2 gamma/v) integral_0^z g(z-s) f(s) ds.
    Local f increments are independent; propagated F(z) values are correlated.
    Returns full complex kernels, not diagonal approximations.
    """
    z = np.asarray(positions, float)
    if z.ndim != 1 or not np.all(np.isfinite(z)) or np.any(z < 0):
        raise ValueError("positions must be a finite nonnegative one-dimensional array")
    na, nr = _occupation(atomic_occupation), _occupation(reservoir_occupation)
    rate, frequency = parameters.gamma/parameters.v, parameters.omega/parameters.v
    lam = -rate+1j*frequency
    g = np.exp(lam*z)
    boundary = g[:, None]*g.conj()[None, :]
    m = np.minimum(z[:, None], z[None, :])
    # Integral of the common history up to min(z,z'); expm1 also handles gamma=0.
    reservoir = (np.exp(lam*(z[:, None]-m)+lam.conjugate()*(z[None, :]-m))
                 * (-np.expm1(-2*rate*m)))
    return {
        "boundary_commutator": boundary,
        "reservoir_commutator": reservoir,
        "total_commutator": boundary+reservoir,
        "boundary_symmetrized": (na+0.5)*boundary,
        "reservoir_symmetrized": (nr+0.5)*reservoir,
        "total_symmetrized": (na+0.5)*boundary+(nr+0.5)*reservoir,
    }


def _max_abs(matrix):
    return float(np.max(np.abs(matrix)))


def _channel_errors(channel):
    t, w = channel.transfer, channel.reservoir
    # The Hermitian part is used ONLY to report real eigenvalues, never as a repair.
    return {
        "full_commutator_error": _max_abs(t @ t.conj().T+w-np.eye(2)),
        "reservoir_hermiticity_error": _max_abs(w-w.conj().T),
        "reservoir_min_eigenvalue": float(np.linalg.eigvalsh((w+w.conj().T)/2).min()),
    }


def _complex_json(value):
    value = np.asarray(value)
    return {"real": value.real.tolist(), "imag": value.imag.tolist()}


def _kernel_control(segment, parameters, channel):
    """Output noise from a double integral of the nonlocal atomic source kernel."""
    optical = optical_diagnostics(channel)
    records = []
    for count in (48, 96, 192, 384):
        nodes, weights = np.polynomial.legendre.leggauss(count)
        z, weights = (nodes+1)*segment.length/2, weights*segment.length/2
        # Exact Volterra resolvent: response to optical forcing at z is U_bb(L,z).
        response = np.array([expm(drift(segment.kappa, parameters)*(segment.length-s))[0, 0]
                             for s in z])
        h = weights*(-1j*segment.kappa)*response
        kernels = atomic_source_covariance(z, parameters)
        boundary = float((h @ kernels["boundary_commutator"] @ h.conj()).real)
        bath = float((h @ kernels["reservoir_commutator"] @ h.conj()).real)
        # ANALYSIS-ONLY negative control: erase cross-z bath correlations, retaining
        # their diagonal entries. This PSD matrix is NOT an alternative physical bath.
        diagonal_bath = float(np.sum(abs(h)**2*np.diag(kernels["reservoir_commutator"]).real))
        records.append({
            "quadrature_nodes": count,
            "atomic_boundary_commutator": boundary,
            "reservoir_commutator": bath,
            "boundary_error": abs(boundary-optical["atomic_boundary_commutator"]),
            "reservoir_error": abs(bath-optical["reservoir_commutator"]),
            "full_commutator_error": abs(optical["optical_input_commutator"]+boundary+bath-1),
            "diagonalized_reservoir_commutator": diagonal_bath,
            "diagonalized_bath_commutator_deficit": 1-optical["optical_input_commutator"]-boundary-diagonal_bath,
        })
    positions = np.array([0.0, 0.2, 0.5, 0.8, 1.0])*segment.length
    kernels = atomic_source_covariance(positions, parameters)
    expected = np.exp(-parameters.gamma/parameters.v*abs(positions[:, None]-positions[None, :])
                      +1j*parameters.omega/parameters.v*(positions[:, None]-positions[None, :]))
    return {
        "positions": positions.tolist(),
        "boundary_commutator_kernel": _complex_json(kernels["boundary_commutator"]),
        "reservoir_commutator_kernel": _complex_json(kernels["reservoir_commutator"]),
        "vacuum_symmetrized_kernel": _complex_json(kernels["total_symmetrized"]),
        "stationary_kernel_identity_error": _max_abs(kernels["total_commutator"]-expected),
        "reservoir_offdiagonal_magnitude": float(abs(kernels["reservoir_commutator"][1, 3])),
        "reservoir_kernel_min_eigenvalue": float(np.linalg.eigvalsh(kernels["reservoir_commutator"]).min()),
        "quadrature_refinements": records,
        "negative_control_status": "analysis only; cross-z reservoir entries discarded, no noise repair",
    }


def build_control():
    """Return a deterministic, strict-JSON-compatible audit; write no files.

    All tolerances are declared validation gates, never noise fitting parameters.
    The smooth-profile ODE is independent of the finite-segment approximation.
    """
    p = Parameters()
    segment = Segment(1.3, 0.9+0.35j)
    channel = segment_channel(segment, p)
    vacuum = optical_diagnostics(channel)
    thermal = optical_diagnostics(channel, atomic_occupation=0.7, reservoir_occupation=0.3)
    direct = direct_profile_ode(lambda z: segment.kappa, segment.length, p)
    constant = {
        "transfer": _complex_json(channel.transfer),
        "reservoir_gramian": _complex_json(channel.reservoir),
        **_channel_errors(channel),
        "ode_transfer_error": _max_abs(channel.transfer-direct.transfer),
        "ode_reservoir_error": _max_abs(channel.reservoir-direct.reservoir),
        "vacuum": vacuum, "thermal": thermal,
    }

    def profile(z):
        u = z/segment.length
        return (0.85+0.20*np.cos(2*np.pi*u))*np.exp(1j*(0.35+1.1*u+0.25*np.sin(2*np.pi*u)))

    direct_profile = direct_profile_ode(profile, segment.length, p)
    refinements = []
    previous_error = None
    for count in (8, 16, 32, 64, 128):
        dz = segment.length/count
        segments = [Segment(dz, profile((i+0.5)*dz)) for i in range(count)]
        current = propagate_segments(segments, p)
        error_t = _max_abs(current.transfer-direct_profile.transfer)
        error_w = _max_abs(current.reservoir-direct_profile.reservoir)
        error = max(error_t, error_w)
        refinements.append({
            "segments": count, "transfer_error": error_t, "reservoir_error": error_w,
            "observed_order": None if previous_error is None else float(np.log2(previous_error/error)),
            **_channel_errors(current), "optical": optical_diagnostics(current),
        })
        previous_error = error
    reversed_channel = propagate_segments(reversed(segments), p)
    spatial = _kernel_control(segment, p, channel)
    boundary = vacuum["atomic_boundary_commutator"]
    reservoir = vacuum["reservoir_commutator"]
    negative = {
        "status": "analysis only; atomic inflow removed without any repair or fit",
        "remaining_added_noise_min_eigenvalue": reservoir/2,
        "commutator_sum_without_inflow": vacuum["optical_input_commutator"]+reservoir,
        "commutator_deficit": 1-vacuum["optical_input_commutator"]-reservoir,
        "expected_deficit_boundary": boundary,
        "deficit_identity_error": abs(1-vacuum["optical_input_commutator"]-reservoir-boundary),
        "gaussian_cp_min_eigenvalue_without_inflow": (reservoir-(1-vacuum["optical_input_commutator"]))/2,
        "apparent_vacuum_variance_without_inflow": (vacuum["optical_input_commutator"]+reservoir)/2,
        "reversed_profile_transfer_error": _max_abs(reversed_channel.transfer-direct_profile.transfer),
    }
    thresholds = {
        "algebraic_absolute": 5e-12, "ode_absolute": 2e-10,
        "smooth_profile_finest_absolute": 2e-5, "smooth_profile_minimum_order": 1.9,
        "spatial_kernel_quadrature_absolute": 5e-6,
        "negative_control_minimum_deficit": 0.05,
    }
    algebraic = thresholds["algebraic_absolute"]
    finest = refinements[-1]
    passed = bool(
        constant["full_commutator_error"] < algebraic
        and constant["reservoir_hermiticity_error"] < algebraic
        and constant["reservoir_min_eigenvalue"] >= -algebraic
        and max(constant["ode_transfer_error"], constant["ode_reservoir_error"]) < thresholds["ode_absolute"]
        and min(vacuum["cp_eigenvalues"]) >= -algebraic
        and min(thermal["cp_eigenvalues"]) >= -algebraic
        and thermal["variance_occupation_identity_error"] < algebraic
        and abs(vacuum["optical_vacuum_input_output_variance"]-0.5) < algebraic
        and max(finest["transfer_error"], finest["reservoir_error"]) < thresholds["smooth_profile_finest_absolute"]
        and min(row["observed_order"] for row in refinements[1:]) > thresholds["smooth_profile_minimum_order"]
        and max(row["full_commutator_error"] for row in refinements) < algebraic
        and max(row["reservoir_hermiticity_error"] for row in refinements) < algebraic
        and min(row["reservoir_min_eigenvalue"] for row in refinements) >= -algebraic
        and min(min(row["optical"]["cp_eigenvalues"]) for row in refinements) >= -algebraic
        and _channel_errors(direct_profile)["full_commutator_error"] < thresholds["ode_absolute"]
        and spatial["stationary_kernel_identity_error"] < algebraic
        and spatial["reservoir_kernel_min_eigenvalue"] >= -algebraic
        and spatial["reservoir_offdiagonal_magnitude"] > thresholds["negative_control_minimum_deficit"]
        and max(row["boundary_error"] for row in spatial["quadrature_refinements"]) < algebraic
        and spatial["quadrature_refinements"][-1]["full_commutator_error"] < thresholds["spatial_kernel_quadrature_absolute"]
        and spatial["quadrature_refinements"][-1]["diagonalized_bath_commutator_deficit"] > thresholds["negative_control_minimum_deficit"]
        and negative["remaining_added_noise_min_eigenvalue"] >= 0
        and negative["commutator_deficit"] > thresholds["negative_control_minimum_deficit"]
        and negative["deficit_identity_error"] < algebraic
        and negative["gaussian_cp_min_eigenvalue_without_inflow"] < -algebraic
        and negative["reversed_profile_transfer_error"] > thresholds["negative_control_minimum_deficit"]
    )
    return {
        "schema_version": 1,
        "source": "analysis/grand_challenge/reference/ballistic_linear_channel.py",
        "frequency": {"omega": p.omega, "units": "rad/t", "fourier_convention": "exp(-i Omega t)"},
        "units": {"length": "ell", "time": "t", "c": "ell/t", "v": "ell/t",
                  "gamma": "1/t", "kappa": "1/ell"},
        "scope": "passive two-mode finite-velocity bosonic testbed only",
        "exclusions": ["finite Rb atom", "FWM squeezing", "general moving-vapor closure"],
        "normalization": "canonical forward flux ports; [x,p]=i; vacuum V=I/2",
        "parameters": {**asdict(p), "length": segment.length, "constant_kappa": _complex_json(segment.kappa),
                       "units": "length ell, time t; kappa 1/ell, omega and gamma 1/t, c and v ell/t",
                       "optical_input_occupation": 0.0, "vacuum_environment_occupations": [0.0, 0.0],
                       "thermal_atomic_occupation": 0.7, "thermal_reservoir_occupation": 0.3},
        "precision_thresholds": thresholds,
        "constant_segment": constant,
        "nonconstant_profile": {
            "formula": "u=z/L; kappa=(0.85+0.20*cos(2*pi*u))*exp(i*(0.35+1.1*u+0.25*sin(2*pi*u)))",
            "sampling": "midpoint piecewise constant; later segments multiply on left",
            "ode_method": "DOP853", "ode_rtol": 2e-12, "ode_atol": 2e-14,
            "direct_ode": {"transfer": _complex_json(direct_profile.transfer),
                           "reservoir": _complex_json(direct_profile.reservoir),
                           **_channel_errors(direct_profile)},
            "refinements": refinements,
        },
        "spatial_atomic_source": spatial,
        "negative_controls": negative,
        "expected_controls_passed": passed,
        "passed": passed,
    }


if __name__ == "__main__":
    print(json.dumps(build_control(), indent=2, allow_nan=False))
