"""Multihyperfine NMOR with velocity-randomized coated-cell ground memory.

The optical excited manifold is eliminated to second order in the optical
coupling using the effective-operator construction of Reiter and Sorensen,
PRA 85, 032111 (2012), Eqs. (3)-(5). All hyperfine/Zeeman dipole amplitudes
and their signs are retained. Ground-F coherences are secularly discarded;
Zeeman coherences, coherent light shifts and open hyperfine pumping remain.
This is a weak-optical-excitation approximation, not a high-power full OBE.

Coated-cell velocity mixing uses one dark reservoir shared by every velocity
class. Its Schur complement is solved exactly for the selected quadrature.
The ``none`` option instead associates an independent dark state with each
class, providing the velocity-preserving limit. All rate arguments ending
in ``_hz`` are exponential rates in s^-1, not rates divided by 2*pi.

Budker et al., PRA 62, 043403 (2000), motivates the wall/transit comparison.
No empirical detuning factor or prescribed off-resonant optimum is used.
"""
from dataclasses import dataclass
from functools import lru_cache
import math

import numpy as np
from scipy.linalg.lapack import get_lapack_funcs
from numpy.polynomial.legendre import leggauss

from . import constants, core, species, zeeman
from .magneto_coordinates import ground_coordinates
from .magneto_optical import assemble_optical_generator


MODEL_SCOPE = (
    "weak-optical-excitation effective-operator model; all F,m dipole channels; "
    "ground-F secular approximation; coherent light shifts; conditional local "
    "uniform-intensity response; optical Zeeman shifts and linear Faraday "
    "background omitted; no spin-exchange, radiation trapping or propagation"
)


@dataclass(frozen=True)
class MultilevelAtom:
    isotope: str
    line: str
    ground: tuple
    excited: tuple
    ground_slices: tuple
    reduced_slices: tuple
    dipoles: dict
    ground_energy_mhz: np.ndarray
    excited_energy_mhz: np.ndarray
    gamma: float
    frequency_hz: float
    dipole_sq: float
    mass: float
    thermal: np.ndarray
    trace: np.ndarray
    ground_g: tuple

    @property
    def dimension(self):
        return self.thermal.size


@lru_cache(maxsize=4)
def build_multilevel_atom(isotope="Rb85", line="D1"):
    """Full angular dipoles plus secular ground-density coordinates.

    Dipole matrices map excited to ground states and satisfy
    sum_q C_q^dagger C_q = I_excited. Thus sqrt(Gamma)*C_q are normalized
    emission jumps and d0*C_q are physical dipoles, with
    d0^2=3*pi*epsilon0*hbar*c^3*Gamma/omega^3.
    """
    if isotope not in ("Rb85", "Rb87") or line not in ("D1", "D2"):
        raise ValueError("multilevel NMOR supports Rb85/Rb87 and D1/D2")
    iso = species.ISOTOPES[isotope]
    Je, nu, width, Ae, Be = iso.line(line)
    fs = tuple(int(f) for f in species.f_values(iso.I, iso.Jg))
    fes = tuple(int(f) for f in species.f_values(iso.I, Je))
    ground = tuple((f, m) for f in fs for m in range(-f, f+1))
    excited = tuple((f, m) for f in fes for m in range(-f, f+1))
    dipoles = {q: np.zeros((len(ground), len(excited))) for q in (-1, 0, 1)}
    for g, (f, m) in enumerate(ground):
        for e, (fe, me) in enumerate(excited):
            q = me-m
            if q not in dipoles:
                continue
            reduced = ((-1)**int(round(Je+iso.I+f+1))
                       * math.sqrt((2*Je+1)*(2*f+1))
                       * species.wigner6j(Je, fe, iso.I, f, iso.Jg, 1))
            dipoles[q][g, e] = reduced * zeeman.clebsch_gordan(f, m, 1, q, fe, me)
    gs, rs, i, j = [], [], 0, 0
    for f in fs:
        n = 2*f+1
        gs.append(slice(i, i+n)); rs.append(slice(j, j+n*n))
        i += n; j += n*n
    trace = np.zeros(j)
    for f, block in zip(fs, rs):
        n = 2*f+1
        trace[block] = np.eye(n).ravel()
    thermal = trace / len(ground)
    gamma = width * species.MHZ
    omega = 2*np.pi*nu
    d2 = 3*np.pi*constants.EPS_0*constants.HBAR*constants.C_LIGHT**3*gamma/omega**3
    return MultilevelAtom(
        isotope, line, ground, excited, tuple(gs), tuple(rs), dipoles,
        np.array([species.hf_energy_mhz(iso.A_S, 0, iso.I, iso.Jg, f) for f in fs]),
        np.array([species.hf_energy_mhz(Ae, Be, iso.I, Je, f) for f, _ in excited]),
        gamma, nu, d2, iso.mass, thermal, trace,
        tuple(species.lande_gf(species.GJ_S12_LEGACY, iso.g_I, iso.I, iso.Jg, f)
              for f in fs))


def component_centroid_mhz(atom, ground_f=None):
    """Optical ground-component line-strength centroid from fine centroid."""
    iso = species.ISOTOPES[atom.isotope]
    Je = iso.line(atom.line)[0]
    fs = [atom.ground[s.start][0] for s in atom.ground_slices]
    f = max(fs) if ground_f is None else int(ground_f)
    if f not in fs:
        raise ValueError("ground_f is not an allowed ground hyperfine level")
    fes = sorted(set(fe for fe, _ in atom.excited))
    strengths = np.array([species.line_strength(f, fe, iso.I, iso.Jg, Je) for fe in fes])
    energies = np.array([atom.excited_energy_mhz[next(i for i, (F, _) in
                                                enumerate(atom.excited) if F == fe)] for fe in fes])
    return float(strengths @ energies / strengths.sum() - atom.ground_energy_mhz[fs.index(f)])


def _polarization(qwp_deg):
    th = math.radians(qwp_deg)
    c, s = math.cos(th), math.sin(th)
    ex, ey = c*c+1j*s*s, (1-1j)*s*c
    return {1: (ex-1j*ey)/math.sqrt(2), -1: (ex+1j*ey)/math.sqrt(2)}


def velocity_quadrature(atom, laser_from_centroid_mhz, temp_k, *, order=6, doppler=True):
    """Resonance-centered composite Gauss-Legendre Maxwell quadrature.

    Integration variable is optical Doppler shift kv/(2*pi), in MHz. Knots
    bracket every allowed hyperfine optical resonance and its Lorentzian
    wings, in addition to a coarse Gaussian grid. No low-order Hermite grid
    is used. The omitted Gaussian tail outside +/-6 sigma is reported.
    ``order`` is nodes per interval, not the resulting velocity count.
    """
    if not doppler:
        return np.array([0.0]), np.array([1.0]), 0.0
    sigma = atom.frequency_hz/constants.C_LIGHT * math.sqrt(constants.KB*temp_k/atom.mass) / 1e6
    hwhm = atom.gamma/(2*species.MHZ)
    lo, hi = -6*sigma, 6*sigma
    knots = list(np.linspace(lo, hi, 13))
    for bi, block in enumerate(atom.ground_slices):
        strength = sum(np.sum(atom.dipoles[q][block]**2, axis=0) for q in (-1, 0, 1))
        for energy in np.unique(atom.excited_energy_mhz[strength > 1e-15]):
            center = laser_from_centroid_mhz - (energy-atom.ground_energy_mhz[bi])
            for scale in (-32, -8, -2, -0.5, 0, 0.5, 2, 8, 32):
                value = center + scale*hwhm
                if lo < value < hi:
                    knots.append(value)
    knots = np.unique(knots)
    nodes, weights = leggauss(int(order))
    half = np.diff(knots)/2
    center = (knots[1:]+knots[:-1])/2
    x = (center[:, None] + half[:, None]*nodes).ravel()
    w = (half[:, None]*weights).ravel() * np.exp(-0.5*(x/sigma)**2)/(math.sqrt(2*np.pi)*sigma)
    # Conditional finite-tail measure; normalization differs by <2e-9.
    w /= w.sum()
    return x, w, sigma


def _left_right(K):
    """Batched row-major left/right matrix-action superoperators."""
    n = K.shape[-1]
    eye = np.eye(n)
    left = np.einsum("vai,bj->vabij", K, eye).reshape(-1, n*n, n*n)
    right = np.einsum("ai,vjb->vabij", eye, K).reshape(-1, n*n, n*n)
    return left, right


def _optical_terms(atom, laser_eff_mhz, intensity_mw_cm2, polarization):
    field = math.sqrt(2*intensity_mw_cm2*10/(constants.EPS_0*constants.C_LIGHT))
    omega_scale = math.sqrt(atom.dipole_sq)*field/constants.HBAR
    excitation = sum(polarization[q]*atom.dipoles[q].T for q in (-1, 1))
    V = -0.5*omega_scale*excitation
    eliminated = []
    resolvents = []
    for bi, (gs, rs) in enumerate(zip(atom.ground_slices, atom.reduced_slices)):
        det = (laser_eff_mhz[:, None] - (atom.excited_energy_mhz[None, :]
               - atom.ground_energy_mhz[bi]))*species.MHZ
        R = 1/(-det-0.5j*atom.gamma)
        Y = R[:, :, None]*V[None, :, gs]
        eliminated.append(Y)
        resolvents.append(R)
    L = assemble_optical_generator(atom, resolvents, polarization, omega_scale**2)
    return L, eliminated, resolvents, excitation


@lru_cache(maxsize=12)
def _magnetic_generator(isotope, line, axis="z"):
    """Ground commutator derivative dL/dB_axis, with B in SI Tesla."""
    atom = build_multilevel_atom(isotope, line)
    result = np.zeros((atom.dimension, atom.dimension), complex)
    ai = {"x": 0, "y": 1, "z": 2}[axis]
    for gs, rs, gf in zip(atom.ground_slices, atom.reduced_slices, atom.ground_g):
        f = atom.ground[gs.start][0]
        angular = zeeman.angular_momentum_matrices(f)[ai]
        hamiltonian = constants.MU_B/constants.HBAR*gf*angular
        n = len(angular)
        result[rs, rs] = -1j*(np.kron(hamiltonian, np.eye(n))
                                    - np.kron(np.eye(n), hamiltonian.T))
    result.flags.writeable = False
    return result


def _ground_generator(atom, b_ut, reset_rate, pure_dephasing_rate,
                      b_x_ut=0.0, b_y_ut=0.0):
    M = atom.dimension
    L = reset_rate*(np.outer(atom.thermal, atom.trace)-np.eye(M)).astype(complex)
    for gs, rs, gf in zip(atom.ground_slices, atom.reduced_slices, atom.ground_g):
        m = np.array([level[1] for level in atom.ground[gs]])
        L[rs, rs] -= np.diag(pure_dephasing_rate*(1-np.eye(len(m))).ravel())
    for axis, field in (("x", b_x_ut), ("y", b_y_ut), ("z", b_ut)):
        if field:
            L += float(field)*1e-6*_magnetic_generator(atom.isotope, atom.line, axis)
    return L


@lru_cache(maxsize=2)
def _normalization_lapack(dtype):
    return get_lapack_funcs(("getrf", "getrs"), dtype=np.dtype(dtype))


def _stationary_factored(generator, trace):
    """Normalized states and reusable LU, with any leading velocity batch."""
    A = np.array(generator, dtype=np.result_type(generator, np.float64), copy=True)
    scale = np.maximum(np.max(np.abs(A), axis=(-2, -1)), 1.)
    A[..., 0, :] = trace*scale[..., None]
    rhs = np.zeros(A.shape[:-1] + (1,), dtype=A.dtype)
    rhs[..., 0, 0] = scale
    # Low-level LAPACK keeps the same LU reusable on SciPy>=1.14, including
    # versions predating scipy.linalg's batched LU wrappers.
    matrices = A.reshape(-1, A.shape[-1], A.shape[-1])
    sources = rhs.reshape(-1, A.shape[-1], 1)
    lu = np.empty_like(matrices)
    pivots = np.empty(matrices.shape[:-1], dtype=np.intc)
    state = np.empty_like(sources)
    getrf, getrs = _normalization_lapack(A.dtype.str)
    for i, (matrix, source) in enumerate(zip(matrices, sources)):
        lu[i], pivots[i], info = getrf(matrix, overwrite_a=True)
        if info:
            raise np.linalg.LinAlgError("singular normalized steady-state generator")
        state[i], info = getrs(lu[i], pivots[i], source, overwrite_b=True)
        if info:
            raise np.linalg.LinAlgError("normalized steady-state solve failed")
    return state.reshape(A.shape[:-1]), (lu, pivots, A.shape[:-2])


def _stationary(L, trace):
    return _stationary_factored(L, trace)[0]


def _response_from_factor(factors, forcing):
    lu, pivots, batch_shape = factors
    rhs = -np.array(forcing, dtype=lu.dtype, copy=True)[..., None]
    rhs[..., 0, 0] = 0.
    sources = rhs.reshape(-1, rhs.shape[-2], 1)
    result = np.empty_like(sources)
    _, getrs = _normalization_lapack(lu.dtype.str)
    for i, source in enumerate(sources):
        result[i], info = getrs(lu[i], pivots[i], source, overwrite_b=True)
        if info:
            raise np.linalg.LinAlgError("normalized magnetic-response solve failed")
    return result.reshape(batch_shape + (rhs.shape[-2],))


def _trace_zero_response(generator, forcing, trace):
    """Differentiate a normalized steady state: L rho' = -forcing, tr rho'=0."""
    _, factors = _stationary_factored(generator, trace)
    return _response_from_factor(factors, forcing)


def _reservoir_states(light, dark, weights, gamma_out, gamma_in, mixing, trace, thermal):
    """Conditional trace-one light states; exact block Schur elimination."""
    states, reservoir, _, _ = _reservoir_solution(
        light, dark, weights, gamma_out, gamma_in, mixing, trace, thermal)
    return states, reservoir


def _reservoir_solution(light, dark, weights, gamma_out, gamma_in, mixing,
                        trace, thermal, magnetic_generator=None):
    """Common-reservoir steady state and exact normalized DC Bz response.

    Derivatives are with respect to Tesla and include the evolution of the
    dark state. They are not derivatives of a frozen-reservoir approximation.
    """
    nv, M, _ = light.shape
    eye = np.eye(M)
    if gamma_out == 0:
        states, factors = _stationary_factored(light, trace)
        response = (None if magnetic_generator is None else
                    _response_from_factor(factors, states@magnetic_generator.T))
        return states, None, response, None
    if gamma_in == 0:
        # Conditional gamma_in -> 0+ limit: the reset dark reservoir is
        # thermal. It supplies a transit-reset source at gamma_out. Absolute
        # light-region occupation is not inferred from this conditional limit.
        refill = gamma_out*(np.outer(thermal, trace)-eye)
        states, factors = _stationary_factored(light+refill, trace)
        response = (None if magnetic_generator is None else
                    _response_from_factor(factors, states@magnetic_generator.T))
        return states, thermal, response, np.zeros(M, dtype=light.dtype) if response is not None else None
    # NumPy handles the velocity batch within its solve instead of one Python
    # scipy.solve call per class. Keep the complete matrix maps: the shared
    # reservoir is unknown until their weighted matrix average is assembled.
    maps = np.linalg.solve(gamma_out*eye-light, np.broadcast_to(gamma_out*eye, light.shape))
    if mixing == "complete":
        mean_map = np.einsum("v,vij->ij", weights, maps)
        dark_generator = dark + gamma_in*(mean_map-eye)
        dark_state, factors = _stationary_factored(dark_generator, trace)
        states = maps @ dark_state
    else:
        dark_generators = dark[None]+gamma_in*(maps-eye)
        dark_state, factors = _stationary_factored(dark_generators, trace)
        states = np.einsum("vij,vj->vi", maps, dark_state)
    if magnetic_generator is None:
        return states, dark_state, None, None
    # A_v^-1(dL/dB rho_v), using the already factored maps g_out A_v^-1.
    forcing = states@magnetic_generator.T
    local_response = np.einsum("vij,vj->vi", maps, forcing)/gamma_out
    if mixing == "complete":
        source = magnetic_generator@dark_state + gamma_in*(weights@local_response)
        dark_response = _response_from_factor(factors, source)
        response = local_response + maps@dark_response
    else:
        source = dark_state@magnetic_generator.T + gamma_in*local_response
        dark_response = _response_from_factor(factors, source)
        response = local_response + np.einsum("vij,vj->vi", maps, dark_response)
    return states, dark_state, response, dark_response


def _readout(atom, states, resolvents, excitation, eliminated, weights, polarization, density,
             check_density_state=True):
    pol = {q: np.zeros(states.shape[0], complex) for q in (-1, 1)}
    excitation_fraction = np.zeros(states.shape[0])
    min_eigenvalue = 1.0
    hermiticity_error = 0.0
    for gs, rs, R, Y in zip(atom.ground_slices, atom.reduced_slices, resolvents, eliminated):
        n = gs.stop-gs.start
        rho = states[:, rs].reshape(-1, n, n)
        if check_density_state:
            hermiticity_error = max(hermiticity_error, float(np.max(np.abs(rho-rho.conj().transpose(0, 2, 1)))))
            min_eigenvalue = min(min_eigenvalue, float(np.min(np.linalg.eigvalsh(0.5*(rho+rho.conj().transpose(0, 2, 1))))))
        optical = np.einsum("ve,ei,vij->vej", R, excitation[:, gs], rho)
        for q in (-1, 1):
            pol[q] += np.einsum("ie,vei->v", atom.dipoles[q][gs], optical)
        if check_density_state:
            excitation_fraction += np.einsum("vei,vij,vej->v", Y, rho, Y.conj()).real
    prefactor = density*atom.dipole_sq/(constants.EPS_0*constants.HBAR)
    averaged = {q: prefactor*np.dot(weights, pol[q]) for q in (-1, 1)}
    probe = sum(np.conj(polarization[q])*averaged[q] for q in (-1, 1))
    components = {q: averaged[q]/polarization[q] if abs(polarization[q]) > 1e-12
                  else complex(np.nan, np.nan) for q in (-1, 1)}
    diagnostics = dict(
        excited_fraction_mean=float(weights@excitation_fraction),
        excited_fraction_max=float(np.max(excitation_fraction)),
        min_ground_eigenvalue=min_eigenvalue,
        hermiticity_error=hermiticity_error,
        trace_error=float(np.max(np.abs(states@atom.trace-1))),
    ) if check_density_state else {}
    return components, probe, diagnostics


def _magnetic_response_error(light, dark, states, reservoir, response,
                             reservoir_response, derivative, weights,
                             gamma_out, gamma_in, mixing, trace):
    """Residual of the differentiated coupled equations, not a finite difference."""
    forcing = states@derivative.T
    error = np.einsum("vij,vj->vi", light, response)+forcing
    scale = max(float(np.max(np.abs(forcing))), 1.0)
    response_scale = max(float(np.max(np.abs(response))), 1.0)
    trace_error = float(np.max(np.abs(response@trace)))/response_scale
    if reservoir is not None:
        error += gamma_out*(reservoir_response-response)
        if gamma_in > 0:
            source = np.atleast_2d(reservoir)@derivative.T
            feed = weights@response if mixing == "complete" else response
            dark_error = (np.atleast_2d(reservoir_response)@dark.T + source
                          + gamma_in*(feed-reservoir_response))
            scale = max(scale, float(np.max(np.abs(source))))
            residual = max(float(np.max(np.abs(error))), float(np.max(np.abs(dark_error))))/scale
        else:
            residual = float(np.max(np.abs(error)))/scale
        response_scale = max(response_scale, float(np.max(np.abs(reservoir_response))))
        trace_error = max(trace_error, float(np.max(np.abs(np.atleast_2d(reservoir_response)@trace)))/response_scale)
    else:
        residual = float(np.max(np.abs(error)))/scale
    return residual, trace_error


def rescale_readout(result, *, density_m3, cell_length_m):
    """Re-evaluate the homogeneous-path readout without another atomic solve.

    Density and path length do not enter this local uniform-intensity model's
    light/dark equations. Keep the solved-state diagnostics and analytic-slope
    validity; only susceptibility scaling and optical readout change. A cached
    positive reference density retains the response even for a requested vacuum.
    """
    density, length = float(density_m3), float(cell_length_m)
    reference_density = float(result["density_m3"])
    reference_length = float(result["cell_length_m"])
    if not np.all(np.isfinite([density, length])) or min(density, length) < 0:
        raise ValueError("finite nonnegative density and cell length required")
    if reference_density == 0 and density > 0:
        raise ValueError("positive reference density needed to restore a nonzero-density readout")
    ratio = density/reference_density if reference_density else 0.
    scaled = dict(result, density_m3=density, cell_length_m=length)
    for name in ("chi_plus", "chi_minus", "chi_probe", "chi_plus_slope_per_t",
                 "chi_minus_slope_per_t", "chi_probe_slope_per_t"):
        if name in result:
            scaled[name] = np.asarray(result[name])*ratio
    if "chi_probe" in scaled:
        k = 2*np.pi*float(result["frequency_hz"])/constants.C_LIGHT
        scaled["optical_depth"] = k*length*scaled["chi_probe"].imag
        scaled["transmission"] = np.exp(-scaled["optical_depth"])
        if "chi_probe_slope_per_t" in scaled:
            scaled["transmission_slope_per_t"] = (
                -scaled["transmission"]*k*length*scaled["chi_probe_slope_per_t"].imag)
        linear = result.get("rotation_status", "").startswith("linear-input")
        if linear:
            scaled["rotation_rad"] = 0.25*k*length*(scaled["chi_plus"]-scaled["chi_minus"]).real
            if "chi_plus_slope_per_t" in scaled:
                scaled["rotation_slope_rad_per_t"] = (
                    0.25*k*length*(scaled["chi_plus_slope_per_t"]
                                      -scaled["chi_minus_slope_per_t"]).real)
        return scaled
    # Compatibility for saved/mock readouts lacking susceptibility arrays.
    # Production caches retain χ, which also handles underflowed transmission.
    if reference_length == 0 and length > 0:
        raise ValueError("positive reference length or susceptibility arrays required")
    path_ratio = ratio*length/reference_length if reference_length else 0.
    old_t = np.asarray(result["transmission"])
    with np.errstate(divide="ignore"):
        old_depth = np.asarray(result.get("optical_depth", -np.log(old_t)))
    scaled["optical_depth"] = old_depth*path_ratio if path_ratio else np.zeros_like(old_depth)
    scaled["transmission"] = np.exp(-scaled["optical_depth"])
    for name in ("rotation_rad", "rotation_slope_rad_per_t"):
        if name in result:
            scaled[name] = np.asarray(result[name])*path_ratio
    if "transmission_slope_per_t" in result:
        normalized_slope = np.divide(np.asarray(result["transmission_slope_per_t"]),
                                     old_t, out=np.zeros_like(old_t, dtype=float), where=old_t > 0)
        scaled["transmission_slope_per_t"] = scaled["transmission"]*path_ratio*normalized_slope
    return scaled


def multilevel_response(detunings_mhz, *, isotope="Rb85", line="D1", ground_f=None,
                        intensity_mw_cm2=1.0, temp_c=20.0, b_ut=0.0,
                        b_x_ut=0.0, b_y_ut=0.0, magnetic_derivative=True,
                        cell_length_m=0.1, gamma_out_hz=5e4, gamma_in_hz=20.0,
                        ground_relaxation_hz=1.0, wall_relaxation_hz=1.0,
                        pure_dephasing_hz=0.0, velocity_mixing="complete",
                        velocity_points=6, doppler=True, qwp_deg=0.0,
                        density_m3=None, density_scale=1.0, solver_basis="auto"):
    """Return local susceptibilities and homogeneous-path transmission/rotation.

    Every numerical output curve has shape (len(b_ut), len(detunings_mhz)),
    including scalar inputs. Optical detuning zero is the addressed ground-F
    component's line-strength centroid. ``density_m3`` denotes this isotope's
    density; an abundance must not be applied again. ``velocity_points`` is
    the Gauss-Legendre order per resonance-centered interval (2..20).

    ``rotation_rad`` uses the small-susceptibility circular-birefringence
    formula only for linear input. Elliptical/circular input returns NaN
    rotation but still supplies the probe susceptibility and transmission.
    ``magnetic_derivative`` differentiates the entire normalized light/dark
    steady-state problem analytically in longitudinal SI Tesla, including
    the common reservoir response. No finite B stencil is needed for slope.

    ``solver_basis='auto'`` uses real Hermitian coordinates and, for exactly
    zero transverse field, the invariant even-Δm sector. ``'complex'`` retains
    every original density coordinate as a numerical reference. Both restore
    the full states for the same readout and physical residual diagnostics.
    """
    delta = np.atleast_1d(np.asarray(detunings_mhz, dtype=float))
    fields = np.atleast_1d(np.asarray(b_ut, dtype=float))
    settings = np.array([intensity_mw_cm2, temp_c, cell_length_m, gamma_out_hz,
                         gamma_in_hz, ground_relaxation_hz, wall_relaxation_hz,
                         pure_dephasing_hz, density_scale, qwp_deg,
                         b_x_ut, b_y_ut], dtype=float)
    if (delta.ndim != 1 or fields.ndim != 1 or not delta.size or not fields.size
            or not np.all(np.isfinite(delta)) or not np.all(np.isfinite(fields))
            or not np.all(np.isfinite(settings)) or temp_c <= -273.15):
        raise ValueError("finite one-dimensional axes and physical finite settings required")
    if min(intensity_mw_cm2, cell_length_m, gamma_out_hz, gamma_in_hz,
           ground_relaxation_hz, wall_relaxation_hz, pure_dephasing_hz, density_scale) < 0:
        raise ValueError("intensity, length, rates and density scale must be nonnegative")
    if ground_relaxation_hz <= 0 or wall_relaxation_hz <= 0:
        raise ValueError("positive population-reset rates are needed for a unique steady state")
    if velocity_mixing not in ("complete", "none"):
        raise ValueError("velocity_mixing must be complete or none")
    if solver_basis not in ("auto", "complex"):
        raise ValueError("solver_basis must be auto or complex")
    if int(velocity_points) != velocity_points or not 2 <= velocity_points <= 20:
        raise ValueError("velocity_points is quadrature order per interval, 2..20")
    atom = build_multilevel_atom(isotope, line)
    even_coherences = float(b_x_ut) == 0. and float(b_y_ut) == 0.
    coordinates = (ground_coordinates(tuple(atom.ground[s.start][0]
                    for s in atom.ground_slices), even_coherences)
                   if solver_basis == "auto" else None)
    solve_trace = atom.trace if coordinates is None else coordinates.pack(atom.trace)
    solve_thermal = atom.thermal if coordinates is None else coordinates.pack(atom.thermal)
    centroid = component_centroid_mhz(atom, ground_f)
    temp_k = temp_c+273.15
    density = species.number_density(species.ISOTOPES[isotope], temp_k) if density_m3 is None else float(density_m3)
    if not np.isfinite(density) or density < 0:
        raise ValueError("density_m3 must be finite and nonnegative")
    density *= density_scale
    polarization = _polarization(qwp_deg)
    shape = (fields.size, delta.size)
    cp, cm, probe = [np.empty(shape, complex) for _ in range(3)]
    dcp, dcm, dprobe = [np.full(shape, complex(np.nan, np.nan)) for _ in range(3)]
    slope_residual = np.full(shape, np.nan)
    slope_valid = np.zeros(shape, bool)
    derivative = (_magnetic_generator(isotope, line, "z")
                  if magnetic_derivative else None)
    solve_derivative = (derivative if coordinates is None or derivative is None
                        else coordinates.generator(derivative))
    diagnostics = []
    counts = []
    with core.blas_single_thread():
        for di, detuning in enumerate(delta):
            shifts, weights, sigma = velocity_quadrature(
                atom, centroid+detuning, temp_k, order=velocity_points, doppler=doppler)
            counts.append(shifts.size)
            optical, eliminated, resolvents, excitation = _optical_terms(
                atom, centroid+detuning-shifts, intensity_mw_cm2, polarization)
            solve_optical = optical if coordinates is None else coordinates.generator(optical)
            virtual_bound = max(float(np.max(np.linalg.eigvalsh(
                np.einsum("vei,vej->vij", Y.conj(), Y)))) for Y in eliminated)
            for bi, field in enumerate(fields):
                ground = _ground_generator(atom, field, ground_relaxation_hz,
                    pure_dephasing_hz, b_x_ut, b_y_ut)
                light = optical + ground
                dark = _ground_generator(atom, field, wall_relaxation_hz,
                    pure_dephasing_hz, b_x_ut, b_y_ut)
                solve_light = (light if coordinates is None else
                               solve_optical+coordinates.generator(ground))
                solve_dark = dark if coordinates is None else coordinates.generator(dark)
                states, dark_state, state_response, dark_response = _reservoir_solution(
                    solve_light, solve_dark, weights, gamma_out_hz, gamma_in_hz,
                    velocity_mixing, solve_trace, solve_thermal, solve_derivative)
                if coordinates is not None:
                    states, dark_state, state_response, dark_response = (
                        coordinates.expand(value) for value in
                        (states, dark_state, state_response, dark_response))
                components, response, diagnostic = _readout(atom, states, resolvents,
                    excitation, eliminated, weights, polarization, density)
                cp[bi, di], cm[bi, di], probe[bi, di] = components[1], components[-1], response
                if magnetic_derivative:
                    dcomponents, dresponse, _ = _readout(atom, state_response, resolvents,
                        excitation, eliminated, weights, polarization, density,
                        check_density_state=False)
                    dcp[bi, di], dcm[bi, di], dprobe[bi, di] = dcomponents[1], dcomponents[-1], dresponse
                    dr, dt = _magnetic_response_error(light, dark, states, dark_state,
                        state_response, dark_response, derivative, weights,
                        gamma_out_hz, gamma_in_hz, velocity_mixing, atom.trace)
                    slope_residual[bi, di] = max(dr, dt)
                    slope_valid[bi, di] = (np.isfinite(slope_residual[bi, di])
                                           and slope_residual[bi, di] <= 1e-8)
                residual = np.einsum("vij,vj->vi", light, states)
                if dark_state is not None:
                    residual += gamma_out_hz*(dark_state-states)
                diagnostic["steady_state_relative_residual"] = float(np.max(np.abs(residual)) / max(1.0, np.max(np.abs(light)), gamma_out_hz))
                diagnostic["dark_state_relative_residual"] = 0.0
                if dark_state is not None:
                    dark_states = np.atleast_2d(dark_state)
                    diagnostic["trace_error"] = max(diagnostic["trace_error"],
                        float(np.max(np.abs(dark_states@atom.trace-1))))
                    for block in atom.reduced_slices:
                        n = round(math.sqrt(block.stop-block.start))
                        matrix = dark_states[:, block].reshape(-1, n, n)
                        adjoint = matrix.conj().transpose(0, 2, 1)
                        diagnostic["hermiticity_error"] = max(diagnostic["hermiticity_error"],
                            float(np.max(np.abs(matrix-adjoint))))
                        diagnostic["min_ground_eigenvalue"] = min(diagnostic["min_ground_eigenvalue"],
                            float(np.min(np.linalg.eigvalsh(0.5*(matrix+adjoint)))))
                    if gamma_in_hz > 0:
                        feed = weights@states if velocity_mixing == "complete" else states
                        dark_residual = dark_states@dark.T + gamma_in_hz*(feed-dark_state)
                        diagnostic["dark_state_relative_residual"] = float(np.max(np.abs(dark_residual))
                            / max(1.0, np.max(np.abs(dark)), gamma_in_hz))
                diagnostic["virtual_excitation_bound"] = virtual_bound
                diagnostics.append(diagnostic)
    k = 2*np.pi*atom.frequency_hz/constants.C_LIGHT
    optical_depth = k*probe.imag*cell_length_m
    transmission = np.exp(-optical_depth)
    linear = math.isclose(abs(polarization[1]), abs(polarization[-1]), abs_tol=1e-12)
    rotation = 0.25*k*cell_length_m*(cp-cm).real if linear else np.full(shape, np.nan)
    rotation_slope = (0.25*k*cell_length_m*(dcp-dcm).real if linear
                      else np.full(shape, np.nan))
    # With exactly linear QWP settings and no transverse field, reversing Bz
    # exchanges the sigma+/- channels. Their incident-power-weighted absorption
    # is therefore even in Bz, so d Im(chi_probe)/dBz is exactly zero at Bz=0.
    # Enforce this proved parity identity rather than converting roundoff into
    # a very large but finite field sensitivity. No small-number cutoff is used.
    symmetry_rows = ((fields == 0.0) & (float(b_x_ut) == 0.0)
                     & (float(b_y_ut) == 0.0) & (float(qwp_deg) % 90.0 == 0.0)
                     & bool(magnetic_derivative))
    dprobe[symmetry_rows] = dprobe[symmetry_rows].real + 0j
    transmission_slope = -transmission*k*cell_length_m*dprobe.imag
    transmission_symmetry_zero = np.broadcast_to(symmetry_rows[:, None], shape).copy()
    health = {name: max(d[name] for d in diagnostics) for name in diagnostics[0] if name != "min_ground_eigenvalue"}
    health["min_ground_eigenvalue"] = min(d["min_ground_eigenvalue"] for d in diagnostics)
    health["optical_excitation_status"] = ("weak excitation" if health["virtual_excitation_bound"] <= 0.05
                                           else "optical elimination requires full-OBE validation")
    health["solver_status"] = ("physical steady state" if
        health["min_ground_eigenvalue"] >= -1e-9 and
        max(health[key] for key in ("trace_error", "hermiticity_error",
            "steady_state_relative_residual", "dark_state_relative_residual")) <= 1e-8
        else "steady-state validation failed")
    health["magnetic_response_status"] = (
        "not requested" if not magnetic_derivative else
        "analytic coupled-reservoir response validated" if np.all(slope_valid)
        else "analytic magnetic response validation failed")
    if health["solver_status"] != "physical steady state":
        slope_valid[:] = False
    return dict(
        detunings_mhz=delta, b_ut=fields, chi_plus=cp, chi_minus=cm, chi_probe=probe,
        transmission=transmission, optical_depth=optical_depth, rotation_rad=rotation,
        chi_plus_slope_per_t=dcp, chi_minus_slope_per_t=dcm, chi_probe_slope_per_t=dprobe,
        rotation_slope_rad_per_t=rotation_slope, transmission_slope_per_t=transmission_slope,
        transmission_slope_symmetry_zero=transmission_symmetry_zero,
        slope_valid=slope_valid, slope_relative_residual=slope_residual,
        slope_method="analytic normalized steady-state response; light and dark states differentiated",
        transverse_field_ut=dict(x=float(b_x_ut), y=float(b_y_ut)),
        density_m3=density, cell_length_m=float(cell_length_m),
        isotope=isotope, line=line, frequency_hz=atom.frequency_hz,
        ground_f=max(f for f, _ in atom.ground) if ground_f is None else int(ground_f),
        reference_centroid_mhz=centroid,
        detuning_reference="line-strength centroid of selected ground-F optical component",
        model_scope=MODEL_SCOPE, diagnostics=health, velocity_mixing=velocity_mixing,
        exchange_status=("no dark exchange" if gamma_out_hz == 0 else
                         "conditional zero-return thermal transit limit" if gamma_in_hz == 0
                         else "conditional light state with coated-cell exchange"),
        velocity_quadrature=dict(method="resonance-centered composite Gauss-Legendre Maxwell",
            interval_order=int(velocity_points), class_counts=counts, sigma_mhz=sigma,
            tail_probability=math.erfc(6/math.sqrt(2)) if doppler else 0.0,
            convergence_status="not assessed; compare interval orders"),
        level_counts=dict(ground=len(atom.ground), excited=len(atom.excited),
                          retained_ground_density=atom.dimension),
        numerical_solver=dict(
            basis="complex density" if coordinates is None else "real Hermitian density",
            coherence_sector=("even delta-m" if coordinates is not None and even_coherences
                              else "all delta-m"),
            full_density_dimension=atom.dimension,
            solve_dimension=atom.dimension if coordinates is None else coordinates.dimension,
            method="batched local maps; normalized state/response share LU"),
        rate_convention="all _hz arguments are exponential rates in s^-1; no 2pi conversion",
        rotation_status="linear-input circular-birefringence approximation" if linear else "unsupported polarization",
    )
