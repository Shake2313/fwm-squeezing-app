"""Single-velocity PSR field response and ordered Langevin noise at finite RF.

The local pump is the stationary, driven Rb87 D1 atom from psr_full_obe.
Fields are photon-flux normalized and share a uniform transverse area A.
For mode dipole d0 and carrier wave number k, the reciprocal atom/field
coupling has g_A² = d0² k/(2 ε0 ℏ A). Both the drift and independent-atom
noise carry the same line density n A, so A cancels at fixed intensity.
"""
from dataclasses import replace
from functools import lru_cache

import numpy as np

from . import constants, core, species
from .psr_conventions import jones_to_circular
from .psr_full_obe import MODEL_ID as PUMP_MODEL_ID
from .quantum.contracts import AnalysisFrequencyAxis, GeneratorFrequencyAxis
from .quantum.diffusion import traceless_hermitian_basis
from .quantum.traveling import LocalNambuGenerator, constant_segment


PSR_NAMBU_MODES = ("x", "y", "x_dagger", "y_dagger")
PSR_RF_FRAME = "psr_pump_rotating_frame"
_SIGNS = (1., 1., -1., -1.)


@lru_cache(maxsize=1)
def _atomic_coordinates():
    atom = species.level_structure("Rb87", "D1")
    ng = len(atom.ground)
    dimension = ng + len(atom.excited)
    basis = traceless_hermitian_basis(dimension).reshape(dimension**2-1, -1).T
    lowering = []
    for polarization in ((1., 0.), (0., 1.)):
        circular = jones_to_circular(polarization)
        operator = np.zeros((dimension, dimension), complex)
        operator[:ng, ng:] = sum(circular[q].conjugate() * atom.dipoles[q] for q in (-1, 1))
        lowering.append(operator)
    operators = np.array([*lowering, *(op.conj().T for op in lowering)])
    return basis, operators


def _jump_gram(jumps, observables, rho):
    """Einstein diffusion projected before constructing its operator Gram matrix.

    Entry ij is the sum of Tr(rho [Ui,J]† [Uj,J]) over the supplied jumps.
    Intermediate storage scales with the four optical ports, independent of
    the number of traceless atomic coordinate pairs.
    """
    left = observables[None, ...] @ jumps[:, None, ...]
    commutators = left - jumps[:, None, ...] @ observables[None, ...]
    return np.einsum("kiab,kjab->ij", commutators.conj(), commutators @ rho)


def _replacement_distribution(jumps):
    """Recognize the declared rank-one replacement family and its actual rate.

    J_ab = sqrt(rate * p_a) |a><b| for every input b. Return None for
    any other family so the optimized path always has a general fallback.
    """
    occupied = jumps != 0
    if not np.all(occupied.sum(axis=(-2, -1)) == 1):
        return None
    rates = np.sum(np.abs(jumps)**2, axis=0)
    if not np.array_equal(rates, np.broadcast_to(rates[:, :1], rates.shape)):
        return None
    rate = float(rates[:, 0].sum())
    if rate <= 0:
        return None
    return rate, np.diag(rates[:, 0]/rate)


def _replacement_gram(replacement, observables, rho):
    """Exact jump Gram for a recognized replacement Lindbladian.

    The adjoint jump map is sum J_dagger A J = rate*Tr(tau A)*I.
    Expanding the two commutators removes the explicit jump dimension.
    """
    rate, target = replacement
    means_rho = np.einsum("iab,ba->i", observables, rho)
    means_target = np.einsum("iab,ba->i", observables, target)
    gram = np.einsum("iab,jab->ij", observables.conj(), observables @ (rho + target))
    return rate*(gram - np.outer(means_target.conj(), means_rho)
                 - np.outer(means_rho.conj(), means_target))


def psr_local_generator(pump_state, rf_hz, *, method="direct"):
    """Return M, D greater and D lesser per metre in the retarded pump frame.

    Supply the full dictionary returned by solve_psr_pump and a strictly
    increasing, signed RF vector in Hz. At f the ports are annihilators x,y
    at +f and creators x,y at -f. DC is the trace-zero static derivative.
    Reservoir terms are grouped by the jump kind declared in the S1 state.
    The default direct method preserves the S2 reference calculation. The
    optimized method contracts recognized replacement jumps algebraically
    and reuses the exact Nambu conjugation symmetry for paired signed RF.
    Both methods retain the same four-port quantum noise and named sources.
    """
    if pump_state.get("model_id") != PUMP_MODEL_ID:
        raise ValueError("A full Rb87 D1 PSR pump state is required")
    if method not in ("direct", "optimized"):
        raise ValueError("method must be direct or optimized")
    rf = AnalysisFrequencyAxis.from_hz(rf_hz)
    axis = GeneratorFrequencyAxis(rf.omega_rad_s, PSR_RF_FRAME)
    atom = species.level_structure("Rb87", "D1")
    density = float(pump_state["conditions"]["density_m3"])
    if not np.isfinite(density) or density < 0:
        raise ValueError("finite nonnegative isotope number density required")
    wave_number = 2*np.pi * atom.frequency_hz/constants.C_LIGHT
    line_coupling = density * atom.dipole_sq * wave_number/(2*constants.EPS_0*constants.HBAR)
    rho = pump_state["rho"]
    channels = pump_state["reservoirs"].channels
    groups = {}
    for channel in channels:
        groups.setdefault(channel.kind, []).append(channel.operator)
    jumps = [np.stack(operators) for operators in groups.values()]
    replacements = ([_replacement_distribution(group) for group in jumps]
                    if method == "optimized" else [None]*len(jumps))
    shape = (len(rf.omega_rad_s), 4, 4)
    field_drift = np.zeros(shape, complex)
    greater = np.zeros((len(jumps),) + shape, complex)
    lesser = np.zeros_like(greater)

    with core.blas_single_thread():
        basis, operators = _atomic_coordinates()
        # Drop the identity (global trace) exactly, including at DC.
        atomic_drift = basis.conj().T @ pump_state["generator_s_inverse"] @ basis
        adjoints = operators.conj().swapaxes(-1, -2)
        # H_int/ℏ = -g_A (b O† + b† O). Drive and Maxwell readout use the
        # same signed Cartesian lowering operators as the classical pump.
        forcing = 1j * (adjoints @ rho - rho @ adjoints)
        drive = basis.conj().T @ forcing.reshape(4, -1).T
        readout = (1j * np.asarray(_SIGNS)[:, None]
                   * (operators.swapaxes(-1, -2).reshape(4, -1) @ basis))
        eye = np.eye(atomic_drift.shape[0])
        completed = {}
        swap = np.array([2, 3, 0, 1])
        for index, omega in enumerate(rf.omega_rad_s):
            mirrored = completed.get(float(-omega)) if method == "optimized" else None
            if mirrored is not None:
                field_drift[index] = field_drift[mirrored].conj()[np.ix_(swap, swap)]
                greater[:, index] = lesser[:, mirrored].conj()[:, swap][:, :, swap]
                lesser[:, index] = greater[:, mirrored].conj()[:, swap][:, :, swap]
                continue
            # Four transpose solves give the rows C (-iΩ - L)^-1, without
            # forming the inverse or storing 255² noise entries per jump.
            resolvent_rows = np.linalg.solve((-1j*omega*eye - atomic_drift).T, readout.T).T
            field_drift[index] = line_coupling * (resolvent_rows @ drive)
            resolved = (resolvent_rows @ basis.T).reshape(4, *rho.shape)
            for source, group in enumerate(jumps):
                # For resolved port Ui, the greater ordering contracts
                # [J†, Ui] [Uj†, J]; the lesser contracts [J†, Uj†] [Ui, J].
                replacement = replacements[source]
                gram = ((lambda observable: _jump_gram(group, observable, rho))
                        if replacement is None else
                        (lambda observable: _replacement_gram(replacement, observable, rho)))
                greater[source, index] = line_coupling * gram(resolved.conj().swapaxes(-1, -2))
                lesser[source, index] = line_coupling * gram(resolved).T
            completed[float(omega)] = index
    local = LocalNambuGenerator(axis, _SIGNS, PSR_NAMBU_MODES, field_drift,
                                greater, lesser, tuple(groups))
    if not local.audit()["passed"]:
        raise ValueError("PSR microscopic response failed the commutator or ordered diffusion check")
    return local


def psr_uniform_transfer(local, pump_state, length_m, *, carrier_phase_rate_rad_per_m=None):
    """Integrate the linearized noise in a uniform, co-phased pump segment.

    Supply the same S1 pump used for the local generator, or the S3 Maxwell
    mean paired with the averaged local generator. The real,
    x-polarized B=0 carrier keeps its intensity and direction
    while accumulating the mean dispersive phase β z. Rotating the optical
    amplitudes into its local phase changes the drift to M - i β diag(signs).
    The noise is already evaluated in those local pump-phase coordinates.
    An explicit carrier_phase_rate_rad_per_m permits an arbitrary Jones pump
    in a prescribed common-phase gauge for one S4 propagation segment.
    """
    if (not isinstance(local, LocalNambuGenerator) or local.mode_labels != PSR_NAMBU_MODES
            or not np.array_equal(local.signs, _SIGNS) or local.frequency_axis.frame != PSR_RF_FRAME):
        raise ValueError("PSR Cartesian pump-frame generator required")
    if pump_state.get("model_id") not in (PUMP_MODEL_ID, "psr_rb87_d1_maxwell_full_obe_v1"):
        raise ValueError("A full Rb87 D1 PSR pump state or its Maxwell mean is required")
    jones = np.asarray(pump_state["conditions"]["jones_xy"])
    if carrier_phase_rate_rad_per_m is None and (abs(jones[1]) > 1e-10 or abs(jones[0].imag) > 1e-10):
        raise ValueError("PSR y-mode homodyne requires an x-polarized pump with real carrier amplitude")
    length = float(length_m)
    if not np.isfinite(length) or length < 0:
        raise ValueError("finite nonnegative cell length required")
    k = 2*np.pi*pump_state["conditions"]["optical_frequency_hz"]/constants.C_LIGHT
    beta = (k*pump_state["chi_probe"].real/2 if carrier_phase_rate_rad_per_m is None
            else float(carrier_phase_rate_rad_per_m))
    if not np.isfinite(beta):
        raise ValueError("finite carrier phase rate required")
    rotating = replace(local, drift=local.drift - 1j*beta*np.diag(local.signs))
    with core.blas_single_thread(), np.errstate(over="raise", invalid="raise"):
        # Power gain exp(2 Re(lambda) L) above 1/epsilon loses the
        # unit commutator to roundoff. Check the necessary eigenvalue
        # bound before exponentiation, then account for nonnormal gain.
        epsilon = np.finfo(float).eps
        if length * np.max(np.linalg.eigvals(rotating.drift).real) >= -0.5*np.log(epsilon):
            raise ValueError("field gain exceeds float64 resolution of one shot-noise unit")
        result = constant_segment(rotating, length)
        gain = np.max(np.linalg.norm(result.transfer, ord=2, axis=(-2, -1)))**2
        noise = max(np.max(np.abs(result.noise_greater)), np.max(np.abs(result.noise_lesser)))
        if epsilon*max(gain, noise) >= 1:
            raise ValueError("field transfer or added noise obscures one shot-noise unit at float64 precision")
        return result
