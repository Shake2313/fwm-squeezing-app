"""
Cluster D — 85Rb D1 double-Λ four-wave mixing.

The physics is ported verbatim from the original fwm_obe.py (Hamiltonians,
χ̄-matrix table, probe scan, Doppler average, Maxwell-Bloch propagation). The
only change is that shared engine pieces now come from gabes.core / .doppler /
.observables and the 4-level structure from gabes.atoms.

`compute_spectrum` / `operating_point` keep their original signatures so the CLI
shim and the regression test can call them directly. `FWMScheme` wraps them for
the generic Streamlit front-end.
"""
import math
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass

import numpy as np

from .. import atoms, constants, doppler, hyperfine, kernels, observables, species
from ..constants import K_VEC, OMEGA_HF, OMEGA_EXCITED_HF, rabi_freq
from ..core import blas_single_thread, build_liouvillian, comm_super, floquet_solve
from ..lineshape import fwhm_interp
from .base import ExtraView, ParamSpec, Preset, Scheme

# =========================================================
# Level indices (this scheme's labelling of the atom model)
# =========================================================
ATOM = atoms.get("double_lambda_rb85")
G1, G2, E2, E3 = 0, 1, 2, 3
GROUND_STATES = (G1, G2)
EXCITED_STATES = (E2, E3)
N_LEVELS = ATOM.n_levels
GROUND_F = {G1: 2, G2: 3}
EXCITED_F = {E2: 2, E3: 3}
TRANSITION_DIPOLE_SCALE = np.zeros((N_LEVELS, N_LEVELS), dtype=float)
for _g in GROUND_STATES:
    for _e in EXCITED_STATES:
        TRANSITION_DIPOLE_SCALE[_g, _e] = np.sqrt(
            3.0 * hyperfine.CF2[(GROUND_F[_g], EXCITED_F[_e])])


def physical_coupling_norm(branch):
    """First-principles macroscopic-coupling normalization for the lumped
    4-level double-Λ model — the factor the old hand-tuned ``line_strength=0.05``
    was standing in for.

    The AutoOD-validated absorption scale (`schemes.absorption._hyperfine_alpha`,
    <0.1 % vs lab) carries an explicit ``p_F / [2(2I+1)]``: the ground-population
    fraction times the ground-sublevel degeneracy. The lumped FWM solve instead
    normalizes Tr(ρ)=1 over just 4 levels and multiplies by the *total* atomic
    density, so this population/degeneracy factor is otherwise omitted. The χ̄
    already carries the Clebsch-Gordan strengths (C² = 3·C_F² via drive×readout),
    so this is the only structural piece missing relative to the validated path.

    The seeded probe couples from F=3 on the (−) Raman branch and F=2 on (+);
    the residual asymmetry between the two ground manifolds folds into the
    dimensionless ``line_strength`` residual knob.
    """
    probe_F = GROUND_F[G2] if branch == -1 else GROUND_F[G1]
    return hyperfine.GROUND_POP[probe_F] / hyperfine.N_GROUND_SUBLEVELS

# =========================================================
# FWM experiment configuration (cell, beams, detection, scan)
# =========================================================
L_CELL = 12.5e-3
W_PUMP = 530e-6
W_PROBE = 330e-6
P_PUMP, P_PROBE = 600e-3, 10e-6
T_CELL = 394.15

# Detection efficiency matched to Sim et al. (Sci. Rep. 15, 7727 (2025)): the
# reported system (detection) loss is 8.0% -> QE 0.92; the 5.5% optical loss is
# the separate `loss_pct` knob. Together QE·(1−0.055)=0.869 reproduces the paper's
# ~13.5% total detection loss. (RESPONSIVITY_AW is display-only: 0.92·795/1240.)
QE_DETECTOR = 0.92
RESPONSIVITY_AW = 0.59
LOSS_FRAC = 0.0
ETA_TOTAL = QE_DETECTOR * (1.0 - LOSS_FRAC)

OMEGA_C_SEED = 0.0
SCAN_MIN_GHZ = -8.0
SCAN_MAX_GHZ = 12.0
SCAN_COARSE_POINTS = 401
RESONANCE_WINDOW_MHZ = 80.0
SCAN_FINE_POINTS = 801
PUMP_OVERLAP_EXCLUSION_MHZ = 1e-3
VELOCITY_STEP_MPS = 1.0
VELOCITY_CUTOFF_SIGMA = 3.0
DELTA_GHZ_LIST = [0.9]
BRANCHES = (-1, +1)
DEFAULT_BRANCH = -1
PHASE_LEGACY = "legacy"
PHASE_BALANCED = "balanced"
PHASE_FINE = "fine"
PHASE_ULTRA = "ultra"
SEEDED_PHASE_ANGLE_DEG = 0.32
ULTRA_PHASE_ITERATIONS = 1   # n(χ) has no Δk feedback; one refractive pass is exact
ULTRA_PROPAGATION_SEGMENTS = 64


# =========================================================
# Generic SFWM / biphoton topology layer
# =========================================================
NM = 1e-9
MHZ_ANG = 2 * np.pi * 1e6

MODE_SEEDED = "Squeezing"
MODE_BIPHOTON = "Biphoton"

# Biphoton source model. "Predictive" solves the Doppler-averaged cascade/double-Λ
# biphoton amplitude from first principles (Chen et al. PRR 4, 023132 (2024)
# Eq. (3-5); Kim et al. QST 9, 045006 (2024) Eq. (2); Du, Wen, Rubin JOSAB 25,
# C98 (2008)) — waveform, FWHM, bandwidth, OD reshaping and the rate scaling are
# computed, with only one residual scalar per topology setting the absolute rate
# (the squeezing-mode `line_strength` philosophy). "Calibrated" is the legacy
# reference-injected estimate kept for comparison.
BIPHOTON_PREDICTIVE = "Predictive (first-principles)"
BIPHOTON_CALIBRATED = "Calibrated (reference)"
BIPHOTON_MODELS = (BIPHOTON_PREDICTIVE, BIPHOTON_CALIBRATED)

# Two-photon (ground) coherence dephasing as a fraction of the intermediate Γ;
# sets the EIT/Raman two-photon linewidth. Chen et al. fit γ ≈ 0.02–0.03 Γ.
GROUND_DEPHASING_FRAC = 0.02
# Regularization clip for the complex longitudinal function ρ̄ at high OD
# (mirrors the dilute-vapor clip in `_safe_refractive_index`).
PRED_RHO_CLIP = 60.0

# Predictive velocity-grid auto-refinement. The nonlinear source |amp(v)| that the
# velocity-class coherent sum integrates is only ~Γ/k — a few m/s — wide, far
# narrower than the Doppler width σ_v the navigate-only `biphoton_velocity_step`
# is sized for. A step coarser than that resonance aliases the (Fourier) sum, so
# the reported absolute BTW width tracks numerical undersampling, not physics
# (e.g. a factor-~20, non-monotonic swing over vstep=1–12 m/s in the 780/1529 nm
# telecom cascade). The predictive path therefore starts from a step that
# oversamples the *measured* resonance, then halves until the |ψ|² FWHM is stable;
# if the point cap is hit first the width is flagged unconverged.
PRED_V_OVERSAMPLE = 16.0      # starting step = probe-measured resonance FWHM / this
PRED_V_FWHM_TOL = 0.03        # relative |ψ|² FWHM change accepted as converged
PRED_V_MAX_REFINE = 10        # max successive halvings of the velocity step
PRED_V_MAX_POINTS = 40000     # velocity-grid point cap (guards runtime)

TOPOLOGY_RB87_TELECOM = "cascade_rb87_telecom"
TOPOLOGY_CS_BTW = "cascade_cs_btw"
TOPOLOGY_DIAMOND = "diamond_generic"
CS_CHANNEL_917 = "6D5/2: 852-917 nm"
CS_CHANNEL_795 = "8S1/2: 852-795 nm"
SIDE_PLUS = "+"
SIDE_MINUS = "-"
SIDE_CHOICES = (SIDE_PLUS, SIDE_MINUS)


@dataclass(frozen=True)
class LevelSpec:
    """Lightweight FWM level metadata used by the generic topology layer."""
    name: str
    energy_hz: float
    gamma_mhz: float = 0.0


@dataclass(frozen=True)
class FieldSpec:
    """A driven or generated optical field in a four-wave-mixing topology."""
    role: str
    lower: int
    upper: int
    wavelength_nm: float
    detuning_mhz: float = 0.0
    rabi_mhz: float = 0.0
    phase_sign: float = 1.0
    direction: float = 1.0
    angle_deg: float = 0.0
    side_sign: float = 0.0

    @property
    def k(self):
        return 2 * np.pi / (self.wavelength_nm * NM)

    @property
    def frequency_hz(self):
        return constants.C_LIGHT / (self.wavelength_nm * NM)


@dataclass(frozen=True)
class TopologySpec:
    """Generic SFWM topology; presets carry the reference-calibrated constants."""
    name: str
    label: str
    family: str
    isotope_name: str
    levels: tuple
    fields: tuple
    signal_role: str
    idler_role: str
    default_temp_c: float
    default_cell_mm: float
    default_pump_uw: float
    default_coupling_mw: float
    pair_rate_cps_per_mw: float
    emission_decay_ns: float
    target_g2_peak: float | None = None
    reference_fwhm_ns: float | None = None
    reference_od: float | None = None
    reference_bandwidth_mhz: float | None = None
    reference_width_ratio: float | None = None
    reference_delta_k: float | None = None
    notes: str = ""

    @property
    def isotope(self):
        return species.ISOTOPES[self.isotope_name]

    @property
    def field_map(self):
        return {f.role: f for f in self.fields}


def _wavevector_nm(wavelength_nm):
    return 2 * np.pi / (float(wavelength_nm) * NM)


def _side_sign(value):
    if isinstance(value, str):
        return 1.0 if value.strip() == SIDE_PLUS else -1.0
    return 1.0 if float(value) >= 0.0 else -1.0


def _side_label(value):
    return SIDE_PLUS if float(value) >= 0.0 else SIDE_MINUS


def transverse_matched_angle_deg(source_wavelength_nm, target_wavelength_nm,
                                 source_angle_deg):
    """Collection angle that cancels transverse k for two generated photons."""
    source_k = _wavevector_nm(source_wavelength_nm)
    target_k = _wavevector_nm(target_wavelength_nm)
    x = source_k / target_k * math.sin(math.radians(float(source_angle_deg)))
    return math.degrees(math.asin(float(np.clip(x, -1.0, 1.0))))


def _field_with(field, *, wavelength_nm=None, detuning_mhz=None, rabi_mhz=None,
                angle_deg=None, side_sign=None):
    return FieldSpec(
        role=field.role,
        lower=field.lower,
        upper=field.upper,
        wavelength_nm=field.wavelength_nm if wavelength_nm is None else wavelength_nm,
        detuning_mhz=field.detuning_mhz if detuning_mhz is None else detuning_mhz,
        rabi_mhz=field.rabi_mhz if rabi_mhz is None else rabi_mhz,
        phase_sign=field.phase_sign,
        direction=field.direction,
        angle_deg=field.angle_deg if angle_deg is None else angle_deg,
        side_sign=field.side_sign if side_sign is None else side_sign,
    )


def phase_mismatch(fields, *, signal_angle_deg=None, idler_angle_deg=None,
                   reference_delta_k=0.0):
    """Longitudinal four-field phase mismatch, with an optional reference offset."""
    total = 0.0
    for field in fields:
        angle = field.angle_deg
        if field.role == "signal" and signal_angle_deg is not None:
            angle = signal_angle_deg
        if field.role == "idler" and idler_angle_deg is not None:
            angle = idler_angle_deg
        total += field.phase_sign * field.k * math.cos(math.radians(angle))
    return total - (reference_delta_k or 0.0)


def phase_mismatch_vector(fields, *, signal_angle_deg=None, idler_angle_deg=None,
                          signal_side_sign=None, idler_side_sign=None,
                          reference_delta_k=0.0):
    """Biphoton vector mismatch: calibrated longitudinal, absolute transverse."""
    delta_k_z_absolute = 0.0
    delta_k_x = 0.0
    for field in fields:
        angle = field.angle_deg
        side = field.side_sign
        if field.role == "signal":
            if signal_angle_deg is not None:
                angle = signal_angle_deg
            if signal_side_sign is not None:
                side = signal_side_sign
        elif field.role == "idler":
            if idler_angle_deg is not None:
                angle = idler_angle_deg
            if idler_side_sign is not None:
                side = idler_side_sign
        angle_rad = math.radians(float(angle))
        delta_k_z_absolute += field.phase_sign * field.k * math.cos(angle_rad)
        delta_k_x += field.phase_sign * field.k * math.sin(angle_rad) * side
    delta_k_z_relative = delta_k_z_absolute - (reference_delta_k or 0.0)
    delta_k_vector = math.hypot(delta_k_z_relative, delta_k_x)
    return {
        "delta_k_z_relative": delta_k_z_relative,
        "delta_k_z_absolute": delta_k_z_absolute,
        "delta_k_x": delta_k_x,
        "delta_k_vector": delta_k_vector,
    }


def phase_matching_weight(delta_k, L):
    """sinc^2(delta_k L / 2), normalized to 1 at perfect phase matching."""
    x = 0.5 * np.asarray(delta_k, dtype=float) * L
    out = np.ones_like(x, dtype=float)
    mask = np.abs(x) > 1e-12
    out[mask] = (np.sin(x[mask]) / x[mask]) ** 2
    return out


def _sinc_complex(x):
    """sinc(x) = sin(x)/x for complex argument (→1 at x→0). The longitudinal
    detuning function Φ = sinc(ρ̄)·e^{iρ̄} (Du et al. Eq. 14, Chen et al. Eq. 3)
    carries a complex ρ̄ when in-cell loss/dispersion (OD) is included."""
    x = np.asarray(x, dtype=complex)
    out = np.ones_like(x)
    mask = np.abs(x) > 1e-9
    out[mask] = np.sin(x[mask]) / x[mask]
    return out


def _bandwidth_from_waveform_mhz(tau_s, psi):
    """Spectral FWHM [MHz] of a biphoton temporal waveform ψ(τ) via its FFT."""
    n = np.asarray(tau_s).size
    if n < 4:
        return float("nan")
    dt = float(tau_s[1] - tau_s[0])
    nfft = 4 * n
    spec = np.fft.fftshift(np.fft.fft(np.asarray(psi), n=nfft))
    freq = np.fft.fftshift(np.fft.fftfreq(nfft, dt))     # Hz
    power = np.abs(spec)**2
    if power.max() <= 0:
        return float("nan")
    above = np.where(power >= 0.5 * power.max())[0]
    if above.size < 2:
        return float("nan")
    return float((freq[above[-1]] - freq[above[0]]) / 1e6)


def phase_mismatch_grid(fields, signal_axis_deg, idler_axis_deg,
                        reference_delta_k=0.0):
    """Vectorized signal/idler longitudinal mismatch grid."""
    signal_axis_deg = np.asarray(signal_axis_deg, dtype=float)
    idler_axis_deg = np.asarray(idler_axis_deg, dtype=float)
    sig, ide = np.meshgrid(signal_axis_deg, idler_axis_deg, indexing="ij")
    total = np.zeros_like(sig, dtype=float)
    for field in fields:
        if field.role == "signal":
            total += field.phase_sign * field.k * np.cos(np.deg2rad(sig))
        elif field.role == "idler":
            total += field.phase_sign * field.k * np.cos(np.deg2rad(ide))
        else:
            total += field.phase_sign * field.k * math.cos(math.radians(field.angle_deg))
    return total - (reference_delta_k or 0.0)


def phase_mismatch_vector_grid(fields, signal_axis_deg, idler_axis_deg,
                               reference_delta_k=0.0):
    """Vectorized signal/idler mismatch magnitude grid for biphoton collection."""
    signal_axis_deg = np.asarray(signal_axis_deg, dtype=float)
    idler_axis_deg = np.asarray(idler_axis_deg, dtype=float)
    sig, ide = np.meshgrid(signal_axis_deg, idler_axis_deg, indexing="ij")
    delta_k_z_absolute = np.zeros_like(sig, dtype=float)
    delta_k_x = np.zeros_like(sig, dtype=float)
    for field in fields:
        if field.role == "signal":
            angle_rad = np.deg2rad(sig)
            delta_k_z_absolute += field.phase_sign * field.k * np.cos(angle_rad)
            delta_k_x += (field.phase_sign * field.k * np.sin(angle_rad)
                          * field.side_sign)
        elif field.role == "idler":
            angle_rad = np.deg2rad(ide)
            delta_k_z_absolute += field.phase_sign * field.k * np.cos(angle_rad)
            delta_k_x += (field.phase_sign * field.k * np.sin(angle_rad)
                          * field.side_sign)
        else:
            angle_rad = math.radians(field.angle_deg)
            delta_k_z_absolute += field.phase_sign * field.k * math.cos(angle_rad)
            delta_k_x += (field.phase_sign * field.k * math.sin(angle_rad)
                          * field.side_sign)
    delta_k_z_relative = delta_k_z_absolute - (reference_delta_k or 0.0)
    return np.hypot(delta_k_z_relative, delta_k_x)


def energy_mismatch_hz(fields):
    signs = {"pump": 1.0, "coupling": 1.0, "signal": -1.0, "idler": -1.0}
    return sum(signs.get(f.role, 0.0) * f.frequency_hz for f in fields)


def _raw_delta_k(fields):
    return phase_mismatch(fields, reference_delta_k=0.0)


def _with_reference_delta_k(spec):
    return TopologySpec(
        name=spec.name,
        label=spec.label,
        family=spec.family,
        isotope_name=spec.isotope_name,
        levels=spec.levels,
        fields=spec.fields,
        signal_role=spec.signal_role,
        idler_role=spec.idler_role,
        default_temp_c=spec.default_temp_c,
        default_cell_mm=spec.default_cell_mm,
        default_pump_uw=spec.default_pump_uw,
        default_coupling_mw=spec.default_coupling_mw,
        pair_rate_cps_per_mw=spec.pair_rate_cps_per_mw,
        emission_decay_ns=spec.emission_decay_ns,
        target_g2_peak=spec.target_g2_peak,
        reference_fwhm_ns=spec.reference_fwhm_ns,
        reference_od=spec.reference_od,
        reference_bandwidth_mhz=spec.reference_bandwidth_mhz,
        reference_width_ratio=spec.reference_width_ratio,
        reference_delta_k=_raw_delta_k(spec.fields),
        notes=spec.notes,
    )


def _rb87_telecom_spec():
    levels = (
        LevelSpec("5S1/2(F=2)", 0.0, 0.0),
        LevelSpec("5P3/2", constants.C_LIGHT / (780.24 * NM), 6.07),
        LevelSpec("4D5/2", constants.C_LIGHT / (780.24 * NM)
                  + constants.C_LIGHT / (1529.37 * NM), 0.66),
        LevelSpec("5P3/2 collection", constants.C_LIGHT / (780.24 * NM), 6.07),
    )
    fields = (
        FieldSpec("pump", 0, 1, 780.24, phase_sign=+1.0),
        FieldSpec("coupling", 1, 2, 1529.37, phase_sign=-1.0, direction=-1.0),
        FieldSpec("signal", 2, 3, 1529.37, phase_sign=-1.0, angle_deg=1.5,
                  side_sign=+1.0),
        FieldSpec("idler", 3, 0, 780.24, phase_sign=+1.0,
                  angle_deg=transverse_matched_angle_deg(1529.37, 780.24, 1.5),
                  side_sign=+1.0),
    )
    return _with_reference_delta_k(TopologySpec(
        name=TOPOLOGY_RB87_TELECOM,
        label="87Rb cascade telecom (5S-5P-4D)",
        family="cascade",
        isotope_name="Rb87",
        levels=levels,
        fields=fields,
        signal_role="signal",
        idler_role="idler",
        default_temp_c=90.0,
        default_cell_mm=12.5,
        default_pump_uw=10.0,
        default_coupling_mw=1.0,
        pair_rate_cps_per_mw=38_000.0,
        emission_decay_ns=0.52,
        target_g2_peak=44.0,
        reference_fwhm_ns=0.56,
        reference_od=112.0,
        reference_bandwidth_mhz=300.0,
        notes=("Reference-calibrated cascade SFWM estimate for the telecom "
               "biphoton source in hot 87Rb."),
    ))


def _cs_btw_spec(channel):
    if channel == CS_CHANNEL_795:
        upper = "8S1/2"
        coupling_nm = 795.0
        upper_gamma_mhz = 1.7
        decay_ns = 1.35
        label = "133Cs cascade BTW (852-795 nm)"
    else:
        upper = "6D5/2"
        coupling_nm = 917.0
        upper_gamma_mhz = 2.6
        decay_ns = 4.1
        label = "133Cs cascade BTW (852-917 nm)"
    pump_nm = 852.35
    levels = (
        LevelSpec("6S1/2(F=4)", 0.0, 0.0),
        LevelSpec("6P3/2(F'=5)", constants.C_LIGHT / (pump_nm * NM), 5.23),
        LevelSpec(upper, constants.C_LIGHT / (pump_nm * NM)
                  + constants.C_LIGHT / (coupling_nm * NM), upper_gamma_mhz),
        LevelSpec("6P3/2 collection", constants.C_LIGHT / (pump_nm * NM), 5.23),
    )
    fields = (
        FieldSpec("pump", 0, 1, pump_nm, phase_sign=+1.0),
        FieldSpec("coupling", 1, 2, coupling_nm, phase_sign=-1.0, direction=-1.0),
        FieldSpec("signal", 2, 3, coupling_nm, phase_sign=-1.0, angle_deg=1.5,
                  side_sign=+1.0),
        FieldSpec("idler", 3, 0, pump_nm, phase_sign=+1.0,
                  angle_deg=transverse_matched_angle_deg(coupling_nm, pump_nm, 1.5),
                  side_sign=+1.0),
    )
    return _with_reference_delta_k(TopologySpec(
        name=TOPOLOGY_CS_BTW,
        label=label,
        family="cascade",
        isotope_name="Cs133",
        levels=levels,
        fields=fields,
        signal_role="signal",
        idler_role="idler",
        default_temp_c=75.0,
        default_cell_mm=12.5,
        default_pump_uw=20.0,
        default_coupling_mw=1.0,
        pair_rate_cps_per_mw=12_000.0,
        emission_decay_ns=decay_ns,
        target_g2_peak=18.0,
        reference_fwhm_ns=decay_ns * 0.42,
        reference_od=10.0,
        reference_width_ratio=3.0,
        notes=("Velocity-class coherent-sum model for the Cs biphoton temporal "
               "waveform comparison."),
    ))


def _diamond_generic_spec(params=None):
    params = params or {}
    pump_nm = float(params.get("diamond_pump_nm", 780.0))
    coupling_nm = float(params.get("diamond_coupling_nm", 776.0))
    signal_nm = float(params.get("diamond_signal_nm", 795.0))
    # Energy conservation: 1/λ_idler = 1/λ_pump + 1/λ_coupling − 1/λ_signal.
    # Guard the reciprocal sum — a zero/negative denominator (e.g. pump=coupling
    # with signal=pump/2) has no physical idler and would otherwise raise
    # ZeroDivisionError. Fall back to NaN; the live UI always supplies an explicit
    # diamond_idler_nm, so this default is only the degenerate-case fallback.
    _inv_idler = 1.0 / pump_nm + 1.0 / coupling_nm - 1.0 / signal_nm
    idler_default = 1.0 / _inv_idler if abs(_inv_idler) > 1e-9 else float("nan")
    idler_nm = float(params.get("diamond_idler_nm", idler_default))
    levels = (
        LevelSpec("g", 0.0, 0.0),
        LevelSpec("e1", constants.C_LIGHT / (pump_nm * NM), 6.0),
        LevelSpec("e2", constants.C_LIGHT / (coupling_nm * NM), 6.0),
        LevelSpec("u", constants.C_LIGHT / (pump_nm * NM)
                  + constants.C_LIGHT / (coupling_nm * NM), 1.0),
    )
    fields = (
        FieldSpec("pump", 0, 1, pump_nm, phase_sign=+1.0),
        FieldSpec("coupling", 0, 2, coupling_nm, phase_sign=+1.0),
        FieldSpec("signal", 3, 1, signal_nm, phase_sign=-1.0, angle_deg=2.0,
                  side_sign=+1.0),
        FieldSpec("idler", 3, 2, idler_nm, phase_sign=-1.0,
                  angle_deg=transverse_matched_angle_deg(signal_nm, idler_nm, 2.0),
                  side_sign=-1.0),
    )
    return _with_reference_delta_k(TopologySpec(
        name=TOPOLOGY_DIAMOND,
        label="Generic diamond four-level SFWM",
        family="diamond",
        isotope_name="Rb87",
        levels=levels,
        fields=fields,
        signal_role="signal",
        idler_role="idler",
        default_temp_c=60.0,
        default_cell_mm=12.5,
        default_pump_uw=20.0,
        default_coupling_mw=1.0,
        pair_rate_cps_per_mw=5_000.0,
        emission_decay_ns=8.0,
        target_g2_peak=None,
        reference_fwhm_ns=None,
        reference_od=None,
        reference_bandwidth_mhz=None,
        reference_width_ratio=None,
        notes=("Generic template only; not tied to a validated diamond reference "
               "preset."),
    ))


def topology_from_params(params):
    topo = params.get("topology", TOPOLOGY_RB87_TELECOM)
    if topo == TOPOLOGY_CS_BTW:
        return _cs_btw_spec(params.get("cs_channel", CS_CHANNEL_917))
    if topo == TOPOLOGY_DIAMOND:
        return _diamond_generic_spec(params)
    return _rb87_telecom_spec()


def _default_biphoton_geometry(params):
    spec = topology_from_params(params)
    signal = spec.field_map["signal"]
    idler = spec.field_map["idler"]
    return dict(
        signal_angle_deg=signal.angle_deg,
        idler_angle_deg=idler.angle_deg,
        signal_side=_side_label(signal.side_sign),
        idler_side=_side_label(idler.side_sign),
    )


class GenericFWMSolver:
    """Reference-calibrated v3 engine for generic SFWM biphoton estimates."""

    def __init__(self, topology):
        self.topology = topology

    def _fields_from_params(self, params):
        pump_rabi = params.get("pump_biphoton_uw", self.topology.default_pump_uw)
        coupling_rabi = params.get("coupling_mw", self.topology.default_coupling_mw)
        out = []
        for field in self.topology.fields:
            rabi = field.rabi_mhz
            detuning = field.detuning_mhz
            if field.role == "pump":
                rabi = math.sqrt(max(pump_rabi, 0.0) / max(self.topology.default_pump_uw, 1e-12))
                rabi *= 2.0
                detuning = params.get("pump_detuning_mhz", 0.0)
            elif field.role == "coupling":
                rabi = math.sqrt(max(coupling_rabi, 0.0) / max(self.topology.default_coupling_mw, 1e-12))
                rabi *= 12.0
                detuning = params.get("coupling_detuning_mhz", 0.0)
            out.append(_field_with(
                field,
                detuning_mhz=detuning,
                rabi_mhz=rabi,
                angle_deg=params.get(f"{field.role}_angle_deg", field.angle_deg),
                side_sign=_side_sign(params[f"{field.role}_side"])
                if field.role in ("signal", "idler")
                and f"{field.role}_side" in params else field.side_sign,
            ))
        return tuple(out)

    def _leg_optical_depth(self, density, L, gamma_mhz, wavelength_nm):
        """Physical on-resonance optical depth of one cascade leg from its natural
        linewidth, via the AutoOD-validated Γ→|d|² route (`species.reduced_dipole_sq`):
        α₀ = 2 N k |d|² /(ε₀ ℏ Γ), OD = α₀ L."""
        gamma_nat = max(float(gamma_mhz), 1e-3) * MHZ_ANG
        lam = float(wavelength_nm) * NM
        k = 2.0 * np.pi / lam
        d2 = species.reduced_dipole_sq(gamma_nat, lam, 0.5, 0.5)
        alpha0 = 2.0 * density * k * d2 / (constants.EPS_0 * constants.HBAR * gamma_nat)
        return float(alpha0 * L)

    def _apply_longitudinal_response(self, kappa_tau, tau_s, v_grid, weights,
                                     residual_k, two_det, Oc, Gamma_e, gamma_g,
                                     od_phys):
        """Convolve the nonlinear response κ̃(τ) with the linear longitudinal
        function Φ̃(τ) (Du Eq. 15), done as a product in the conjugate domain:
        ψ = FFT⁻¹[FFT(κ̃)·Φ(δ)], with Φ(δ)=sinc(ρ̄)·e^{iρ̄} (Du Eq. 14) and ρ̄(δ)
        the OD-weighted EIT/slow-light phase (Chen Eq. 5). At OD→0, ρ̄→0, Φ→1 and
        ψ→κ̃ exactly (no reshaping)."""
        if od_phys <= 0:
            return np.asarray(kappa_tau)
        n = tau_s.size
        dt = float(tau_s[1] - tau_s[0])
        nfft = 4 * n
        K = np.fft.fft(np.asarray(kappa_tau), n=nfft)
        omega = 2.0 * np.pi * np.fft.fftfreq(nfft, dt)        # δ axis [rad/s]
        d = omega[:, None]
        v = v_grid[None, :]
        w = weights[None, :]
        twoden = Oc**2 - 4.0 * (d + 1j * gamma_g) * (
            d + two_det - residual_k * v + 0.5j * Gamma_e)
        rho = (0.5 * od_phys * Gamma_e) * ((d + 1j * gamma_g) / twoden * w).sum(axis=1)
        rho = (np.clip(rho.real, -PRED_RHO_CLIP, PRED_RHO_CLIP)
               + 1j * np.clip(rho.imag, 0.0, PRED_RHO_CLIP))
        phi = _sinc_complex(rho) * np.exp(1j * rho)
        return np.fft.ifft(K * phi, n=nfft)[:n]

    def _predictive_waveform(self, params, fields, v_grid, weights, pump, coupling,
                             signal, idler, residual_k, density, L, pm_weight,
                             tau_axis_ns, od_value):
        """First-principles Doppler-averaged biphoton amplitude.

        Frequency-domain form of Chen et al. (Phys. Rev. Research 4, 023132 (2024))
        Eq. (3-5), equivalent to Kim et al. (Quantum Sci. Technol. 9, 045006 (2024))
        Eq. (2) and Du, Wen, Rubin (J. Opt. Soc. Am. B 25, C98 (2008)) Eq. (13-18):
        for each signal detuning δ the Maxwell velocity classes are coherently
        summed into the nonlinear coupling κ̃(δ) and the linear longitudinal
        function ρ̄(δ); the joint amplitude A(δ)=κ̃·sinc(ρ̄)·e^{iρ̄} is
        inverse-Fourier-transformed to relative time τ. The two-photon denominator
        carries the Ω_c² Autler-Townes term (no weak-coupling approximation), the
        decay envelope emerges from the transform (no injected lifetime), and the
        optical depth α enters κ̃/ρ̄ (group-delay / Sommerfeld-precursor reshaping).
        """
        Gamma_i = max(self.topology.levels[1].gamma_mhz, 0.1) * MHZ_ANG   # intermediate (idler leg)
        Gamma_e = max(self.topology.levels[2].gamma_mhz, 0.1) * MHZ_ANG   # excited Γ₃
        gamma_g = GROUND_DEPHASING_FRAC * Gamma_i                         # two-photon (ground) dephasing
        Op = pump.rabi_mhz * MHZ_ANG
        Oc = max(coupling.rabi_mhz, 1e-6) * MHZ_ANG
        dp = params.get("pump_detuning_mhz", 0.0) * MHZ_ANG
        two_det = (params.get("pump_detuning_mhz", 0.0)
                   + params.get("coupling_detuning_mhz", 0.0)) * MHZ_ANG

        # optical depth seen by the near-resonant idler leg → Chen's α. OD is a
        # measured quantity in these sources (like cell temperature), so use the
        # reference-anchored, density/L-scaled value where available; the in-cell
        # reshaping (ρ̄) it drives is still computed from first principles.
        od_phys = float(od_value)

        tau_s = tau_axis_ns * 1e-9

        # ---- Nonlinear response κ̃(τ): time-domain Kim et al. Eq. (2) ----
        # Per velocity class amp(v) = Ω_p Ω_c / [4·f₁·f₂ + Ω_c²] with the Ω_c²
        # Autler-Townes term in the two-photon denominator (no weak-coupling
        # approximation). The collective two-photon coherence is the coherent sum
        # over Maxwell velocity classes carrying the single-photon phase
        # e^{i k_P v τ}, ×natural decay e^{−Γ τ/2}, ×H(τ) (τ≥0 implicit). The
        # velocity-sum dephasing — not an injected lifetime — sets the BTW width.
        f1 = 0.5 * Gamma_i + 1j * (dp - pump.direction * pump.k * v_grid)
        f2 = 0.5 * Gamma_e + 1j * (two_det - residual_k * v_grid)
        amp_v = weights * (Op * Oc) / (4.0 * f1 * f2 + Oc**2)
        amp_v = amp_v * od_phys                          # κ ∝ α (OD) → rate ∝ OD²
        coherent = np.exp(1j * pump.k * v_grid[:, None] * tau_s[None, :])
        kappa_tau = (amp_v[:, None] * coherent).sum(axis=0)
        kappa_tau = kappa_tau * np.exp(-0.5 * Gamma_i * tau_s)

        # ---- Linear longitudinal response Φ̃(τ): Du Eq. (15) convolution ----
        # ρ̄(δ) (Chen Eq. 5) is the OD-weighted EIT / slow-light phase; the
        # longitudinal function Φ(δ)=sinc(ρ̄)·e^{iρ̄} (Du Eq. 14) reshapes the
        # waveform (group delay / Sommerfeld precursor) at high OD. Convolved with
        # κ̃(τ); at low OD ρ̄→0, Φ→1, ψ→κ̃ (no reshaping). OFF by default: the
        # lumped 4-level model overestimates the high-OD group-delay broadening
        # (it would smear the validated narrow telecom BTW), so the reshaping is a
        # diagnostic opt-in (`biphoton_od_reshaping`) rather than the default path.
        if params.get("biphoton_od_reshaping", False):
            psi_tau = self._apply_longitudinal_response(
                kappa_tau, tau_s, v_grid, weights, residual_k, two_det,
                Oc, Gamma_e, gamma_g, od_phys)
        else:
            psi_tau = kappa_tau
        psi_tau = psi_tau * math.sqrt(max(pm_weight, 0.0))   # transverse collection

        # Source spectral width from the waveform itself (predictive).
        bandwidth_mhz = _bandwidth_from_waveform_mhz(tau_s, psi_tau)

        # Du regime split: group-delay time τ_g≈(2γ/Ω_c²)·OD·Γ vs Rabi time 2π/Ω_c.
        tau_group = (2.0 * gamma_g / max(Oc**2, 1e-30)) * od_phys * Gamma_e
        tau_rabi = 2.0 * np.pi / max(Oc, 1e-30)
        regime = "group-delay" if tau_group > tau_rabi else "damped-Rabi"

        source_v = amp_v   # velocity-class source for the existing diagnostic plot

        return {
            "psi_tau": psi_tau,
            "source_v": source_v,
            "od_phys": float(od_phys),
            "bandwidth_mhz": float(bandwidth_mhz),
            "regime": regime,
            # Chen et al. ultimate spectral-brightness ceiling ≈ (π/2)·10⁶ pairs/s/MHz
            "brightness_limit_cps_per_mhz": float(0.5 * np.pi * 1e6),
        }

    def _predictive_velocity_step(self, params, pump, coupling, residual_k, T, iso):
        """Velocity-grid step that oversamples the velocity-space resonance of the
        nonlinear source |amp(v)| — the integrand of the velocity-class coherent
        sum in `_predictive_waveform`.

        The resonance is measured directly on a fine probe grid (so it follows the
        Autler-Townes / detuning broadening of the denominator, not just the bare
        linewidth), and the step is set to its FWHM / `PRED_V_OVERSAMPLE`. Falls
        back to the natural-linewidth width Γ/2k if the probe feature is
        ill-defined. This is only the *starting* step; `_converged_predictive_
        waveform` halves it further until the waveform width converges.
        """
        sigma_v = math.sqrt(constants.KB * T / iso.mass)
        Gamma_i = max(self.topology.levels[1].gamma_mhz, 0.1) * MHZ_ANG
        Gamma_e = max(self.topology.levels[2].gamma_mhz, 0.1) * MHZ_ANG
        Oc = max(coupling.rabi_mhz, 1e-6) * MHZ_ANG
        Op = pump.rabi_mhz * MHZ_ANG
        dp = params.get("pump_detuning_mhz", 0.0) * MHZ_ANG
        two_det = (params.get("pump_detuning_mhz", 0.0)
                   + params.get("coupling_detuning_mhz", 0.0)) * MHZ_ANG
        kp = max(abs(pump.k), 1e-30)
        rk = max(abs(residual_k), 1e-30)
        # probe spacing resolves the narrowest bare-linewidth velocity scale
        narrow = min(Gamma_i / (2.0 * kp), Gamma_e / (2.0 * rk))
        dv_probe = max(narrow / 6.0, 1e-3)
        n_probe = int(min(2.0 * 3.2 * sigma_v / dv_probe + 1.0, 200000))
        v = np.linspace(-3.2 * sigma_v, 3.2 * sigma_v, max(n_probe, 64))
        w = np.exp(-v**2 / (2.0 * sigma_v**2))
        f1 = 0.5 * Gamma_i + 1j * (dp - pump.direction * pump.k * v)
        f2 = 0.5 * Gamma_e + 1j * (two_det - residual_k * v)
        integ = np.abs(w * (Op * Oc) / (4.0 * f1 * f2 + Oc**2))
        res_fwhm = fwhm_interp(v, integ)
        if not np.isfinite(res_fwhm) or res_fwhm <= 0:
            res_fwhm = Gamma_i / (2.0 * kp)
        return max(res_fwhm / PRED_V_OVERSAMPLE, 1e-3)

    def _converged_predictive_waveform(self, params, fields, pump, coupling, signal,
                                       idler, residual_k, density, L, pm_weight,
                                       tau_axis_ns, od_value, T, iso, user_step):
        """Predictive biphoton waveform on a velocity grid auto-refined to
        convergence.

        The coherent sum over Maxwell velocity classes is a discretized Fourier
        integral of the narrow source |amp(v)|; a step coarser than that resonance
        aliases it, so the navigate-only `biphoton_velocity_step` (sized for the
        σ_v-wide Doppler profile) is far too coarse and the absolute BTW width
        swings with it. Start from a step that oversamples the measured resonance —
        never coarser than the user step — then halve until the |ψ|² FWHM is stable
        within `PRED_V_FWHM_TOL`. Returns the converged waveform, its velocity grid
        and weights, the step used, and a convergence flag (False if the point cap
        `PRED_V_MAX_POINTS` is reached first → width is qualitative only).
        """
        def solve_on(v_grid, weights):
            wf = self._predictive_waveform(
                params, fields, v_grid, weights, pump, coupling, signal, idler,
                residual_k=residual_k, density=density, L=L, pm_weight=pm_weight,
                tau_axis_ns=tau_axis_ns, od_value=od_value)
            fw = float(fwhm_interp(tau_axis_ns, np.abs(wf["psi_tau"])**2))
            return wf, fw

        dv = min(float(user_step),
                 self._predictive_velocity_step(params, pump, coupling,
                                                residual_k, T, iso))
        v_grid, weights = doppler.velocity_grid(T, dv=dv, cutoff_sigma=3.0,
                                                mass=iso.mass)
        wf, fw = solve_on(v_grid, weights)
        converged = False
        for _ in range(PRED_V_MAX_REFINE):
            step2 = 0.5 * dv
            v2, w2 = doppler.velocity_grid(T, dv=step2, cutoff_sigma=3.0,
                                           mass=iso.mass)
            if v2.size > PRED_V_MAX_POINTS:
                break
            wf2, fw2 = solve_on(v2, w2)
            rel = abs(fw2 - fw) / max(fw2, 1e-12)
            v_grid, weights, wf, dv, fw = v2, w2, wf2, step2, fw2
            if np.isfinite(rel) and rel < PRED_V_FWHM_TOL:
                converged = True
                break
        return wf, v_grid, weights, float(dv), bool(converged)

    def compute_biphoton(self, params):
        T = params.get("biphoton_temp_c", self.topology.default_temp_c) + 273.15
        L = params.get("biphoton_cell_mm", self.topology.default_cell_mm) * 1e-3
        fields = self._fields_from_params(params)
        fmap = {f.role: f for f in fields}
        detail = params.get("phase_detail", "Balanced")
        fine_phase = str(detail).lower() == "fine"
        iso = self.topology.isotope
        v_step = params.get("biphoton_velocity_step", 2.0)
        v_grid, weights = doppler.velocity_grid(T, dv=v_step, cutoff_sigma=3.0,
                                                mass=iso.mass)
        pump = fmap["pump"]
        coupling = fmap["coupling"]
        signal = fmap["signal"]
        idler = fmap["idler"]
        delta_k = phase_mismatch(
            fields,
            signal_angle_deg=signal.angle_deg,
            idler_angle_deg=idler.angle_deg,
            reference_delta_k=self.topology.reference_delta_k,
        )
        pm_weight_longitudinal = float(phase_matching_weight(np.array([delta_k]), L)[0])
        delta_k_absolute = phase_mismatch(
            fields,
            signal_angle_deg=signal.angle_deg,
            idler_angle_deg=idler.angle_deg,
            reference_delta_k=0.0,
        )
        pm_weight_absolute = float(phase_matching_weight(
            np.array([delta_k_absolute]), L)[0])
        vector_pm = phase_mismatch_vector(
            fields,
            signal_angle_deg=signal.angle_deg,
            idler_angle_deg=idler.angle_deg,
            reference_delta_k=self.topology.reference_delta_k,
        )
        pm_weight = float(phase_matching_weight(
            np.array([vector_pm["delta_k_vector"]]), L)[0])
        density = species.number_density(iso, T)
        default_density = species.number_density(iso, self.topology.default_temp_c + 273.15)
        residual_k = pump.direction * pump.k + coupling.direction * coupling.k
        if self.topology.reference_od is not None:
            od_estimate = (self.topology.reference_od
                           * density / max(default_density, 1e-30)
                           * L / max(self.topology.default_cell_mm * 1e-3, 1e-30))
        else:
            od_estimate = np.nan

        # Absolute pair rate is reference-anchored with physical scaling (pump
        # power, density/OD, vector phase matching). The lumped 4-level model does
        # not pin the absolute collection coefficient, so — exactly like the
        # squeezing-mode `line_strength` residual — the rate magnitude is anchored
        # to the reference while its dependences are physical. Shared by both modes.
        pump_mw = params.get("pump_biphoton_uw", self.topology.default_pump_uw) * 1e-3
        coupling_scale = max(params.get("coupling_mw", self.topology.default_coupling_mw),
                             0.0) / max(self.topology.default_coupling_mw, 1e-12)
        pair_rate = (self.topology.pair_rate_cps_per_mw * pump_mw
                     * math.sqrt(max(coupling_scale, 0.0)) * pm_weight
                     * math.sqrt(max(density, 1e-30) / max(default_density, 1e-30)))

        model = params.get("biphoton_model", BIPHOTON_PREDICTIVE)
        predictive = model == BIPHOTON_PREDICTIVE
        tau_axis_ns = np.linspace(0.0, params.get("tau_max_ns", 12.0), 481)

        if predictive:
            od_value = (od_estimate if np.isfinite(od_estimate)
                        else self._leg_optical_depth(
                            density, L, self.topology.levels[1].gamma_mhz,
                            idler.wavelength_nm))
            # The velocity-class coherent sum aliases on a step coarser than the
            # narrow source resonance, so the predictive grid is auto-refined to
            # convergence (the user `biphoton_velocity_step` acts as an upper bound).
            wf, v_grid, weights, v_step, velocity_converged = (
                self._converged_predictive_waveform(
                    params, fields, pump, coupling, signal, idler,
                    residual_k=residual_k, density=density, L=L,
                    pm_weight=pm_weight, tau_axis_ns=tau_axis_ns,
                    od_value=od_value, T=T, iso=iso, user_step=v_step))
            psi_tau = wf["psi_tau"]
            source_v = wf["source_v"]
            od_estimate = wf["od_phys"]
            source_bandwidth_mhz = wf["bandwidth_mhz"]
            regime = wf["regime"]
            brightness_limit = wf["brightness_limit_cps_per_mhz"]
        else:
            lower_gamma = max(self.topology.levels[1].gamma_mhz, 0.1) * MHZ_ANG
            upper_gamma = max(self.topology.levels[2].gamma_mhz, 0.1) * MHZ_ANG
            pump_det = params.get("pump_detuning_mhz", 0.0) * MHZ_ANG
            two_det = (params.get("pump_detuning_mhz", 0.0)
                       + params.get("coupling_detuning_mhz", 0.0)) * MHZ_ANG
            lower = lower_gamma / 2.0 + 1j * (pump_det - pump.direction * pump.k * v_grid)
            upper = upper_gamma / 2.0 + 1j * (two_det - residual_k * v_grid)
            drive = (pump.rabi_mhz * coupling.rabi_mhz) * (MHZ_ANG ** 2)
            source_v = weights * drive / (lower * upper)
            source_v *= math.sqrt(max(pm_weight, 0.0))
            tau_s = tau_axis_ns * 1e-9
            phase_k = abs(idler.k)
            coherent = np.exp(1j * phase_k * v_grid[:, None] * tau_s[None, :])
            psi_tau = (source_v[:, None] * coherent).sum(axis=0)
            psi_tau *= np.exp(-tau_axis_ns / max(self.topology.emission_decay_ns, 1e-12))
            if np.max(np.abs(psi_tau)) > 0:
                psi_tau = psi_tau / np.max(np.abs(psi_tau))
            source_bandwidth_mhz = float(self.topology.reference_bandwidth_mhz or 300.0)
            regime = "calibrated"
            brightness_limit = float("nan")
            # Calibrated width is set by the injected emission lifetime, not the
            # velocity-sum dephasing, so it is stable on the user grid (no refine).
            velocity_converged = True

        angle_axis = np.linspace(max(idler.angle_deg - 4.0, 0.0),
                                 idler.angle_deg + 4.0, 181)
        angle_dk = np.array([
            phase_mismatch_vector(
                fields, idler_angle_deg=a,
                reference_delta_k=self.topology.reference_delta_k
            )["delta_k_vector"]
            for a in angle_axis
        ])
        acceptance = phase_matching_weight(angle_dk, L)
        signal_axis = phase_matching_2d = None
        idler_axis_2d = None
        if fine_phase:
            signal_axis = np.linspace(max(signal.angle_deg - 3.0, 0.0),
                                      signal.angle_deg + 3.0, 121)
            idler_axis_2d = np.linspace(max(idler.angle_deg - 3.0, 0.0),
                                        idler.angle_deg + 3.0, 121)
            dk_2d = phase_mismatch_vector_grid(
                fields, signal_axis, idler_axis_2d,
                reference_delta_k=self.topology.reference_delta_k)
            phase_matching_2d = phase_matching_weight(dk_2d, L)

        return {
            "kind": "biphoton",
            "topology": self.topology,
            "fields": fields,
            "tau_axis_ns": tau_axis_ns,
            "psi_tau": psi_tau,
            "v_grid": v_grid,
            "velocity_weights": weights,
            "source_v": source_v,
            "angle_axis_deg": angle_axis,
            "phase_matching": acceptance,
            "signal_angle_axis_deg": signal_axis,
            "idler_angle_axis_2d_deg": idler_axis_2d,
            "phase_matching_2d": phase_matching_2d,
            "delta_k": float(delta_k),
            "delta_k_absolute": float(delta_k_absolute),
            "delta_k_z_relative": float(vector_pm["delta_k_z_relative"]),
            "delta_k_z_absolute": float(vector_pm["delta_k_z_absolute"]),
            "delta_k_x": float(vector_pm["delta_k_x"]),
            "delta_k_vector": float(vector_pm["delta_k_vector"]),
            "phase_match_weight": pm_weight,
            "phase_match_weight_vector": pm_weight,
            "phase_match_weight_longitudinal": pm_weight_longitudinal,
            "phase_match_weight_absolute": pm_weight_absolute,
            "phase_detail": "Fine" if fine_phase else "Balanced",
            "energy_mismatch_hz": float(energy_mismatch_hz(fields)),
            "pair_rate_cps": float(pair_rate),
            "density": float(density),
            "od_estimate": float(od_estimate),
            "source_bandwidth_mhz": float(source_bandwidth_mhz),
            "temperature_K": float(T),
            "cell_length_m": float(L),
            "residual_two_photon_k": float(residual_k),
            "biphoton_model": model,
            "predictive": bool(predictive),
            "regime": regime,
            "brightness_limit_cps_per_mhz": float(brightness_limit),
            "velocity_step_used": float(v_step),
            "velocity_converged": bool(velocity_converged),
            "notes": self.topology.notes,
        }

# =========================================================
# Hamiltonians
# =========================================================
def _add_static_drive(H, ground, omega):
    for excited in EXCITED_STATES:
        omega_ge = omega * TRANSITION_DIPOLE_SCALE[ground, excited]
        H[ground, excited] += omega_ge / 2
        H[excited, ground] += omega_ge / 2


def _add_sideband_drive(H, ground, omega):
    for excited in EXCITED_STATES:
        H[excited, ground] += omega * TRANSITION_DIPOLE_SCALE[ground, excited] / 2


def _polarization_coherence(rho, ground):
    return sum(TRANSITION_DIPOLE_SCALE[ground, e] * rho[:, e, ground]
               for e in EXCITED_STATES)


def static_hamiltonian_at_Deff_zero(Op_A, Op_B, Os, delta, branch):
    """H₀ with Δ_eff = 0, so the only v / Δ_eff dependence is added later."""
    H0 = np.zeros((N_LEVELS, N_LEVELS), dtype=complex)
    H0[G2, G2] = delta
    H0[E2, E2] = -OMEGA_EXCITED_HF
    H0[E3, E3] = 0.0
    if branch == -1:
        _add_static_drive(H0, G1, Op_A)
        _add_static_drive(H0, G2, Os)
        return H0
    if branch == +1:
        _add_static_drive(H0, G1, Os)
        _add_static_drive(H0, G2, Op_B)
        return H0
    raise ValueError(f"branch must be one of {BRANCHES}, got {branch}")


def sideband_hamiltonian(Op_A, Op_B, Oc, branch):
    Hp = np.zeros((N_LEVELS, N_LEVELS), dtype=complex)
    if branch == -1:
        _add_sideband_drive(Hp, G1, Oc)
        _add_sideband_drive(Hp, G2, Op_B)
        return Hp
    if branch == +1:
        _add_sideband_drive(Hp, G1, Op_A)
        _add_sideband_drive(Hp, G2, Oc)
        return Hp
    raise ValueError(f"branch must be one of {BRANCHES}, got {branch}")


def sideband_template(Op_A, Op_B, Oc, branch):
    Hp = sideband_hamiltonian(Op_A, Op_B, Oc, branch)
    Cp = comm_super(Hp)
    Cm = comm_super(Hp.conj().T)
    return Cp, Cm


# =========================================================
# χ-matrix table   (T- and Δ-independent)
# =========================================================
def _coherence_weights(ground):
    """w such that Σ_k w[k]·vec(ρ)[k] = Σ_e scale[g,e]·ρ[e,g] (vec row-major)."""
    w = np.zeros(N_LEVELS * N_LEVELS, dtype=complex)
    for e in EXCITED_STATES:
        w[e * N_LEVELS + ground] = TRANSITION_DIPOLE_SCALE[ground, e]
    return w


def chi_matrix_table(Op_A, Op_B, Os_ref, Oc_ref, delta_axis, Delta_eff_axis, branch):
    """
    Two solves per probe-detuning point to extract (χ̄_ss, χ̄_cs, χ̄_sc, χ̄_cc)
    on a 2-D (δ, Δ_eff) grid. Returns each array (n_delta, n_deff), complex.

    With numba available the whole grid runs in one fused compiled kernel
    (`kernels.floquet_chi_grid`); the NumPy δ-loop below is the fallback and
    the reference implementation (tests/test_kernels.py pins them together).
    """
    probe_ground = G2 if branch == -1 else G1
    conj_ground = G1 if branch == -1 else G2
    n_d = delta_axis.size
    n_de = Delta_eff_axis.size

    Cp_no_c, Cm_no_c = sideband_template(Op_A, Op_B, 0.0, branch)      # solve 1
    Cp_c, Cm_c = sideband_template(Op_A, Op_B, Oc_ref, branch)        # solve 2

    if kernels.available():
        # H₀(δ) is affine in δ (only H₀[G2,G2] = δ), so L₀(δ) = L₀(0) + δ·C_δ.
        E_g2 = np.zeros((N_LEVELS, N_LEVELS), dtype=complex)
        E_g2[G2, G2] = 1.0
        C_delta = comm_super(E_g2)
        w_probe = _coherence_weights(probe_ground)
        w_conj = _coherence_weights(conj_ground)
        delta_axis = np.ascontiguousarray(delta_axis, dtype=float)
        deff_axis = np.ascontiguousarray(Delta_eff_axis, dtype=float)

        L0_1 = build_liouvillian(
            static_hamiltonian_at_Deff_zero(Op_A, Op_B, Os_ref, 0.0, branch), ATOM)
        probe_a, conj_a = kernels.floquet_chi_grid(
            L0_1, C_delta, ATOM.S_v, Cp_no_c, Cm_no_c, delta_axis, deff_axis,
            OMEGA_HF, float(branch), w_probe, w_conj, N_LEVELS)

        L0_2 = build_liouvillian(
            static_hamiltonian_at_Deff_zero(Op_A, Op_B, 0.0, 0.0, branch), ATOM)
        probe_b, conj_b = kernels.floquet_chi_grid(
            L0_2, C_delta, ATOM.S_v, Cp_c, Cm_c, delta_axis, deff_axis,
            OMEGA_HF, float(branch), w_probe, w_conj, N_LEVELS)

        return probe_a / Os_ref, conj_a / Os_ref, probe_b / Oc_ref, conj_b / Oc_ref

    chi_ss = np.zeros((n_d, n_de), dtype=complex)
    chi_cs = np.zeros((n_d, n_de), dtype=complex)
    chi_sc = np.zeros((n_d, n_de), dtype=complex)
    chi_cc = np.zeros((n_d, n_de), dtype=complex)

    for i, delta in enumerate(delta_axis):
        Omega_beat = OMEGA_HF - branch * delta

        # ---- Solve 1: probe drive only ----
        H0_1 = static_hamiltonian_at_Deff_zero(Op_A, Op_B, Os_ref, delta, branch)
        L0_1 = build_liouvillian(H0_1, ATOM)
        rho0_a, rhop_a = floquet_solve(
            L0_1, Cp_no_c, Cm_no_c, Omega_beat, Delta_eff_axis, ATOM.S_v, N_LEVELS)
        probe_a = _polarization_coherence(rho0_a, probe_ground)
        conj_a = _polarization_coherence(rhop_a, conj_ground)
        chi_ss[i] = probe_a / Os_ref
        chi_cs[i] = conj_a / Os_ref

        # ---- Solve 2: conjugate seed only ----
        H0_2 = static_hamiltonian_at_Deff_zero(Op_A, Op_B, 0.0, delta, branch)
        L0_2 = build_liouvillian(H0_2, ATOM)
        rho0_b, rhop_b = floquet_solve(
            L0_2, Cp_c, Cm_c, Omega_beat, Delta_eff_axis, ATOM.S_v, N_LEVELS)
        probe_b = _polarization_coherence(rho0_b, probe_ground)
        conj_b = _polarization_coherence(rhop_b, conj_ground)
        chi_sc[i] = probe_b / Oc_ref
        chi_cc[i] = conj_b / Oc_ref

    return chi_ss, chi_cs, chi_sc, chi_cc


# =========================================================
# Probe-detuning axis
# =========================================================
def branch_center_GHz(Delta_GHz, branch):
    if branch not in BRANCHES:
        raise ValueError(f"branch must be one of {BRANCHES}, got {branch}")
    return Delta_GHz + branch * constants.NU_HF / 1e9


def probe_scan_axis_GHz(Delta_GHz, coarse_points=None, fine_points=None,
                        window_mhz=None, scan_min=None, scan_max=None,
                        branches=BRANCHES):
    """Probe-frequency scan axis [GHz], referenced to the 85Rb F=2 → F'=3 line."""
    coarse_points = SCAN_COARSE_POINTS if coarse_points is None else coarse_points
    fine_points = SCAN_FINE_POINTS if fine_points is None else fine_points
    window_mhz = RESONANCE_WINDOW_MHZ if window_mhz is None else window_mhz
    scan_min = SCAN_MIN_GHZ if scan_min is None else scan_min
    scan_max = SCAN_MAX_GHZ if scan_max is None else scan_max

    coarse = np.linspace(scan_min, scan_max, coarse_points)
    half_window = window_mhz * 1e-3
    parts = [coarse]
    for branch in branches:
        center = branch_center_GHz(Delta_GHz, branch)
        fmin = max(scan_min, center - half_window)
        fmax = min(scan_max, center + half_window)
        if fmin < fmax:
            parts.append(np.linspace(fmin, fmax, fine_points))
    axis = np.unique(np.concatenate(parts))
    exclusion_GHz = PUMP_OVERLAP_EXCLUSION_MHZ * 1e-3
    return axis[np.abs(axis - Delta_GHz) > exclusion_GHz]


def two_photon_detuning_from_probe_scan(probe_GHz, Delta_GHz, branch):
    delta_Hz = (probe_GHz - branch_center_GHz(Delta_GHz, branch)) * 1e9
    return 2 * np.pi * delta_Hz


def _optical_k_from_offset(offset_rad_s):
    return (constants.OMEGA_D1 + np.asarray(offset_rad_s, dtype=float)) / constants.C_LIGHT


def seeded_phase_mismatch_z(D_GHz, probe_axis_GHz, angle_deg=SEEDED_PHASE_ANGLE_DEG,
                            n_seed=None, n_conj=None):
    """Longitudinal seeded-FWM mismatch: Delta k_z = 2k_p - k_s - k_c."""
    theta = math.radians(float(angle_deg))
    pump_offset = 2 * np.pi * float(D_GHz) * 1e9
    seed_offset = 2 * np.pi * np.asarray(probe_axis_GHz, dtype=float) * 1e9
    conj_offset = 2.0 * pump_offset - seed_offset
    k_pump = _optical_k_from_offset(pump_offset)
    k_seed = _optical_k_from_offset(seed_offset)
    k_conj = _optical_k_from_offset(conj_offset)
    if n_seed is not None:
        k_seed = k_seed * np.asarray(n_seed, dtype=float)
    if n_conj is not None:
        k_conj = k_conj * np.asarray(n_conj, dtype=float)
    return 2.0 * k_pump - (k_seed + k_conj) * math.cos(theta)


def _safe_refractive_index(chi_bar, N_atoms, line_strength):
    chi = observables.chi_phys(chi_bar, N_atoms, line_strength=line_strength)
    # Keep the cheap refractive correction in the dilute-vapor regime. Very large
    # dispersive excursions usually mean the simplified propagation model is being
    # pushed outside its calibrated range, so clip rather than letting k explode.
    return 1.0 + np.clip(0.5 * np.real(chi), -1e-5, 1e-5)


def _segment_profile_from_absorption(chi_bar, N_atoms, line_strength, nseg, L=L_CELL):
    chi = observables.chi_phys(chi_bar, N_atoms, line_strength=line_strength)
    alpha = np.maximum(K_VEC * np.imag(chi), 0.0)
    od = float(np.nanmedian(alpha) * L) if alpha.size else 0.0
    od = float(np.clip(od, 0.0, 2.0))
    z_frac = (np.arange(nseg, dtype=float) + 0.5) / nseg
    return np.exp(-0.5 * od * z_frac), od


def _gaussian_overlap_profile(nseg, L, w_pump, w_probe, angle_deg):
    """Amplitude overlap of crossed Gaussian beams along the cell."""
    theta = math.radians(float(angle_deg))
    z = ((np.arange(nseg, dtype=float) + 0.5) / nseg - 0.5) * L
    separation = np.abs(z * math.tan(theta))
    waist_sq = max(float(w_pump) ** 2 + float(w_probe) ** 2, 1e-30)
    profile = np.exp(-(separation ** 2) / waist_sq)
    return profile / max(float(np.nanmax(profile)), 1e-30)


def _ultra_phase_mismatch(D_GHz, probe_axis_GHz, chi_ss_avg, chi_cc_avg,
                          N_atoms, line_strength, angle_deg):
    """Single-pass refractive phase-mismatch correction for Ultra.

    The refractive indices depend only on χ (not on Δk), so there is no fixed
    point to iterate — one pass is exact. ``max_change`` reports the
    vacuum→refractive shift in Δk (a meaningful diagnostic), rather than the
    last-iteration delta of a no-op loop.
    """
    delta_k = seeded_phase_mismatch_z(D_GHz, probe_axis_GHz, angle_deg=angle_deg)
    max_change = 0.0
    n_seed = np.ones_like(probe_axis_GHz, dtype=float)
    n_conj = np.ones_like(probe_axis_GHz, dtype=float)
    for _ in range(ULTRA_PHASE_ITERATIONS):
        n_seed = _safe_refractive_index(chi_ss_avg, N_atoms, line_strength)
        n_conj = _safe_refractive_index(chi_cc_avg, N_atoms, line_strength)
        new_delta_k = seeded_phase_mismatch_z(
            D_GHz, probe_axis_GHz, angle_deg=angle_deg,
            n_seed=n_seed, n_conj=n_conj)
        max_change = float(np.nanmax(np.abs(new_delta_k - delta_k)))
        delta_k = new_delta_k
    return delta_k, n_seed, n_conj, max_change


def _ultra_segmented_gain(chi_ss_avg, chi_sc_avg, chi_cs_avg, chi_cc_avg,
                          k_probe, k_conj, L, N_atoms, line_strength,
                          delta_k_z, segment_profile, spatial_profile,
                          P_pump, P_seed):
    """Segmented propagation with approximate dynamic pump depletion.

    In-cell propagation loss is applied once, downstream, as the
    ``segmented_loss_noise_squeezing_dB`` vacuum admixture (the codebase's
    beamsplitter loss convention); it is deliberately not re-applied here as
    field attenuation, which double-counted the same ``segment_od``.
    """
    nseg = int(segment_profile.size)
    dz = L / max(nseg, 1)
    M = observables._gain_matrix_from_chi(
        chi_ss_avg, chi_sc_avg, chi_cs_avg, chi_cc_avg,
        k_probe, k_conj, N_atoms, constants.DIPOLE_D1, line_strength,
        delta_k_z=delta_k_z)
    n = M.shape[0]
    amp = np.zeros((n, 2), dtype=complex)
    amp[:, 0] = math.sqrt(max(float(P_seed), 1e-30))
    T_total = np.broadcast_to(np.eye(2, dtype=complex), (n, 2, 2)).copy()
    pump_remaining = np.full(n, max(float(P_pump), 1e-30), dtype=float)

    for coupling_scale, spatial_scale in zip(segment_profile, spatial_profile):
        pump_scale = np.sqrt(np.clip(pump_remaining / max(float(P_pump), 1e-30),
                                     0.0, 1.0))
        Mz = M.copy()
        scale = float(coupling_scale) * float(spatial_scale) * pump_scale
        Mz[:, 0, 1] *= scale
        Mz[:, 1, 0] *= scale
        Tseg = observables.matrix_exp_2x2(Mz, dz)
        amp = np.einsum("nij,nj->ni", Tseg, amp)
        T_total = Tseg @ T_total
        seed_added = np.maximum(np.abs(amp[:, 0]) ** 2 - float(P_seed), 0.0)
        conj_power = np.maximum(np.abs(amp[:, 1]) ** 2, 0.0)
        pump_remaining = np.maximum(float(P_pump) - seed_added - conj_power, 0.0)

    G_s = np.abs(amp[:, 0]) ** 2 / max(float(P_seed), 1e-30)
    G_c = np.abs(amp[:, 1]) ** 2 / max(float(P_seed), 1e-30)
    return G_s, G_c, T_total, pump_remaining


def _single_branch(branch, branches):
    """Resolve one physical FWM channel; do not merge distinct Raman branches."""
    if branches is None:
        if branch not in BRANCHES:
            raise ValueError(f"branch must be one of {BRANCHES}, got {branch}")
        return branch

    branches = tuple(branches)
    if len(branches) != 1:
        raise ValueError(
            "FWM Raman branches are separate probe/conjugate mode pairs; "
            "compute them one at a time instead of summing susceptibilities."
        )
    only = branches[0]
    if only not in BRANCHES:
        raise ValueError(f"branch must be one of {BRANCHES}, got {only}")
    return only


# =========================================================
# High-level spectrum  (one call → gain + squeezing curves)
# =========================================================
def compute_spectrum(D_GHz, *,
                     T=T_CELL, P_pump=P_PUMP, P_probe=P_PROBE,
                     w_pump=W_PUMP, w_probe=W_PROBE,
                     line_strength=None, loss_frac=LOSS_FRAC, qe=QE_DETECTOR,
                     L=L_CELL,
                     coarse_points=None, fine_points=None, window_mhz=None,
                     scan_min=None, scan_max=None,
                     velocity_step=None, velocity_cutoff=None,
                     branch=DEFAULT_BRANCH, branches=None,
                     phase_detail=PHASE_LEGACY,
                     pump_probe_angle_deg=SEEDED_PHASE_ANGLE_DEG,
                     model_fidelity=None):
    """Full pipeline for one one-photon detuning Δ = 2π·D_GHz·1e9 (see README)."""
    branch = _single_branch(branch, branches)
    if line_strength is None:
        line_strength = constants.LINE_STRENGTH_FACTOR
    # `line_strength` is now a dimensionless residual (≈1.0); the physical
    # macroscopic-coupling scale is computed from first principles and applied
    # on top of it. `coupling_ls` is what enters every χ_phys / gain call below.
    coupling_norm = physical_coupling_norm(branch)
    coupling_ls = line_strength * coupling_norm
    eta = qe * (1.0 - loss_frac)

    Op_A = rabi_freq(P_pump, w_pump)
    Op_B = Op_A
    Os = rabi_freq(P_probe, w_probe)
    Os_ref = Os
    Oc_ref = Os                              # χ̄ is independent of |Ω_ref|

    # Pure-85Rb CRC vapor density, consistent with the AutoOD-validated
    # absorption path (`hyperfine.number_density`). The other (natural-abundance,
    # Steck) `atoms.rb85_density` understated N for the enriched FWM cell.
    N_atoms = hyperfine.number_density(T)
    Delta = 2 * np.pi * D_GHz * 1e9

    probe_axis_GHz = probe_scan_axis_GHz(
        D_GHz, coarse_points, fine_points, window_mhz, scan_min, scan_max,
        branches=(branch,))
    velocity_step = VELOCITY_STEP_MPS if velocity_step is None else velocity_step
    velocity_cutoff = VELOCITY_CUTOFF_SIGMA if velocity_cutoff is None else velocity_cutoff
    v_grid, weights = doppler.velocity_grid(
        T, dv=velocity_step, cutoff_sigma=velocity_cutoff)
    Delta_eff_axis = doppler.build_Delta_eff_axis(Delta, Delta, v_grid)

    chi_ss_avg = np.zeros(probe_axis_GHz.size, dtype=complex)
    chi_cs_avg = np.zeros(probe_axis_GHz.size, dtype=complex)
    chi_sc_avg = np.zeros(probe_axis_GHz.size, dtype=complex)
    chi_cc_avg = np.zeros(probe_axis_GHz.size, dtype=complex)

    delta_axis = two_photon_detuning_from_probe_scan(probe_axis_GHz, D_GHz, branch)
    ch_ss, ch_cs, ch_sc, ch_cc = chi_matrix_table(
        Op_A, Op_B, Os_ref, Oc_ref, delta_axis, Delta_eff_axis, branch)
    chi_ss_avg += doppler.doppler_average(ch_ss, Delta_eff_axis, Delta, v_grid, weights)
    chi_cs_avg += doppler.doppler_average(ch_cs, Delta_eff_axis, Delta, v_grid, weights)
    chi_sc_avg += doppler.doppler_average(ch_sc, Delta_eff_axis, Delta, v_grid, weights)
    chi_cc_avg += doppler.doppler_average(ch_cc, Delta_eff_axis, Delta, v_grid, weights)

    phase_detail = (phase_detail or PHASE_LEGACY).lower()
    delta_k_z = None
    delta_k_z_vacuum = None
    k_probe_prop = K_VEC
    k_conj_prop = K_VEC
    propagation_segments = 1
    segment_profile = None
    segment_od = 0.0
    spatial_profile = None
    ultra_phase_iterations = 0
    ultra_phase_max_change = 0.0
    ultra_pump_remaining_min = float(P_pump)
    zeeman_status = "inactive"
    zeeman_correction = 1.0
    if phase_detail != PHASE_LEGACY:
        delta_k_z_vacuum = seeded_phase_mismatch_z(
            D_GHz, probe_axis_GHz, angle_deg=pump_probe_angle_deg)
        delta_k_z = delta_k_z_vacuum
        if phase_detail == PHASE_FINE:
            n_seed = _safe_refractive_index(chi_ss_avg, N_atoms, coupling_ls)
            n_conj = _safe_refractive_index(chi_cc_avg, N_atoms, coupling_ls)
            seed_offset = 2 * np.pi * probe_axis_GHz * 1e9
            pump_offset = 2 * np.pi * float(D_GHz) * 1e9
            conj_offset = 2.0 * pump_offset - seed_offset
            k_probe_prop = _optical_k_from_offset(seed_offset) * n_seed
            k_conj_prop = _optical_k_from_offset(conj_offset) * n_conj
            delta_k_z = seeded_phase_mismatch_z(
                D_GHz, probe_axis_GHz, angle_deg=pump_probe_angle_deg,
                n_seed=n_seed, n_conj=n_conj)
            propagation_segments = 16
            segment_profile, segment_od = _segment_profile_from_absorption(
                chi_ss_avg, N_atoms, coupling_ls, propagation_segments, L=L)
        elif phase_detail == PHASE_ULTRA:
            delta_k_z, n_seed, n_conj, ultra_phase_max_change = _ultra_phase_mismatch(
                D_GHz, probe_axis_GHz, chi_ss_avg, chi_cc_avg, N_atoms,
                coupling_ls, pump_probe_angle_deg)
            ultra_phase_iterations = ULTRA_PHASE_ITERATIONS
            seed_offset = 2 * np.pi * probe_axis_GHz * 1e9
            pump_offset = 2 * np.pi * float(D_GHz) * 1e9
            conj_offset = 2.0 * pump_offset - seed_offset
            k_probe_prop = _optical_k_from_offset(seed_offset) * n_seed
            k_conj_prop = _optical_k_from_offset(conj_offset) * n_conj
            propagation_segments = ULTRA_PROPAGATION_SEGMENTS
            segment_profile, segment_od = _segment_profile_from_absorption(
                chi_ss_avg, N_atoms, coupling_ls, propagation_segments, L=L)
            spatial_profile = _gaussian_overlap_profile(
                propagation_segments, L, w_pump, w_probe, pump_probe_angle_deg)
            try:
                from .. import zeeman as _zeeman
                z_atom = _zeeman.rb85_d1_double_lambda_zeeman()
                zeeman_correction = float(
                    getattr(z_atom, "lumped_strength_correction", 1.0))
                zeeman_status = (
                    f"24-level CG-sum consistency check = {zeeman_correction:.4f} "
                    "(lumped 3·C_F² reproduced); diagnostic only — full Floquet "
                    "scan not run in Ultra v1")
            except Exception as exc:                  # pragma: no cover
                zeeman_status = f"unavailable: {exc}"

    if phase_detail == PHASE_ULTRA:
        if segment_profile is None:
            segment_profile = np.ones(propagation_segments, dtype=float)
        if spatial_profile is None:
            spatial_profile = np.ones(propagation_segments, dtype=float)
        # zeeman_correction is a CG-sum consistency diagnostic (≡1.0 by
        # construction), not an active correction, so it is no longer multiplied
        # into the coupling.
        G_s, G_c, _T, pump_remaining = _ultra_segmented_gain(
            chi_ss_avg, chi_sc_avg, chi_cs_avg, chi_cc_avg,
            k_probe_prop, k_conj_prop, L, N_atoms,
            coupling_ls, delta_k_z,
            segment_profile, spatial_profile, P_pump, P_probe)
        ultra_pump_remaining_min = float(np.nanmin(pump_remaining))
    else:
        G_s, G_c, _ = observables.gain_from_chi(
            chi_ss_avg, chi_sc_avg, chi_cs_avg, chi_cc_avg,
            k_probe_prop, k_conj_prop, L, N_atoms, line_strength=coupling_ls,
            delta_k_z=delta_k_z, propagation_segments=propagation_segments,
            segment_profile=segment_profile)

    # The propagation above is linear in the (undepleted) pump. At high density
    # it overshoots the energy the pump can supply, so apply Manley-Rowe pump-
    # depletion saturation. This is negligible in the validated regime and only
    # caps the runaway when (G_s−1)·P_seed approaches P_pump/2.
    G_s_smallsignal = G_s
    G_s, G_c = observables.pump_depletion_saturation(G_s, G_c, P_pump, P_probe)
    if phase_detail == PHASE_ULTRA:
        in_cell_loss_frac = float(np.clip(1.0 - np.exp(-segment_od), 0.0, 1.0))
        S_dB = observables.segmented_loss_noise_squeezing_dB(
            G_s, G_c, eta, in_cell_loss_frac=in_cell_loss_frac,
            seed_excess_noise=0.0, pump_scatter_noise=0.0,
            eom_residual_noise=0.0)
    else:
        S_dB = observables.intensity_difference_squeezing_dB(G_s, G_c, eta)

    return {
        "D_GHz": D_GHz,
        "probe_axis_GHz": probe_axis_GHz,
        "G_s": G_s,
        "G_c": G_c,
        "S_dB": S_dB,
        "G_s_smallsignal_peak": float(np.nanmax(G_s_smallsignal)),
        "pump_depletion_cap": 1.0 + 0.5 * P_pump / max(P_probe, 1e-30),
        "phase_detail": phase_detail,
        "model_fidelity": model_fidelity or phase_detail,
        "pump_probe_angle_deg": pump_probe_angle_deg,
        "delta_k_z": delta_k_z,
        "delta_k_z_vacuum": delta_k_z_vacuum,
        "phase_segments": propagation_segments,
        "segment_absorption_od": segment_od,
        "ultra_phase_iterations": ultra_phase_iterations,
        "ultra_phase_max_change": ultra_phase_max_change,
        "ultra_dynamic_depletion": phase_detail == PHASE_ULTRA,
        "ultra_in_cell_loss_noise": phase_detail == PHASE_ULTRA,
        "ultra_pump_remaining_min": ultra_pump_remaining_min,
        "ultra_spatial_overlap_min": (float(np.nanmin(spatial_profile))
                                      if spatial_profile is not None else 1.0),
        "zeeman_status": zeeman_status,
        "zeeman_correction": zeeman_correction,
        "eta": eta,
        "cell_length_m": L,
        "branch": branch,
        "coupling_norm": coupling_norm,
        "line_strength_residual": line_strength,
        "effective_line_strength": coupling_ls,
        "N_atoms": N_atoms,
        "sigma_v": np.sqrt(constants.KB * T / constants.MASS_85RB),
        "n_velocity": v_grid.size,
        "Op_A_2pi_GHz": Op_A / (2 * np.pi) / 1e9,
        "Os_2pi_MHz": Os / (2 * np.pi) / 1e6,
        "raman_center_minus_GHz": branch_center_GHz(D_GHz, -1),
        "raman_center_plus_GHz": branch_center_GHz(D_GHz, +1),
    }


def operating_point(spectrum, delta_mhz, branch=-1):
    """Read G_s / G_c / squeezing at a chosen δ (MHz) on the selected branch."""
    probe_GHz = (spectrum["raman_center_minus_GHz"] if branch == -1
                 else spectrum["raman_center_plus_GHz"]) + delta_mhz * 1e-3
    x = spectrum["probe_axis_GHz"]
    out = {
        "probe_GHz": probe_GHz,
        "G_s": float(np.interp(probe_GHz, x, spectrum["G_s"])),
        "G_c": float(np.interp(probe_GHz, x, spectrum["G_c"])),
        "S_dB": float(np.interp(probe_GHz, x, spectrum["S_dB"])),
    }
    if spectrum.get("delta_k_z") is not None:
        out["delta_k_z"] = float(np.interp(probe_GHz, x, spectrum["delta_k_z"]))
    if spectrum.get("delta_k_z_vacuum") is not None:
        out["delta_k_z_vacuum"] = float(np.interp(
            probe_GHz, x, spectrum["delta_k_z_vacuum"]))
    return out


# =========================================================
# Scheme wrapper for the generic front-end
# =========================================================
WINDOW_GHZ = 0.55          # half-width of the focused probe window around (−) Raman
TPD_LIMIT_MHZ = 500.0
FIDELITY_FAST = "Fast  (~3 s)"
FIDELITY_BALANCED = "Balanced  (~6 s)"
FIDELITY_HIGH = "High fidelity  (~20 s)"
FIDELITY_ULTRA = "Ultra  (slow)"
FIDELITY_LEGACY_FINE = "Fine  (~20 s)"
FWM_FIDELITY = {
    FIDELITY_FAST:     dict(coarse_points=121, velocity_step=5.0,
                            velocity_cutoff=3.0, phase_detail=PHASE_BALANCED),
    FIDELITY_BALANCED: dict(coarse_points=181, velocity_step=4.0,
                            velocity_cutoff=3.0, phase_detail=PHASE_BALANCED),
    FIDELITY_HIGH:     dict(coarse_points=301, velocity_step=2.0,
                            velocity_cutoff=3.0, phase_detail=PHASE_FINE),
    FIDELITY_ULTRA:    dict(coarse_points=401, velocity_step=1.0,
                            velocity_cutoff=4.0, phase_detail=PHASE_ULTRA),
}
RESOLUTION = FWM_FIDELITY


def normalize_fidelity(value):
    """Map old saved labels onto current user-facing fidelity labels."""
    if value == FIDELITY_LEGACY_FINE:
        return FIDELITY_HIGH
    return value if value in FWM_FIDELITY else FIDELITY_BALANCED


class FWMScheme(Scheme):
    name = "fwm"
    cluster = "D — Wave mixing"
    title = "Four-wave mixing (Squeezing / Biphoton)"
    cache_version = "biphoton-slim-ui-v1"
    defaults_version = "biphoton-slim-ui-defaults-v1"
    cache_observables = True
    caption = ("Squeezing keeps the legacy 85Rb double-Lambda gain model; "
               "Biphoton shows generic SFWM source estimates for cascade and "
               "diamond topologies.")

    def param_schema(self):
        seeded = {"mode": MODE_SEEDED}
        biphoton = {"mode": MODE_BIPHOTON}
        cs_btw = {"mode": MODE_BIPHOTON, "topology": TOPOLOGY_CS_BTW}
        diamond = {"mode": MODE_BIPHOTON, "topology": TOPOLOGY_DIAMOND}
        return [
            ParamSpec("mode", "Mode", "Mode", MODE_SEEDED,
                      choices=(MODE_SEEDED, MODE_BIPHOTON),
                      control="segmented", applies_defaults=True,
                      help="Pick the readout. Selecting a mode resets every knob "
                           "to that mode's recommended default values."),
            ParamSpec("topology", "Topology", "Model", TOPOLOGY_RB87_TELECOM,
                      choices=(TOPOLOGY_RB87_TELECOM, TOPOLOGY_CS_BTW, TOPOLOGY_DIAMOND),
                      visible_if=biphoton, applies_defaults=True,
                      help="Level topology for the spontaneous biphoton source model."),
            ParamSpec("cs_channel", "Cs BTW channel", "Model", CS_CHANNEL_917,
                      choices=(CS_CHANNEL_917, CS_CHANNEL_795), visible_if=cs_btw,
                      applies_defaults=True,
                      help="Selects the Cs cascade channel used for BTW comparison."),
            ParamSpec("biphoton_model", "Source model", "Model", BIPHOTON_PREDICTIVE,
                      choices=BIPHOTON_MODELS, control="segmented",
                      visible_if=biphoton, advanced=True,
                      help="Predictive solves the Doppler-averaged cascade biphoton "
                           "amplitude from first principles (Ω_c² Autler-Townes term, "
                           "velocity-class coherent sum, natural-linewidth decay, OD "
                           "reshaping): the BTW shape, decay, bandwidth and the "
                           "wavelength-dependent width ordering emerge from physics. "
                           "Absolute widths are approximate and the pair rate stays "
                           "reference-anchored (the collection coefficient is not "
                           "derivable in the lumped model — like the squeezing "
                           "line-strength residual). Calibrated is the legacy "
                           "reference-injected estimate."),
            ParamSpec("opd", "OPD — one-photon detuning Δ", "Detunings", 0.9,
                      -1.0, 3.0, 0.1, "GHz",
                      visible_if=seeded,
                      help="ω_pump = ω(F=2→F'=3) + Δ. Sets where the pump sits; recomputes."),
            ParamSpec("tpd", "TPD — two-photon detuning δ", "Detunings", -8.0,
                      -TPD_LIMIT_MHZ, TPD_LIMIT_MHZ, 1.0, "MHz", recompute=False,
                      visible_if=seeded,
                      help="ω_seed = ω_pump − ν_HF + δ. Navigates the curve instantly (no recompute)."),
            ParamSpec("temp_c", "Temperature", "Cell", 121.0,
                      60.0, 150.0, 1.0, "°C", visible_if=seeded),
            ParamSpec("cell_mm", "Cell length", "Cell", 12.5,
                      1.0, 100.0, 0.5, "mm", visible_if=seeded,
                      help="Vapor-cell length L. Enters the Maxwell-Bloch "
                           "propagation exp(M·L), so it recomputes the gain."),
            ParamSpec("pump_mw", "Pump power", "Beams", 600.0,
                      50.0, 1200.0, 10.0, "mW", visible_if=seeded),
            ParamSpec("probe_uw", "Seed / probe power", "Beams", 8.0,
                      1.0, 200.0, 1.0, "µW", visible_if=seeded),
            ParamSpec("seeded_angle_deg", "Pump-probe angle", "Beams",
                      SEEDED_PHASE_ANGLE_DEG, 0.0, 2.0, 0.05, "deg",
                      visible_if=seeded,
                      help="Pump-seed crossing angle θ. Real beam geometry: sets "
                           "the seeded-FWM longitudinal phase mismatch Δk_z (active "
                           "from Balanced fidelity up), so it recomputes the gain."),
            ParamSpec("loss_pct", "Loss after cell", "Detection & scaling", 5.5,
                      0.0, 50.0, 0.5, "%", visible_if=seeded,
                      help="Folds into eta = QE x (1 - loss)."),
            ParamSpec("line_strength", "Line-strength factor", "Detection & scaling", 0.74,
                      0.2, 5.0, 0.01, "×",
                      visible_if=seeded, advanced=True,
                      help="Dimensionless residual coupling calibration. The physical "
                           "macroscopic normalization — Rb85 D1 hyperfine Clebsch-Gordan "
                           "strengths × p_F/[2(2I+1)] (ground population × sublevel "
                           "degeneracy) — is computed from first principles in code; this "
                           "residual (0.74) is anchored to Sim et al. (Sci. Rep. 15, 7727 "
                           "(2025)): at Ultra fidelity it reproduces the measured gain "
                           "G_s≈14 and −7.8 dB squeezing at the paper's operating point. "
                           "The ~0.74 (vs 1.0) is the m-sublevel participation / geometry "
                           "reduction not yet derived from first principles. NOT an "
                           "experimentally tunable variable — lives under Advanced."),
            ParamSpec("biphoton_temp_c", "Temperature", "Cell & beams", 90.0,
                      30.0, 160.0, 1.0, "°C", visible_if=biphoton),
            ParamSpec("biphoton_cell_mm", "Cell length", "Cell & beams", 12.5,
                      1.0, 100.0, 0.5, "mm", visible_if=biphoton,
                      advanced=True),
            ParamSpec("pump_biphoton_uw", "Pump power", "Fields", 10.0,
                      0.1, 200.0, 0.1, "µW", visible_if=biphoton),
            ParamSpec("coupling_mw", "Coupling drive scale", "Fields", 1.0,
                      0.01, 50.0, 0.01, "×", visible_if=biphoton,
                      help="Relative coupling drive, √-scaled against the topology's "
                           "reference coupling power — dimensionless, not an absolute mW."),
            ParamSpec("pump_detuning_mhz", "Pump detuning", "Detunings", 0.0,
                      -2000.0, 2000.0, 10.0, "MHz", visible_if=biphoton),
            ParamSpec("two_photon_detuning_mhz", "Two-photon detuning", "Detunings", 0.0,
                      -2000.0, 2000.0, 10.0, "MHz", visible_if=biphoton,
                      help="Delta_p + Delta_c. Internally this sets the coupling "
                           "detuning relative to the pump detuning."),
            ParamSpec("coupling_detuning_mhz", "Coupling detuning", "Detunings", 0.0,
                      -2000.0, 2000.0, 10.0, "MHz", visible_if=biphoton,
                      hidden=True, recompute=False),
            ParamSpec("signal_angle_deg", "Signal collection angle", "Phase matching", 1.5,
                      0.0, 10.0, 0.1, "deg", visible_if=biphoton),
            ParamSpec("idler_angle_offset_deg", "Idler angle offset", "Phase matching",
                      0.0, -5.0, 5.0, 0.05, "deg", visible_if=biphoton,
                      help="Offset from the transverse phase-matched idler angle "
                           "derived from the selected topology and signal angle."),
            ParamSpec("idler_angle_deg", "Idler angle", "Phase matching",
                      transverse_matched_angle_deg(1529.37, 780.24, 1.5),
                      0.0, 10.0, 0.1, "deg", visible_if=biphoton,
                      hidden=True, recompute=False),
            ParamSpec("signal_side", "Signal side", "Phase matching", SIDE_PLUS,
                      choices=SIDE_CHOICES, visible_if=biphoton,
                      hidden=True,
                      help="Transverse collection side used in vector phase matching."),
            ParamSpec("idler_side", "Idler side", "Phase matching", SIDE_PLUS,
                      choices=SIDE_CHOICES, visible_if=biphoton,
                      hidden=True,
                      help="Opposite side flips the idler transverse wavevector."),
            ParamSpec("diamond_pump_nm", "Diamond pump wavelength", "Fields", 780.0,
                      300.0, 2000.0, 1.0, "nm", visible_if=diamond,
                      advanced=True),
            ParamSpec("diamond_coupling_nm", "Diamond coupling wavelength", "Fields", 776.0,
                      300.0, 2000.0, 1.0, "nm", visible_if=diamond,
                      advanced=True),
            ParamSpec("diamond_signal_nm", "Diamond signal wavelength", "Fields", 795.0,
                      300.0, 2000.0, 1.0, "nm", visible_if=diamond,
                      advanced=True),
            ParamSpec("diamond_idler_nm", "Diamond idler wavelength", "Fields", 761.702,
                      300.0, 2500.0, 0.001, "nm", visible_if=diamond,
                      advanced=True),
            ParamSpec("signal_eff_pct", "Signal efficiency", "Detection & scaling", 10.0,
                      0.1, 95.0, 0.1, "%", visible_if=biphoton,
                      advanced=True, recompute=False),
            ParamSpec("idler_eff_pct", "Idler efficiency", "Detection & scaling", 10.0,
                      0.1, 95.0, 0.1, "%", visible_if=biphoton,
                      advanced=True, recompute=False),
            ParamSpec("dark_signal_cps", "Signal background", "Detection & scaling", 2000.0,
                      0.0, 100000.0, 100.0, "cps", visible_if=biphoton,
                      advanced=True, recompute=False),
            ParamSpec("dark_idler_cps", "Idler background", "Detection & scaling", 2000.0,
                      0.0, 100000.0, 100.0, "cps", visible_if=biphoton,
                      advanced=True, recompute=False),
            ParamSpec("coincidence_window_ns", "Coincidence window", "Detection & scaling", 1.0,
                      0.01, 100.0, 0.01, "ns", visible_if=biphoton,
                      recompute=False),
            ParamSpec("timing_jitter_ns", "Timing jitter FWHM", "Detection & scaling", 0.55,
                      0.0, 5.0, 0.01, "ns", visible_if=biphoton,
                      advanced=True, recompute=False),
            ParamSpec("filter_bandwidth_mhz", "Filter bandwidth", "Detection & scaling", 300.0,
                      1.0, 5000.0, 1.0, "MHz", visible_if=biphoton,
                      recompute=False),
            ParamSpec("tau_max_ns", "Temporal window", "Numerics", 12.0,
                      1.0, 100.0, 1.0, "ns", visible_if=biphoton, advanced=True),
            ParamSpec("biphoton_velocity_step", "Velocity step", "Numerics", 2.0,
                      0.5, 20.0, 0.5, "m/s", visible_if=biphoton, advanced=True,
                      help="Maxwell velocity-grid step. The calibrated source model "
                           "uses it directly; the predictive model treats it as an "
                           "upper bound and auto-refines finer until the biphoton "
                           "width converges (a coarse step aliases the velocity-"
                           "class coherent sum)."),
            ParamSpec("resolution", "Model fidelity", "Numerics", FIDELITY_BALANCED,
                      choices=tuple(FWM_FIDELITY.keys()), advanced=True,
                      visible_if=seeded),
            ParamSpec("phase_detail", "Phase detail", "Phase matching", "Balanced",
                      choices=("Balanced", "Fine"), visible_if=biphoton,
                      advanced=True,
                      help="Balanced uses calibrated 1D phase matching; Fine adds "
                           "absolute diagnostics and a 2D signal-idler map."),
        ]

    def presets(self):
        return []

    def recommended_defaults(self, params):
        return {
            MODE_SEEDED: self._squeezing_defaults(),
            MODE_BIPHOTON: self._biphoton_defaults(params),
        }

    def _squeezing_defaults(self):
        return dict(mode=MODE_SEEDED, opd=0.9, tpd=-8.0, temp_c=121.0,
                    cell_mm=12.5, pump_mw=600.0, probe_uw=8.0, loss_pct=5.5,
                    line_strength=0.74, resolution=FIDELITY_BALANCED,
                    seeded_angle_deg=SEEDED_PHASE_ANGLE_DEG)

    def _biphoton_defaults(self, params):
        topology = params.get("topology", TOPOLOGY_RB87_TELECOM)
        base = dict(
            mode=MODE_BIPHOTON,
            topology=topology,
            biphoton_model=params.get("biphoton_model", BIPHOTON_PREDICTIVE),
            biphoton_cell_mm=12.5,
            pump_detuning_mhz=0.0,
            two_photon_detuning_mhz=0.0,
            coupling_detuning_mhz=0.0,
            signal_angle_deg=1.5,
            idler_angle_offset_deg=0.0,
            idler_angle_deg=transverse_matched_angle_deg(1529.37, 780.24, 1.5),
            signal_side=SIDE_PLUS,
            idler_side=SIDE_PLUS,
            signal_eff_pct=10.0,
            idler_eff_pct=10.0,
            dark_signal_cps=2000.0,
            dark_idler_cps=2000.0,
            coincidence_window_ns=1.0,
            filter_bandwidth_mhz=300.0,
            timing_jitter_ns=0.55,
            tau_max_ns=12.0,
            biphoton_velocity_step=2.0,
            phase_detail="Balanced",
        )
        if topology == TOPOLOGY_CS_BTW:
            base.update(cs_channel=params.get("cs_channel", CS_CHANNEL_917),
                        biphoton_temp_c=75.0, pump_biphoton_uw=20.0,
                        coupling_mw=1.0)
        elif topology == TOPOLOGY_DIAMOND:
            base.update(biphoton_temp_c=60.0, pump_biphoton_uw=20.0,
                        coupling_mw=1.0, diamond_pump_nm=780.0,
                        diamond_coupling_nm=776.0, diamond_signal_nm=795.0,
                        diamond_idler_nm=761.702)
        else:
            base.update(topology=TOPOLOGY_RB87_TELECOM, biphoton_temp_c=90.0,
                        pump_biphoton_uw=10.0, coupling_mw=1.0)
        base.update(_default_biphoton_geometry(base))
        return base

    def _biphoton_runtime_params(self, params):
        """Map the compact lab-facing controls onto the backend parameters."""
        out = dict(params)
        if out.get("mode", MODE_SEEDED) != MODE_BIPHOTON:
            return out

        topology = topology_from_params(out)
        signal = topology.field_map["signal"]
        idler = topology.field_map["idler"]
        signal_angle = float(out.get("signal_angle_deg", signal.angle_deg))
        out["signal_angle_deg"] = signal_angle

        if "two_photon_detuning_mhz" in out:
            two_det = float(out.get("two_photon_detuning_mhz", 0.0))
            pump_det = float(out.get("pump_detuning_mhz", 0.0))
            out["coupling_detuning_mhz"] = two_det - pump_det
        else:
            out["two_photon_detuning_mhz"] = (
                float(out.get("pump_detuning_mhz", 0.0))
                + float(out.get("coupling_detuning_mhz", 0.0))
            )

        if "idler_angle_offset_deg" in out:
            matched_idler = transverse_matched_angle_deg(
                signal.wavelength_nm, idler.wavelength_nm, signal_angle)
            out["idler_angle_deg"] = float(np.clip(
                matched_idler + float(out.get("idler_angle_offset_deg", 0.0)),
                0.0, 10.0))
        return out

    def info(self):
        return (
            "**85Rb D1 double-Lambda four-wave mixing.** The legacy seeded gain "
            "and intensity-difference squeezing remain regression-anchored to Sim "
            "*et al.*\n\n"
            "The propagation is linear in the undepleted pump, so the bare gain "
            "grows exponentially with density and would exceed the pump's energy "
            "budget at high T. A Manley-Rowe pump-depletion saturation "
            "((G_s−1)·P_seed, G_c·P_seed → P_pump/2) caps the gain at the energy-"
            "conservation bound; it is negligible in the validated regime. "
            "Balanced fidelity includes the reference 0.32 deg seeded phase "
            "mismatch; High fidelity also applies a chi-reused refractive "
            "correction and segmented propagation profile. Ultra adds a slow "
            "self-consistent phase refinement, dynamic segmented depletion, "
            "in-cell loss/noise, Gaussian overlap, and a visible Zeeman "
            "diagnostic correction.\n\n"
            "**Source model toggle.** *Predictive* solves the Doppler-averaged "
            "cascade biphoton amplitude from first principles: the two-photon "
            "denominator carries the Ω_c² Autler-Townes term (not a weak-coupling "
            "drive), the BTW is the collective velocity-class coherent sum with "
            "natural-linewidth decay (no injected lifetime), the source bandwidth "
            "comes from the waveform, and the wavelength-dependent width *ordering* "
            "(852-917 nm narrower than 852-795 nm) emerges. **Honest limits:** "
            "absolute ns-widths are approximate — the Cs cascade channels land "
            "within ~30 % but the extreme-wavelength-ratio 780/1529 nm telecom case "
            "over-weights the natural-decay tail; and the absolute pair rate stays "
            "*reference-anchored* with physical scaling (pump power, OD, phase "
            "matching), because the lumped 4-level model does not pin the absolute "
            "collection coefficient (the same reason the squeezing mode keeps a "
            "line-strength residual). An OD waveform-reshaping path (Du/Chen ρ̄, "
            "group-delay/precursor) is implemented but off by default — the lumped "
            "model overestimates it at high OD. *Calibrated* is the legacy "
            "reference-injected estimate. Full quantum-Langevin noise remains future "
            "work (see docs/checklist.json).\n\n"
            "**References**\n"
            "- G. Sim, H. Kim, H. S. Moon, *Sci. Rep.* **15**, 7727 (2025) "
            "(seeded 85Rb D1 double-Lambda gain & squeezing, regression anchor).\n"
            "- H. Kim, H. Jeong, H. S. Moon, [*Quantum Sci. Technol.* **9**, "
            "045006 (2024)](https://arxiv.org/abs/2402.06872) (Cs cascade BTW, Eq. 2).\n"
            "- H. Jeong, H. Kim, H. S. Moon, *Adv. Quantum Technol.* **7**, "
            "2300108 (2024) (87Rb telecom).\n"
            "- S. Du, J. Wen, M. H. Rubin, [*J. Opt. Soc. Am. B* **25**, C98 "
            "(2008)](https://arxiv.org/abs/0804.3981) (biphoton = nonlinear ⊛ "
            "linear response, Eq. 15).\n"
            "- Chen *et al.*, [*Phys. Rev. Research* **4**, 023132 (2022)]"
            "(https://arxiv.org/abs/2109.09062) (Doppler-averaged hot-vapor SFWM, Eq. 3-5).\n"
            "- J. Park, T. Jeong, H. S. Moon, [*Sci. Rep.* **10**, 16413 "
            "(2020)](https://www.nature.com/articles/s41598-020-73610-2) "
            "(cascade-type warm-atom biphoton waveform)."
        )

    def compute(self, params):
        if params.get("mode", MODE_SEEDED) == MODE_BIPHOTON:
            params = self._biphoton_runtime_params(params)
            topology = topology_from_params(params)
            return GenericFWMSolver(topology).compute_biphoton(params)
        center = branch_center_GHz(params["opd"], -1)
        fidelity = normalize_fidelity(params["resolution"])
        res = FWM_FIDELITY[fidelity]
        return compute_spectrum(
            params["opd"],
            T=params["temp_c"] + 273.15,
            P_pump=params["pump_mw"] * 1e-3,
            P_probe=params["probe_uw"] * 1e-6,
            line_strength=params["line_strength"],
            L=params["cell_mm"] * 1e-3,
            loss_frac=params["loss_pct"] / 100.0,
            coarse_points=res["coarse_points"], fine_points=0,
            scan_min=center - WINDOW_GHZ, scan_max=center + WINDOW_GHZ,
            velocity_step=res["velocity_step"],
            velocity_cutoff=res.get("velocity_cutoff", 3.0),
            phase_detail=res["phase_detail"],
            pump_probe_angle_deg=params.get("seeded_angle_deg", SEEDED_PHASE_ANGLE_DEG),
            model_fidelity=fidelity,
            branch=DEFAULT_BRANCH,
        )

    def observables(self, raw, params):
        if raw.get("kind") == "biphoton":
            params = self._biphoton_runtime_params(params)
            return self._biphoton_observables(raw, params)
        return self._seeded_observables(raw, params)

    def _seeded_observables(self, raw, params):
        import matplotlib.pyplot as plt

        tpd = params["tpd"]
        op = operating_point(raw, tpd, branch=-1)
        d_axis = (raw["probe_axis_GHz"] - raw["raman_center_minus_GHz"]) * 1e3

        fig, (axG, axS) = plt.subplots(2, 1, figsize=(8.5, 6.4), sharex=True)
        for ax in (axG, axS):
            ax.grid(alpha=0.3)
        axG.plot(d_axis, raw["G_s"], color="#1f77b4", lw=1.8)
        axG.axvline(tpd, color="crimson", ls="--", lw=1.2)
        axG.axhline(1.0, color="black", lw=0.6)
        axG.scatter([tpd], [op["G_s"]], color="crimson", zorder=5)
        axG.set_ylabel("Seed / probe gain G_s")
        axG.set_title(f"Delta = {params['opd']:.1f} GHz,  T = {params['temp_c']:.0f} C,  "
                      f"eta = {raw['eta']:.3f}")
        if np.nanmax(raw["G_s"]) > 50:
            axG.set_yscale("log")
        axS.plot(d_axis, raw["S_dB"], color="#2ca02c", lw=1.8)
        axS.axvline(tpd, color="crimson", ls="--", lw=1.2)
        axS.axhline(0.0, color="black", lw=0.6)
        axS.scatter([tpd], [op["S_dB"]], color="crimson", zorder=5)
        axS.set_ylabel("Intensity-difference\nsqueezing  [dB]")
        axS.set_xlabel("Two-photon detuning delta [MHz]   (probe on the - Raman branch)")
        axS.set_xlim(-TPD_LIMIT_MHZ, TPD_LIMIT_MHZ)
        fig.tight_layout()

        metrics = [
            dict(label="Seed / probe gain  G_s", value=f"{op['G_s']:.2f}",
                 help="Power gain of the seeded probe through the cell. Absolute "
                      "scale uses the first-principles macroscopic normalization "
                      "(Clebsch-Gordan strengths × p_F/[2(2I+1)]); the Advanced "
                      "Line-strength factor is a residual ×1.0 calibration."),
            dict(label="Squeezing", value=f"{op['S_dB']:.2f} dB",
                 delta="below shot noise" if op["S_dB"] < 0 else "above shot noise",
                 delta_color="inverse"),
            dict(label="Conjugate gain  G_c", value=f"{op['G_c']:.2f}",
                 help="Generated conjugate power gain (drives the twin-beam squeezing)."),
        ]
        cap = raw.get("pump_depletion_cap", float("inf"))
        small_signal = raw.get("G_s_smallsignal_peak", op["G_s"])
        depletion_limited = small_signal > 1.1 * cap
        phase_rows = ""
        if raw.get("phase_detail", PHASE_LEGACY) != PHASE_LEGACY:
            phase_rows = (
                f"| Model fidelity | {raw.get('model_fidelity', raw['phase_detail'])} |\n"
                f"| Phase detail | {raw['phase_detail']} |\n"
                f"| Pump-probe angle | {raw['pump_probe_angle_deg']:.3f} deg |\n"
                f"| Operating Delta k_z | {op.get('delta_k_z', np.nan):.3e} 1/m |\n"
                f"| Vacuum Delta k_z | {op.get('delta_k_z_vacuum', np.nan):.3e} 1/m |\n"
                f"| Propagation segments | {raw.get('phase_segments', 1)} |\n"
                f"| Segment absorption OD estimate | {raw.get('segment_absorption_od', 0.0):.3f} |\n"
            )
            if raw.get("phase_detail") == PHASE_ULTRA:
                phase_rows += (
                    f"| Ultra fixed-point iterations | {raw.get('ultra_phase_iterations', 0)} |\n"
                    f"| Ultra final Delta-k change | {raw.get('ultra_phase_max_change', 0.0):.3e} 1/m |\n"
                    f"| Dynamic depletion | {raw.get('ultra_dynamic_depletion', False)} |\n"
                    f"| Min pump remaining | {raw.get('ultra_pump_remaining_min', np.nan):.3e} W |\n"
                    f"| In-cell loss/noise | {raw.get('ultra_in_cell_loss_noise', False)} |\n"
                    f"| Min Gaussian overlap | {raw.get('ultra_spatial_overlap_min', 1.0):.4f} |\n"
                    f"| Zeeman status | {raw.get('zeeman_status', 'inactive')} |\n"
                    f"| Zeeman correction | {raw.get('zeeman_correction', 1.0):.4f} |\n"
                )
        derived = (
            f"| Quantity | Value |\n|---|---|\n"
            f"| N(85Rb) | {raw['N_atoms']:.3e} /m³ |\n"
            f"| σ_v (1-D thermal) | {raw['sigma_v']:.1f} m/s |\n"
            f"| Velocity classes | {raw['n_velocity']} |\n"
            f"| Ω_pump / 2π | {raw['Op_A_2pi_GHz']:.3f} GHz |\n"
            f"| Ω_seed / 2π | {raw['Os_2pi_MHz']:.3f} MHz |\n"
            f"| (−) Raman line (probe axis) | {raw['raman_center_minus_GHz']:.3f} GHz |\n"
            f"| Detection η = QE·(1−loss) | {raw['eta']:.4f} |\n"
            f"| Coupling norm p_F/[2(2I+1)] (first-principles) | "
            f"{raw.get('coupling_norm', float('nan')):.4f} |\n"
            f"| Line-strength residual (Advanced, ≈1.0) | {params['line_strength']:.3f}× |\n"
            f"| Effective coupling scale | {raw.get('effective_line_strength', float('nan')):.4f} |\n"
            f"| Pump-depletion cap on G_s (Manley-Rowe) | {cap:.3e} |\n"
            f"| Small-signal peak G_s (pre-saturation) | {small_signal:.3e} |\n"
            f"| Operating probe detuning | {op['probe_GHz']:.4f} GHz |\n"
            + phase_rows
            + "\n"
            + ("⚠️ **Pump-depletion limited:** the small-signal gain exceeds what "
               "the pump can supply (Manley-Rowe), so the shown gain is capped by "
               "energy conservation. Lower T / raise seed power to stay in the "
               "linear regime.\n\n" if depletion_limited else "")
            + f"Cell L = {raw.get('cell_length_m', L_CELL)*1e3:.1f} mm · "
            f"Fixed: pump w₀ {W_PUMP*1e6:.0f} µm · "
            f"seed w₀ {W_PROBE*1e6:.0f} µm · QE {QE_DETECTOR*100:.2f}% · "
            f"responsivity {RESPONSIVITY_AW} A/W @ 795 nm · pump⊥probe at PBS."
        )
        return {
            "metrics": metrics,
            "figure": fig,
            "tables": [
                {"title": "Derived quantities", "markdown": derived},
            ],
        }

    def _biphoton_observables(self, raw, params):
        import matplotlib.pyplot as plt

        topo = raw["topology"]
        predictive = raw.get("predictive", False)
        velocity_converged = raw.get("velocity_converged", True)
        # Predictive: g²_SI(τ) comes from the computed waveform |ψ|² and the
        # physical accidentals (no target-g² forcing). Calibrated: legacy
        # added-accidental anchoring to the reference g² peak.
        stats = observables.biphoton_stats(
            raw["tau_axis_ns"], raw["psi_tau"], raw["pair_rate_cps"],
            signal_eff=params["signal_eff_pct"] / 100.0,
            idler_eff=params["idler_eff_pct"] / 100.0,
            dark_signal_cps=params["dark_signal_cps"],
            dark_idler_cps=params["dark_idler_cps"],
            coincidence_window_ns=params["coincidence_window_ns"],
            timing_jitter_ns=params["timing_jitter_ns"],
            filter_bandwidth_mhz=params["filter_bandwidth_mhz"],
            source_bandwidth_mhz=raw["source_bandwidth_mhz"],
            target_g2_peak=None if predictive else topo.target_g2_peak,
        )

        fig, (axG2, axPM) = plt.subplots(2, 1, figsize=(8.5, 6.4))
        axG2.plot(stats["tau_axis_ns"], stats["g2_SI_tau"], color="#1f77b4", lw=1.8)
        axG2.axhline(2.0, color="black", lw=0.7, ls=":")
        axG2.set_ylabel(r"$g^{(2)}_{SI}(\tau)$")
        title_tag = ("predictive (waveform + anchored rate)" if predictive
                     else "calibrated spontaneous SFWM estimate")
        axG2.set_title(f"{topo.label}: {title_tag}")
        axG2.grid(alpha=0.3)

        axPM.plot(raw["angle_axis_deg"], raw["phase_matching"], color="#2ca02c", lw=1.8)
        axPM.axvline(params["idler_angle_deg"], color="crimson", lw=1.1, ls="--")
        axPM.set_xlabel("Idler collection angle [deg]")
        axPM.set_ylabel(r"$\mathrm{sinc}^2(|\Delta\mathbf{k}| L / 2)$")
        axPM.grid(alpha=0.3)
        fig.tight_layout()

        amp = np.abs(raw["source_v"])
        amp = amp / np.nanmax(amp) if np.nanmax(amp) > 0 else amp
        phase = np.unwrap(np.angle(raw["source_v"]))
        figV, (axA, axP) = plt.subplots(2, 1, figsize=(8.5, 6.0), sharex=True)
        axA.plot(raw["v_grid"], amp, color="#ff7f0e", lw=1.5)
        axA.set_ylabel("Velocity source amplitude")
        axA.grid(alpha=0.3)
        axP.plot(raw["v_grid"], phase, color="#9467bd", lw=1.3)
        axP.set_xlabel("Atomic velocity [m/s]")
        axP.set_ylabel("Velocity source phase [rad]")
        axP.grid(alpha=0.3)
        figV.tight_layout()

        extra_figures = [("Velocity-class coherent source", figV)]
        pm2d = raw.get("phase_matching_2d")
        if pm2d is not None:
            figPM2, ax2 = plt.subplots(1, 1, figsize=(7.2, 5.2))
            im = ax2.imshow(
                pm2d.T, origin="lower", aspect="auto",
                extent=[
                    raw["signal_angle_axis_deg"][0],
                    raw["signal_angle_axis_deg"][-1],
                    raw["idler_angle_axis_2d_deg"][0],
                    raw["idler_angle_axis_2d_deg"][-1],
                ],
                cmap="viridis", vmin=0.0, vmax=1.0)
            ax2.scatter([params["signal_angle_deg"]], [params["idler_angle_deg"]],
                        color="crimson", s=28, zorder=5)
            ax2.set_xlabel("Signal angle [deg]")
            ax2.set_ylabel("Idler angle [deg]")
            ax2.set_title("2D phase-matching acceptance")
            figPM2.colorbar(im, ax=ax2,
                            label=r"$\mathrm{sinc}^2(|\Delta\mathbf{k}| L / 2)$")
            figPM2.tight_layout()
            extra_figures.append(("2D phase matching", figPM2))

        if predictive:
            g2_help = ("Peak signal-idler cross-correlation from the computed "
                       "waveform |ψ|² and physical accidentals (not forced).")
            rate_help = ("Pair rate is reference-anchored with physical scaling "
                         "(pump power, OD, phase matching) — the absolute collection "
                         "coefficient is not derivable in the lumped 4-level model.")
            fwhm_help = (
                "FWHM of the predicted biphoton temporal waveform, computed on a "
                "Maxwell velocity grid auto-refined until the width converges (the "
                "velocity-class coherent sum aliases on a coarse step). Width "
                "ordering vs wavelength is physical; the absolute ns-scale still "
                "carries the lumped-model per-source calibration uncertainty."
                if velocity_converged else
                "Predicted waveform width did NOT converge: the velocity-grid "
                "refinement hit its point cap, so the value is qualitative only "
                "(shape/ordering still indicative). Raise the temporal window or "
                "report it as unconverged.")
            status_value = (f"predictive · {raw.get('regime', '—')}"
                            if velocity_converged
                            else f"predictive · {raw.get('regime', '—')} · "
                                 "width unconverged")
            status = dict(
                label="Model status", value=status_value,
                help="Waveform shape, decay, bandwidth and the wavelength-dependent "
                     "width ordering are solved from first principles (Ω_c² Autler-"
                     "Townes term, velocity-class coherent sum, natural-linewidth "
                     "decay, OD reshaping) on a velocity grid auto-refined to FWHM "
                     "convergence; pair rate is reference-anchored. Du regime: "
                     "damped-Rabi vs group-delay.")
        else:
            g2_help = ("Peak normalized signal-idler cross-correlation (calibrated to "
                       "the reference target, not predicted).")
            rate_help = "Reference-calibrated generated pair-rate estimate."
            fwhm_help = "FWHM of the modeled biphoton temporal waveform."
            status = dict(
                label="Model status", value="calibrated · non-predictive",
                help="Reference-calibrated source estimate: pair rate and g2 are "
                     "anchored to the reference numbers by construction. Switch the "
                     "Source model to Predictive for the first-principles waveform.")
        metrics = [
            status,
            dict(label="g2_SI peak", value=f"{stats['g2_peak']:.2f}", help=g2_help),
            dict(label="CAR", value=f"{stats['CAR']:.1f}",
                 help="True coincidence divided by accidental coincidence."),
            dict(label="Pair rate", value=f"{stats['pair_rate_cps']:.1f} cps",
                 help=rate_help),
            dict(label="BTW FWHM",
                 value=(f"{stats['fwhm_ns']:.2f} ns"
                        if (velocity_converged or not predictive) else "unconverged"),
                 help=fwhm_help),
            dict(label="Source bandwidth", value=f"{raw['source_bandwidth_mhz']:.0f} MHz",
                 help=("Spectral FWHM from the waveform (predictive)." if predictive
                       else "Reference source bandwidth.")),
            dict(label="Phase match", value=f"{raw['phase_match_weight']:.3f}",
                 help="Vector sinc^2 phase-matching collection weight."),
        ]

        field_rows = "".join(
            f"| {f.role} | {raw['topology'].levels[f.lower].name} -> "
            f"{raw['topology'].levels[f.upper].name} | {f.wavelength_nm:.2f} nm | "
            f"{f.angle_deg:.2f} deg | {_side_label(f.side_sign) if f.side_sign else '0'} |\n"
            for f in raw["fields"]
        )
        pm_warning = (
            "| Warning | Vector phase match < 1e-3; collection geometry is not "
            "physically phase matched. |\n"
            if raw["phase_match_weight"] < 1e-3 else ""
        )
        topology_table = (
            f"| Quantity | Value |\n|---|---|\n"
            f"| Topology | {topo.label} |\n"
            f"| Family | {topo.family} |\n"
            f"| Density | {raw['density']:.3e} /m^3 |\n"
            f"| Cell length | {raw['cell_length_m']*1e3:.2f} mm |\n"
            f"| Delta k z relative | {raw['delta_k_z_relative']:.3e} 1/m |\n"
            f"| Delta k z absolute | {raw['delta_k_z_absolute']:.3e} 1/m |\n"
            f"| Delta k x transverse | {raw['delta_k_x']:.3e} 1/m |\n"
            f"| Delta k vector | {raw['delta_k_vector']:.3e} 1/m |\n"
            f"| Vector phase match | {raw['phase_match_weight']:.3f} |\n"
            f"| Longitudinal phase match | {raw['phase_match_weight_longitudinal']:.3f} |\n"
            f"| Vacuum phase match | {raw['phase_match_weight_absolute']:.3f} |\n"
            + pm_warning +
            f"| Phase detail | {raw['phase_detail']} |\n"
            f"| Energy mismatch | {raw['energy_mismatch_hz']/1e6:.3f} MHz |\n"
            f"| Residual two-photon Doppler k | {raw['residual_two_photon_k']:.3e} 1/m |\n"
            f"| Velocity step{' (auto-refined)' if predictive else ''} | "
            f"{raw.get('velocity_step_used', float('nan')):.3f} m/s"
            f"{'' if velocity_converged else ' · UNCONVERGED'} |\n\n"
            f"| Field | Transition | Wavelength | Angle | Side |\n|---|---|---:|---:|---:|\n"
            + field_rows
        )
        detection_table = (
            "Calibrated source estimate with detector/background model:\n\n"
            f"| Quantity | Value |\n|---|---|\n"
            f"| Signal singles | {stats['singles_signal_cps']:.2f} cps |\n"
            f"| Idler singles | {stats['singles_idler_cps']:.2f} cps |\n"
            f"| True coincidence | {stats['coincidence_cps']:.3f} cps |\n"
            f"| Accidental coincidence | {stats['accidental_cps']:.3e} cps |\n"
            f"| Raw accidental before reference calibration | {stats['raw_accidental_cps']:.3e} cps |\n"
            f"| Added unmodelled accidental/background | {stats['added_accidental_cps']:.3e} cps |\n"
            f"| Raw g2 peak before reference calibration | {stats['raw_g2_peak']:.2f} |\n"
            f"| Heralding signal | {stats['heralding_signal']:.3e} |\n"
            f"| Heralding idler | {stats['heralding_idler']:.3e} |\n"
            f"| Cauchy-Schwarz R | {stats['cauchy_schwarz_R']:.2f} |\n"
            f"| Filter transmission estimate | {stats['filter_transmission']:.3f} |\n\n"
            f"{raw['notes']}"
        )
        validation_table = self._reference_validation_table(raw, params, stats)
        val_title = ("Reference comparison (predictive waveform · anchored rate)"
                     if predictive else
                     "Reference reproduction (calibrated · non-predictive)")
        return {
            "metrics": metrics,
            "figure": fig,
            "figures": extra_figures,
            "tables": [
                {"title": "Generic SFWM topology", "markdown": topology_table},
                {"title": "Biphoton detection model", "markdown": detection_table},
                {"title": val_title, "markdown": validation_table},
            ],
        }

    def _reference_validation_table(self, raw, params, stats):
        topo = raw["topology"]
        predictive = raw.get("predictive", False)
        velocity_converged = raw.get("velocity_converged", True)

        def verdict(ok, kind="physical"):
            # Calibration anchors are matched to the reference *by construction*
            # (the reference number is injected into the model), so a "PASS" there
            # would be circular — label it honestly instead. Predictive waveform
            # quantities are computed but absolute-approximate, so they are flagged
            # "predicted" rather than given a pass/fail they would often fail.
            if kind == "calibrated":
                return "by construction"
            if kind == "predicted":
                return "predicted (approx)"
            return "PASS" if ok else "CHECK"

        def row(name, calc, ref, ok, note="", kind="physical"):
            return f"| {name} | {calc} | {ref} | {verdict(ok, kind)} | {note} |\n"

        rows = []
        if topo.name == TOPOLOGY_RB87_TELECOM:
            pump_mw = max(params["pump_biphoton_uw"] * 1e-3, 1e-30)
            rate_per_mw = stats["pair_rate_cps"] / pump_mw
            rows.append(row(
                "Pair rate / pump", f"{rate_per_mw:.0f} cps/mW",
                "38000 cps/mW", abs(rate_per_mw / 38000.0 - 1.0) < 0.15,
                "anchored: pair rate is scaled from this reference number", kind="calibrated"))
            rows.append(row(
                "g2 peak", f"{stats['g2_peak']:.2f} (raw {stats['raw_g2_peak']:.1f})",
                "44(3)", abs(stats["g2_peak"] - 44.0) <= 3.0,
                ("predicted from waveform |ψ|² + physical accidentals"
                 if predictive else
                 "anchored: forced to the target via added-accidental calibration"),
                kind="predicted" if predictive else "calibrated"))
            if predictive and not velocity_converged:
                rows.append(row(
                    "BTW FWHM", "unconverged",
                    "0.56(4) ns", False,
                    "velocity-grid auto-refine hit the point cap — absolute width "
                    "not converged, qualitative only (shape/ordering still indicative)",
                    kind="predicted"))
            else:
                rows.append(row(
                    "BTW FWHM", f"{stats['fwhm_ns']:.3f} ns",
                    "0.56(4) ns", abs(stats["fwhm_ns"] - 0.56) <= 0.04,
                    ("predicted waveform (velocity-converged); absolute ns carries "
                     "per-source calibration uncertainty, ordering physical"
                     if predictive else "model waveform + detector jitter"),
                    kind="predicted" if predictive else "physical"))
            rows.append(row(
                "OD estimate", f"{raw['od_estimate']:.1f}",
                "112(3)", abs(raw["od_estimate"] - 112.0) <= 3.0,
                "anchored: density/cell scaling of the reference OD", kind="calibrated"))
            rows.append(row(
                "Bandwidth setting", f"{params['filter_bandwidth_mhz']:.0f} MHz",
                "about 300 MHz", abs(params["filter_bandwidth_mhz"] - 300.0) <= 40.0,
                "user filter-bandwidth setting check"))
        elif topo.name == TOPOLOGY_CS_BTW:
            other_channel = CS_CHANNEL_795 if params.get("cs_channel") == CS_CHANNEL_917 else CS_CHANNEL_917
            other_params = dict(params)
            other_params["cs_channel"] = other_channel
            other_raw = GenericFWMSolver(topology_from_params(other_params)).compute_biphoton(other_params)
            other_stats = observables.biphoton_stats(
                other_raw["tau_axis_ns"], other_raw["psi_tau"], other_raw["pair_rate_cps"],
                signal_eff=params["signal_eff_pct"] / 100.0,
                idler_eff=params["idler_eff_pct"] / 100.0,
                dark_signal_cps=params["dark_signal_cps"],
                dark_idler_cps=params["dark_idler_cps"],
                coincidence_window_ns=params["coincidence_window_ns"],
                timing_jitter_ns=params["timing_jitter_ns"],
                filter_bandwidth_mhz=params["filter_bandwidth_mhz"],
                source_bandwidth_mhz=other_raw["source_bandwidth_mhz"],
                target_g2_peak=other_raw["topology"].target_g2_peak,
            )
            ratio = max(stats["fwhm_ns"], other_stats["fwhm_ns"]) / max(
                min(stats["fwhm_ns"], other_stats["fwhm_ns"]), 1e-30)
            rows.append(row(
                "BTW width ratio", f"{ratio:.2f}",
                "about 3", abs(ratio - 3.0) <= 0.5,
                ("predicted ordering (917 narrower than 795); absolute ratio approximate"
                 if predictive else
                 "medium model only; full Cs BTW theory is not yet included"),
                kind="predicted" if predictive else "physical"))
            rows.append(row(
                "OD estimate", f"{raw['od_estimate']:.1f}",
                "about 10", abs(raw["od_estimate"] - 10.0) <= 2.0,
                "anchored: density/cell scaling of the reference OD", kind="calibrated"))
        else:
            rows.append(row(
                "Reference validation", "generic diamond template",
                "no paper anchor", False,
                "configure wavelengths manually; no validated default"))

        rows.append(row(
            "Phase matching", f"{raw['phase_match_weight']:.3f}",
            "> 0.90", raw["phase_match_weight"] > 0.90,
            "vector sinc^2(|Delta k| L / 2)"))
        rows.append(row(
            "Energy conservation", f"{raw['energy_mismatch_hz']/1e6:.3f} MHz",
            "near 0 MHz", abs(raw["energy_mismatch_hz"]) < 1e6,
            "wavelength bookkeeping"))

        if predictive:
            intro = (
                "**Predictive waveform · reference-anchored rate.** The biphoton "
                "amplitude (BTW shape, decay, bandwidth, OD reshaping, and the "
                "wavelength-dependent width *ordering*) is solved from first "
                "principles — rows tagged *predicted (approx)* are computed, but the "
                "absolute ns-widths are approximate (closed-form over-sensitivity in "
                "this regime, as the source papers fit Rabi/dephasing per source). "
                "The pair rate stays *by construction* (reference-anchored: the "
                "absolute collection coefficient is not derivable in the lumped "
                "4-level model, like the squeezing line-strength residual). The "
                "*physical* rows — phase matching and energy conservation — are "
                "genuine PASS/CHECK. Full quantum-Langevin noise is still future "
                "work (see docs/checklist.json).\n\n")
        else:
            intro = (
                "**Reference reproduction — calibrated, non-predictive.** Rows tagged "
                "*by construction* are anchored to the reference number (it is injected "
                "into the model), so agreement there is not an independent validation. "
                "Only the *physical* rows — phase matching and energy conservation, "
                "computed from the geometry/wavelength bookkeeping — are genuine "
                "PASS/CHECK tests. Switch the Source model to Predictive for the "
                "first-principles waveform. Full quantum-Langevin noise and full Cs "
                "BTW theory are not modelled (see docs/checklist.json).\n\n")
        return (
            intro
            + "| Check | Calculated | Reference | Verdict | Note |\n|---|---:|---:|---|---|\n"
            + "".join(rows)
        )

    def extra_views(self):
        def _compute_full(params):
            return full_spectrum(
                params["opd"], params["temp_c"] + 273.15,
                params["pump_mw"], params["probe_uw"], params["line_strength"],
                params["loss_pct"], L=params["cell_mm"] * 1e-3)

        def _render_full(full):
            import matplotlib.pyplot as plt
            figF, (aG, aS) = plt.subplots(2, 1, figsize=(8.5, 6.4), sharex=True)
            for ax in (aG, aS):
                ax.grid(alpha=0.3)
            styles = {
                "minus": dict(color="#1f77b4", label="minus Raman branch"),
                "plus": dict(color="#ff7f0e", label="plus Raman branch"),
            }
            for key, style in styles.items():
                spec = full[key]
                aG.plot(spec["probe_axis_GHz"], spec["G_s"], lw=1.4, **style)
                aS.plot(spec["probe_axis_GHz"], spec["S_dB"], lw=1.4, **style)
            aG.axhline(1.0, color="black", lw=0.6)
            aG.set_ylabel("Seed / probe gain G_s")
            if max(np.nanmax(full[key]["G_s"]) for key in styles) > 50:
                aG.set_yscale("log")
            aS.axhline(0.0, color="black", lw=0.6)
            aS.set_ylabel("Squeezing [dB]")
            aS.set_xlabel("Probe detuning from F=2 -> F'=3 [GHz]")
            for a in (aG, aS):
                a.axvline(full["D_GHz"], color="gray", ls=":", lw=0.8)
            aG.legend(fontsize=9)
            aS.legend(fontsize=9)
            figF.tight_layout()
            return figF

        return [ExtraView(
            key="Full −8…12 GHz probe scan (slow, both Raman branches)",
            description="The focused view zooms on the (−) Raman line. This runs the "
                        "original wide scan showing both branches.",
            compute=_compute_full, render=_render_full,
        )]


def full_spectrum(D_GHz, T_K, P_pump_mW, P_probe_uW, line_strength, loss_pct,
                  L=L_CELL):
    """Wide scan with the two Raman channels calculated independently.

    The ∓ branches are independent pure solves, so run them concurrently. With
    BLAS pinned to one thread per branch (see `core.blas_single_thread`) the two
    threads occupy two cores instead of contending — roughly halves this view.
    """
    common = dict(
        T=T_K, P_pump=P_pump_mW * 1e-3, P_probe=P_probe_uW * 1e-6,
        line_strength=line_strength, L=L, loss_frac=loss_pct / 100.0,
        coarse_points=301, fine_points=401, velocity_step=2.0)

    def _branch(b):
        with blas_single_thread():
            return compute_spectrum(D_GHz, branch=b, **common)

    with ThreadPoolExecutor(max_workers=2) as ex:
        fut = {b: ex.submit(_branch, b) for b in (-1, +1)}
        minus, plus = fut[-1].result(), fut[+1].result()
    return {"D_GHz": D_GHz, "minus": minus, "plus": plus}
