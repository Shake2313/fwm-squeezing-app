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

from .. import atoms, constants, doppler, hyperfine, observables, species
from ..constants import K_VEC, OMEGA_HF, OMEGA_EXCITED_HF, rabi_freq
from ..core import blas_single_thread, build_liouvillian, comm_super, floquet_solve
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

# =========================================================
# FWM experiment configuration (cell, beams, detection, scan)
# =========================================================
L_CELL = 12.5e-3
W_PUMP = 530e-6
W_PROBE = 330e-6
P_PUMP, P_PROBE = 600e-3, 10e-6
T_CELL = 394.15

QE_DETECTOR = 0.9047
RESPONSIVITY_AW = 0.58
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


# =========================================================
# Generic SFWM / biphoton topology layer
# =========================================================
NM = 1e-9
MHZ_ANG = 2 * np.pi * 1e6

MODE_SEEDED = "Seeded gain / squeezing"
MODE_BIPHOTON = "Spontaneous biphoton"

TOPOLOGY_RB87_TELECOM = "cascade_rb87_telecom"
TOPOLOGY_CS_BTW = "cascade_cs_btw"
TOPOLOGY_DIAMOND = "diamond_generic"
CS_CHANNEL_917 = "6D5/2: 852-917 nm"
CS_CHANNEL_795 = "8S1/2: 852-795 nm"


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


def _field_with(field, *, wavelength_nm=None, detuning_mhz=None, rabi_mhz=None,
                angle_deg=None):
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


def phase_matching_weight(delta_k, L):
    """sinc^2(delta_k L / 2), normalized to 1 at perfect phase matching."""
    x = 0.5 * np.asarray(delta_k, dtype=float) * L
    out = np.ones_like(x, dtype=float)
    mask = np.abs(x) > 1e-12
    out[mask] = (np.sin(x[mask]) / x[mask]) ** 2
    return out


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
        FieldSpec("signal", 2, 3, 1529.37, phase_sign=-1.0, angle_deg=1.5),
        FieldSpec("idler", 3, 0, 780.24, phase_sign=+1.0, angle_deg=1.5),
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
        FieldSpec("signal", 2, 3, coupling_nm, phase_sign=-1.0, angle_deg=1.5),
        FieldSpec("idler", 3, 0, pump_nm, phase_sign=+1.0, angle_deg=1.5),
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
    idler_default = 1.0 / (1.0 / pump_nm + 1.0 / coupling_nm - 1.0 / signal_nm)
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
        FieldSpec("signal", 3, 1, signal_nm, phase_sign=-1.0, angle_deg=2.0),
        FieldSpec("idler", 3, 2, idler_nm, phase_sign=-1.0, angle_deg=2.0),
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


class GenericFWMSolver:
    """Reference-calibrated V1 engine for generic SFWM biphoton estimates."""

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
            ))
        return tuple(out)

    def compute_biphoton(self, params):
        T = params.get("biphoton_temp_c", self.topology.default_temp_c) + 273.15
        L = params.get("biphoton_cell_mm", self.topology.default_cell_mm) * 1e-3
        fields = self._fields_from_params(params)
        fmap = {f.role: f for f in fields}
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
        pm_weight = float(phase_matching_weight(np.array([delta_k]), L)[0])
        density = species.number_density(iso, T)
        default_density = species.number_density(iso, self.topology.default_temp_c + 273.15)
        if self.topology.reference_od is not None:
            od_estimate = (self.topology.reference_od
                           * density / max(default_density, 1e-30)
                           * L / max(self.topology.default_cell_mm * 1e-3, 1e-30))
        else:
            od_estimate = np.nan
        pump_mw = params.get("pump_biphoton_uw", self.topology.default_pump_uw) * 1e-3
        coupling_scale = max(params.get("coupling_mw", self.topology.default_coupling_mw),
                             0.0) / max(self.topology.default_coupling_mw, 1e-12)
        pair_rate = (self.topology.pair_rate_cps_per_mw * pump_mw
                     * math.sqrt(max(coupling_scale, 0.0)) * pm_weight)
        pair_rate *= math.sqrt(max(density, 1e-30)
                               / max(species.number_density(iso, self.topology.default_temp_c + 273.15), 1e-30))

        lower_gamma = max(self.topology.levels[1].gamma_mhz, 0.1) * MHZ_ANG
        upper_gamma = max(self.topology.levels[2].gamma_mhz, 0.1) * MHZ_ANG
        pump_det = params.get("pump_detuning_mhz", 0.0) * MHZ_ANG
        two_det = (params.get("pump_detuning_mhz", 0.0)
                   + params.get("coupling_detuning_mhz", 0.0)) * MHZ_ANG
        lower = lower_gamma / 2.0 + 1j * (pump_det - pump.direction * pump.k * v_grid)
        residual_k = pump.direction * pump.k + coupling.direction * coupling.k
        upper = upper_gamma / 2.0 + 1j * (two_det - residual_k * v_grid)
        drive = (pump.rabi_mhz * coupling.rabi_mhz) * (MHZ_ANG ** 2)
        source_v = weights * drive / (lower * upper)
        source_v *= math.sqrt(max(pm_weight, 0.0))

        tau_axis_ns = np.linspace(0.0, params.get("tau_max_ns", 12.0), 481)
        tau_s = tau_axis_ns * 1e-9
        phase_k = abs(idler.k)
        coherent = np.exp(1j * phase_k * v_grid[:, None] * tau_s[None, :])
        psi_tau = (source_v[:, None] * coherent).sum(axis=0)
        psi_tau *= np.exp(-tau_axis_ns / max(self.topology.emission_decay_ns, 1e-12))
        if np.max(np.abs(psi_tau)) > 0:
            psi_tau = psi_tau / np.max(np.abs(psi_tau))

        angle_axis = np.linspace(max(idler.angle_deg - 4.0, 0.0),
                                 idler.angle_deg + 4.0, 181)
        angle_dk = np.array([
            phase_mismatch(fields, idler_angle_deg=a,
                           reference_delta_k=self.topology.reference_delta_k)
            for a in angle_axis
        ])
        acceptance = phase_matching_weight(angle_dk, L)

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
            "delta_k": float(delta_k),
            "phase_match_weight": pm_weight,
            "energy_mismatch_hz": float(energy_mismatch_hz(fields)),
            "pair_rate_cps": float(pair_rate),
            "density": float(density),
            "od_estimate": float(od_estimate),
            "source_bandwidth_mhz": float(self.topology.reference_bandwidth_mhz or 300.0),
            "temperature_K": float(T),
            "cell_length_m": float(L),
            "residual_two_photon_k": float(residual_k),
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
def chi_matrix_table(Op_A, Op_B, Os_ref, Oc_ref, delta_axis, Delta_eff_axis, branch):
    """
    Two solves per probe-detuning point to extract (χ̄_ss, χ̄_cs, χ̄_sc, χ̄_cc)
    on a 2-D (δ, Δ_eff) grid. Returns each array (n_delta, n_deff), complex.
    """
    probe_ground = G2 if branch == -1 else G1
    conj_ground = G1 if branch == -1 else G2
    n_d = delta_axis.size
    n_de = Delta_eff_axis.size

    chi_ss = np.zeros((n_d, n_de), dtype=complex)
    chi_cs = np.zeros((n_d, n_de), dtype=complex)
    chi_sc = np.zeros((n_d, n_de), dtype=complex)
    chi_cc = np.zeros((n_d, n_de), dtype=complex)

    Cp_no_c, Cm_no_c = sideband_template(Op_A, Op_B, 0.0, branch)      # solve 1
    Cp_c, Cm_c = sideband_template(Op_A, Op_B, Oc_ref, branch)        # solve 2

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
                     coarse_points=None, fine_points=None, window_mhz=None,
                     scan_min=None, scan_max=None,
                     velocity_step=None, velocity_cutoff=None,
                     branch=DEFAULT_BRANCH, branches=None):
    """Full pipeline for one one-photon detuning Δ = 2π·D_GHz·1e9 (see README)."""
    branch = _single_branch(branch, branches)
    if line_strength is None:
        line_strength = constants.LINE_STRENGTH_FACTOR
    eta = qe * (1.0 - loss_frac)

    Op_A = rabi_freq(P_pump, w_pump)
    Op_B = Op_A
    Os = rabi_freq(P_probe, w_probe)
    Os_ref = Os
    Oc_ref = Os                              # χ̄ is independent of |Ω_ref|

    N_atoms = atoms.rb85_density(T)
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

    G_s, G_c, _ = observables.gain_from_chi(
        chi_ss_avg, chi_sc_avg, chi_cs_avg, chi_cc_avg,
        K_VEC, K_VEC, L_CELL, N_atoms, line_strength=line_strength)
    S_dB = observables.intensity_difference_squeezing_dB(G_s, G_c, eta)

    return {
        "D_GHz": D_GHz,
        "probe_axis_GHz": probe_axis_GHz,
        "G_s": G_s,
        "G_c": G_c,
        "S_dB": S_dB,
        "eta": eta,
        "branch": branch,
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
    return {
        "probe_GHz": probe_GHz,
        "G_s": float(np.interp(probe_GHz, x, spectrum["G_s"])),
        "G_c": float(np.interp(probe_GHz, x, spectrum["G_c"])),
        "S_dB": float(np.interp(probe_GHz, x, spectrum["S_dB"])),
    }


# =========================================================
# Scheme wrapper for the generic front-end
# =========================================================
WINDOW_GHZ = 0.55          # half-width of the focused probe window around (−) Raman
TPD_LIMIT_MHZ = 500.0
RESOLUTION = {
    "Fast  (~3 s)":     dict(coarse_points=121, velocity_step=5.0),
    "Balanced  (~6 s)": dict(coarse_points=181, velocity_step=4.0),
    "Fine  (~20 s)":    dict(coarse_points=301, velocity_step=2.0),
}


class FWMScheme(Scheme):
    name = "fwm"
    cluster = "D — Wave mixing"
    title = "Four-wave mixing: gain and biphotons"
    cache_version = "generic-sfwm-biphoton-v1"
    defaults_version = "generic-sfwm-defaults-v1"
    cache_observables = True
    caption = ("Seeded gain/squeezing for the legacy 85Rb double-Lambda system, "
               "or spontaneous SFWM biphoton estimates for cascade and diamond "
               "topologies.")

    def param_schema(self):
        seeded = {"mode": MODE_SEEDED}
        biphoton = {"mode": MODE_BIPHOTON}
        cs_btw = {"mode": MODE_BIPHOTON, "topology": TOPOLOGY_CS_BTW}
        diamond = {"mode": MODE_BIPHOTON, "topology": TOPOLOGY_DIAMOND}
        return [
            ParamSpec("mode", "Mode", "Mode", MODE_SEEDED,
                      choices=(MODE_SEEDED, MODE_BIPHOTON),
                      control="segmented",
                      help="Seeded legacy gain/squeezing or spontaneous SFWM biphoton readout."),
            ParamSpec("topology", "Topology", "Model", TOPOLOGY_RB87_TELECOM,
                      choices=(TOPOLOGY_RB87_TELECOM, TOPOLOGY_CS_BTW, TOPOLOGY_DIAMOND),
                      visible_if=biphoton,
                      help="Level topology for the spontaneous biphoton source model."),
            ParamSpec("cs_channel", "Cs BTW channel", "Model", CS_CHANNEL_917,
                      choices=(CS_CHANNEL_917, CS_CHANNEL_795), visible_if=cs_btw,
                      help="Selects the Cs cascade channel used for BTW comparison."),
            ParamSpec("opd", "OPD — one-photon detuning Δ", "Detunings", 0.9,
                      -1.0, 3.0, 0.1, "GHz",
                      visible_if=seeded,
                      help="ω_pump = ω(F=2→F'=3) + Δ. Sets where the pump sits; recomputes."),
            ParamSpec("tpd", "TPD — two-photon detuning δ", "Detunings", -8.0,
                      -TPD_LIMIT_MHZ, TPD_LIMIT_MHZ, 1.0, "MHz", recompute=False,
                      visible_if=seeded,
                      help="ω_seed = ω_pump − ν_HF + δ. Navigates the curve instantly (no recompute)."),
            ParamSpec("temp_c", "Temperature", "Cell & beams", 121.0,
                      60.0, 150.0, 1.0, "°C", visible_if=seeded),
            ParamSpec("pump_mw", "Pump power", "Cell & beams", 600.0,
                      50.0, 1200.0, 10.0, "mW", visible_if=seeded),
            ParamSpec("probe_uw", "Seed / probe power", "Cell & beams", 8.0,
                      1.0, 200.0, 1.0, "µW", visible_if=seeded),
            ParamSpec("loss_pct", "Loss after cell", "Detection & scaling", 5.5,
                      0.0, 50.0, 0.5, "%", visible_if=seeded,
                      help="Folds into eta = QE x (1 - loss)."),
            ParamSpec("ls", "FWM coupling scale", "Detection & scaling", 0.05,
                      0.01, 1.0, 0.01, "",
                      visible_if=seeded,
                      help="Residual propagation-coupling scale after applying "
                           "Rb85 D1 hyperfine Clebsch-Gordan strengths."),
            ParamSpec("biphoton_temp_c", "Temperature", "Cell & beams", 90.0,
                      30.0, 160.0, 1.0, "°C", visible_if=biphoton),
            ParamSpec("biphoton_cell_mm", "Cell length", "Cell & beams", 12.5,
                      1.0, 100.0, 0.5, "mm", visible_if=biphoton),
            ParamSpec("pump_biphoton_uw", "Pump power", "Fields", 10.0,
                      0.1, 200.0, 0.1, "µW", visible_if=biphoton),
            ParamSpec("coupling_mw", "Coupling power scale", "Fields", 1.0,
                      0.01, 50.0, 0.01, "mW", visible_if=biphoton),
            ParamSpec("pump_detuning_mhz", "Pump detuning", "Detunings", 0.0,
                      -2000.0, 2000.0, 10.0, "MHz", visible_if=biphoton),
            ParamSpec("coupling_detuning_mhz", "Coupling detuning", "Detunings", 0.0,
                      -2000.0, 2000.0, 10.0, "MHz", visible_if=biphoton),
            ParamSpec("signal_angle_deg", "Signal angle", "Phase matching", 1.5,
                      0.0, 10.0, 0.1, "deg", visible_if=biphoton),
            ParamSpec("idler_angle_deg", "Idler angle", "Phase matching", 1.5,
                      0.0, 10.0, 0.1, "deg", visible_if=biphoton),
            ParamSpec("diamond_pump_nm", "Diamond pump wavelength", "Fields", 780.0,
                      300.0, 2000.0, 1.0, "nm", visible_if=diamond),
            ParamSpec("diamond_coupling_nm", "Diamond coupling wavelength", "Fields", 776.0,
                      300.0, 2000.0, 1.0, "nm", visible_if=diamond),
            ParamSpec("diamond_signal_nm", "Diamond signal wavelength", "Fields", 795.0,
                      300.0, 2000.0, 1.0, "nm", visible_if=diamond),
            ParamSpec("diamond_idler_nm", "Diamond idler wavelength", "Fields", 761.702,
                      300.0, 2500.0, 0.001, "nm", visible_if=diamond),
            ParamSpec("signal_eff_pct", "Signal efficiency", "Detection & scaling", 10.0,
                      0.1, 95.0, 0.1, "%", visible_if=biphoton),
            ParamSpec("idler_eff_pct", "Idler efficiency", "Detection & scaling", 10.0,
                      0.1, 95.0, 0.1, "%", visible_if=biphoton),
            ParamSpec("dark_signal_cps", "Signal background", "Detection & scaling", 2000.0,
                      0.0, 100000.0, 100.0, "cps", visible_if=biphoton),
            ParamSpec("dark_idler_cps", "Idler background", "Detection & scaling", 2000.0,
                      0.0, 100000.0, 100.0, "cps", visible_if=biphoton),
            ParamSpec("coincidence_window_ns", "Coincidence window", "Detection & scaling", 1.0,
                      0.01, 100.0, 0.01, "ns", visible_if=biphoton),
            ParamSpec("timing_jitter_ns", "Timing jitter FWHM", "Detection & scaling", 0.55,
                      0.0, 5.0, 0.01, "ns", visible_if=biphoton),
            ParamSpec("filter_bandwidth_mhz", "Filter bandwidth", "Detection & scaling", 300.0,
                      1.0, 5000.0, 1.0, "MHz", visible_if=biphoton),
            ParamSpec("tau_max_ns", "Temporal window", "Numerics", 12.0,
                      1.0, 100.0, 1.0, "ns", visible_if=biphoton, advanced=True),
            ParamSpec("biphoton_velocity_step", "Velocity step", "Numerics", 2.0,
                      0.5, 20.0, 0.5, "m/s", visible_if=biphoton, advanced=True),
            ParamSpec("resolution", "Resolution", "Numerics", "Balanced  (~6 s)",
                      choices=tuple(RESOLUTION.keys()), advanced=True,
                      visible_if=seeded),
        ]

    def presets(self):
        return [
            Preset(
                "Sim et al. 85Rb optimum",
                values=dict(mode=MODE_SEEDED, opd=0.9, tpd=-8.0, temp_c=121.0,
                            pump_mw=600.0, probe_uw=8.0, loss_pct=5.5, ls=0.05),
                icon="FWM",
                help="Seeded double-Lambda gain/squeezing conditions.",
            ),
            Preset(
                "87Rb telecom biphoton",
                values=dict(mode=MODE_BIPHOTON, topology=TOPOLOGY_RB87_TELECOM,
                            biphoton_temp_c=90.0, biphoton_cell_mm=12.5,
                            pump_biphoton_uw=10.0, coupling_mw=1.0,
                            pump_detuning_mhz=0.0, coupling_detuning_mhz=0.0,
                            signal_angle_deg=1.5, idler_angle_deg=1.5,
                            signal_eff_pct=10.0, idler_eff_pct=10.0,
                            dark_signal_cps=2000.0, dark_idler_cps=2000.0,
                            coincidence_window_ns=1.0, filter_bandwidth_mhz=300.0,
                            timing_jitter_ns=0.55, tau_max_ns=12.0,
                            biphoton_velocity_step=2.0),
                icon="Rb",
                help="Cascade 5S-5P-4D telecom-wavelength biphoton reference.",
            ),
            Preset(
                "Cs BTW comparison",
                values=dict(mode=MODE_BIPHOTON, topology=TOPOLOGY_CS_BTW,
                            cs_channel=CS_CHANNEL_917, biphoton_temp_c=75.0,
                            biphoton_cell_mm=12.5, pump_biphoton_uw=20.0,
                            coupling_mw=1.0,
                            pump_detuning_mhz=0.0, coupling_detuning_mhz=0.0,
                            signal_angle_deg=1.5, idler_angle_deg=1.5,
                            signal_eff_pct=10.0, idler_eff_pct=10.0,
                            dark_signal_cps=2000.0, dark_idler_cps=2000.0,
                            coincidence_window_ns=1.0, filter_bandwidth_mhz=300.0,
                            timing_jitter_ns=0.55, tau_max_ns=12.0,
                            biphoton_velocity_step=2.0),
                icon="Cs",
                help="Cascade Cs velocity-class biphoton temporal waveform reference.",
            ),
        ]

    def info(self):
        return (
            "**Seeded reference:** G. Sim, H. Kim, H. S. Moon, *Sci. Rep.* "
            "**15**, 7727 (2025). Legacy 85Rb D1 double-Lambda gain and "
            "intensity-difference squeezing remain regression-anchored.\n\n"
            "**Biphoton references:** Heewoo Kim, Hansol Jeong and Han Seb Moon, "
            "[*Quantum Sci. Technol.* 9, 045006 (2024)]"
            "(https://arxiv.org/abs/2402.06872); Hansol Jeong, Heewoo Kim and "
            "Han Seb Moon, [*Advanced Quantum Technologies* 7, 2300108 (2024)]"
            "(https://www.citedrive.com/en/discovery/highperformance-telecomwavelength-biphoton-source-from-a-hot-atomic-vapor-cell/).\n\n"
            "| Biphoton reference quantity | V1 treatment |\n|---|---|\n"
            "| Velocity-class coherent BTW | coherent sum over Maxwell velocity classes |\n"
            "| Phase matching | generic longitudinal Delta-k with sinc^2(Delta-k L / 2) |\n"
            "| 87Rb telecom source | calibrated to order 38,000 cps/mW and g2 peak ~44 |\n"
            "| Cs BTW channels | separate 852-917 nm and 852-795 nm temporal-width presets |\n\n"
            "The spontaneous mode is a calibrated source estimate, not a full "
            "quantum-Langevin propagation model."
        )

    def compute(self, params):
        if params.get("mode", MODE_SEEDED) == MODE_BIPHOTON:
            topology = topology_from_params(params)
            return GenericFWMSolver(topology).compute_biphoton(params)
        center = branch_center_GHz(params["opd"], -1)
        res = RESOLUTION[params["resolution"]]
        return compute_spectrum(
            params["opd"],
            T=params["temp_c"] + 273.15,
            P_pump=params["pump_mw"] * 1e-3,
            P_probe=params["probe_uw"] * 1e-6,
            line_strength=params["ls"],
            loss_frac=params["loss_pct"] / 100.0,
            coarse_points=res["coarse_points"], fine_points=0,
            scan_min=center - WINDOW_GHZ, scan_max=center + WINDOW_GHZ,
            velocity_step=res["velocity_step"], velocity_cutoff=3.0,
            branch=DEFAULT_BRANCH,
        )

    def observables(self, raw, params):
        if raw.get("kind") == "biphoton":
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
        axG.set_ylabel("Seed / probe gain  $G_s$")
        axG.set_title(f"Δ = {params['opd']:.1f} GHz,  T = {params['temp_c']:.0f} °C,  "
                      f"η = {raw['eta']:.3f}")
        if np.nanmax(raw["G_s"]) > 50:
            axG.set_yscale("log")
        axS.plot(d_axis, raw["S_dB"], color="#2ca02c", lw=1.8)
        axS.axvline(tpd, color="crimson", ls="--", lw=1.2)
        axS.axhline(0.0, color="black", lw=0.6)
        axS.scatter([tpd], [op["S_dB"]], color="crimson", zorder=5)
        axS.set_ylabel("Intensity-difference\nsqueezing  [dB]")
        axS.set_xlabel("Two-photon detuning δ  [MHz]   (probe on the − Raman branch)")
        axS.set_xlim(-TPD_LIMIT_MHZ, TPD_LIMIT_MHZ)
        fig.tight_layout()

        metrics = [
            dict(label="Seed / probe gain  G_s", value=f"{op['G_s']:.2f}",
                 help="Power gain of the seeded probe through the cell."),
            dict(label="Squeezing", value=f"{op['S_dB']:.2f} dB",
                 delta="below shot noise" if op["S_dB"] < 0 else "above shot noise",
                 delta_color="inverse"),
            dict(label="Conjugate gain  G_c", value=f"{op['G_c']:.2f}",
                 help="Generated conjugate power gain (drives the twin-beam squeezing)."),
        ]
        derived = (
            f"| Quantity | Value |\n|---|---|\n"
            f"| N(85Rb) | {raw['N_atoms']:.3e} /m³ |\n"
            f"| σ_v (1-D thermal) | {raw['sigma_v']:.1f} m/s |\n"
            f"| Velocity classes | {raw['n_velocity']} |\n"
            f"| Ω_pump / 2π | {raw['Op_A_2pi_GHz']:.3f} GHz |\n"
            f"| Ω_seed / 2π | {raw['Os_2pi_MHz']:.3f} MHz |\n"
            f"| (−) Raman line (probe axis) | {raw['raman_center_minus_GHz']:.3f} GHz |\n"
            f"| Detection η = QE·(1−loss) | {raw['eta']:.4f} |\n"
            f"| Operating probe detuning | {op['probe_GHz']:.4f} GHz |\n\n"
            f"Fixed: cell L = {L_CELL*1e3:.1f} mm · pump w₀ {W_PUMP*1e6:.0f} µm · "
            f"seed w₀ {W_PROBE*1e6:.0f} µm · QE {QE_DETECTOR*100:.2f}% · "
            f"responsivity {RESPONSIVITY_AW} A/W @ 795 nm · pump⊥probe at PBS."
        )
        # ---- Twin-beam coincidence / correlations (ideal parametric) ----
        coinc = observables.coincidence_stats(raw["G_s"], raw["G_c"])
        ci = int(np.argmin(np.abs(d_axis - tpd)))      # nearest operating-point index
        g2sc_op = coinc["g2_sc"][ci]
        R_op = coinc["cauchy_schwarz"][ci]
        ns_op = coinc["n_s"][ci]
        nc_op = coinc["n_c"][ci]
        has_gain = bool(np.isfinite(g2sc_op))

        figC, (axN, axG2) = plt.subplots(2, 1, figsize=(8.5, 6.4), sharex=True)
        for ax in (axN, axG2):
            ax.grid(alpha=0.3)
        axN.plot(d_axis, coinc["n_s"], color="#1f77b4", lw=1.6, label="signal $n_s = G_s-1$")
        axN.plot(d_axis, coinc["n_c"], color="#ff7f0e", lw=1.2, ls="--", label="conjugate $n_c = G_c$")
        axN.axvline(tpd, color="crimson", ls="--", lw=1.2)
        axN.set_ylabel("Generated photons / mode")
        axN.set_title("Twin-beam photon pairs (ideal parametric, gain region)")
        axN.legend(fontsize=9)
        axG2.plot(d_axis, coinc["g2_sc"], color="#2ca02c", lw=1.6)
        axG2.axhline(2.0, color="black", lw=0.7, ls=":")
        axG2.text(d_axis.min(), 2.05, "classical bound g²=2", fontsize=8, va="bottom")
        axG2.axvline(tpd, color="crimson", ls="--", lw=1.2)
        axG2.set_ylabel("Cross-correlation $g^{(2)}_{sc}(0)$")
        axG2.set_xlabel("Two-photon detuning δ  [MHz]   (probe on the − Raman branch)")
        axG2.set_xlim(-TPD_LIMIT_MHZ, TPD_LIMIT_MHZ)
        axG2.set_ylim(1.8, 12)
        figC.tight_layout()

        if has_gain:
            coinc_rows = (
                f"| Generated photons/mode n_s = G_s−1 | {ns_op:.3f} |\n"
                f"| Conjugate photons n_c = G_c | {nc_op:.3f} |\n"
                f"| Cross-correlation g²_sc(0) = 2 + 1/n_s | {g2sc_op:.3f} |\n"
                f"| Auto-correlation g²_ss = g²_cc | 2.000 (thermal each arm) |\n"
                f"| Cauchy-Schwarz R = g²_sc²/(g²_ss g²_cc) | {R_op:.3f} |\n"
                f"| Nonclassical (R > 1) | {'yes' if R_op > 1 else 'no'} |\n"
            )
        else:
            coinc_rows = "| (no net gain at this δ — no photon pairs) | — |\n"
        coinc_table = (
            f"Ideal (lossless) twin-beam photon-pair statistics at the operating "
            f"point (δ = {tpd:.0f} MHz), from the parametric gain — the spontaneous "
            f"/ coincidence-counting regime:\n\n"
            f"| Quantity | Value |\n|---|---|\n"
            + coinc_rows +
            "\ng²_sc > 2 (and R > 1) violates the classical Cauchy-Schwarz "
            "inequality — the hallmark of FWM photon-pair correlations. g²_sc → 2 "
            "at high gain, → large near threshold. Propagation loss is not modelled "
            "here (as in the squeezing panel)."
        )

        return {
            "metrics": metrics,
            "figure": fig,
            "figures": [("Twin-beam coincidence / correlations", figC)],
            "tables": [
                {"title": "Derived quantities", "markdown": derived},
                {"title": "Twin-beam coincidence (spontaneous)", "markdown": coinc_table},
            ],
        }

    def _biphoton_observables(self, raw, params):
        import matplotlib.pyplot as plt

        topo = raw["topology"]
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
            target_g2_peak=topo.target_g2_peak,
        )

        fig, (axG2, axPM) = plt.subplots(2, 1, figsize=(8.5, 6.4))
        axG2.plot(stats["tau_axis_ns"], stats["g2_SI_tau"], color="#1f77b4", lw=1.8)
        axG2.axhline(2.0, color="black", lw=0.7, ls=":")
        axG2.set_ylabel(r"$g^{(2)}_{SI}(\tau)$")
        axG2.set_title(f"{topo.label}: calibrated spontaneous SFWM estimate")
        axG2.grid(alpha=0.3)

        axPM.plot(raw["angle_axis_deg"], raw["phase_matching"], color="#2ca02c", lw=1.8)
        axPM.axvline(params["idler_angle_deg"], color="crimson", lw=1.1, ls="--")
        axPM.set_xlabel("Idler collection angle [deg]")
        axPM.set_ylabel(r"$\mathrm{sinc}^2(\Delta k L / 2)$")
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

        metrics = [
            dict(label="g2_SI peak", value=f"{stats['g2_peak']:.2f}",
                 help="Peak normalized signal-idler cross-correlation."),
            dict(label="CAR", value=f"{stats['CAR']:.1f}",
                 help="True coincidence divided by accidental coincidence."),
            dict(label="Pair rate", value=f"{stats['pair_rate_cps']:.1f} cps",
                 help="Reference-calibrated generated pair-rate estimate."),
            dict(label="BTW FWHM", value=f"{stats['fwhm_ns']:.2f} ns",
                 help="FWHM of the modeled biphoton temporal waveform."),
            dict(label="Phase match", value=f"{raw['phase_match_weight']:.3f}",
                 help="sinc^2 phase-matching collection weight."),
        ]

        field_rows = "".join(
            f"| {f.role} | {raw['topology'].levels[f.lower].name} -> "
            f"{raw['topology'].levels[f.upper].name} | {f.wavelength_nm:.2f} nm | "
            f"{f.angle_deg:.2f} deg |\n"
            for f in raw["fields"]
        )
        topology_table = (
            f"| Quantity | Value |\n|---|---|\n"
            f"| Topology | {topo.label} |\n"
            f"| Family | {topo.family} |\n"
            f"| Density | {raw['density']:.3e} /m^3 |\n"
            f"| Cell length | {raw['cell_length_m']*1e3:.2f} mm |\n"
            f"| Delta k | {raw['delta_k']:.3e} 1/m |\n"
            f"| Energy mismatch | {raw['energy_mismatch_hz']/1e6:.3f} MHz |\n"
            f"| Residual two-photon Doppler k | {raw['residual_two_photon_k']:.3e} 1/m |\n\n"
            f"| Field | Transition | Wavelength | Angle |\n|---|---|---:|---:|\n"
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
        return {
            "metrics": metrics,
            "figure": fig,
            "figures": [("Velocity-class coherent source", figV)],
            "tables": [
                {"title": "Generic SFWM topology", "markdown": topology_table},
                {"title": "Biphoton detection model", "markdown": detection_table},
                {"title": "Reference validation (medium model)", "markdown": validation_table},
            ],
        }

    def _reference_validation_table(self, raw, params, stats):
        topo = raw["topology"]

        def verdict(ok):
            return "PASS" if ok else "CHECK"

        def row(name, calc, ref, ok, note=""):
            return f"| {name} | {calc} | {ref} | {verdict(ok)} | {note} |\n"

        rows = []
        if topo.name == TOPOLOGY_RB87_TELECOM:
            pump_mw = max(params["pump_biphoton_uw"] * 1e-3, 1e-30)
            rate_per_mw = stats["pair_rate_cps"] / pump_mw
            rows.append(row(
                "Pair rate / pump", f"{rate_per_mw:.0f} cps/mW",
                "38000 cps/mW", abs(rate_per_mw / 38000.0 - 1.0) < 0.15,
                "calibrated source-rate anchor"))
            rows.append(row(
                "g2 peak", f"{stats['g2_peak']:.2f} (raw {stats['raw_g2_peak']:.1f})",
                "44(3)", abs(stats["g2_peak"] - 44.0) <= 3.0,
                "uses explicit added-accidental calibration, not a noise first-principles result"))
            rows.append(row(
                "BTW FWHM", f"{stats['fwhm_ns']:.3f} ns",
                "0.56(4) ns", abs(stats["fwhm_ns"] - 0.56) <= 0.04,
                "detector jitter is part of this medium model"))
            rows.append(row(
                "OD estimate", f"{raw['od_estimate']:.1f}",
                "112(3)", abs(raw["od_estimate"] - 112.0) <= 3.0,
                "density/cell scaling from reference OD"))
            rows.append(row(
                "Bandwidth setting", f"{params['filter_bandwidth_mhz']:.0f} MHz",
                "about 300 MHz", abs(params["filter_bandwidth_mhz"] - 300.0) <= 40.0,
                "filter/source bandwidth check"))
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
                "medium model only; full Cs BTW theory is not yet included"))
            rows.append(row(
                "OD estimate", f"{raw['od_estimate']:.1f}",
                "about 10", abs(raw["od_estimate"] - 10.0) <= 2.0,
                "density/cell scaling from reference note"))
        else:
            rows.append(row(
                "Reference validation", "generic diamond template",
                "no paper anchor", False,
                "configure wavelengths manually; no validated default"))

        rows.append(row(
            "Phase matching", f"{raw['phase_match_weight']:.3f}",
            "> 0.90", raw["phase_match_weight"] > 0.90,
            "sinc^2(Delta k L / 2)"))
        rows.append(row(
            "Energy conservation", f"{raw['energy_mismatch_hz']/1e6:.3f} MHz",
            "near 0 MHz", abs(raw["energy_mismatch_hz"]) < 1e6,
            "wavelength bookkeeping"))

        return (
            "This table is the medium-complexity reference check. It verifies "
            "bookkeeping, calibrated rate/OD anchors, detector-level g2 calibration, "
            "and phase matching. It does not claim full quantum-Langevin noise or "
            "full Cs BTW theory.\n\n"
            "| Check | Calculated | Reference | Verdict | Note |\n|---|---:|---:|---|---|\n"
            + "".join(rows)
        )

    def extra_views(self):
        def _compute_full(params):
            return full_spectrum(
                params["opd"], params["temp_c"] + 273.15,
                params["pump_mw"], params["probe_uw"], params["ls"], params["loss_pct"])

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
            aG.set_ylabel("Seed / probe gain  $G_s$")
            if max(np.nanmax(full[key]["G_s"]) for key in styles) > 50:
                aG.set_yscale("log")
            aS.axhline(0.0, color="black", lw=0.6)
            aS.set_ylabel("Squeezing [dB]")
            aS.set_xlabel(r"Probe detuning from $F=2\to F'=3$  [GHz]")
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


def full_spectrum(D_GHz, T_K, P_pump_mW, P_probe_uW, line_strength, loss_pct):
    """Wide scan with the two Raman channels calculated independently.

    The ∓ branches are independent pure solves, so run them concurrently. With
    BLAS pinned to one thread per branch (see `core.blas_single_thread`) the two
    threads occupy two cores instead of contending — roughly halves this view.
    """
    common = dict(
        T=T_K, P_pump=P_pump_mW * 1e-3, P_probe=P_probe_uW * 1e-6,
        line_strength=line_strength, loss_frac=loss_pct / 100.0,
        coarse_points=301, fine_points=401, velocity_step=2.0)

    def _branch(b):
        with blas_single_thread():
            return compute_spectrum(D_GHz, branch=b, **common)

    with ThreadPoolExecutor(max_workers=2) as ex:
        fut = {b: ex.submit(_branch, b) for b in (-1, +1)}
        minus, plus = fut[-1].result(), fut[+1].result()
    return {"D_GHz": D_GHz, "minus": minus, "plus": plus}
