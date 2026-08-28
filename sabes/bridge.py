"""
Lab settings -> GABES parameters -> gain diagnostics, and back out to the detector.

The inversion this whole package exists for happens in `to_gabes_params`. GABES
asks for a one-photon detuning, a two-photon detuning, powers, waists, a crossing
angle and a detection efficiency. None of those is a thing anyone turns. Here each
one is derived:

    pump_mw          <- amplifier current, split half-wave plate, coupling losses
    probe_uw         <- modulator drive, etalon chain, trim waveplates
    pump/probe waist <- collimator diameter x telescope magnification
    seeded_angle_deg <- millimetres of separation at the D-shaped mirrors
    tpd              <- signal-generator frequency (exactly, no coefficient)
    loss_pct         <- iris clipping and post-cell optics
    qe_pct           <- photodiode responsivity

Pipeline order matters and is not arbitrary: detection *geometry* is computed
before the solve because `loss_pct` is a solve input, while detection *readout*
needs the gains the solve produces. `run()` sequences it correctly.

Fidelity note: the pump waist has an interior optimum at Ultra fidelity but looks
monotonic at Fast, because Fast carries only the longitudinal phase mismatch and
not the Gaussian crossing overlap. Any geometry optimisation must run at Ultra.
"""
from dataclasses import dataclass, replace
from typing import Any, Dict, Optional

import numpy as np

from gabes.schemes import fwm

from . import detection as detection_module
from .beamline import SetupSettings, build_source_chain
from .calibration import default_calibration
from .detection import DetectionSettings

#: GABES knobs SABES is responsible for. Everything else keeps its scheme default,
#: so a SABES run stays comparable with a hand-driven GABES run.
DERIVED_KEYS = (
    "opd", "tpd", "temp_c", "cell_mm",
    "pump_mw", "probe_uw", "pump_waist_um", "probe_waist_um",
    "seeded_angle_deg", "loss_pct", "qe_pct",
    "seed_wanted_sideband_uw", "eom_residual_carrier_uw",
    "eom_other_sidebands_uw", "eom_residual_carrier_to_wanted_ratio",
    "eom_other_sidebands_to_wanted_ratio",
)

# Non-numeric audit metadata changes the returned trust ledger even though it
# does not change the Bloch equations.  It must therefore survive the SABES
# cached-solve reconstruction and participate in that result's cache key.
SOLVE_AUDIT_KEYS = (
    "eom_seed_spectrum_provenance",
    "eom_seed_spectrum_status",
    "eom_seed_spectrum_application",
)

#: Rounding grid for derived parameters, in each parameter's own unit.
#:
#: Derived values carry the full float noise of the chain -- nudging a waveplate
#: by a thousandth of a degree moves `pump_mw` in its seventh digit. Left alone
#: that busts the solve cache on knob movements far below any physical
#: significance. Quantising the parameters themselves (rather than only the cache
#: key) keeps the displayed number identical to the solved one; a separate key
#: rounding would let the two drift apart, which is worse than a slow UI.
PARAM_GRID = {
    "pump_mw": 0.01,
    "probe_uw": 0.001,
    "seed_wanted_sideband_uw": 0.001,
    "eom_residual_carrier_uw": 0.001,
    "eom_other_sidebands_uw": 0.001,
    "pump_waist_um": 0.01,
    "probe_waist_um": 0.01,
    "seeded_angle_deg": 1e-4,
    "loss_pct": 1e-3,
    "qe_pct": 1e-3,
    "tpd": 1e-3,
}


def quantize_params(params, grid=None):
    """Snap derived parameters onto `PARAM_GRID`. Settings pass through."""
    grid = PARAM_GRID if grid is None else grid
    out = dict(params)
    for key, step in grid.items():
        if key in out and isinstance(out[key], (int, float)):
            out[key] = round(float(out[key]) / step) * step
    return out


@dataclass(frozen=True)
class SabesResult:
    """One complete run: settings in, diagnostics and a full audit trail out."""
    settings: SetupSettings
    chain: Any
    geometry: Any
    params: Dict[str, Any]
    raw: Optional[Dict[str, Any]] = None
    observables: Optional[Dict[str, Any]] = None
    readout: Optional[Any] = None

    @property
    def gains(self):
        if self.raw is None:
            return None
        index = _operating_index(self.raw, self.params)
        return (float(np.asarray(self.raw["G_s"])[index]),
                float(np.asarray(self.raw["G_c"])[index]))

    @property
    def gain_referred_noise_db(self):
        """Gain-referred algebraic diagnostic; not a squeezing prediction."""
        if self.raw is None:
            return None
        values = self.raw.get("gain_referred_noise_dB")
        if values is None:
            values = self.raw["S_dB"]
        return float(np.asarray(values)[
            _operating_index(self.raw, self.params)])

    @property
    def squeezing_db(self):
        """Backward-compatible alias for :attr:`gain_referred_noise_db`."""
        return self.gain_referred_noise_db

    @property
    def physical_squeezing_db(self):
        """Physical squeezing prediction, intentionally unavailable for now."""
        if self.raw is None:
            return None
        return self.raw.get("physical_squeezing_dB")

    @property
    def validation_level(self):
        """Machine-readable claim level emitted by the seeded-FWM solver."""
        if self.raw is None:
            return None
        return self.raw.get("validation_level")

    @property
    def claim_gate(self):
        """Machine-readable permissions and reasons for physical claims."""
        if self.raw is None:
            return None
        return self.raw.get("claim_gate")

    @property
    def eta(self):
        """Detection efficiency the solve actually used.

        Read from `params`, not from `geometry`: an explicit override (a test
        pinning the paper's numbers, say) changes what was solved without
        changing what the optics predict, and the reported floor must follow the
        former.
        """
        return self.params["qe_pct"] / 100.0 * (1.0 - self.params["loss_pct"] / 100.0)

    @property
    def gain_referred_detection_floor_db(self):
        """Symmetric-loss floor of the algebraic diagnostic, not a noise model."""
        return 10.0 * np.log10(max(1.0 - self.eta, 1e-300))

    @property
    def squeezing_floor_db(self):
        """Backward-compatible alias for the gain-referred detection floor."""
        return self.gain_referred_detection_floor_db

    @property
    def warnings(self):
        parts = list(self.chain.warnings) + list(self.geometry.warnings)
        if self.readout is not None:
            parts += list(self.readout.warnings)
        return tuple(parts)


def _operating_index(raw, params):
    """Index of the two-photon detuning the modulator is actually set to."""
    axis = np.asarray(raw["probe_axis_GHz"])
    center = fwm.branch_center_GHz(params["opd"], -1)
    return int(np.argmin(np.abs((axis - center) * 1e3 - params["tpd"])))


def to_gabes_params(chain, settings, geometry, *, base=None, quantize=True,
                    **overrides):
    """Build the GABES FWM parameter dict from a source chain and geometry.

    ``probe_uw`` remains the wanted-sideband power for backward compatibility.
    The additional EOM fields preserve the rest of the cell-plane seed spectrum
    as model inputs only. No carrier or sideband power is converted into atomic
    or detector noise without an independently calibrated transfer model.
    """
    params = dict(base if base is not None else fwm.FWMScheme().defaults())
    wanted_seed_uw = chain.wanted_seed_sideband_power_w * 1e6
    params.update({
        "mode": fwm.MODE_SEEDED,
        "opd": settings.opd_ghz,
        "tpd": settings.tpd_mhz,
        "temp_c": settings.cell_temp_c,
        "pump_mw": chain.pump_power_w * 1e3,
        "probe_uw": wanted_seed_uw,
        "seed_wanted_sideband_uw": wanted_seed_uw,
        "eom_residual_carrier_uw": chain.eom_residual_carrier_power_w * 1e6,
        "eom_other_sidebands_uw": chain.eom_other_sidebands_power_w * 1e6,
        "eom_residual_carrier_to_wanted_ratio": (
            chain.eom_residual_carrier_to_wanted_ratio),
        "eom_other_sidebands_to_wanted_ratio": (
            chain.eom_other_sidebands_to_wanted_ratio),
        "eom_seed_spectrum_provenance": chain.seed_spectrum_provenance,
        "eom_seed_spectrum_status": "unsupported",
        "eom_seed_spectrum_application": "unapplied",
        "pump_waist_um": chain.pump_waist_m * 1e6,
        "probe_waist_um": chain.seed_waist_m * 1e6,
        "seeded_angle_deg": settings.crossing_angle_deg,
        "loss_pct": geometry.loss_pct,
        "qe_pct": geometry.qe_pct,
    })
    if ("probe_uw" in overrides and "seed_wanted_sideband_uw" in overrides
            and not np.isclose(
                float(overrides["probe_uw"]),
                float(overrides["seed_wanted_sideband_uw"]),
                rtol=0.0, atol=1e-12)):
        raise ValueError(
            "probe_uw and seed_wanted_sideband_uw must describe the same "
            "wanted EOM sideband power")
    if str(overrides.get(
            "eom_seed_spectrum_status", "unsupported")).strip().lower() \
            != "unsupported":
        raise ValueError(
            "SABES EOM spectrum remains unsupported until a calibrated "
            "transfer model is implemented")
    if str(overrides.get(
            "eom_seed_spectrum_application", "unapplied")).strip().lower() \
            != "unapplied":
        raise ValueError(
            "SABES EOM residual components must remain unapplied")
    if "probe_uw" in overrides:
        overrides["seed_wanted_sideband_uw"] = overrides["probe_uw"]
    elif "seed_wanted_sideband_uw" in overrides:
        overrides["probe_uw"] = overrides["seed_wanted_sideband_uw"]
    params.update(overrides)
    if quantize:
        params = quantize_params(params)

    # Ratios are derived from the exact powers handed to the solver, after any
    # cache-grid quantization, so the bridge ledger and FWM raw audit cannot
    # silently disagree.
    wanted = float(params["probe_uw"])
    params["seed_wanted_sideband_uw"] = wanted
    denominator = wanted if wanted > 0.0 else float("nan")
    params["eom_residual_carrier_to_wanted_ratio"] = (
        float(params["eom_residual_carrier_uw"]) / denominator
        if wanted > 0.0 else float("inf"))
    params["eom_other_sidebands_to_wanted_ratio"] = (
        float(params["eom_other_sidebands_uw"]) / denominator
        if wanted > 0.0 else float("inf"))
    return params


def solve_key(params, scheme=None):
    """Cache key for the heavy solve: the recompute knobs and nothing else.

    Keyed on the *derived* parameters rather than the lab settings, so a knob
    that does not change what reaches the cell -- a detection lens, an iris, a
    trim angle whose effect the grid rounds away -- is a cache hit and the
    readout updates instantly.
    """
    scheme = scheme or fwm.FWMScheme()
    keys = set(scheme.recompute_keys()) | {"mode"} | set(SOLVE_AUDIT_KEYS)
    return tuple(sorted((k, v) for k, v in params.items() if k in keys))


def run(settings=None, calibration=None, detection=None, *,
        layout=None, solve=True, scheme=None, **param_overrides):
    """Full pipeline: source chain -> geometry -> GABES solve -> detector readout.

    `solve=False` stops after the parameter dict, which is the cheap path the UI
    uses when only optics knobs moved: beam sizes, purity and the power budget
    all update without touching the Bloch solve.
    """
    settings = settings or SetupSettings()
    calibration = calibration or default_calibration()
    detection = detection or DetectionSettings()
    scheme = scheme or fwm.FWMScheme()

    chain = build_source_chain(settings, calibration)
    geom = detection_module.geometry(chain, settings, detection, calibration,
                                     layout=layout)
    params = to_gabes_params(chain, settings, geom, **param_overrides)

    if not solve:
        return SabesResult(settings=settings, chain=chain, geometry=geom,
                           params=params)

    raw = scheme.compute(params)
    observables = scheme.headless_observables(raw, params)
    result = SabesResult(settings=settings, chain=chain, geometry=geom,
                         params=params, raw=raw, observables=observables)
    read = detection_module.readout(geom, chain, result.gains, settings,
                                    detection, calibration)
    return replace(result, readout=read)


def derived_table(result):
    """Primary setting -> derived GABES quantity, for the UI's trust panel.

    Showing both columns side by side is what lets a SABES number be checked
    against a hand-driven GABES run instead of being taken on faith.
    """
    settings, params, chain = result.settings, result.params, result.chain
    return (
        ("Signal generator", f"{settings.eom_frequency_hz / 1e9:.6f} GHz",
         "TPD δ", f"{params['tpd']:+.3f} MHz"),
        ("Split HWP", f"{settings.hwp_split_deg:.3f}°",
         "Pump power", f"{params['pump_mw']:.1f} mW"),
        ("Seed trim HWP", f"{settings.seed_trim_hwp_deg:.3f}°",
         "Seed power", f"{params['probe_uw']:.2f} µW"),
        ("Pump telescope",
         f"f{settings.pump_telescope_f1_mm:.0f}→f{settings.pump_telescope_f2_mm:.0f}",
         "Pump waist", f"{params['pump_waist_um']:.0f} µm"),
        ("Seed telescope",
         f"f{settings.seed_telescope_f1_mm:.0f}→f{settings.seed_telescope_f2_mm:.0f}",
         "Seed waist", f"{params['probe_waist_um']:.0f} µm"),
        ("D-mirror separation", f"{settings.dmirror_separation_mm:.2f} mm",
         "Crossing angle", f"{params['seeded_angle_deg']:.3f}°"),
        ("Post-cell optics",
         f"{result.geometry.optics_transmission * 100:.2f} % transmitted",
         "Loss after cell", f"{params['loss_pct']:.2f} %"),
        ("Photodiode", "Hamamatsu S3883",
         "Detector QE", f"{params['qe_pct']:.2f} %"),
        ("Etalon chain",
         f"{sum(1 for s in chain.stages if s.name.startswith('etalon'))} stages",
         "Carrier : seed (solver grid)",
         f"{params['eom_residual_carrier_to_wanted_ratio'] * 100:.4f} %"),
        ("EOM spectrum provenance", params["eom_seed_spectrum_provenance"],
         "Carrier/noise application",
         f"{params['eom_seed_spectrum_application']} "
         f"({params['eom_seed_spectrum_status']})"),
    )


def paper_reference_params(scheme=None):
    """The seeded-FWM defaults, for checking the bridge plumbing."""
    return dict((scheme or fwm.FWMScheme()).defaults())
