"""
Lab settings -> GABES parameters -> squeezing, and back out to the detector.

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
    """One complete run: settings in, squeezing and a full audit trail out."""
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
    def squeezing_db(self):
        if self.raw is None:
            return None
        return float(np.asarray(self.raw["S_dB"])[
            _operating_index(self.raw, self.params)])

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
    def squeezing_floor_db(self):
        return 10.0 * np.log10(max(1.0 - self.eta, 1e-300))

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
    """Build the GABES FWM parameter dict from a source chain and geometry."""
    params = dict(base if base is not None else fwm.FWMScheme().defaults())
    params.update({
        "mode": fwm.MODE_SEEDED,
        "opd": settings.opd_ghz,
        "tpd": settings.tpd_mhz,
        "temp_c": settings.cell_temp_c,
        "pump_mw": chain.pump_power_w * 1e3,
        "probe_uw": chain.seed_power_w * 1e6,
        "pump_waist_um": chain.pump_waist_m * 1e6,
        "probe_waist_um": chain.seed_waist_m * 1e6,
        "seeded_angle_deg": settings.crossing_angle_deg,
        "loss_pct": geometry.loss_pct,
        "qe_pct": geometry.qe_pct,
    })
    params.update(overrides)
    return quantize_params(params) if quantize else params


def solve_key(params, scheme=None):
    """Cache key for the heavy solve: the recompute knobs and nothing else.

    Keyed on the *derived* parameters rather than the lab settings, so a knob
    that does not change what reaches the cell -- a detection lens, an iris, a
    trim angle whose effect the grid rounds away -- is a cache hit and the
    readout updates instantly.
    """
    scheme = scheme or fwm.FWMScheme()
    keys = set(scheme.recompute_keys()) | {"mode"}
    return tuple(sorted((k, v) for k, v in params.items() if k in keys))


def run(settings=None, calibration=None, detection=None, *,
        solve=True, scheme=None, **param_overrides):
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
    geom = detection_module.geometry(chain, settings, detection, calibration)
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
        ("Iris", f"{result.geometry.arms[0].iris_transmission * 100:.2f} % of each twin kept",
         "Loss after cell", f"{params['loss_pct']:.2f} %"),
        ("Photodiode", "Hamamatsu S3883",
         "Detector QE", f"{params['qe_pct']:.2f} %"),
        ("Etalon chain",
         f"{sum(1 for s in chain.stages if s.name.startswith('etalon'))} stages",
         "Carrier : seed", f"{chain.carrier_ratio * 100:.4f} %"),
    )


def paper_reference_params(scheme=None):
    """The GABES Squeezing defaults, for checking the bridge plumbing."""
    return dict((scheme or fwm.FWMScheme()).defaults())
