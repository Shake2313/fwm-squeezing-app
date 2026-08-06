"""
The Sim et al. source chain, assembled.

    ECDL -> TA -> PBS ---- pump ---> fiber -> collimator -> telescope -> cell
                       \
                        -- seed ---> fiber EOM -> etalon x3 -> fiber
                                     -> collimator -> telescope -> trim -> cell

Two paths with different jobs, which is why they are modelled at different
depths. The pump path is a *power* budget: nothing about its spectrum is in
question, and the pump/seed PBS half-wave plate is the knob that sets it. The
seed path is a *purity* budget: it has 30x more power than it needs and is
deliberately attenuated, so what matters is the residual carrier the etalons let
through, which is the quantity Sim et al. show degrading the squeezing.

What is deliberately NOT modelled
  - ECDL frequency from temperature / current / PZT. Mode hops make it
    non-monotonic and hysteretic, so a first-principles model would be fiction.
    `SetupSettings.opd_ghz` is the detuning; the three ECDL knobs are recorded as
    annotations so a working point can be written down and returned to.
  - Amplifier gain competition between spectral components. The modulator sits
    *after* the split, so the amplifier only ever sees one carrier.
  - Anything time-dependent: locks, drift, phase noise.
"""
from dataclasses import dataclass, replace
from typing import Tuple

import math

from gabes import constants

from .beamstate import BeamState, monochromatic
from .calibration import default_calibration
from .layout.parts import CELL_TO_DMIRROR_M
from .devices import (EtalonFilter, FiberCollimator, FiberCoupling, FreeSpace,
                      HalfWavePlate, LossElement, PhaseModulator, Polarizer,
                      PolarizingBeamSplitter, QuarterWavePlate,
                      TaperedAmplifier, ThinLens, WaistRescale)

WAVELENGTH_M = constants.WAVELENGTH_D1_85RB
NU_HF_HZ = constants.NU_GROUND_HF


@dataclass(frozen=True)
class SetupSettings:
    """Primary lab variables -- what someone actually turns.

    Everything here is a dial, a current, a temperature or a lens choice. Derived
    physics quantities (Rabi frequencies, detunings in linewidths, waists in the
    cell) are *outputs*, computed from these. That inversion is the whole
    difference between SABES and driving the GABES scheme directly.
    """
    # --- source ---
    ecdl_power_mw: float = 40.0
    ta_current_a: float = 2.5
    hwp_ta_deg: float = 0.0

    # --- pump / seed split ---
    # 12.3106 deg puts 600 mW on the pump with the shipped calibration. If a
    # coefficient changes, this angle no longer hits 600 mW -- that is the point.
    # `nominal_settings()` re-solves all three power knobs against whatever
    # calibration is in force.
    hwp_split_deg: float = 12.3106

    # --- seed modulation ---
    # The fiber EOM stops modulating cleanly above ~10 mW, well below what the
    # split alone leaves in the seed arm once the pump is set. A linear polarizer
    # in front of the modulator is what actually holds the input inside its
    # rating, so its angle is a primary variable, not a fixed loss.
    seed_polarizer_deg: float = 12.1500
    hwp_eom_deg: float = 0.0
    eom_frequency_hz: float = NU_HF_HZ + 8.0e6
    eom_rf_dbm: float = 18.0
    etalon_detune_ghz: Tuple[float, float, float] = (0.0, 0.0, 0.0)

    # --- seed power trim (QWP + HWP + GTP) ---
    seed_trim_qwp_deg: float = 0.0
    seed_trim_hwp_deg: float = 38.6033
    seed_gtp_deg: float = 0.0

    # --- beam shaping (lenses are swappable, so they are settings) ---
    pump_telescope_f1_mm: float = 300.0
    pump_telescope_f2_mm: float = 150.0
    seed_telescope_f1_mm: float = 400.0
    seed_telescope_f2_mm: float = 150.0

    # --- geometry / cell ---
    opd_ghz: float = 0.9
    cell_temp_c: float = 121.0
    dmirror_separation_mm: float = 3.2628

    # --- recorded but not modelled (see module docstring) ---
    ecdl_temp_c: float = float("nan")
    ecdl_current_ma: float = float("nan")
    ecdl_pzt_v: float = float("nan")
    ta_temp_c: float = float("nan")

    @property
    def tpd_mhz(self):
        """Two-photon detuning implied by the modulator setting.

        The -1 sideband sits at -f_EOM, and the GABES convention is
        omega_seed = omega_pump - nu_HF + delta, so delta = nu_HF - f_EOM.
        Exact, with no calibration coefficient: the signal generator sets this
        to sub-Hz precision, which makes it the one primary variable that needs
        no tuning factor at all -- and the most sensitive one.
        """
        return (NU_HF_HZ - self.eom_frequency_hz) / 1e6

    @property
    def crossing_angle_deg(self):
        """Angle from the separation the beams are actually set to.

        Alignment is done by overlapping at the cell and separating by a stated
        number of mm at the D-shaped mirrors, so the millimetre reading is the
        primary variable and the angle is derived.
        """
        return math.degrees(math.atan(
            self.dmirror_separation_mm * 1e-3 / CELL_TO_DMIRROR_M))


@dataclass(frozen=True)
class Stage:
    """One row of the power budget."""
    path: str
    name: str
    beam: BeamState

    @property
    def power_w(self):
        return self.beam.total_power_w

    @property
    def radius_m(self):
        return self.beam.mode.radius_m

    @property
    def waist_m(self):
        """The quantity GABES wants: over a 12.5 mm cell the beam is collimated
        to within 0.01 %, so the waist is the effective radius in the vapour."""
        return self.beam.mode.waist_m


@dataclass(frozen=True)
class SourceChain:
    """Result of running the chain: the two beams at the cell, plus the audit."""
    pump: BeamState
    seed: BeamState
    seed_offset_hz: float
    stages: Tuple[Stage, ...] = ()
    warnings: Tuple[str, ...] = ()

    @property
    def pump_power_w(self):
        return self.pump.total_power_w

    @property
    def seed_power_w(self):
        """Power in the wanted sideband only -- the residual carrier is not seed."""
        return self.seed.power_at(self.seed_offset_hz)

    @property
    def carrier_ratio(self):
        return self.seed.carrier_ratio(self.seed_offset_hz)

    @property
    def seed_purity(self):
        return self.seed.purity(self.seed_offset_hz)

    @property
    def pump_waist_m(self):
        return self.pump.mode.waist_m

    @property
    def seed_waist_m(self):
        return self.seed.mode.waist_m

    def stage_power(self, path, name):
        for stage in self.stages:
            if stage.path == path and stage.name == name:
                return stage.power_w
        raise KeyError(f"no stage {path}/{name!r} in this chain")

    def budget_rows(self, path=None):
        """(path, stage, power_W, waist_m) for the budget table."""
        return tuple(
            (s.path, s.name, s.power_w, s.waist_m)
            for s in self.stages if path is None or s.path in ("common", path))


def _telescope_recorded(beam, f1_m, f2_m, record, path):
    """Afocal pair, but recording the state between the lenses.

    `Telescope.apply` collapses the whole thing into one step, which is fine for
    a power budget and wrong for a table view: a probe dropped between L1 and L2
    would otherwise be handed the state from *after* the telescope. Recording the
    intermediate makes the layout's labels true. The optics are identical, so the
    delivered waist is unchanged.
    """
    stage = FreeSpace(f1_m).apply(beam)
    stage = ThinLens(f1_m).apply(stage)
    record(path, "telescope L1", stage)
    stage = FreeSpace(f1_m + f2_m).apply(stage)
    stage = ThinLens(f2_m).apply(stage)
    return FreeSpace(f2_m).apply(stage)


def _etalon(calibration, seed_offset_hz, detune_ghz):
    """One filter stage, tuned so the wanted sideband sits on the passband."""
    return EtalonFilter(
        fsr_hz=calibration.value("etalon_fsr_hz"),
        fwhm_hz=calibration.value("etalon_fwhm_hz"),
        peak_transmission=calibration.value("etalon_peak_transmission"),
        center_offset_hz=seed_offset_hz + float(detune_ghz) * 1e9,
        leak=calibration.value("etalon_leak"),
    )


def build_source_chain(settings=None, calibration=None):
    """Run the full source chain and return both beams at the cell face."""
    settings = settings or SetupSettings()
    calibration = calibration or default_calibration()
    c = calibration.value
    stages = []
    warnings = []

    def record(path, name, beam):
        stages.append(Stage(path, name, beam))
        return beam

    # ---------------- common: ECDL -> TA -> split ----------------
    beam = monochromatic(
        settings.ecdl_power_mw * 1e-3, WAVELENGTH_M,
        waist_m=0.5e-3, m2=1.0)
    record("common", "ECDL output", beam)

    beam = LossElement(c("loss_ecdl_to_ta")).apply(beam)
    record("common", "to amplifier", beam)

    beam = HalfWavePlate(settings.hwp_ta_deg).apply(beam)
    amplifier = TaperedAmplifier(
        current_a=settings.ta_current_a,
        threshold_a=c("ta_threshold_a"),
        slope_w_per_a=c("ta_slope_w_per_a"),
        saturation_w=c("ta_saturation_w"),
        m2_out=c("ta_m2_out"),
        ase_fraction=c("ta_ase_fraction"),
    )
    beam = amplifier.apply(beam)
    record("common", "tapered amplifier", beam)

    beam = LossElement(c("loss_ta_to_split")).apply(beam)
    beam = HalfWavePlate(settings.hwp_split_deg).apply(beam)
    pump, seed = PolarizingBeamSplitter().split(beam)
    record("pump", "PBS transmitted", pump)
    record("seed", "PBS reflected", seed)

    # ---------------- pump path ----------------
    pump = LossElement(c("loss_pump_path")).apply(pump)
    pump = FiberCoupling(c("pump_fiber_efficiency")).apply(pump)
    record("pump", "delivery fiber", pump)

    pump = FiberCollimator(c("pump_collimator_diameter_mm") * 1e-3).apply(pump)
    record("pump", "collimator", pump)
    pump = _telescope_recorded(pump, settings.pump_telescope_f1_mm * 1e-3,
                               settings.pump_telescope_f2_mm * 1e-3,
                               record, "pump")
    pump = WaistRescale(c("waist_scale_factor")).apply(pump)
    record("pump", "collimator + telescope", pump)

    pump = LossElement(c("loss_cell_windows")).apply(pump)
    record("pump", "at cell", pump)

    # ---------------- seed path ----------------
    seed_offset_hz = -settings.eom_frequency_hz

    seed = LossElement(c("loss_seed_path")).apply(seed)
    seed = Polarizer(settings.seed_polarizer_deg).apply(seed)
    record("seed", "input polarizer", seed)

    seed = HalfWavePlate(settings.hwp_eom_deg).apply(seed)
    seed = FiberCoupling(c("seed_fiber_in_efficiency")).apply(seed)
    record("seed", "into fiber EOM", seed)

    modulator = PhaseModulator(
        frequency_hz=settings.eom_frequency_hz,
        v_pi_v=c("eom_v_pi_v"),
        rf_power_dbm=settings.eom_rf_dbm,
        insertion_loss_db=c("eom_insertion_loss_db"),
        max_input_w=c("eom_max_input_w"),
        orthogonal_leakage=c("eom_orthogonal_leakage"),
    )
    if modulator.overdriven(seed):
        warnings.append(
            f"fiber EOM input {seed.total_power_w * 1e3:.1f} mW exceeds its "
            f"{c('eom_max_input_w') * 1e3:.0f} mW rating; modulation depth is "
            f"not trustworthy above it")
    seed = modulator.apply(seed)
    record("seed", "EOM sidebands", seed)

    stage_count = int(round(c("etalon_stages")))
    detunes = list(settings.etalon_detune_ghz)
    detunes += [0.0] * max(stage_count - len(detunes), 0)
    for index in range(stage_count):
        seed = _etalon(calibration, seed_offset_hz, detunes[index]).apply(seed)
        record("seed", f"etalon {index + 1}", seed)

    seed = FiberCoupling(c("seed_fiber_out_efficiency")).apply(seed)
    seed = FiberCollimator(c("seed_collimator_diameter_mm") * 1e-3).apply(seed)
    record("seed", "delivery fiber + collimator", seed)

    seed = QuarterWavePlate(settings.seed_trim_qwp_deg).apply(seed)
    seed = HalfWavePlate(settings.seed_trim_hwp_deg).apply(seed)
    record("seed", "collimator", seed)
    seed = _telescope_recorded(seed, settings.seed_telescope_f1_mm * 1e-3,
                               settings.seed_telescope_f2_mm * 1e-3,
                               record, "seed")
    seed = WaistRescale(c("waist_scale_factor")).apply(seed)
    seed = Polarizer(settings.seed_gtp_deg).apply(seed)
    record("seed", "trim + telescope", seed)

    seed = LossElement(c("loss_cell_windows")).apply(seed)
    record("seed", "at cell", seed)

    if seed.power_at(seed_offset_hz) <= 0.0:
        warnings.append("no power left in the wanted sideband at the cell")

    return SourceChain(pump=pump, seed=seed, seed_offset_hz=seed_offset_hz,
                       stages=tuple(stages), warnings=tuple(warnings))


# ---------------------------------------------------------------------------
# Knob inversion: "what angle gives me X?"
#
# Solved numerically rather than from the cos^2 law, because the beam reaching a
# waveplate is never exactly the polarization the ideal law assumes -- PBS
# extinction leaves a couple of degrees of the wrong component, which is small
# for power but large enough to miss a target by 30 %. A coarse scan brackets the
# root on a monotone branch and bisection finishes it; the chain is closed-form,
# so this costs microseconds.
# ---------------------------------------------------------------------------


def _solve_angle(field, target, readout, settings, calibration,
                 lo=0.0, hi=90.0, samples=181, iterations=60):
    def value(angle):
        chain = build_source_chain(replace(settings, **{field: angle}),
                                   calibration)
        return readout(chain)

    step = (hi - lo) / (samples - 1)
    grid = [lo + step * i for i in range(samples)]
    values = [value(a) for a in grid]
    bracket = None
    for i in range(samples - 1):
        if (values[i] - target) * (values[i + 1] - target) <= 0.0:
            bracket = (grid[i], grid[i + 1])
            break
    if bracket is None:
        best = min(range(samples), key=lambda i: abs(values[i] - target))
        raise ValueError(
            f"{field}: target {target:.6g} is out of reach; the closest "
            f"achievable value is {values[best]:.6g} at {grid[best]:.3f} deg")

    a, b = bracket
    fa = value(a) - target
    for _ in range(iterations):
        mid = 0.5 * (a + b)
        fm = value(mid) - target
        if fa * fm <= 0.0:
            b = mid
        else:
            a, fa = mid, fm
    return 0.5 * (a + b)


def solve_split_angle_deg(target_pump_w, settings=None, calibration=None):
    """Half-wave-plate angle at the pump/seed PBS giving a target pump power."""
    return _solve_angle(
        "hwp_split_deg", float(target_pump_w), lambda ch: ch.pump_power_w,
        settings or SetupSettings(), calibration or default_calibration(),
        lo=0.0, hi=45.0)


def solve_seed_polarizer_deg(target_eom_input_w, settings=None, calibration=None):
    """Seed polarizer angle giving a target optical power at the fiber EOM.

    Separate from the split angle on purpose. Once the pump is set, the seed arm
    carries far more than the modulator's ~10 mW rating, so these are two
    independent knobs solving two independent constraints -- which is exactly why
    the setup needs a polarizer there at all.
    """
    return _solve_angle(
        "seed_polarizer_deg", float(target_eom_input_w),
        lambda ch: ch.stage_power("seed", "into fiber EOM"),
        settings or SetupSettings(), calibration or default_calibration(),
        lo=0.0, hi=88.0)


def nominal_settings(calibration=None, *, pump_w=0.600, eom_input_w=7.0e-3,
                     seed_w=8.0e-6, base=None):
    """Settings whose three power knobs reproduce the paper operating point.

    Re-solved against whatever calibration is in force, so it stays honest when a
    coefficient is replaced by a measured one. The dataclass defaults are the
    answer for the shipped calibration; this is how to get it for any other.
    """
    calibration = calibration or default_calibration()
    settings = base or SetupSettings()
    settings = replace(settings, hwp_split_deg=solve_split_angle_deg(
        pump_w, settings, calibration))
    settings = replace(settings, seed_polarizer_deg=solve_seed_polarizer_deg(
        eom_input_w, settings, calibration))
    return replace(settings, seed_trim_hwp_deg=solve_seed_trim_deg(
        seed_w, settings, calibration))


def solve_seed_trim_deg(target_seed_w, settings=None, calibration=None):
    """Trim half-wave-plate angle giving a target seed power at the cell.

    The seed arm arrives with tens of times more power than the FWM wants, so
    this knob exists purely to throw the excess away -- which is also why the
    seed budget is not a constraint on the design.
    """
    return _solve_angle(
        "seed_trim_hwp_deg", float(target_seed_w), lambda ch: ch.seed_power_w,
        settings or SetupSettings(), calibration or default_calibration(),
        lo=0.0, hi=45.0)
