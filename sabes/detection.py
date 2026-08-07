"""
Everything after the cell: separate the twins from the pump, get them onto the
photodiodes, and decide whether the measurement can see squeezing at all.

Split in two on purpose.

  `geometry()`  runs BEFORE the atomic solve. Beam sizes, crossing separation
                and the post-cell transmission depend only on where the optics
                are, so they can produce the `loss_pct` GABES needs as a solve
                input.
  `readout()`   runs AFTER, because photodiode power density, amplifier
                saturation and shot-noise clearance all need the parametric
                gains that only the solve knows.

That ordering is what keeps the pipeline acyclic.

What this module adds over GABES: GABES takes a detection efficiency as a number.
Here it is a consequence -- of how far the D-shaped mirrors are, the per-optic
transmissions along the way, which lens is in the mount, and how far the diode
sits from focus. Every one of those can be wrong in a way GABES cannot express.

Pump rejection is NOT derived. A Gaussian tail at five beam radii says 1e-9,
which is an answer so good it is useless and one the paper contradicts by
naming pump scatter as a real limit. `pump_leakage_dbm` is an observed number
instead: measure it at the detector, and the spectrum analyser uses it.
"""
from dataclasses import dataclass
from typing import Tuple

import math

from .calibration import default_calibration
from .layout import parts as layout_parts

ELEMENTARY_CHARGE = 1.602176634e-19

# Separation over beam radius below which a Gaussian tail is no longer a fair
# description of "these beams are resolved": at 3 the overlap integral is already
# ~1e-4 and real scatter dominates long before the Gaussian does.
SEPARATION_MARGIN_WARN = 3.0

# Measuring S dB of squeezing needs the squeezed trace to stay clear of the
# amplifier floor, so the shot-noise clearance has to exceed |S| by a usable
# margin -- roughly |S| + 6 dB. At the -8 dB operating point that is ~14 dB, so
# 15 dB is the threshold below which the readout is called into question.
CLEARANCE_WARN_DB = 15.0


@dataclass(frozen=True)
class DetectionSettings:
    """Post-cell knobs. Every one of these is something adjustable on the table.

    The iris is deliberately absent. It is opened past the twin beams in
    practice, and modelling pump rejection as a Gaussian tail claimed 1e-9
    suppression -- an answer so good it was useless, and one the paper
    contradicts. Pump leakage is now an OBSERVED quantity, `pump_leakage_dbm`,
    measured at the detector and fed to the spectrum analyser.
    """
    pump_leakage_dbm: float = -60.0
    probe_lens_focal_mm: float = 75.0
    conjugate_lens_focal_mm: float = 100.0
    pd_defocus_mm: float = 0.0
    bpd_gain_v_per_a: float = 1e5
    analysis_bandwidth_hz: float = 1.0e5


@dataclass(frozen=True)
class ArmGeometry:
    """One twin-beam arm from the cell to its photodiode."""
    name: str
    radius_at_dmirror_m: float
    spot_radius_m: float
    lens_focal_m: float


@dataclass(frozen=True)
class DetectionGeometry:
    """Pre-solve result: what the optics do to power, independent of gain."""
    arms: Tuple[ArmGeometry, ...]
    twin_separation_m: float
    pump_separation_m: float
    pump_radius_at_dmirror_m: float
    twin_radius_at_dmirror_m: float
    optics_transmission: float
    loss_pct: float
    qe_pct: float
    warnings: Tuple[str, ...] = ()

    @property
    def separation_margin(self):
        """Pump-to-twin centre separation in twin beam radii.

        It rises with distance and saturates at `theta * pi * w0 / lambda`, the
        crossing angle divided by the divergence angle. Near the cell separation
        outruns the beam's growth; far away the two grow together and the ratio
        locks onto that ceiling. So moving the D-shaped mirrors back buys a
        bounded amount -- at the actual 23 inch the margin is already ~82 % of
        the ceiling -- and raising the ceiling itself needs a larger waist or a
        larger angle.
        """
        return self.pump_separation_m / self.twin_radius_at_dmirror_m

    @property
    def eta(self):
        return self.qe_pct / 100.0 * (1.0 - self.loss_pct / 100.0)

    @property
    def squeezing_floor_db(self):
        """10 log10(1 - eta): no operating point can beat this."""
        return 10.0 * math.log10(max(1.0 - self.eta, 1e-300))


@dataclass(frozen=True)
class ArmReadout:
    name: str
    power_w: float
    residual_pump_w: float
    power_density_w_per_cm2: float
    spot_radius_m: float
    fills_active_area: bool


@dataclass(frozen=True)
class DetectionReadout:
    """Post-solve result: can this measurement actually see the squeezing?"""
    arms: Tuple[ArmReadout, ...]
    total_power_w: float
    shot_noise_a_per_rthz: float
    electronic_noise_a_per_rthz: float
    clearance_db: float
    warnings: Tuple[str, ...] = ()

    @property
    def measurable(self):
        return self.clearance_db >= CLEARANCE_WARN_DB

    def margin_above_electronic_db(self, squeezing_db):
        """How far a squeezed trace sits above the amplifier noise floor.

        The squeezed level is |S| below shot noise, and shot noise is
        `clearance_db` above the electronics, so what is left is the difference.
        Under a few dB the electronic-noise subtraction dominates the result.
        """
        return self.clearance_db - abs(float(squeezing_db))


def _propagate(mode, distance_m):
    return mode.propagated(distance_m)


def geometry(chain, settings, detection=None, calibration=None, layout=None):
    """Beam sizes, separation, and the post-cell loss GABES should solve with.

    Distances come from `layout`, not from calibration coefficients: the
    table geometry and the physics are then the same edit, and each segment
    carries its own provenance instead of hiding behind a single number.
    """
    detection = detection or DetectionSettings()
    calibration = calibration or default_calibration()
    layout = layout or layout_parts.PART2
    c = calibration.value

    distance = layout_parts.cell_to_dmirror_m(layout)
    angle_rad = math.radians(settings.crossing_angle_deg)

    # The conjugate is born in the cell; phase matching ties its transverse mode
    # to the probe's, so both twins are propagated as the seed mode.
    twin_mode = chain.seed.mode
    pump_mode = chain.pump.mode

    twin_at_dmirror = _propagate(twin_mode, distance)
    pump_at_dmirror = _propagate(pump_mode, distance)

    pump_separation = distance * math.tan(angle_rad)
    twin_separation = 2.0 * pump_separation

    lens_distance = layout_parts.cell_to_lens_m(layout)
    arms = []
    for name, focal_mm in (("probe", detection.probe_lens_focal_mm),
                           ("conjugate", detection.conjugate_lens_focal_mm)):
        focal_m = focal_mm * 1e-3
        at_lens = _propagate(twin_mode, lens_distance)
        focused = at_lens.through_thin_lens(focal_m)
        # Place the diode at the post-lens waist, then apply the defocus knob.
        at_pd = focused.propagated(-focused.z_m + detection.pd_defocus_mm * 1e-3)
        arms.append(ArmGeometry(
            name=name,
            radius_at_dmirror_m=twin_at_dmirror.radius_m,
            spot_radius_m=at_pd.radius_m,
            lens_focal_m=focal_m,
        ))

    # Per-optic transmissions along the post-cell route, times the residual that
    # reconciles those nominal guesses with the one measured total.
    transmission = (layout_parts.transmission_of(layout_parts.ROUTE_POST_CELL,
                                                 calibration)
                    * c("loss_post_cell_residual"))
    qe_pct = _quantum_efficiency_pct(calibration)

    warnings = []
    margin = pump_separation / twin_at_dmirror.radius_m
    if margin < SEPARATION_MARGIN_WARN:
        warnings.append(
            f"pump sits {margin:.1f} twin radii away at the D-shaped mirrors; "
            f"below {SEPARATION_MARGIN_WARN:.0f} the beams are not cleanly "
            f"resolved, and the assumption that the pump is cleanly "
            f"rejected stops holding")

    return DetectionGeometry(
        arms=tuple(arms),
        twin_separation_m=twin_separation,
        pump_separation_m=pump_separation,
        pump_radius_at_dmirror_m=pump_at_dmirror.radius_m,
        twin_radius_at_dmirror_m=twin_at_dmirror.radius_m,
        optics_transmission=transmission,
        loss_pct=100.0 * (1.0 - transmission),
        qe_pct=qe_pct,
        warnings=tuple(warnings),
    )


def _quantum_efficiency_pct(calibration):
    """QE implied by the photodiode responsivity: QE = R * hc / (e * lambda)."""
    responsivity = calibration.value("pd_responsivity_a_per_w")
    return 100.0 * responsivity * 1239.841984 / 795.0


def readout(geom, chain, gains, settings, detection=None, calibration=None):
    """Photodiode powers, linearity, and shot-noise clearance.

    `gains` is (G_s, G_c) from the GABES solve -- the twins' power is the seed
    power times the parametric gain, so nothing here can be known before it.
    """
    detection = detection or DetectionSettings()
    calibration = calibration or default_calibration()
    c = calibration.value

    seed_w = chain.seed_power_w
    pump_out_w = chain.pump_power_w
    pbs_leak = 1.0 / c("pbs_extinction_ratio")
    density_limit = c("pd_max_power_density_w_per_cm2")
    active_radius = 0.5 * c("pd_active_diameter_mm") * 1e-3
    responsivity = c("pd_responsivity_a_per_w")
    nep = c("bpd_nep_w_per_rthz")
    saturation = c("bpd_cw_saturation_w")

    warnings = []
    arms = []
    total = 0.0
    # Pump leakage is measured at the detector, not derived: a Gaussian tail
    # claims 1e-9 rejection, which the paper contradicts, and the part that
    # actually matters is diffuse scatter no spatial filter removes.
    leakage_w = 1e-3 * 10.0 ** (detection.pump_leakage_dbm / 10.0)
    for arm, gain in zip(geom.arms, gains):
        power = seed_w * float(gain) * geom.optics_transmission
        residual = leakage_w
        density = 2.0 * power / (math.pi * arm.spot_radius_m ** 2) / 1e4
        total += power
        arms.append(ArmReadout(
            name=arm.name, power_w=power, residual_pump_w=residual,
            power_density_w_per_cm2=density, spot_radius_m=arm.spot_radius_m,
            fills_active_area=arm.spot_radius_m > 0.25 * active_radius,
        ))
        if density > density_limit:
            warnings.append(
                f"{arm.name}: {density:.1f} W/cm^2 on the diode exceeds the "
                f"{density_limit:.0f} W/cm^2 linearity limit -- use a longer "
                f"focal length or move the diode off focus")
        if power > saturation:
            warnings.append(
                f"{arm.name}: {power * 1e3:.1f} mW exceeds the "
                f"{saturation * 1e3:.0f} mW CW saturation of the amplifier")
        if arm.spot_radius_m > active_radius:
            warnings.append(
                f"{arm.name}: spot radius {arm.spot_radius_m * 1e6:.0f} um is "
                f"larger than the {active_radius * 1e6:.0f} um active radius")
        if residual > 0.05 * power:
            warnings.append(
                f"{arm.name}: residual pump is {100 * residual / power:.1f} % of "
                f"the signal; pump scatter will dominate the noise floor")

    bandwidth = 1.0
    shot = math.sqrt(2.0 * ELEMENTARY_CHARGE * responsivity * total * bandwidth)
    electronic = responsivity * nep
    clearance = 20.0 * math.log10(shot / electronic) if electronic > 0 else math.inf
    if clearance < CLEARANCE_WARN_DB:
        warnings.append(
            f"shot noise is only {clearance:.1f} dB above the amplifier noise; "
            f"a squeezed trace would sit close to the electronic floor. More "
            f"optical power fixes this and the atoms do not care about seed power")

    return DetectionReadout(
        arms=tuple(arms), total_power_w=total,
        shot_noise_a_per_rthz=shot, electronic_noise_a_per_rthz=electronic,
        clearance_db=clearance, warnings=tuple(warnings))


def power_for_clearance(target_db, calibration=None):
    """Total optical power giving a target shot-noise-to-electronic clearance."""
    calibration = calibration or default_calibration()
    responsivity = calibration.value("pd_responsivity_a_per_w")
    nep = calibration.value("bpd_nep_w_per_rthz")
    ratio_sq = 10.0 ** (float(target_db) / 10.0)
    return ratio_sq * responsivity * nep ** 2 / (2.0 * ELEMENTARY_CHARGE)


def focal_length_for_density(target_w_per_cm2, power_w, radius_at_lens_m,
                             wavelength_m):
    """Lens focal length whose focal spot keeps a diode inside its linearity limit.

    w_f = lambda f / (pi w_L), and peak density is 2P / (pi w_f^2), so the
    required focal length grows as sqrt(P / I_max).
    """
    density = float(target_w_per_cm2) * 1e4
    waist = math.sqrt(2.0 * float(power_w) / (math.pi * density))
    return math.pi * waist * float(radius_at_lens_m) / float(wavelength_m)
