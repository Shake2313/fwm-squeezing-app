"""
Optic models. Each device maps a `BeamState` to a `BeamState`.

Everything here is a closed-form device law with a small number of coefficients.
Coefficients that cannot be predicted (alignment efficiency, etalon leakage, the
tapered-amplifier P-I curve) are inputs, not guesses baked into the code -- see
`sabes.calibration` for how they carry their provenance.

Scope note: these are steady-state power/mode/polarization laws. No thermal
transients, no feedback loops, no phase noise. SABES answers "what arrives at the
cell, and what arrives at the detector", not "how does the lock behave".
"""
from dataclasses import dataclass

import math

from .beamstate import SpectralLine, apply_jones_matrix, linear_jones

# ---------------------------------------------------------------------------
# Bessel functions of the first kind, integer order.
# scipy is not a dependency of this repo, and the phase-modulator sideband
# ladder needs integer orders up to about |x| + 10 over the UI's full RF range;
# the ascending series still converges rapidly there.
# ---------------------------------------------------------------------------


def bessel_j(order, x, terms=60):
    """J_n(x) for integer n by the ascending series, via stable term recursion."""
    n = int(order)
    x = float(x)
    if n < 0:
        return (-1.0) ** (-n) * bessel_j(-n, x, terms)
    if x == 0.0:
        return 1.0 if n == 0 else 0.0
    half = 0.5 * x
    term = half ** n / math.factorial(n)
    total = term
    for k in range(terms):
        term *= -(half * half) / ((k + 1.0) * (k + 1.0 + n))
        total += term
        if abs(term) < 1e-18 * max(abs(total), 1e-30):
            break
    return total


def dbm_to_watt(dbm):
    return 10.0 ** (float(dbm) / 10.0) * 1e-3


def rf_peak_volts(dbm, impedance_ohm=50.0):
    """Peak voltage of a sinusoid delivering `dbm` into `impedance_ohm`."""
    return math.sqrt(2.0 * float(impedance_ohm) * dbm_to_watt(dbm))


def db_to_fraction(db):
    """Transmission fraction of a positive insertion loss in dB."""
    return 10.0 ** (-abs(float(db)) / 10.0)


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class FreeSpace:
    """Propagation. The only device that makes divergence visible."""
    distance_m: float

    def apply(self, beam):
        return beam.with_mode(beam.mode.propagated(self.distance_m))


@dataclass(frozen=True)
class ThinLens:
    focal_m: float

    def apply(self, beam):
        return beam.with_mode(beam.mode.through_thin_lens(self.focal_m))


@dataclass(frozen=True)
class Telescope:
    """Afocal Keplerian pair, waist to waist: magnification exactly f2 / f1.

    This is the device that sets the waist in the cell. In the Sim et al. setup
    the pump is a 2.0 mm collimator through f300 -> f150 and the seed a 1.6 mm
    collimator through f400 -> f150, which is where 530 / 330 um come from.

    The input waist is taken to sit at L1's front focal plane and the output
    waist lands on L2's back focal plane, which is the arrangement that makes the
    magnification purely geometric. Putting the waist at the lens instead would
    leave the output slightly converging, so the waist depends on the cell
    position. Use the designed arrangement until those distances are measured.
    """
    f1_m: float
    f2_m: float

    def apply(self, beam):
        stage = FreeSpace(self.f1_m).apply(beam)
        stage = ThinLens(self.f1_m).apply(stage)
        stage = FreeSpace(self.f1_m + self.f2_m).apply(stage)
        stage = ThinLens(self.f2_m).apply(stage)
        return FreeSpace(self.f2_m).apply(stage)

    @property
    def magnification(self):
        return self.f2_m / self.f1_m


@dataclass(frozen=True)
class FiberCollimator:
    """Fiber output collimated to a stated 1/e^2 diameter.

    Also resets the mode origin: the collimator's output waist is the reference
    plane for everything downstream.
    """
    output_diameter_m: float
    m2: float = 1.0

    def apply(self, beam):
        mode = beam.mode
        return beam.with_mode(type(mode)(
            waist_m=0.5 * float(self.output_diameter_m),
            wavelength_m=mode.wavelength_m,
            z_m=0.0,
            m2=float(self.m2)))


@dataclass(frozen=True)
class WaistRescale:
    """Multiply the waist by a calibration factor, holding z/z_R.

    Not an optic. It carries the one shared factor by which the ideal
    collimator-times-magnification chain misses the measured cell waists, so that
    discrepancy stays visible as a named coefficient instead of being absorbed
    into a lens focal length that does not match the one on the table.
    """
    factor: float = 1.0

    def apply(self, beam):
        mode = beam.mode
        ratio = mode.z_m / mode.rayleigh_m
        waist = mode.waist_m * float(self.factor)
        rayleigh = math.pi * waist ** 2 / mode.wavelength_eff_m
        return beam.with_mode(type(mode)(
            waist_m=waist, wavelength_m=mode.wavelength_m,
            z_m=ratio * rayleigh, m2=mode.m2))


@dataclass(frozen=True)
class Aperture:
    """Iris or knife edge. `offset_m` is the beam centre's distance from centre."""
    radius_m: float
    offset_m: float = 0.0

    def apply(self, beam):
        return beam.scaled(
            beam.mode.clipping_transmission(self.radius_m, self.offset_m))


@dataclass(frozen=True)
class LossElement:
    """Lumped alignment / surface / handling loss for one stage.

    Every stated loss in SABES is one of these, so the budget table can show
    where the power went instead of hiding it in a device's efficiency.
    """
    transmission: float = 1.0
    label: str = ""

    @classmethod
    def from_db(cls, db, label=""):
        return cls(db_to_fraction(db), label)

    def apply(self, beam):
        return beam.scaled(self.transmission)


# ---------------------------------------------------------------------------
# Polarization
# ---------------------------------------------------------------------------


def _conjugate_by_rotation(diagonal, angle_rad):
    """R(theta) . diag(a, b) . R(-theta) as a nested tuple."""
    a, b = diagonal
    c, s = math.cos(angle_rad), math.sin(angle_rad)
    return ((a * c * c + b * s * s, (a - b) * c * s),
            ((a - b) * c * s, a * s * s + b * c * c))


@dataclass(frozen=True)
class HalfWavePlate:
    """Rotates linear polarization by 2*(angle - zero). `zero_deg` is mount zero."""
    angle_deg: float
    zero_deg: float = 0.0

    def apply(self, beam):
        a = math.radians(self.angle_deg - self.zero_deg)
        c, s = math.cos(2.0 * a), math.sin(2.0 * a)
        fraction, jones = apply_jones_matrix(beam.jones, ((c, s), (s, -c)))
        return beam.scaled(fraction).with_jones(jones)


@dataclass(frozen=True)
class QuarterWavePlate:
    angle_deg: float
    zero_deg: float = 0.0

    def apply(self, beam):
        a = math.radians(self.angle_deg - self.zero_deg)
        matrix = _conjugate_by_rotation((1.0 + 0j, 1j), a)
        fraction, jones = apply_jones_matrix(beam.jones, matrix)
        return beam.scaled(fraction).with_jones(jones)


@dataclass(frozen=True)
class Polarizer:
    """Glan-Taylor / sheet polarizer. `extinction_ratio` is power, e.g. 1e5."""
    angle_deg: float = 0.0
    extinction_ratio: float = 1e5
    insertion_loss: float = 0.0

    def apply(self, beam):
        leak_amp = math.sqrt(1.0 / max(float(self.extinction_ratio), 1.0))
        matrix = _conjugate_by_rotation((1.0 + 0j, complex(leak_amp)),
                                        math.radians(self.angle_deg))
        fraction, jones = apply_jones_matrix(beam.jones, matrix)
        return beam.scaled(fraction * (1.0 - self.insertion_loss)).with_jones(jones)


@dataclass(frozen=True)
class PolarizingBeamSplitter:
    """Two-port PBS. `split` returns (transmitted, reflected).

    The finite-extinction leakage is transferred to the wrong port, so the two
    output powers sum to the input power after insertion loss.
    """
    extinction_ratio: float = 3000.0
    insertion_loss: float = 0.0

    def split(self, beam):
        loss = float(self.insertion_loss)
        if not 0.0 <= loss <= 1.0:
            raise ValueError("insertion_loss must lie between zero and one")
        ratio = max(float(self.extinction_ratio), 1.0)
        main = ratio / (ratio + 1.0)
        leak = 1.0 / (ratio + 1.0)
        main_amp, leak_amp = math.sqrt(main), math.sqrt(leak)
        t_fraction, t_jones = apply_jones_matrix(
            beam.jones, ((main_amp, 0j), (0j, leak_amp)))
        r_fraction, r_jones = apply_jones_matrix(
            beam.jones, ((leak_amp, 0j), (0j, main_amp)))
        transmitted = beam.scaled(t_fraction * (1.0 - loss)).with_jones(t_jones)
        reflected = beam.scaled(r_fraction * (1.0 - loss)).with_jones(r_jones)
        return transmitted, reflected

    def apply(self, beam):
        return self.split(beam)[0]


# ---------------------------------------------------------------------------
# Active devices
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class TaperedAmplifier:
    """Saturable MOPA stage.

        P_out = P_max(I) * P_in / (P_in + P_sat),   P_max = slope * (I - I_th)

    Three coefficients, all fixable from one measured P-I curve at a known seed
    power. Saturation is what makes the output insensitive to seed drift, which
    is why a bare linear gain would misreport the pump budget.

    The amplifier also degrades beam quality (`m2_out`) and adds broadband ASE.
    ASE is removed from the coherent output rather than added as a fake spectral
    line: it is not at any FWM-relevant frequency, and in the seed path the
    etalons erase it. `ase_fraction` therefore only reduces useful power.
    """
    current_a: float
    threshold_a: float
    slope_w_per_a: float
    saturation_w: float
    m2_out: float = 1.0
    ase_fraction: float = 0.0
    input_axis_deg: float = 0.0
    output_axis_deg: float = 0.0
    input_extinction_ratio: float = 1e4

    @property
    def max_power_w(self):
        return self.slope_w_per_a * max(self.current_a - self.threshold_a, 0.0)

    def gain_output_w(self, input_w):
        input_w = max(float(input_w), 0.0)
        if input_w <= 0.0:
            return 0.0
        return self.max_power_w * input_w / (input_w + self.saturation_w)

    def apply(self, beam):
        # Only the component along the gain axis is amplified; the HWP in front
        # of the amplifier exists to put the light there.
        coupled = Polarizer(self.input_axis_deg,
                            self.input_extinction_ratio).apply(beam)
        p_in = coupled.total_power_w
        p_out = self.gain_output_w(p_in) * (1.0 - self.ase_fraction)
        scale = 0.0 if p_in <= 0.0 else p_out / p_in
        out = coupled.scaled(scale)
        out = out.with_jones(linear_jones(self.output_axis_deg))
        return out.with_mode(out.mode.with_m2(self.m2_out))


@dataclass(frozen=True)
class PhaseModulator:
    """Fiber phase EOM: carrier -> Bessel sideband ladder at n * frequency_hz.

        beta = pi * V_peak / V_pi,   P_n = P_in * J_n(beta)^2

    The -1 sideband is the FWM seed; the residual carrier (J_0) is the undesired
    mode the etalons must kill. Because J_0 and J_1 are comparable at the drive
    that maximizes J_1, the carrier leaving the modulator is *not* small -- the
    filter chain, not the EOM, is what makes the seed pure.

    Polarization: a single-polarization LiNbO3 waveguide modulates only light on
    its axis. `orthogonal_leakage` is the fraction of the mismatched component
    that still reaches the output as *unmodulated carrier*; the EOSPACE default
    of 0 means the waveguide rejects it entirely (it becomes coupling loss, not
    carrier contamination).
    """
    frequency_hz: float
    v_pi_v: float
    rf_power_dbm: float
    impedance_ohm: float = 50.0
    insertion_loss_db: float = 0.0
    max_input_w: float = float("inf")
    axis_deg: float = 0.0
    orthogonal_leakage: float = 0.0
    # ``None`` sums Bessel sideband power adaptively.  The old fixed order
    # of six silently discarded 17 % of the optical power at the UI's 30 dBm
    # limit (beta ~= 7.1).
    max_order: int | None = None
    sideband_tail_tolerance: float = 1.0e-12

    def __post_init__(self):
        if self.frequency_hz <= 0.0:
            raise ValueError("frequency_hz must be positive")
        if self.v_pi_v <= 0.0 or self.impedance_ohm <= 0.0:
            raise ValueError("v_pi_v and impedance_ohm must be positive")
        if self.max_order is not None and self.max_order < 0:
            raise ValueError("max_order must be non-negative or None")
        if not 0.0 < self.sideband_tail_tolerance < 1.0:
            raise ValueError("sideband_tail_tolerance must lie between 0 and 1")

    @property
    def modulation_index(self):
        return math.pi * rf_peak_volts(self.rf_power_dbm,
                                       self.impedance_ohm) / self.v_pi_v

    @property
    def resolved_max_order(self):
        """Sideband order that closes ``sum J_n(beta)^2 = 1`` to tolerance."""
        if self.max_order is not None:
            return int(self.max_order)
        beta = self.modulation_index
        represented = bessel_j(0, beta) ** 2
        for order in range(1, 65):
            represented += 2.0 * bessel_j(order, beta) ** 2
            if max(1.0 - represented, 0.0) <= self.sideband_tail_tolerance:
                return order
        raise RuntimeError(
            "EOM sideband ladder did not converge by order 64; reduce the RF "
            "drive or supply an explicit max_order")

    def sideband_fractions(self):
        """``{order: J_n(beta)^2}`` with a power-closed adaptive ladder."""
        beta = self.modulation_index
        max_order = self.resolved_max_order
        return {n: bessel_j(n, beta) ** 2
                for n in range(-max_order, max_order + 1)}

    def apply(self, beam):
        axis = linear_jones(self.axis_deg)
        aligned = abs(beam.jones[0] * axis[0].conjugate()
                      + beam.jones[1] * axis[1].conjugate()) ** 2
        p_in = beam.total_power_w
        transmission = db_to_fraction(self.insertion_loss_db)

        modulated_w = p_in * aligned * transmission
        leaked_w = p_in * (1.0 - aligned) * transmission * self.orthogonal_leakage

        lines = []
        for order, fraction in sorted(self.sideband_fractions().items()):
            power = modulated_w * fraction
            if order == 0:
                power += leaked_w
            lines.append(SpectralLine(order * self.frequency_hz, power,
                                      _sideband_label(order)))
        return beam.with_lines(lines).with_jones(axis)

    def overdriven(self, beam):
        """True when the input exceeds the modulator's rated optical power."""
        return beam.total_power_w > self.max_input_w


def _sideband_label(order):
    if order == 0:
        return "carrier"
    return f"sideband{order:+d}"


@dataclass(frozen=True)
class EtalonFilter:
    """Fabry-Perot filter, Airy transmission plus a lumped leakage floor.

        T(nu) = (T_peak - leak) / (1 + F sin^2(pi (nu - nu_0) / FSR)) + leak

    written so that T(peak) = T_peak exactly and T(stopband) -> leak: the floor
    is a stopband pedestal, not extra throughput.

    `leak` is the part an ideal Airy cascade cannot explain: surface scatter,
    the finite angular spread of a real beam washing out the resonance dip,
    alignment, and the noise floor of whatever measured the suppression. It
    matters because without it three cascaded ideal etalons predict a carrier
    suppression far below what Sim et al. report, which would overestimate the
    achievable squeezing. See `sabes.calibration` for the fitted default and its
    provenance.

    Tuning: the passband centre is moved by tilt and by temperature. Both are
    lumped into `center_offset_hz`, which `tune_offset_hz` computes from lab
    knob settings when calibration coefficients are available.
    """
    fsr_hz: float
    fwhm_hz: float
    peak_transmission: float = 0.85
    center_offset_hz: float = 0.0
    leak: float = 0.0

    def __post_init__(self):
        if self.fsr_hz <= 0.0 or not 0.0 < self.fwhm_hz < self.fsr_hz:
            raise ValueError("require 0 < fwhm_hz < fsr_hz")
        if not 0.0 <= self.leak <= self.peak_transmission <= 1.0:
            raise ValueError("require 0 <= leak <= peak_transmission <= 1")

    @property
    def finesse(self):
        return self.fsr_hz / self.fwhm_hz

    @property
    def coefficient_of_finesse(self):
        # Exact finite-FSR coefficient: at +/-FWHM/2 the Airy contribution is
        # one half.  (2 F/pi)^2 is only its narrow-line approximation.
        angle = math.pi * self.fwhm_hz / (2.0 * self.fsr_hz)
        return 1.0 / math.sin(angle) ** 2

    def airy_shape(self, offset_hz):
        """Normalized Airy profile: 1 on resonance, ~1/F deep in the stopband."""
        detune = (float(offset_hz) - self.center_offset_hz) / self.fsr_hz
        return 1.0 / (1.0 + self.coefficient_of_finesse
                      * math.sin(math.pi * detune) ** 2)

    def transmission(self, offset_hz):
        return ((self.peak_transmission - self.leak) * self.airy_shape(offset_hz)
                + self.leak)

    def field_amplitude_transmission(self, offset_hz):
        """Calibration has no phase data; return ``sqrt(T)`` only."""
        return math.sqrt(max(self.transmission(offset_hz), 0.0))

    def apply(self, beam):
        return beam.with_lines(
            line.scaled(self.transmission(line.offset_hz)) for line in beam.lines)


def leak_from_carrier_ratio(etalon, stages, offset_hz, input_ratio, target_ratio):
    """Leakage floor per stage that reproduces a measured carrier:seed ratio.

    Solves `input_ratio * (T(offset) / T_peak)^stages = target_ratio` for the
    stopband floor, assuming the seed sits on the passband peak. This is the hook
    that turns a single measured suppression number into a model coefficient.
    """
    if stages <= 0 or target_ratio <= 0.0 or input_ratio <= 0.0:
        raise ValueError("stages, input_ratio and target_ratio must be positive")
    peak = etalon.peak_transmission
    # wanted = (peak - leak) * shape + leak  ->  solve the linear equation.
    wanted = (target_ratio / input_ratio) ** (1.0 / stages) * peak
    shape = etalon.airy_shape(offset_hz)
    denominator = 1.0 - shape
    if denominator <= 0.0:
        raise ValueError("offset_hz sits on the passband; no floor can be solved")
    return max((wanted - peak * shape) / denominator, 0.0)


@dataclass(frozen=True)
class FiberCoupling:
    """Free space -> single-mode fiber.

    Two separate effects, kept separate because only one of them is calibratable:

      alignment  a measured efficiency, the dominant term in practice;
      mode       a beam-quality penalty 4 M^2 / (1 + M^2)^2, the standard
                 Gaussian-overlap form (1 at M^2 = 1, falling monotonically).

    The fiber is also a spatial filter: whatever couples in emerges with
    M^2 ~ 1. That is why the tapered amplifier's degraded beam quality costs
    power here but does not propagate downstream.
    """
    alignment_efficiency: float = 1.0
    apply_mode_penalty: bool = True
    output_m2: float = 1.0

    @staticmethod
    def mode_penalty(m2):
        m2 = max(float(m2), 1.0)
        return 4.0 * m2 / (1.0 + m2) ** 2

    def apply(self, beam):
        efficiency = float(self.alignment_efficiency)
        if self.apply_mode_penalty:
            efficiency *= self.mode_penalty(beam.mode.m2)
        out = beam.scaled(efficiency)
        return out.with_mode(out.mode.with_m2(self.output_m2))
