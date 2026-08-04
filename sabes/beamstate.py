"""
The state a beam carries between optics: spectrum, transverse mode, polarization.

SABES threads one immutable `BeamState` through an ordered list of devices. Three
things travel together because every optic in the Sim et al. setup touches at
least one of them:

  spectrum      EOM sidebands and the residual carrier. Etalons act here.
  mode          Gaussian w0 / M^2. Collimators, telescopes, free space act here.
  polarization  Jones vector. HWP / QWP / PBS / GTP act here; the fiber EOM's
                effective drive depends on it too.

Conventions
  - SI everywhere inside this module (W, m, Hz, rad). Lab units are converted at
    the settings boundary, not here.
  - Waists are 1/e^2 *radii*, matching the paper and `gabes.schemes.fwm`.
  - Spectral offsets are signed frequencies relative to the pump carrier, so the
    EOM -1 sideband sits at -f_EOM.
  - One Jones vector for the whole beam. Valid here because every line in a given
    path is born from the same input polarization; a device that would split the
    spectrum by polarization is not part of this setup.

M^2 handling uses the embedded-Gaussian trick: a real beam of quality M^2
propagates exactly like an ideal Gaussian at wavelength M^2*lambda, so every
propagation law below is the textbook one with `lambda_eff`.
"""
from dataclasses import dataclass, replace
from typing import Tuple

import math

# Jones vectors for the two linear basis states. "h" is the PBS transmission
# axis throughout SABES, so `Polarizer(0 deg)` and a PBS transmit port agree.
JONES_H = (1.0 + 0j, 0.0 + 0j)
JONES_V = (0.0 + 0j, 1.0 + 0j)


def normalize_jones(jones):
    """Unit-norm Jones vector. Power lives in the spectrum, not in this vector."""
    ex, ey = complex(jones[0]), complex(jones[1])
    norm = math.hypot(abs(ex), abs(ey))
    if norm <= 0.0:
        raise ValueError("Jones vector has zero norm")
    return (ex / norm, ey / norm)


def linear_jones(angle_deg):
    """Linear polarization at `angle_deg` from the PBS transmission axis."""
    a = math.radians(float(angle_deg))
    return (complex(math.cos(a)), complex(math.sin(a)))


def apply_jones_matrix(jones, matrix):
    """(transmitted power fraction, renormalized Jones vector) for a 2x2 matrix.

    Returning the fraction separately keeps `BeamState.jones` a pure direction:
    the loss it implies is applied to the spectrum, where power is tracked.
    """
    (m00, m01), (m10, m11) = matrix
    ex, ey = jones
    out = (m00 * ex + m01 * ey, m10 * ex + m11 * ey)
    fraction = abs(out[0]) ** 2 + abs(out[1]) ** 2
    if fraction <= 0.0:
        # Fully extinguished: keep the incoming direction so the state stays
        # well-formed; the caller multiplies the power by 0.
        return 0.0, jones
    return fraction, normalize_jones(out)


@dataclass(frozen=True)
class SpectralLine:
    """One optical frequency component.

    `offset_hz` is relative to the carrier (0.0 = the un-modulated laser line),
    so the seed the FWM wants is the -1 sideband at `-f_EOM`.
    """
    offset_hz: float
    power_w: float
    label: str = ""

    def scaled(self, factor):
        return replace(self, power_w=self.power_w * float(factor))


@dataclass(frozen=True)
class GaussianMode:
    """Gaussian beam with a beam-quality factor, at some plane.

    `z_m` is the signed distance from the waist to the current plane: positive
    means the waist is upstream (the beam is diverging).
    """
    waist_m: float
    wavelength_m: float
    z_m: float = 0.0
    m2: float = 1.0

    @property
    def wavelength_eff_m(self):
        """Embedded-Gaussian wavelength: M^2 * lambda."""
        return self.m2 * self.wavelength_m

    @property
    def rayleigh_m(self):
        return math.pi * self.waist_m ** 2 / self.wavelength_eff_m

    @property
    def radius_m(self):
        """1/e^2 radius at the current plane."""
        return self.waist_m * math.sqrt(1.0 + (self.z_m / self.rayleigh_m) ** 2)

    @property
    def divergence_rad(self):
        """Far-field half-angle M^2 * lambda / (pi * w0)."""
        return self.wavelength_eff_m / (math.pi * self.waist_m)

    @property
    def q(self):
        """Complex beam parameter q = z + i*z_R at the current plane."""
        return complex(self.z_m, self.rayleigh_m)

    def _from_q(self, q):
        rayleigh = q.imag
        if rayleigh <= 0.0:
            raise ValueError("beam parameter lost its Rayleigh range")
        waist = math.sqrt(rayleigh * self.wavelength_eff_m / math.pi)
        return replace(self, waist_m=waist, z_m=q.real)

    def propagated(self, distance_m):
        return replace(self, z_m=self.z_m + float(distance_m))

    def through_thin_lens(self, focal_m):
        """1/q_out = 1/q_in - 1/f. f <= 0 is allowed (diverging lens)."""
        focal_m = float(focal_m)
        if focal_m == 0.0:
            raise ValueError("focal length must be non-zero")
        return self._from_q(1.0 / (1.0 / self.q - 1.0 / focal_m))

    def with_m2(self, m2):
        """Degrade (or restore) beam quality while holding the physical radius.

        An amplifier does not change how wide the beam is at its output face; it
        changes how fast that width grows afterwards. Holding `radius_m` fixed
        and re-solving the waist is what expresses that.
        """
        m2 = float(m2)
        if m2 < 1.0:
            raise ValueError("M^2 cannot be below 1")
        radius = self.radius_m
        z_ratio = self.z_m / self.rayleigh_m
        # w = w0 sqrt(1 + r^2) with r = z/z_R held fixed -> same shape, new lambda_eff.
        waist = radius / math.sqrt(1.0 + z_ratio ** 2)
        rayleigh = math.pi * waist ** 2 / (m2 * self.wavelength_m)
        return replace(self, m2=m2, waist_m=waist, z_m=z_ratio * rayleigh)

    def clipping_transmission(self, aperture_radius_m, offset_m=0.0):
        """Power through a centred circular aperture, or past a knife edge.

        With `offset_m` = 0 this is the exact circular-aperture result
        1 - exp(-2 a^2 / w^2). With a non-zero offset it degrades to the 1-D
        knife-edge form, which is what an iris cutting a beam whose centre is
        `offset_m` away actually approximates. Good enough to answer "does the
        pump leak past the iris"; it is not a diffraction calculation.
        """
        w = self.radius_m
        a = max(float(aperture_radius_m), 0.0)
        offset = abs(float(offset_m))
        if offset == 0.0:
            return 1.0 - math.exp(-2.0 * a * a / (w * w))
        # Fraction of a Gaussian centred at `offset` that falls inside |x| < a.
        lo = (offset - a) * math.sqrt(2.0) / w
        hi = (offset + a) * math.sqrt(2.0) / w
        return 0.5 * (math.erfc(lo) - math.erfc(hi))


@dataclass(frozen=True)
class BeamState:
    """Spectrum + transverse mode + polarization at one plane of the setup."""
    lines: Tuple[SpectralLine, ...]
    mode: GaussianMode
    jones: Tuple[complex, complex] = JONES_H
    carrier_hz: float = 0.0

    # ---- power ----
    @property
    def total_power_w(self):
        return sum(line.power_w for line in self.lines)

    def power_at(self, offset_hz, tolerance_hz=1.0):
        """Total power within `tolerance_hz` of a spectral offset."""
        return sum(line.power_w for line in self.lines
                   if abs(line.offset_hz - float(offset_hz)) <= tolerance_hz)

    def line(self, label):
        for candidate in self.lines:
            if candidate.label == label:
                return candidate
        raise KeyError(f"no spectral line labelled {label!r}")

    def scaled(self, factor):
        """Uniform (wavelength-independent) attenuation or gain."""
        return replace(self, lines=tuple(l.scaled(factor) for l in self.lines))

    def with_lines(self, lines):
        return replace(self, lines=tuple(lines))

    def with_mode(self, mode):
        return replace(self, mode=mode)

    def with_jones(self, jones):
        return replace(self, jones=normalize_jones(jones))

    # ---- spectral purity ----
    def carrier_ratio(self, seed_offset_hz):
        """P(carrier) / P(seed line) -- the quantity Sim et al. Fig. 5 tracks.

        The paper reports << 0.5 % with three etalons and ~0.5 % with two, and
        shows squeezing degrading as it rises (still visible at 25 %).
        """
        seed = self.power_at(seed_offset_hz)
        if seed <= 0.0:
            return float("inf")
        return self.power_at(0.0) / seed

    def purity(self, seed_offset_hz):
        """Fraction of the total power that sits in the wanted seed line."""
        total = self.total_power_w
        if total <= 0.0:
            return 0.0
        return self.power_at(seed_offset_hz) / total

    def dominant_line(self):
        return max(self.lines, key=lambda l: l.power_w)


def monochromatic(power_w, wavelength_m, waist_m, *, m2=1.0, z_m=0.0,
                  jones=JONES_H, carrier_hz=0.0, label="carrier"):
    """A single-line beam -- what a laser emits before any modulation."""
    return BeamState(
        lines=(SpectralLine(0.0, float(power_w), label),),
        mode=GaussianMode(waist_m=float(waist_m), wavelength_m=float(wavelength_m),
                          z_m=float(z_m), m2=float(m2)),
        jones=normalize_jones(jones),
        carrier_hz=float(carrier_hz),
    )
