"""
The instruments whose head is a beam: power meter, beam profiler, wavemeter.

All three read a `BeamState` and nothing else, so any point the layout can
resolve is a place they can be put. They are `COMPUTED` throughout -- these are
the model's own quantities reported in lab units, with no reconstruction.

A wavemeter reports the dominant optical frequency or wavelength. Spectral-line
powers remain available from the power-meter view.
"""
from dataclasses import dataclass

import numpy as np

from .base import COMPUTED, HEAD_OPTICAL, Quantity, Reading, Trace


@dataclass(frozen=True)
class PowerMeter:
    """Total optical power, and where it sits in the spectrum.

    A real power meter integrates everything that lands on it, which is exactly
    why the seed's residual carrier is invisible to one: the number it reports
    is `total`, not `wanted`. Reporting both is the point.
    """
    head = HEAD_OPTICAL
    name: str = "Power meter"
    wanted_offset_hz: float = None

    def measure(self, beam) -> Reading:
        total = beam.total_power_w
        quantities = [Quantity("Total power", total * 1e3, "mW")]

        wanted = None
        if self.wanted_offset_hz is not None:
            wanted = beam.power_at(self.wanted_offset_hz)
            quantities += [
                Quantity("Wanted line", wanted * 1e6, "µW"),
                Quantity("In the wanted line", 100.0 * wanted / total if total
                         else float("nan"), "%",
                         "A power meter cannot tell these apart; it reads the "
                         "total."),
            ]

        lines = sorted(beam.lines, key=lambda l: l.offset_hz)
        trace = Trace(
            x=np.array([l.offset_hz / 1e9 for l in lines]),
            series={"power": np.array([l.power_w * 1e6 for l in lines])},
            x_label="Offset from carrier", x_unit="GHz",
            y_label="Power", y_unit="µW", kind="stem",
        )
        warnings = []
        if wanted is not None and total > 0 and wanted / total < 0.9:
            warnings.append(
                f"only {100 * wanted / total:.1f} % of the power here is in the "
                f"wanted line; a power-meter reading overstates the useful beam")
        return Reading(self.name, tuple(quantities), trace, tuple(warnings),
                       COMPUTED)


@dataclass(frozen=True)
class BeamProfiler:
    """Transverse size and quality at the probed plane.

    Reports the 1/e^2 radius where it is put, not the waist: those differ by a
    factor of two over a metre of table, and confusing them is how a beam gets
    clipped on something that looked comfortable in a design note.
    """
    head = HEAD_OPTICAL
    name: str = "Beam profiler"
    samples: int = 121

    def measure(self, beam) -> Reading:
        mode = beam.mode
        radius = mode.radius_m
        power = beam.total_power_w
        peak_intensity = (2.0 * power / (np.pi * radius ** 2)) if radius > 0 else 0.0

        quantities = (
            Quantity("1/e² radius here", radius * 1e6, "µm"),
            Quantity("Waist w₀", mode.waist_m * 1e6, "µm"),
            Quantity("Distance from waist", mode.z_m * 1e3, "mm"),
            Quantity("Rayleigh range", mode.rayleigh_m * 1e3, "mm"),
            Quantity("Divergence (half-angle)", mode.divergence_rad * 1e3, "mrad"),
            Quantity("M²", mode.m2, ""),
            Quantity("Peak intensity", peak_intensity / 1e4, "W/cm²"),
        )

        span = 2.5 * radius
        x = np.linspace(-span, span, self.samples)
        profile = np.exp(-2.0 * (x / radius) ** 2) if radius > 0 else np.zeros_like(x)
        trace = Trace(
            x=x * 1e6, series={"intensity": profile},
            x_label="Transverse position", x_unit="µm",
            y_label="Normalised intensity", y_unit="",
        )
        return Reading(self.name, quantities, trace, (), COMPUTED)


@dataclass(frozen=True)
class Wavemeter:
    """Dominant optical frequency and vacuum wavelength."""
    head = HEAD_OPTICAL
    name: str = "Wavemeter"
    wanted_offset_hz: float = None

    def measure(self, beam) -> Reading:
        powered = [line for line in beam.lines if line.power_w > 0.0]
        if not powered:
            quantities = (
                Quantity("Optical frequency", float("nan"), "THz"),
                Quantity("Vacuum wavelength", float("nan"), "nm"),
            )
            return Reading(self.name, quantities, None,
                           ("no optical power reaches the wavemeter",), COMPUTED)

        speed_of_light = 299_792_458.0
        largest_power = max(line.power_w for line in powered)
        dominant_lines = [
            line for line in powered
            if np.isclose(line.power_w, largest_power, rtol=1.0e-9, atol=0.0)
        ]
        if len(dominant_lines) != 1:
            quantities = (
                Quantity("Optical frequency", float("nan"), "THz"),
                Quantity("Vacuum wavelength", float("nan"), "nm"),
                Quantity("Dominant-line offset", float("nan"), "GHz"),
            )
            return Reading(
                self.name, quantities, None,
                ("multiple optical lines have the same largest power",), COMPUTED)

        dominant = dominant_lines[0]
        carrier_hz = (beam.carrier_hz if beam.carrier_hz > 0.0 else
                      speed_of_light / beam.mode.wavelength_m)
        frequency_hz = carrier_hz + dominant.offset_hz
        quantities = [
            Quantity("Optical frequency", frequency_hz / 1e12, "THz"),
            Quantity("Vacuum wavelength", speed_of_light / frequency_hz * 1e9,
                     "nm"),
            Quantity("Dominant-line offset", dominant.offset_hz / 1e9, "GHz"),
        ]
        warnings = []
        if self.wanted_offset_hz is not None:
            error_hz = dominant.offset_hz - float(self.wanted_offset_hz)
            quantities.append(Quantity("Target offset error", error_hz / 1e6,
                                       "MHz"))
            if abs(error_hz) > 1.0:
                warnings.append("the dominant line is not the selected seed line")
        return Reading(self.name, tuple(quantities), None, tuple(warnings),
                       COMPUTED)
