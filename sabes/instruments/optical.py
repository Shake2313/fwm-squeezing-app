"""
The instruments whose head is a beam: power meter, beam profiler, wavemeter.

All three read a `BeamState` and nothing else, so any point the layout can
resolve is a place they can be put. They are `COMPUTED` throughout -- these are
the model's own quantities reported in lab units, with no reconstruction.

The wavemeter is the one to reach for first when explaining SABES. Dropped either
side of the etalon chain it shows the carrier and the unwanted sidebands being
killed, which is the entire reason the seed path has three filters in it.
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
    """The spectrum as the etalons see it.

    Frequencies are reported relative to the carrier, which is the axis the whole
    seed path is designed around: the wanted sideband sits at -f_EOM, and every
    other line is something the filters have to remove.
    """
    head = HEAD_OPTICAL
    name: str = "Wavemeter"
    wanted_offset_hz: float = None
    floor_db: float = -90.0

    def measure(self, beam) -> Reading:
        lines = sorted(beam.lines, key=lambda l: l.offset_hz)
        total = beam.total_power_w
        strongest = max((l.power_w for l in lines), default=0.0)

        quantities = [
            Quantity("Lines present", float(len(lines)), ""),
            Quantity("Total power", total * 1e3, "mW"),
        ]
        warnings = []
        if self.wanted_offset_hz is not None:
            wanted = beam.power_at(self.wanted_offset_hz)
            carrier = beam.power_at(0.0)
            ratio = (carrier / wanted) if wanted > 0 else float("inf")
            quantities += [
                Quantity("Wanted sideband", wanted * 1e6, "µW"),
                Quantity("Residual carrier", carrier * 1e6, "µW"),
                Quantity("Carrier : wanted", 100.0 * ratio, "%",
                         "Sim et al. report << 0.5 % with three etalons; "
                         "squeezing degrades as it rises and survives to 25 %."),
            ]
            if ratio > 0.25:
                warnings.append(
                    f"carrier is {100 * ratio:.1f} % of the wanted sideband — "
                    f"past the point where squeezing was still observed")

        # dB relative to the strongest line is how a spectrum analyser trace is
        # read, and it is the only way the suppressed lines are visible at all.
        relative = np.array([
            10.0 * np.log10(max(l.power_w, 1e-30) / max(strongest, 1e-30))
            for l in lines])
        relative = np.maximum(relative, self.floor_db)
        trace = Trace(
            x=np.array([l.offset_hz / 1e9 for l in lines]),
            series={"relative": relative},
            x_label="Offset from carrier", x_unit="GHz",
            y_label="Relative power", y_unit="dB", kind="stem",
        )
        return Reading(self.name, tuple(quantities), trace, tuple(warnings),
                       COMPUTED)
