"""
What a virtual instrument is, and what it is allowed to claim.

Instruments take a `BeamState` (or, for the electrical ones, a
`PhotocurrentSignal`) and return a `Reading`. Nothing here knows about the table
layout or about Streamlit, so the package is reusable outside SABES: the only
import is the beam state itself.

The load-bearing field is `Reading.provenance`.

    COMPUTED      the numbers follow from the model's own physics. A power
                  meter adding up spectral lines is computed; so is a beam
                  profiler propagating a Gaussian.
    SYNTHESISED   the trace is a re-expression of a power spectral density the
                  model already knows, not a simulation that produced new
                  information. A noise time series drawn from a PSD looks like
                  instrument output and is statistically right, but it tells you
                  nothing the PSD did not already say.

The distinction matters because a plot that looks like an oscilloscope invites
being read as evidence. Anything SYNTHESISED must carry that label on screen.
"""
from dataclasses import dataclass
from typing import Dict, Optional, Tuple

import numpy as np

#: A reading whose numbers come from the model's physics.
COMPUTED = "computed"
#: A reading constructed from a PSD the model already knows. Not new evidence.
SYNTHESISED = "synthesised"

#: What an instrument's head accepts.
HEAD_OPTICAL = "optical"
HEAD_PHOTOCURRENT = "photocurrent"

ELEMENTARY_CHARGE = 1.602176634e-19
PLANCK_EV_NM = 1239.841984


@dataclass(frozen=True)
class Quantity:
    """One scalar readout, with the unit it is quoted in."""
    label: str
    value: float
    unit: str = ""
    note: str = ""

    def formatted(self, digits=4):
        if not np.isfinite(self.value):
            return "—"
        return f"{self.value:.{digits}g}{(' ' + self.unit) if self.unit else ''}"


@dataclass(frozen=True)
class Trace:
    """One plot: a shared x axis and any number of named series."""
    x: np.ndarray
    series: Dict[str, np.ndarray]
    x_label: str = ""
    y_label: str = ""
    x_unit: str = ""
    y_unit: str = ""
    kind: str = "line"          # line | stem

    def __post_init__(self):
        for name, values in self.series.items():
            if len(values) != len(self.x):
                raise ValueError(
                    f"series {name!r} has {len(values)} points against "
                    f"{len(self.x)} on the x axis")


@dataclass(frozen=True)
class Reading:
    """What an instrument reports."""
    instrument: str
    quantities: Tuple[Quantity, ...] = ()
    trace: Optional[Trace] = None
    warnings: Tuple[str, ...] = ()
    provenance: str = COMPUTED
    note: str = ""

    def __post_init__(self):
        if self.provenance not in (COMPUTED, SYNTHESISED):
            raise ValueError(f"unknown provenance {self.provenance!r}")

    @property
    def synthesised(self):
        return self.provenance == SYNTHESISED

    def quantity(self, label) -> Quantity:
        for item in self.quantities:
            if item.label == label:
                return item
        raise KeyError(f"{self.instrument}: no quantity {label!r}")


@dataclass(frozen=True)
class PhotocurrentSignal:
    """A photodiode's output: the head the electrical instruments plug into.

    Carries the noise terms rather than a single current, because everything the
    oscilloscope and spectrum analyser do downstream is about noise, and
    recomputing it from a bare current would need the responsivity again.
    """
    optical_power_w: float
    responsivity_a_per_w: float
    dc_current_a: float
    shot_noise_a_per_rthz: float
    electronic_noise_a_per_rthz: float
    transimpedance_v_per_a: float
    bandwidth_hz: float
    label: str = "photodiode"

    @property
    def dc_voltage_v(self):
        return self.dc_current_a * self.transimpedance_v_per_a

    @property
    def clearance_db(self):
        """Shot noise above the amplifier floor, in noise power."""
        if self.electronic_noise_a_per_rthz <= 0:
            return float("inf")
        return 20.0 * np.log10(self.shot_noise_a_per_rthz
                               / self.electronic_noise_a_per_rthz)

    def total_noise_a_per_rthz(self, excess_a_per_rthz=0.0):
        """Quadrature sum of shot, electronic and any excess term."""
        return float(np.sqrt(self.shot_noise_a_per_rthz ** 2
                             + self.electronic_noise_a_per_rthz ** 2
                             + float(excess_a_per_rthz) ** 2))


def shot_noise_a_per_rthz(responsivity, optical_power_w):
    """sqrt(2 e I) for a photocurrent I = R * P, per root hertz."""
    current = max(float(responsivity) * float(optical_power_w), 0.0)
    return float(np.sqrt(2.0 * ELEMENTARY_CHARGE * current))


def quantum_efficiency(responsivity_a_per_w, wavelength_nm):
    """QE = R * hc / (e * lambda). The number that sets the squeezing floor."""
    return (float(responsivity_a_per_w) * PLANCK_EV_NM
            / float(wavelength_nm))


def to_dbm(power_w):
    """Power in dBm, floored so an empty trace does not become -inf."""
    watts = np.maximum(np.asarray(power_w, dtype=float), 1e-30)
    return 10.0 * np.log10(watts / 1e-3)
