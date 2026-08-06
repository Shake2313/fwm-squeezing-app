"""
Virtual lab instruments, driven by the real model.

Each takes a head and returns a `Reading`. Optical heads (power meter, beam
profiler, wavemeter, photodiode) accept a `BeamState`, which is whatever
`sabes.layout.beam_at` resolves at a point on the table. Electrical ones
(oscilloscope, spectrum analyser) take a `PhotocurrentSignal` from a photodiode,
exactly as on the bench.

Nothing here imports the layout or Streamlit, so the package stands alone: the
only dependency is the beam state and numpy.

`Reading.provenance` says whether a number follows from the model's physics
(`COMPUTED`) or re-expresses a PSD the model already knows (`SYNTHESISED`). Only
the oscilloscope's time-series mode and the spectrum analyser are the latter, and
they must be labelled wherever they are drawn.
"""
from .base import (COMPUTED, HEAD_OPTICAL, HEAD_PHOTOCURRENT,
                   PhotocurrentSignal, Quantity, Reading, SYNTHESISED, Trace,
                   quantum_efficiency, shot_noise_a_per_rthz)
from .electrical import Oscilloscope, Photodiode, SpectrumAnalyzer
from .optical import BeamProfiler, PowerMeter, Wavemeter

#: Instruments whose head is a point on a beam.
OPTICAL_INSTRUMENTS = (PowerMeter, BeamProfiler, Wavemeter, Photodiode)
#: Instruments that plug into a photodiode.
ELECTRICAL_INSTRUMENTS = (Oscilloscope, SpectrumAnalyzer)

__all__ = [
    "Reading", "Quantity", "Trace", "PhotocurrentSignal",
    "COMPUTED", "SYNTHESISED", "HEAD_OPTICAL", "HEAD_PHOTOCURRENT",
    "PowerMeter", "BeamProfiler", "Wavemeter",
    "Photodiode", "Oscilloscope", "SpectrumAnalyzer",
    "OPTICAL_INSTRUMENTS", "ELECTRICAL_INSTRUMENTS",
    "shot_noise_a_per_rthz", "quantum_efficiency",
]
