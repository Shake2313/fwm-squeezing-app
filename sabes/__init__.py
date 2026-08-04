"""
SABES -- Specific Atomic Bloch Equation Solver.

A lab-facing layer over GABES for one experiment: the EOM-based single-diode-laser
intensity-difference squeezing setup of G. Sim, H. Kim and H. S. Moon,
Sci. Rep. 15, 7727 (2025).

Where GABES takes already-processed physics parameters (one-photon detuning,
Rabi frequency, cell temperature), SABES takes the settings someone turns --
signal-generator frequency, half-wave-plate angles, amplifier current, lens
choices -- and derives the GABES parameters from them through explicit device
models. Every coefficient in that derivation carries a provenance tag, so a
result can honestly be labelled predicted or reproduced.

Dependency direction is one-way and must stay that way: `sabes` imports `gabes`,
never the reverse.
"""
from .beamstate import BeamState, GaussianMode, SpectralLine, monochromatic
from .calibration import Calibration, Coefficient, default_calibration

__all__ = [
    "BeamState", "GaussianMode", "SpectralLine", "monochromatic",
    "Calibration", "Coefficient", "default_calibration",
]
